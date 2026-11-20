#pragma once

#include "mu_minus_m1_conversion.hxx"

namespace stopped_muon {

// F-block (secular) Hamiltonian plus isotropic one-way M1 conversion.
// Explicitly requires an F-dephased initial state; no Cu parameters or RNG.
class MuMinusSecularTrajectory {
public:
    MuMinusSecularTrajectory(const MuNuclearDensity& initial,
                            const MuMinusNuclearParameters& parameters,
                            double conversion_rate_ns)
        : initial_(initial), rate_(conversion_rate_ns) {
        mu_minus_detail::validate_density(initial_);
        (void)mu_minus_detail::hamiltonian(parameters);
        require_nonnegative(rate_, "M1 conversion rate must be finite and >= 0");
        if (parameters.hyperfine_rad_ns == 0)
            throw std::invalid_argument("secular trajectory requires nonzero hyperfine splitting");
        upper_F_ = parameters.hyperfine_rad_ns > 0 ? 2 : 1;
        const double field = parameters.field_T.R();
        const Vec3 axis = field == 0 ? Vec3(0, 0, 1) : parameters.field_T / field;
        const double g1 = -parameters.gamma_mu_rad_ns_T / 4 +
                          5 * parameters.gamma_nucleus_rad_ns_T / 4;
        const double g2 = parameters.gamma_mu_rad_ns_T / 4 +
                          3 * parameters.gamma_nucleus_rad_ns_T / 4;
        omega_upper_ = (upper_F_ == 2 ? g2 : g1) * field;
        omega_lower_ = (upper_F_ == 2 ? g1 : g2) * field;
        weak_field_ratio_ = (std::abs(parameters.gamma_mu_rad_ns_T) / 2 +
                            1.5 * std::abs(parameters.gamma_nucleus_rad_ns_T)) *
                           field / (2 * std::abs(parameters.hyperfine_rad_ns));
        if (!std::isfinite(omega_upper_) || !std::isfinite(omega_lower_) ||
            !std::isfinite(weak_field_ratio_))
            throw std::overflow_error("secular trajectory frequency overflow");
        using mu_minus_m1_detail::multiply;
        const auto coupling = mu_minus_detail::hamiltonian({{0, 0, 0}, 1, -1, 0});
        const auto along = mu_minus_detail::hamiltonian({axis, 0, -1, -1});
        MuNuclearDensity identity{}, p1{}, pu{}, pl{};
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
            identity[i][j] = double(i == j);
            p1[i][j] = (.75 * identity[i][j] - coupling[i][j]) / 2.;
            pu[i][j] = upper_F_ == 1 ? p1[i][j] : identity[i][j] - p1[i][j];
            pl[i][j] = identity[i][j] - pu[i][j];
        }
        const auto off_block = multiply(multiply(pu, initial_), pl);
        for (const auto& row : off_block) for (const auto& value : row)
            if (std::abs(value) > 1e-12)
                throw std::invalid_argument("secular trajectory requires F-dephased density");
        // Spectral projectors of F.axis, separately in each F manifold.
        std::array<MuNuclearDensity, 8> projectors{};
        std::array<int, 8> f{}, m{};
        int count = 0;
        for (int ff : {1, 2}) for (int mm = -ff; mm <= ff; ++mm) {
            f[count] = ff;
            m[count] = mm;
            auto projector = ff == 1 ? p1 : (upper_F_ == 1 ? pl : pu);
            for (int other = -ff; other <= ff; ++other) if (other != mm) {
                MuNuclearDensity factor{};
                for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                    factor[i][j] = (along[i][j] - double(other) * identity[i][j]) /
                                   double(mm - other);
                projector = multiply(projector, factor);
            }
            projectors[count++] = projector;
        }
        for (int a = 0; a < 8; ++a) for (int b = 0; b < 8; ++b) if (f[a] == f[b]) {
            const auto part = multiply(multiply(projectors[a], initial_), projectors[b]);
            auto& target = f[a] == upper_F_ ? upper_[m[a] - m[b] + 4] : lower_[m[a] - m[b] + 4];
            for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) target[i][j] += part[i][j];
        }
        const double strength = upper_F_ == 2 ? 3. / 8 : 5. / 8;
        for (const Vec3& direction : {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)}) {
            const auto spin = mu_minus_detail::hamiltonian({direction, 0, -1, 0});
            const auto jump = multiply(multiply(pl, spin), pu);
            for (int q = 0; q < 9; ++q) {
                const auto gain = multiply(multiply(jump, upper_[q]), mu_minus_m1_detail::adjoint(jump));
                for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                    transferred_[q][i][j] += gain[i][j] / strength;
            }
        }
    }

    MuNuclearDensity density(double elapsed_ns) const {
        require_nonnegative(elapsed_ns, "trajectory elapsed time must be finite and >= 0");
        if (elapsed_ns == 0) return initial_;
        MuNuclearDensity out{};
        const double decay = std::exp(-rate_ * elapsed_ns);
        for (int q = -4; q <= 4; ++q) {
            const double upper_phase = (q * omega_upper_) * elapsed_ns;
            const double lower_phase = (q * omega_lower_) * elapsed_ns;
            if (!std::isfinite(upper_phase) || !std::isfinite(lower_phase))
                throw std::overflow_error("secular trajectory phase overflow");
            const auto upper_phase_factor = std::polar(decay, upper_phase);
            const auto lower_phase_factor = std::polar(1., lower_phase);
            const auto weight = conversion_weight(q * (omega_upper_ - omega_lower_), elapsed_ns);
            for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                out[i][j] += upper_phase_factor * upper_[q + 4][i][j] +
                            lower_phase_factor * (lower_[q + 4][i][j] +
                                                  weight * transferred_[q + 4][i][j]);
        }
        (void)mu_minus_detail::real_embedding(out);
        return out;
    }

    Vec3 polarization(double elapsed_ns) const { return mu_minus_polarization(density(elapsed_ns)); }
    const MuNuclearDensity& initial_density() const { return initial_; }
    int upper_F() const { return upper_F_; }
    double zeeman_to_hyperfine_ratio() const { return weak_field_ratio_; }

private:
    std::complex<double> conversion_weight(double frequency, double time) const {
        if (rate_ == 0) return 0.;
        if (frequency == 0) return -std::expm1(-rate_ * time);
        const double a = rate_ * time, phase = frequency * time;
        if (!std::isfinite(phase)) throw std::overflow_error("conversion phase overflow");
        const std::complex<double> z(a, -phase);
        if (std::abs(z) < 1e-3) {
            std::complex<double> term(1), sum(1);
            for (int k = 1; k <= 6; ++k) { term *= -z / double(k + 1); sum += term; }
            return a * sum;
        }
        return (rate_ / std::complex<double>(rate_, -frequency)) *
               (1. - std::polar(std::exp(-a), phase));
    }

    MuNuclearDensity initial_;
    std::array<MuNuclearDensity, 9> upper_{}, lower_{}, transferred_{};
    double rate_, omega_upper_{}, omega_lower_{}, weak_field_ratio_{};
    int upper_F_{};
};

} // namespace stopped_muon
