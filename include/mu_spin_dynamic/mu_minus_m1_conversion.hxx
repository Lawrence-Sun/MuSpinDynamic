#pragma once

#include "mu_minus_nuclear.hxx"

namespace stopped_muon {
namespace mu_minus_m1_detail {
inline MuNuclearDensity multiply(const MuNuclearDensity& a, const MuNuclearDensity& b) {
    MuNuclearDensity out{};
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
        for (int k = 0; k < 8; ++k) out[i][j] += a[i][k] * b[k][j];
    return out;
}
inline MuNuclearDensity adjoint(const MuNuclearDensity& a) {
    MuNuclearDensity out{};
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
        out[i][j] = std::conj(a[j][i]);
    return out;
}
} // namespace mu_minus_m1_detail

// Isotropic, one-way, rank-1 M1 conversion only; no Hamiltonian or Cu defaults.
// For F-dephased initial states this is exact zero-field spin evolution within
// the stated stationary-bath model. General F coherence is retained and damped,
// but its free hyperfine phase is NOT supplied by this dissipative channel.
class MuMinusM1Conversion {
public:
    MuMinusM1Conversion(const MuNuclearDensity& initial, int upper_F, double rate_ns)
        : initial_(initial), rate_(rate_ns) {
        mu_minus_detail::validate_density(initial_);
        require_nonnegative(rate_, "M1 conversion rate must be finite and >= 0");
        if (upper_F != 1 && upper_F != 2)
            throw std::invalid_argument("upper F must be explicitly 1 or 2");
        const auto coupling = mu_minus_detail::hamiltonian({{0, 0, 0}, 1, -1, 0});
        MuNuclearDensity pu{}, pl{};
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
            const double identity = i == j ? 1 : 0;
            const auto pi1 = (.75 * identity - coupling[i][j]) / 2.;
            pu[i][j] = upper_F == 1 ? pi1 : identity - pi1;
            pl[i][j] = identity - pu[i][j];
        }
        using mu_minus_m1_detail::multiply;
        upper_ = multiply(multiply(pu, initial_), pu);
        lower_ = multiply(multiply(pl, initial_), pl);
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
            cross_[i][j] = initial_[i][j] - upper_[i][j] - lower_[i][j];
        const double strength = upper_F == 2 ? 3. / 8 : 5. / 8;
        for (const Vec3& axis : {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)}) {
            const auto spin = mu_minus_detail::hamiltonian({axis, 0, -1, 0});
            const auto jump = multiply(multiply(pl, spin), pu);
            const auto gain = multiply(multiply(jump, initial_), mu_minus_m1_detail::adjoint(jump));
            for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
                transferred_[i][j] += gain[i][j] / strength;
        }
    }

    MuNuclearDensity density(double elapsed_ns) const {
        require_nonnegative(elapsed_ns, "M1 elapsed time must be finite and >= 0");
        if (elapsed_ns == 0 || rate_ == 0) return initial_;
        const double s = std::exp(-rate_ * elapsed_ns);
        const double amplitude = std::exp(-(.5 * rate_) * elapsed_ns);
        const double converted = -std::expm1(-rate_ * elapsed_ns);
        MuNuclearDensity out{};
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
            out[i][j] = lower_[i][j] + s * upper_[i][j] + amplitude * cross_[i][j] +
                        converted * transferred_[i][j];
        return out;
    }

    Vec3 polarization(double elapsed_ns) const { return mu_minus_polarization(density(elapsed_ns)); }

private:
    MuNuclearDensity initial_, upper_{}, lower_{}, cross_{}, transferred_{};
    double rate_;
};

} // namespace stopped_muon
