#pragma once

#include "stopped_muon.hxx"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace stopped_muon {

using DensityMatrix = std::array<std::array<std::complex<double>, 4>, 4>;

namespace muonium_detail {
inline DensityMatrix initial_density(const Vec3& p) {
    (void)SpinEvolution(SpinState::Muonium, p, {});
    DensityMatrix rho{};
    for (int e = 0; e < 2; ++e) {
        rho[e][e] = (1 + p.Z()) / 4;
        rho[e + 2][e + 2] = (1 - p.Z()) / 4;
        rho[e][e + 2] = std::complex<double>(p.X(), -p.Y()) / 4.;
        rho[e + 2][e] = std::conj(rho[e][e + 2]);
    }
    return rho;
}

inline DensityMatrix hamiltonian(const Vec3& field, double hyperfine) {
    DensityMatrix h{};
    for (int i = 0; i < 4; ++i) {
        const double m = i < 2 ? 0.5 : -0.5;
        const double e = i % 2 == 0 ? 0.5 : -0.5;
        h[i][i] = hyperfine * m * e + field.Z() * (-gamma_muon * m + gamma_electron * e);
        h[i][i ^ 2] = -gamma_muon / 2 * std::complex<double>(field.X(), -2 * m * field.Y());
        h[i][i ^ 1] = gamma_electron / 2 * std::complex<double>(field.X(), -2 * e * field.Y());
    }
    h[1][2] = h[2][1] = hyperfine / 2;
    return h;
}

// A Hermitian 4x4 matrix becomes the real symmetric block matrix [Re H,-Im H; Im H,Re H].
inline TMatrixDSym real_embedding(const DensityMatrix& h) {
    TMatrixDSym real(8);
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
        if (!std::isfinite(h[i][j].real()) || !std::isfinite(h[i][j].imag())) {
            throw std::overflow_error("non-finite muonium matrix");
        }
        real(i, j) = real(i + 4, j + 4) = h[i][j].real();
        real(i + 4, j) = h[i][j].imag();
        real(i, j + 4) = -h[i][j].imag();
    }
    return real;
}
} // namespace muonium_detail

inline Vec3 muon_polarization(const DensityMatrix& rho) {
    const auto coherence = rho[2][0] + rho[3][1];
    return {2 * coherence.real(), 2 * coherence.imag(),
            (rho[0][0] + rho[1][1] - rho[2][2] - rho[3][3]).real()};
}

struct MuoniumParameters {
    Vec3 external_field{0, 0, 0}; // T; the same magnetic field acts on both spins.
    double local_field_sigma_T = 0; // Standard deviation of each Cartesian component, T.
    double hop_rate = 0; // ns^-1.
    double hyperfine = omega_muonium_hfs; // Isotropic A in rad/ns; zero allows a decoupled test.
};

class MuoniumHoppingTrajectory;

// Coherent evolution only: initially unpolarized electron, no T1/T2 multiplier.
class MuoniumEvolution {
public:
    MuoniumEvolution(const Vec3& p0, const Vec3& field, double hyperfine = omega_muonium_hfs)
        : MuoniumEvolution(muonium_detail::initial_density(p0), field, hyperfine) {}

    DensityMatrix density(double time) const {
        require_nonnegative(time, "muonium time must be finite and >= 0");
        if (time == 0) return initial_;
        DensityMatrix u{};
        for (int k = 0; k < 8; ++k) {
            const double phase = energies_[k] * time;
            if (!std::isfinite(phase)) throw std::overflow_error("muonium phase overflow");
            const double c = std::cos(phase), s = std::sin(phase);
            for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) {
                const double a = eigenvectors_[i][k] * eigenvectors_[j][k];
                const double b = eigenvectors_[i + 4][k] * eigenvectors_[j][k];
                // Recover cos(Ht)-i sin(Ht) from the real block matrix functions.
                u[i][j] += std::complex<double>(a * c + b * s, b * c - a * s);
            }
        }
        DensityMatrix temporary{}, result{};
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) temporary[i][j] += u[i][k] * initial_[k][j];
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) result[i][j] += temporary[i][k] * std::conj(u[j][k]);
        (void)muonium_detail::real_embedding(result);
        return result;
    }

    Vec3 polarization(double time) const { return muon_polarization(density(time)); }
    const DensityMatrix& initial_density() const { return initial_; }

private:
    friend class MuoniumHoppingTrajectory;
    MuoniumEvolution(const DensityMatrix& initial, const Vec3& field, double hyperfine)
        : initial_(initial) {
        require_nonnegative(hyperfine, "hyperfine frequency must be finite and >= 0");
        if (!finite(field) || !std::isfinite(field.R())) {
            throw std::invalid_argument("muonium magnetic field must be finite");
        }
        const TMatrixDSymEigen spectrum(muonium_detail::real_embedding(
            muonium_detail::hamiltonian(field, hyperfine)));
        const auto& values = spectrum.GetEigenValues();
        const auto& vectors = spectrum.GetEigenVectors();
        for (int k = 0; k < 8; ++k) {
            energies_[k] = values[k];
            if (!std::isfinite(energies_[k])) throw std::overflow_error("muonium spectrum overflow");
            for (int i = 0; i < 8; ++i) {
                eigenvectors_[i][k] = vectors(i, k);
                if (!std::isfinite(eigenvectors_[i][k])) {
                    throw std::overflow_error("muonium eigenvector overflow");
                }
            }
        }
    }

    DensityMatrix initial_;
    std::array<double, 8> energies_{};
    std::array<std::array<double, 8>, 8> eigenvectors_{};
};

class MuoniumHoppingTrajectory {
public:
    struct Segment {
        double start;
        Vec3 local_field;
        MuoniumEvolution evolution;
    };

    template <typename Rng>
    MuoniumHoppingTrajectory(const Vec3& p0, const MuoniumParameters& parameters,
                            double duration, Rng& rng, std::size_t max_hops = 100000)
        : duration_(duration) {
        // Validate deterministic inputs before drawing a field history.
        const MuoniumEvolution validated(p0, parameters.external_field, parameters.hyperfine);
        require_nonnegative(parameters.local_field_sigma_T, "local field sigma must be finite and >= 0");
        const double delta = parameters.local_field_sigma_T * gamma_muon;
        if (!std::isfinite(delta)) throw std::invalid_argument("local field sigma conversion overflow");
        // Reuse the tested Poisson/Gaussian scheduler; its dummy spin is identically zero.
        const HoppingTrajectory history(SpinState::MuPlus, Vec3{},
            {parameters.external_field, delta, parameters.hop_rate}, duration, rng, max_hops);
        DensityMatrix rho = validated.initial_density();
        for (const auto& field : history.segments()) {
            if (!segments_.empty()) {
                const auto& previous = segments_.back();
                rho = previous.evolution.density(field.start - previous.start);
            }
            segments_.push_back({field.start, field.local_field,
                MuoniumEvolution(rho, parameters.external_field + field.local_field, parameters.hyperfine)});
        }
    }

    DensityMatrix density(double time) const {
        require_nonnegative(time, "muonium query time must be finite and >= 0");
        if (time > duration_) throw std::invalid_argument("query exceeds prepared muonium duration");
        const auto after = std::upper_bound(segments_.begin(), segments_.end(), time,
            [](double t, const Segment& s) { return t < s.start; });
        const auto& segment = *(after - 1);
        return segment.evolution.density(time - segment.start);
    }
    Vec3 polarization(double time) const { return muon_polarization(density(time)); }
    const std::vector<Segment>& segments() const { return segments_; }
    std::size_t hop_count() const { return segments_.size() - 1; }
    double duration() const { return duration_; }

private:
    double duration_;
    std::vector<Segment> segments_;
};

} // namespace stopped_muon
