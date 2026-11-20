#pragma once

#include "stopped_muon.hxx"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include <array>
#include <complex>

namespace stopped_muon {

using MuNuclearDensity = std::array<std::array<std::complex<double>, 8>, 8>;

// Explicit coherent model for mu-(1/2) x nucleus(3/2), not a Cu parameter set.
struct MuMinusNuclearParameters {
    Vec3 field_T;
    double hyperfine_rad_ns;       // Signed A in H/hbar = A S_mu.I.
    double gamma_mu_rad_ns_T;      // Physical signed gyromagnetic ratio, negative.
    double gamma_nucleus_rad_ns_T; // Physical signed nuclear gyromagnetic ratio.
};

// Only an explicitly chosen uncorrelated POST-cascade test state, not a cascade model.
inline MuNuclearDensity mu_minus_product_density(const Vec3& p_post) {
    if (!finite(p_post) || !std::isfinite(p_post.R()) || p_post.R() > 1) {
        throw std::invalid_argument("post-cascade polarization must be finite with magnitude <= 1");
    }
    MuNuclearDensity rho{};
    for (int n = 0; n < 4; ++n) {
        rho[n][n] = (1 + p_post.Z()) / 8;
        rho[n + 4][n + 4] = (1 - p_post.Z()) / 8;
        rho[n][n + 4] = std::complex<double>(p_post.X(), -p_post.Y()) / 8.;
        rho[n + 4][n] = std::conj(rho[n][n + 4]);
    }
    return rho;
}

inline Vec3 mu_minus_polarization(const MuNuclearDensity& rho) {
    std::complex<double> coherence{}, longitudinal{};
    for (int n = 0; n < 4; ++n) {
        coherence += rho[n + 4][n];
        longitudinal += rho[n][n] - rho[n + 4][n + 4];
    }
    return {2 * coherence.real(), 2 * coherence.imag(), longitudinal.real()};
}

namespace mu_minus_detail {
// Same ROOT real-symmetric embedding technique as the existing Mu core.
inline TMatrixDSym real_embedding(const MuNuclearDensity& h) {
    TMatrixDSym out(16);
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
        if (!std::isfinite(h[i][j].real()) || !std::isfinite(h[i][j].imag())) {
            throw std::overflow_error("non-finite mu-nuclear matrix");
        }
        out(i, j) = out(i + 8, j + 8) = h[i][j].real();
        out(i + 8, j) = h[i][j].imag();
        out(i, j + 8) = -h[i][j].imag();
    }
    return out;
}

inline void validate_density(const MuNuclearDensity& rho) {
    std::complex<double> trace{};
    for (int i = 0; i < 8; ++i) {
        trace += rho[i][i];
        for (int j = 0; j < 8; ++j) {
            if (!std::isfinite(rho[i][j].real()) || !std::isfinite(rho[i][j].imag()) ||
                std::abs(rho[i][j] - std::conj(rho[j][i])) > 1e-12) {
                throw std::invalid_argument("post-cascade density must be finite and Hermitian");
            }
        }
    }
    if (std::abs(trace - std::complex<double>(1, 0)) > 1e-12) {
        throw std::invalid_argument("post-cascade density must have trace one");
    }
    const TMatrixDSymEigen spectrum(real_embedding(rho));
    for (int k = 0; k < 16; ++k) {
        const double value = spectrum.GetEigenValues()[k];
        if (!std::isfinite(value) || value < -1e-12) {
            throw std::invalid_argument("post-cascade density must be positive semidefinite");
        }
    }
    // Tolerance permits roundoff only; never normalize, symmetrize or clip the state.
}

inline MuNuclearDensity hamiltonian(const MuMinusNuclearParameters& p) {
    if (!finite(p.field_T) || !std::isfinite(p.field_T.R()) ||
        !std::isfinite(p.hyperfine_rad_ns) || !std::isfinite(p.gamma_mu_rad_ns_T) ||
        p.gamma_mu_rad_ns_T >= 0 || !std::isfinite(p.gamma_nucleus_rad_ns_T)) {
        throw std::invalid_argument("mu-nuclear parameters must be finite with negative mu- gamma");
    }
    MuNuclearDensity h{};
    for (int i = 0; i < 8; ++i) {
        const double s = i < 4 ? .5 : -.5;
        const int n = i % 4;
        const double m = 1.5 - n;
        h[i][i] = p.hyperfine_rad_ns * s * m -
            p.field_T.Z() * (p.gamma_mu_rad_ns_T * s + p.gamma_nucleus_rad_ns_T * m);
        h[i][i ^ 4] = -p.gamma_mu_rad_ns_T / 2 *
            std::complex<double>(p.field_T.X(), -2 * s * p.field_T.Y());
        if (n > 0) {
            const double ladder = std::sqrt(double(n * (4 - n))) / 2;
            h[i - 1][i] = -p.gamma_nucleus_rad_ns_T * ladder *
                std::complex<double>(p.field_T.X(), -p.field_T.Y());
            h[i][i - 1] = std::conj(h[i - 1][i]);
            if (i < 4) h[i][i + 3] = h[i + 3][i] = p.hyperfine_rad_ns * ladder;
        }
    }
    return h;
}
} // namespace mu_minus_detail

// t=0 is the supplied POST-cascade joint state, not a redefinition of P_stop.
// Constant field, coherent evolution only; no RNG, cascade, relaxation or termination.
class MuMinusNuclearEvolution {
public:
    MuMinusNuclearEvolution(const MuNuclearDensity& rho_post,
                           const MuMinusNuclearParameters& parameters)
        : initial_(rho_post) {
        mu_minus_detail::validate_density(initial_);
        const TMatrixDSymEigen spectrum(mu_minus_detail::real_embedding(
            mu_minus_detail::hamiltonian(parameters)));
        const auto& values = spectrum.GetEigenValues();
        const auto& vectors = spectrum.GetEigenVectors();
        for (int k = 0; k < 16; ++k) {
            energies_[k] = values[k];
            if (!std::isfinite(energies_[k])) throw std::overflow_error("mu-nuclear spectrum overflow");
            for (int i = 0; i < 16; ++i) {
                eigenvectors_[i][k] = vectors(i, k);
                if (!std::isfinite(eigenvectors_[i][k])) {
                    throw std::overflow_error("mu-nuclear eigenvector overflow");
                }
            }
        }
    }

    MuNuclearDensity density(double elapsed_ns) const {
        require_nonnegative(elapsed_ns, "post-cascade elapsed time must be finite and >= 0");
        if (elapsed_ns == 0) return initial_;
        MuNuclearDensity u{}, temporary{}, result{};
        for (int k = 0; k < 16; ++k) {
            const double phase = energies_[k] * elapsed_ns;
            if (!std::isfinite(phase)) throw std::overflow_error("mu-nuclear phase overflow");
            const double c = std::cos(phase), s = std::sin(phase);
            for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
                const double a = eigenvectors_[i][k] * eigenvectors_[j][k];
                const double b = eigenvectors_[i + 8][k] * eigenvectors_[j][k];
                u[i][j] += std::complex<double>(a * c + b * s, b * c - a * s);
            }
        }
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 8; ++k) temporary[i][j] += u[i][k] * initial_[k][j];
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 8; ++k) result[i][j] += temporary[i][k] * std::conj(u[j][k]);
        (void)mu_minus_detail::real_embedding(result);
        return result;
    }

    Vec3 polarization(double elapsed_ns) const { return mu_minus_polarization(density(elapsed_ns)); }
    const MuNuclearDensity& initial_density() const { return initial_; }

private:
    MuNuclearDensity initial_;
    std::array<double, 16> energies_{};
    std::array<std::array<double, 16>, 16> eigenvectors_{};
};

} // namespace stopped_muon
