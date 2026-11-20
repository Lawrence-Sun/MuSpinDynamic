#pragma once

#include "mu_minus_nuclear.hxx"

namespace stopped_muon {

// Isotropic, F-dephased spin-1/2 -> (spin-1/2 x spin-3/2) channel.
// Coefficients multiply the laboratory P_stop, not its unit vector.
// No Cu defaults: atomic cascade coefficients must be supplied explicitly.
struct MuMinusIsotropicInitialMap {
    double population_F1;
    double orientation_F1;
    double orientation_F2;
};

inline MuNuclearDensity mu_minus_isotropic_post_density(
        const Vec3& p_stop, const MuMinusIsotropicInitialMap& map) {
    const double p1 = map.population_F1, p2 = 1 - p1;
    const double c1 = map.orientation_F1, c2 = map.orientation_F2;
    if (!finite(p_stop) || !std::isfinite(p_stop.R()) || p_stop.R() > 1 ||
        !std::isfinite(p1) || p1 < 0 || p1 > 1 ||
        !std::isfinite(c1) || !std::isfinite(c2) ||
        c1 < -p1 / 6 || c1 > p1 / 3 || c2 < -p2 / 15 || c2 > p2 / 10) {
        throw std::invalid_argument("invalid P_stop or non-completely-positive isotropic initial map");
    }
    const auto coupling = mu_minus_detail::hamiltonian({{0, 0, 0}, 1, -1, 0});
    MuNuclearDensity f1{}, f2{}, fp{}, rho{};
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
        const double identity = i == j ? 1 : 0;
        f1[i][j] = (.75 * identity - coupling[i][j]) / 2.;
        f2[i][j] = identity - f1[i][j];
        rho[i][j] = p1 / 3 * f1[i][j] + p2 / 5 * f2[i][j];
    }
    for (int i = 0; i < 8; ++i) {
        const double s = i < 4 ? .5 : -.5;
        const int n = i % 4;
        fp[i][i] = (s + 1.5 - n) * p_stop.Z();
        fp[i][i ^ 4] = .5 * std::complex<double>(p_stop.X(), -2 * s * p_stop.Y());
        if (n > 0) {
            fp[i - 1][i] = std::sqrt(double(n * (4 - n))) / 2 *
                std::complex<double>(p_stop.X(), -p_stop.Y());
            fp[i][i - 1] = std::conj(fp[i - 1][i]);
        }
    }
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j)
        for (int k = 0; k < 8; ++k)
            rho[i][j] += (c1 * f1[i][k] + c2 * f2[i][k]) * fp[k][j];
    return rho;
}

} // namespace stopped_muon
