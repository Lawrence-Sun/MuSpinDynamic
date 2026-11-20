#pragma once

#include "stopped_muon.hxx"
#include "bound_dio_data.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace stopped_muon {

enum class DioNuclide { Cu63 = 0, Cu65 = 1, Fe56 = 2, S32 = 3 };

struct DioElectron {
    double energy_MeV;
    Vec3 momentum_MeV_c;
    int isotope;
    const char* table_sha256;
    std::size_t energy_trials, angular_trials;
    int pdg = 11;
    const char* model = "Cu_static_finite_charge_Watanabe_conditional";
};

inline std::size_t dio_nuclide_index(DioNuclide nuclide) {
    switch (nuclide) {
        case DioNuclide::Cu63: return 0;
        case DioNuclide::Cu65: return 1;
        case DioNuclide::Fe56: return 2;
        case DioNuclide::S32:  return 3;
    }
    throw std::invalid_argument("unknown DioNuclide");
}

inline const char* dio_nuclide_spin_category(DioNuclide nuclide) {
    switch (nuclide) {
        case DioNuclide::Cu63: return "cu63_mu_minus";
        case DioNuclide::Cu65: return "cu65_mu_minus";
        case DioNuclide::Fe56: return "fe56_mu_minus";
        case DioNuclide::S32:  return "s32_mu_minus";
    }
    throw std::invalid_argument("unknown DioNuclide");
}

// Hash-certified static finite-charge Watanabe DIO spectrum sampler (Cu63, Cu65, Fe56, S32).
// Certified shape only: Q_shape never alters caller-supplied competing rates.
class BoundDioSpectrum {
public:
    explicit BoundDioSpectrum(DioNuclide nuclide) : nuclide_(nuclide) {
        table_ = &bound_dio_detail::tables[dio_nuclide_index(nuclide)];
        for (std::size_t i = 0; i < 14; ++i) {
            const auto& r = table_->regions[i];
            envelope_integral_ += r.upper_N * (r.right - r.left);
            cumulative_[i] = envelope_integral_;
        }
    }

    DioNuclide nuclide() const { return nuclide_; }
    const char* name() const { return table_->name; }
    int atomic_number() const { return table_->Z; }
    int isotope() const { return table_->isotope; }
    const char* model() const { return table_->model; }
    double endpoint_MeV() const { return table_->endpoint; }
    double Q_shape() const { return table_->Q; }
    const char* table_sha256() const { return table_->sha256; }
    double envelope_integral() const { return envelope_integral_; }

    // Returns N(E), beta(E); total electron energy in MeV, N in MeV^-1.
    std::pair<double, double> density(double energy) const {
        if (!std::isfinite(energy) || energy < bound_dio_detail::electron_mass_MeV ||
            energy > table_->endpoint) {
            throw std::invalid_argument("DIO total energy outside [m_e, W]");
        }
        std::size_t i = 0;
        while (i < 13 && energy >= table_->regions[i].right) ++i;
        const auto& r = table_->regions[i];
        const double x = (2 * energy - r.left - r.right) / (r.right - r.left);
        const double gap = table_->endpoint - energy;
        const double N = std::exp(chebyshev(r.log_N, x)) * std::pow(gap, 5);
        const double beta = std::sqrt((energy - bound_dio_detail::electron_mass_MeV) *
                                      (energy + bound_dio_detail::electron_mass_MeV)) *
                            chebyshev(r.beta_over_p, x);
        if (!std::isfinite(N) || !std::isfinite(beta) || N < 0 || std::abs(beta) > 1) {
            throw std::runtime_error("invalid DIO interpolation; no clipping");
        }
        return {N, beta};
    }

    template <typename Rng>
    DioElectron sample(
            const Vec3& polarization, Rng& rng, std::size_t max_trials = 100000) const {
        const double magnitude = polarization.R();
        if (!finite(polarization) || !std::isfinite(magnitude) || magnitude > 1) {
            throw std::invalid_argument("DIO polarization must be finite with |P|<=1");
        }
        if (max_trials == 0) {
            throw std::invalid_argument("DIO trial budget must be positive");
        }
        double energy = 0, beta = 0;
        std::size_t energy_trials = 0, angular_trials = 0;
        for (;;) {
            if (energy_trials == max_trials) {
                throw std::runtime_error("DIO energy trial budget exhausted");
            }
            ++energy_trials;
            const double choice = rng.Rndm() * envelope_integral_;
            std::size_t i = 0;
            while (i < 13 && choice >= cumulative_[i]) ++i;
            const auto& r = table_->regions[i];
            energy = r.left + (r.right - r.left) * rng.Rndm();
            if (energy <= r.left || energy >= r.right) continue;
            const auto value = density(energy);
            // Bounds are generated from Bernstein enclosures, not fitted maxima.
            if (value.first > r.upper_N) {
                throw std::runtime_error("DIO energy envelope violated; no repair");
            }
            if (rng.Rndm() * r.upper_N < value.first) {
                beta = value.second;
                break;
            }
        }
        Vec3 direction;
        for (;;) {
            if (angular_trials == max_trials) {
                throw std::runtime_error("DIO angular trial budget exhausted");
            }
            ++angular_trials;
            const double z = 2 * rng.Rndm() - 1;
            const double phi = 2 * pi * rng.Rndm();
            const double transverse = std::sqrt(std::max(0., 1 - z * z));
            direction = {transverse * std::cos(phi), transverse * std::sin(phi), z};
            const double weight = 1 + beta * polarization.Dot(direction);
            const double upper = 1 + std::abs(beta) * magnitude + 1e-14;
            if (weight < 0 || weight > upper) {
                throw std::runtime_error("DIO angular envelope violated; no repair");
            }
            if (rng.Rndm() * upper < weight) break;
        }
        const double momentum = std::sqrt((energy - bound_dio_detail::electron_mass_MeV) *
                                          (energy + bound_dio_detail::electron_mass_MeV));
        return {energy, momentum * direction, isotope(), table_sha256(),
                energy_trials, angular_trials, 11, model()};
    }

private:
    static double chebyshev(const std::array<double, 8>& c, double x) {
        double b1 = 0, b2 = 0;
        for (int j = 7; j >= 1; --j) {
            const double b = c[j] + 2 * x * b1 - b2;
            b2 = b1;
            b1 = b;
        }
        return c[0] + x * b1 - b2;
    }

    DioNuclide nuclide_;
    const bound_dio_detail::Table* table_;
    std::array<double, 14> cumulative_{};
    double envelope_integral_ = 0;
};

} // namespace stopped_muon
