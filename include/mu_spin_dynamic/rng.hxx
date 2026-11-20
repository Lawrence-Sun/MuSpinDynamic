#pragma once

#include "types.hxx"

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

namespace stopped_muon {

// Adapts any callable returning uniform doubles in [0, 1) (for example
// `[] { return G4UniformRand(); }` in Geant4 or a std::mt19937_64 wrapper)
// to the Rndm() / Exp(tau) / Gaus(mean, sigma) interface expected by the core
// physics routines. Existing callers passing ROOT's TRandom3& directly work
// without this adapter.
template <typename UniformFn>
class UniformRngAdapter {
public:
    explicit UniformRngAdapter(UniformFn uniform) : uniform_(std::move(uniform)) {}

    double Rndm() { return static_cast<double>(uniform_()); }

    // Exponential variate with mean `tau` (matches TRandom::Exp(tau)).
    double Exp(double tau) {
        if (tau <= 0.0) return 0.0;
        return -tau * std::log(next_positive_uniform());
    }

    // Normal variate with given mean and standard deviation (matches TRandom::Gaus).
    double Gaus(double mean = 0.0, double sigma = 1.0) {
        if (sigma == 0.0) return mean;
        if (has_spare_) {
            has_spare_ = false;
            return mean + sigma * spare_;
        }
        const double u1 = next_positive_uniform();
        const double u2 = Rndm();
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * pi * u2;
        spare_ = r * std::sin(theta);
        has_spare_ = true;
        return mean + sigma * r * std::cos(theta);
    }

private:
    double next_positive_uniform() {
        double u = 0.0;
        for (int attempt = 0; attempt < 64; ++attempt) {
            u = Rndm();
            if (u > 0.0) return u;
        }
        return std::numeric_limits<double>::min();
    }

    UniformFn uniform_;
    bool has_spare_ = false;
    double spare_ = 0.0;
};

template <typename UniformFn>
inline UniformRngAdapter<std::decay_t<UniformFn>> make_uniform_rng(UniformFn&& uniform) {
    return UniformRngAdapter<std::decay_t<UniformFn>>(std::forward<UniformFn>(uniform));
}

} // namespace stopped_muon
