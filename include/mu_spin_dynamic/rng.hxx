#pragma once

#include "types.hxx"

#include <array>
#include <cmath>
#include <cstdint>
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

// Deterministic 64-bit SplitMix64 mixer for deriving independent per-(seed, event_id, replica_id) streams.
inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline std::uint64_t derive_event_replica_seed(
        std::uint64_t base_seed,
        std::int64_t event_id,
        std::int32_t replica_id = 0) noexcept {
    std::uint64_t s = splitmix64(base_seed ^ 0x243f6a8885a308d3ULL);
    s = splitmix64(s ^ static_cast<std::uint64_t>(event_id));
    s = splitmix64(s ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(replica_id)) +
                        0x9e3779b97f4a7c15ULL));
    return s;
}

// Fast, compact Xoshiro256** uniform [0, 1) generator (32-byte state, cheap per-event seeding).
class Xoshiro256StarStarUniform {
public:
    explicit Xoshiro256StarStarUniform(std::uint64_t seed) noexcept {
        std::uint64_t state = seed;
        for (std::uint64_t& word : s_) {
            state = splitmix64(state);
            word = state;
        }
    }

    double operator()() noexcept {
        // Exact 53-bit IEEE-754 uniform double in [0, 1).
        return static_cast<double>(next_u64() >> 11) * (1.0 / 9007199254740992.0);
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) noexcept {
        return (x << k) | (x >> (64 - k));
    }

    std::uint64_t next_u64() noexcept {
        const std::uint64_t result = rotl(s_[1] * 5ULL, 7) * 9ULL;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    std::array<std::uint64_t, 4> s_{};
};

using EventReplicaRng = UniformRngAdapter<Xoshiro256StarStarUniform>;

inline EventReplicaRng make_event_replica_rng(
        std::uint64_t base_seed,
        std::int64_t event_id,
        std::int32_t replica_id = 0) noexcept {
    return EventReplicaRng(Xoshiro256StarStarUniform(
        derive_event_replica_seed(base_seed, event_id, replica_id)));
}

} // namespace stopped_muon
