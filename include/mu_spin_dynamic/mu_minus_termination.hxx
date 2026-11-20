#pragma once

#include "mu_minus_secular_trajectory.hxx"

namespace stopped_muon {

// Spin-/F-independent constant competing rates; ns^-1, no material defaults.
struct MuMinusTerminationRates {
    double dio_ns;
    double capture_ns;
};

enum class MuMinusTerminationChannel { DIO, NuclearCapture };

// A termination record, NOT a complete DIO electron event.
// Capture has no decay electron; a DIO record still requires the Cu DIO sampler.
struct MuMinusTermination {
    MuMinusTerminationChannel channel;
    double time_ns;
    MuNuclearDensity density;

    Vec3 polarization() const { return mu_minus_polarization(density); }
};

template <typename Rng>
inline MuMinusTermination sample_mu_minus_termination(
        const MuMinusSecularTrajectory& trajectory,
        const MuMinusTerminationRates& rates, Rng& rng) {
    require_nonnegative(rates.dio_ns, "DIO rate must be finite and >= 0");
    require_nonnegative(rates.capture_ns, "capture rate must be finite and >= 0");
    const double total = rates.dio_ns + rates.capture_ns;
    if (!std::isfinite(total) || total <= 0 || !std::isfinite(1 / total))
        throw std::invalid_argument("total termination rate must be positive with finite mean time");
    // Exactly one time draw and one channel draw, including zero-channel limits.
    const double time = rng.Exp(1 / total);
    if (!std::isfinite(time)) throw std::overflow_error("nonfinite mu- termination time");
    const auto channel = rng.Rndm() < rates.dio_ns / total
        ? MuMinusTerminationChannel::DIO : MuMinusTerminationChannel::NuclearCapture;
    return {channel, time, trajectory.density(time)};
}

} // namespace stopped_muon
