#pragma once

#include "copper_dio.hxx"
#include "mu_minus_initial_state.hxx"
#include "mu_minus_termination.hxx"
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

namespace stopped_muon {

// One explicitly selected isotope; no automatic abundance/stopping/DIO mixture.
struct MuMinusMaterial {
    int isotope;
    MuMinusIsotropicInitialMap preparation;
    MuMinusNuclearParameters spin;
    double conversion_rate_ns;
    MuMinusTerminationRates rates;
    std::size_t max_dio_trials = 100000;
};

struct MuMinusEvent {
    int isotope;
    Vec3 p_stop;
    MuMinusTermination termination;
    std::optional<DioElectron> electron; // Empty for nuclear capture, never a zero-momentum placeholder.
    MuMinusSecularTrajectory trajectory;
    int charge = -1;
};

template <typename Rng>
inline MuMinusEvent simulate_stopped_mu_minus(
        const MuMinusMaterial& material, const Vec3& p_stop, Rng& rng) {
    const CopperDioSpectrum spectrum(material.isotope);
    if (material.max_dio_trials == 0)
        throw std::invalid_argument("Cu DIO trial budget must be positive");
    const auto initial = mu_minus_isotropic_post_density(p_stop, material.preparation);
    MuMinusSecularTrajectory trajectory(initial, material.spin, material.conversion_rate_ns);
    auto termination = sample_mu_minus_termination(trajectory, material.rates, rng);
    std::optional<DioElectron> electron;
    if (termination.channel == MuMinusTerminationChannel::DIO)
        electron = spectrum.sample(termination.polarization(), rng, material.max_dio_trials);
    return {material.isotope, p_stop, std::move(termination), std::move(electron), std::move(trajectory)};
}

} // namespace stopped_muon
