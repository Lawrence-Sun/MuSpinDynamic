#pragma once

#include "bound_dio.hxx"
#include "mu_minus_termination.hxx"

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace stopped_muon {

// Bound 1s negative muon on a zero-nuclear-spin (I = 0) target nuclide (e.g., 56Fe, 32S).
// No nuclear hyperfine splitting or M1 conversion; explicit post-cascade residual polarization
// factor f_res in [0, 1] is applied once to p_stop (theoretical I=0 cascade limit ~1/6).
struct MuMinusSpin0Material {
    DioNuclide nuclide = DioNuclide::Fe56;
    double residual_polarization = 1.0 / 6.0; // Post-cascade factor f_res in [0, 1].
    Vec3 external_field{0, 0, 0};             // Laboratory external field [T].
    DomainFieldModel domain_field{};          // Static ferromagnetic internal field (None for non-magnetic S).
    double local_sigma_T = 0;                 // Per-component Gaussian local-field width [T].
    double hop_rate = 0;                      // Local-field fluctuation rate [ns^-1].
    double gyromagnetic_ratio = gamma_muon;   // Positive magnitude |gamma_mu-| [rad/(ns*T)].
    MuMinusTerminationRates rates{0, 0};      // Constant competing DIO and nuclear-capture rates [ns^-1].
    std::size_t max_hops = 100000;
    std::size_t max_dio_trials = 100000;
};

struct MuMinusSpin0Event {
    DioNuclide nuclide;
    int Z;
    int isotope;
    Vec3 p_stop;
    Vec3 p_post_cascade;
    Vec3 domain_field_T;
    MuMinusTerminationChannel channel;
    double time_ns;
    Vec3 decay_polarization;
    std::optional<DioElectron> electron; // Empty on nuclear capture; never a zero-momentum placeholder.
    HoppingTrajectory trajectory;
    int charge = -1;

    Vec3 polarization(double t_ns) const { return trajectory.polarization(t_ns); }
    std::size_t hop_count() const { return trajectory.hop_count(); }
};

inline MuMinusSpin0Event simulate_stopped_mu_minus_spin0(
        const MuMinusSpin0Material& material, const Vec3& p_stop, TRandom3& rng) {
    // Validate all deterministic inputs before consuming any RNG draws.
    if (!finite(p_stop) || !std::isfinite(p_stop.R()) || p_stop.R() > 1) {
        throw std::invalid_argument("mu- spin-0 p_stop must be finite with |p_stop| <= 1");
    }
    if (!std::isfinite(material.residual_polarization) ||
        material.residual_polarization < 0 || material.residual_polarization > 1) {
        throw std::invalid_argument("mu- spin-0 residual_polarization must be in [0, 1]");
    }
    if (!finite(material.external_field) || !std::isfinite(material.external_field.R())) {
        throw std::invalid_argument("mu- spin-0 external_field must be finite");
    }
    validate_domain_field(material.domain_field);
    require_nonnegative(material.local_sigma_T, "mu- spin-0 local_sigma_T must be finite and >= 0");
    require_nonnegative(material.hop_rate, "mu- spin-0 hop_rate must be finite and >= 0");
    if (!std::isfinite(material.gyromagnetic_ratio) || material.gyromagnetic_ratio <= 0) {
        throw std::invalid_argument("mu- spin-0 gyromagnetic_ratio magnitude must be finite and > 0");
    }
    require_nonnegative(material.rates.dio_ns, "DIO rate must be finite and >= 0");
    require_nonnegative(material.rates.capture_ns, "capture rate must be finite and >= 0");
    const double total_rate = material.rates.dio_ns + material.rates.capture_ns;
    if (!std::isfinite(total_rate) || total_rate <= 0 || !std::isfinite(1.0 / total_rate)) {
        throw std::invalid_argument("mu- spin-0 total termination rate must be positive with finite mean time");
    }
    if (material.max_dio_trials == 0) {
        throw std::invalid_argument("mu- spin-0 DIO trial budget must be positive");
    }
    const BoundDioSpectrum spectrum(material.nuclide);
    const double gamma_scale = material.gyromagnetic_ratio / gamma_muon;
    const double effective_delta = material.local_sigma_T * material.gyromagnetic_ratio;
    if (!std::isfinite(gamma_scale) || !std::isfinite(effective_delta)) {
        throw std::invalid_argument("mu- spin-0 gyromagnetic scaling overflow");
    }

    // Stochastic draws: single termination time, single channel, optional domain field, trajectory, optional DIO.
    const double time_ns = rng.Exp(1.0 / total_rate);
    if (!std::isfinite(time_ns)) {
        throw std::overflow_error("nonfinite mu- spin-0 termination time");
    }
    const auto channel = (rng.Rndm() < material.rates.dio_ns / total_rate)
        ? MuMinusTerminationChannel::DIO
        : MuMinusTerminationChannel::NuclearCapture;

    const Vec3 p_post = material.residual_polarization * p_stop;
    const Vec3 dom_field = sample_domain_field(material.domain_field, rng);
    const Vec3 static_field = material.external_field + dom_field;
    if (!finite(static_field) || !std::isfinite(static_field.R())) {
        throw std::overflow_error("mu- spin-0 static field overflow");
    }
    const Vec3 scaled_static_field = gamma_scale * static_field;
    if (!finite(scaled_static_field) || !std::isfinite(scaled_static_field.R())) {
        throw std::overflow_error("mu- spin-0 scaled static field overflow");
    }

    HoppingTrajectory trajectory(
        SpinState::MuMinus, p_post,
        {scaled_static_field, effective_delta, material.hop_rate},
        time_ns, rng, material.max_hops);
    const Vec3 decay_p = trajectory.polarization(time_ns);

    std::optional<DioElectron> electron;
    if (channel == MuMinusTerminationChannel::DIO) {
        electron = spectrum.sample(decay_p, rng, material.max_dio_trials);
    }

    return {
        material.nuclide,
        spectrum.atomic_number(),
        spectrum.isotope(),
        p_stop,
        p_post,
        dom_field,
        channel,
        time_ns,
        decay_p,
        std::move(electron),
        std::move(trajectory),
    };
}

} // namespace stopped_muon
