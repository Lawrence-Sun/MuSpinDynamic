#pragma once

#include "anisotropic_muonium.hxx"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <variant>

namespace stopped_muon {

// Adopted muon rest mass in MeV; massless LO daughter endpoint is half this value.
inline constexpr double muon_mass_MeV = 105.6583755;
inline constexpr const char* free_michel_model = "free_michel_LO_massless";

struct MichelPositron {
    double energy_MeV;
    Vec3 momentum_MeV_c;
    double polarization_scale; // 1 unless roundoff made the supplied |P| slightly exceed 1.
    int pdg = -11;
    const char* model = free_michel_model;
};

template <typename Rng>
inline MichelPositron sample_muplus_michel(const Vec3& polarization, Rng& rng) {
    const double magnitude = polarization.R();
    if (!finite(polarization) || !std::isfinite(magnitude) || magnitude > 1 + 2e-9) {
        throw std::invalid_argument("Michel polarization must be finite with |P|<=1 (2e-9 roundoff tolerance)");
    }
    // Only the angular sampler bounds numerical norm drift; the stored spin state is never changed.
    const double scale = 1 / std::max(1., magnitude);
    const Vec3 p = scale * polarization;
    double x;
    do {
        x = rng.Rndm();
    } while (x == 0 || rng.Rndm() >= x * x * (3 - 2 * x)); // f(x)/2, envelope height 2.
    const double a = (2 * x - 1) / (3 - 2 * x);
    const double envelope = 1 + std::abs(a) * std::min(1., magnitude);
    Vec3 direction;
    double weight;
    do {
        const double z = 2 * rng.Rndm() - 1;
        const double phi = 2 * pi * rng.Rndm();
        const double transverse = std::sqrt(std::max(0., 1 - z * z));
        direction = {transverse * std::cos(phi), transverse * std::sin(phi), z};
        weight = 1 + a * p.Dot(direction); // mu+ sign: high-energy positrons favor +P.
    } while (rng.Rndm() * envelope >= weight);
    const double energy = .5 * muon_mass_MeV * x;
    return {energy, energy * direction, scale};
}

struct MuPlusMaterial {
    Vec3 external_field{0, 0, 0}; // Same external field for both spin categories, T.
    double muonium_fraction = 0;
    double ordinary_sigma_T = 0; // Per Cartesian component, T.
    double ordinary_hop_rate = 0; // ns^-1.
    double muonium_sigma_T = 0;
    double muonium_hop_rate = 0;
    double hyperfine = omega_muonium_hfs; // Signed isotropic A_iso in rad/ns.
    std::size_t max_hops = 100000;
    DomainFieldModel domain_field{}; // Static ferromagnetic domain internal field (default None).
    double muonium_dipolar_D = 0; // Signed axial dipolar D = A_parallel - A_perp, rad/ns.
    double muonium_spin_exchange_rate = 0; // Electron spin randomization rate, ns^-1.
    double muonium_reorientation_rate = 0; // Molecular axis reorientation rate, ns^-1.
    Vec3 muonium_fixed_axis{0, 0, 0}; // Zero => isotropic powder average.
    // Two-state chemical exchange (Cox et al. 2011 JPCM 23 315801):
    // State 1 uses (hyperfine, muonium_dipolar_D); State 2 uses (muonium_state2_hyperfine, muonium_state2_dipolar_D).
    // Enabled when muonium_exchange_w12 > 0 || muonium_exchange_w21 > 0 || non-zero state2/initial_p1.
    double muonium_state2_hyperfine = 0; // State 2 signed isotropic A2, rad/ns.
    double muonium_state2_dipolar_D = 0; // State 2 signed axial dipolar D2, rad/ns.
    double muonium_exchange_w12 = 0;     // State 1 -> State 2 exit rate W12, ns^-1.
    double muonium_exchange_w21 = 0;     // State 2 -> State 1 return rate W21, ns^-1.
    std::optional<double> muonium_initial_p1{}; // Initial State 1 fraction in [0, 1]; empty => equilibrium.
};

using MuPlusTrajectory = std::variant<
    HoppingTrajectory,
    MuoniumHoppingTrajectory,
    AnisotropicMuoniumTrajectory,
    TwoStateExchangeTrajectory>;

struct MuPlusEvent {
    int charge;
    SpinState spin_state;
    double decay_time_ns;
    Vec3 initial_polarization;
    Vec3 decay_polarization;
    MichelPositron positron;
    MuPlusTrajectory trajectory;
    const char* termination = "decay";
    bool has_positron = true;

    Vec3 polarization(double time_ns) const {
        return std::visit([&](const auto& track) { return track.polarization(time_ns); }, trajectory);
    }
    std::size_t hop_count() const {
        return std::visit([](const auto& track) { return track.hop_count(); }, trajectory);
    }
};

template <typename Rng>
inline MuPlusEvent simulate_stopped_muon(const MuPlusMaterial& material, int charge,
                                         const Vec3& p0, Rng& rng) {
    if (charge != 1) throw std::invalid_argument("stage B supports charge=+1 only; mu- capture/DIO is not implemented");
    // Deterministic input failures do not consume RNG; stochastic trajectory-budget failures can.
    (void)SpinEvolution(SpinState::MuPlus, p0, {material.external_field, 0, 0});
    if (!std::isfinite(material.muonium_fraction) || material.muonium_fraction < 0 || material.muonium_fraction > 1) {
        throw std::invalid_argument("muonium fraction must be in [0,1]");
    }
    require_nonnegative(material.ordinary_sigma_T, "ordinary field sigma must be finite and >=0");
    require_nonnegative(material.ordinary_hop_rate, "ordinary hop rate must be finite and >=0");
    require_nonnegative(material.muonium_sigma_T, "muonium field sigma must be finite and >=0");
    require_nonnegative(material.muonium_hop_rate, "muonium hop rate must be finite and >=0");
    if (!std::isfinite(material.hyperfine)) {
        throw std::invalid_argument("hyperfine must be finite");
    }
    if (!std::isfinite(material.muonium_dipolar_D)) {
        throw std::invalid_argument("muonium dipolar D must be finite");
    }
    if (!std::isfinite(material.muonium_state2_hyperfine)) {
        throw std::invalid_argument("muonium state2 hyperfine must be finite");
    }
    if (!std::isfinite(material.muonium_state2_dipolar_D)) {
        throw std::invalid_argument("muonium state2 dipolar D must be finite");
    }
    require_nonnegative(material.muonium_exchange_w12, "muonium exchange W12 must be finite and >=0");
    require_nonnegative(material.muonium_exchange_w21, "muonium exchange W21 must be finite and >=0");
    if (material.muonium_initial_p1.has_value()) {
        if (!std::isfinite(*material.muonium_initial_p1) ||
            *material.muonium_initial_p1 < 0.0 || *material.muonium_initial_p1 > 1.0) {
            throw std::invalid_argument("muonium initial_p1 must be finite and in [0, 1]");
        }
    }
    require_nonnegative(material.muonium_spin_exchange_rate, "muonium spin exchange rate must be finite and >=0");
    require_nonnegative(material.muonium_reorientation_rate, "muonium reorientation rate must be finite and >=0");
    if (!finite(material.muonium_fixed_axis) || !std::isfinite(material.muonium_fixed_axis.R())) {
        throw std::invalid_argument("muonium fixed axis must be finite");
    }
    validate_domain_field(material.domain_field);
    if (!std::isfinite(material.ordinary_sigma_T * gamma_muon) ||
        !std::isfinite(material.muonium_sigma_T * gamma_muon)) {
        throw std::invalid_argument("field width conversion overflow");
    }
    const bool use_two_state_exchange =
        material.muonium_exchange_w12 > 0 ||
        material.muonium_exchange_w21 > 0 ||
        material.muonium_state2_hyperfine != 0 ||
        material.muonium_state2_dipolar_D != 0 ||
        material.muonium_initial_p1.has_value();
    const bool use_anisotropic_radical =
        !use_two_state_exchange && (
            material.hyperfine < 0 ||
            material.muonium_dipolar_D != 0 ||
            material.muonium_spin_exchange_rate > 0 ||
            material.muonium_reorientation_rate > 0 ||
            material.muonium_fixed_axis.R() > 0);
    if (material.muonium_fraction > 0) {
        if (use_two_state_exchange) {
            const Vec3 probe_axis = material.muonium_fixed_axis.R() > 0
                ? material.muonium_fixed_axis : Vec3{0, 0, 1};
            (void)TwoStateExchangeEvolution(
                p0, material.external_field,
                material.hyperfine, material.muonium_dipolar_D,
                material.muonium_state2_hyperfine, material.muonium_state2_dipolar_D,
                material.muonium_exchange_w12, material.muonium_exchange_w21,
                probe_axis, material.muonium_spin_exchange_rate, material.muonium_initial_p1);
        } else if (use_anisotropic_radical) {
            const Vec3 probe_axis = material.muonium_fixed_axis.R() > 0
                ? material.muonium_fixed_axis : Vec3{0, 0, 1};
            (void)AnisotropicMuoniumEvolution(
                p0, material.external_field, material.hyperfine, material.muonium_dipolar_D, probe_axis);
        } else {
            (void)MuoniumEvolution(p0, material.external_field, material.hyperfine);
        }
    }
    const bool muonium = material.muonium_fraction == 1 ||
        (material.muonium_fraction > 0 && rng.Rndm() < material.muonium_fraction);
    const double time = rng.Exp(muon_lifetime);
    if (!std::isfinite(time)) throw std::overflow_error("nonfinite muon decay time");
    const Vec3 static_field = material.external_field + sample_domain_field(material.domain_field, rng);
    if (!finite(static_field) || !std::isfinite(static_field.R())) {
        throw std::overflow_error("static magnetic field overflow");
    }
    MuPlusTrajectory trajectory = [&]() -> MuPlusTrajectory {
        if (!muonium) {
            return HoppingTrajectory(SpinState::MuPlus, p0,
                {static_field, material.ordinary_sigma_T * gamma_muon, material.ordinary_hop_rate},
                time, rng, material.max_hops);
        }
        if (use_two_state_exchange) {
            TwoStateExchangeParameters tp{};
            tp.external_field = static_field;
            tp.local_field_sigma_T = material.muonium_sigma_T;
            tp.hop_rate = material.muonium_hop_rate;
            tp.state1_hyperfine_iso = material.hyperfine;
            tp.state1_dipolar_D = material.muonium_dipolar_D;
            tp.state2_hyperfine_iso = material.muonium_state2_hyperfine;
            tp.state2_dipolar_D = material.muonium_state2_dipolar_D;
            tp.exchange_w12 = material.muonium_exchange_w12;
            tp.exchange_w21 = material.muonium_exchange_w21;
            tp.spin_exchange_rate = material.muonium_spin_exchange_rate;
            tp.reorientation_rate = material.muonium_reorientation_rate;
            tp.fixed_axis = material.muonium_fixed_axis;
            tp.initial_p1 = material.muonium_initial_p1;
            return TwoStateExchangeTrajectory(p0, tp, time, rng, material.max_hops);
        }
        if (use_anisotropic_radical) {
            AnisotropicMuoniumParameters ap{};
            ap.external_field = static_field;
            ap.local_field_sigma_T = material.muonium_sigma_T;
            ap.hop_rate = material.muonium_hop_rate;
            ap.hyperfine_iso = material.hyperfine;
            ap.dipolar_D = material.muonium_dipolar_D;
            ap.spin_exchange_rate = material.muonium_spin_exchange_rate;
            ap.reorientation_rate = material.muonium_reorientation_rate;
            ap.fixed_axis = material.muonium_fixed_axis;
            return AnisotropicMuoniumTrajectory(p0, ap, time, rng, material.max_hops);
        }
        return MuoniumHoppingTrajectory(p0,
            {static_field, material.muonium_sigma_T, material.muonium_hop_rate, material.hyperfine},
            time, rng, material.max_hops);
    }();
    const Vec3 decay_p = std::visit([&](const auto& track) { return track.polarization(time); }, trajectory);
    const auto positron = sample_muplus_michel(decay_p, rng);
    return {1, muonium ? SpinState::Muonium : SpinState::MuPlus, time, p0, decay_p, positron, std::move(trajectory)};
}

} // namespace stopped_muon
