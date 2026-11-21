#pragma once

#include "spin_dynamic_decay.hxx"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace stopped_muon {

enum class CopperIsotopeMode {
    Natural,
    Cu63,
    Cu65
};

inline CopperIsotopeMode parse_copper_isotope_mode(const std::string& text) {
    if (text == "natural" || text == "nat") return CopperIsotopeMode::Natural;
    if (text == "63" || text == "cu63" || text == "Cu63") return CopperIsotopeMode::Cu63;
    if (text == "65" || text == "cu65" || text == "Cu65") return CopperIsotopeMode::Cu65;
    throw std::invalid_argument(
        "未知铜同位素模式: " + text + "（仅支持 natural, 63 或 65）");
}

// Shared declarative material configuration used by both Geant4 Mode 1 integration and CLI Mode 3 batch processing.
struct MaterialConfig {
    std::string preset = "copper";
    bool field_explicit = false;
    Vec3 external_field{-1762.0e-9, 37638.2e-9, -26120.2e-9};

    // mu+ overrides
    std::optional<double> mu_plus_sigma_T;
    std::optional<double> mu_plus_hop_rate;
    std::optional<double> muonium_fraction;
    std::optional<double> muonium_sigma_T;
    std::optional<double> muonium_hop_rate;
    std::optional<double> muonium_hyperfine;
    std::optional<double> muonium_dipolar_D;
    std::optional<double> muonium_state2_hyperfine;
    std::optional<double> muonium_state2_dipolar_D;
    std::optional<double> muonium_exchange_w12;
    std::optional<double> muonium_exchange_w21;
    std::optional<double> muonium_spin_exchange_rate;
    std::optional<double> muonium_reorientation_rate;
    std::optional<std::size_t> max_hops;

    // Ferromagnetic domain overrides (domain mode/dir shared by mu+ and I=0 mu- in the same target)
    std::optional<DomainFieldMode> domain_mode;
    std::optional<Vec3> domain_dir;
    std::optional<double> mu_plus_internal_field_T;

    // mu- (Copper I=3/2 and Fe56/S32 I=0) overrides
    CopperIsotopeMode copper_isotope = CopperIsotopeMode::Natural;
    bool cu_map_B = false;
    std::optional<double> f1_pop;
    std::optional<double> f1_orient;
    std::optional<double> f2_orient;
    std::optional<double> cu_hyperfine;
    std::optional<double> m1_rate_ns;
    std::optional<double> residual_pol;
    std::optional<double> mu_minus_internal_field_T;
    std::optional<double> mu_minus_sigma_T;
    std::optional<double> mu_minus_hop_rate;
    std::optional<double> dio_rate_ns;
    std::optional<double> capture_rate_ns;
};

struct MaterialPair {
    MuPlusMaterial mu_plus{};
    std::optional<AnyMuMinusMaterial> mu_minus{};
    Vec3 effective_field{0, 0, 0};
};

inline AnyMuMinusMaterial make_copper_mu_minus_for_mode(
        CopperIsotopeMode mode,
        const Vec3& field,
        bool map_B,
        double m1_rate_ns) {
    if (mode == CopperIsotopeMode::Natural) {
        return make_natural_copper_mu_minus_conditional(field, map_B, m1_rate_ns);
    }
    const int isotope = (mode == CopperIsotopeMode::Cu63) ? 63 : 65;
    return make_copper_mu_minus_conditional(isotope, field, map_B, m1_rate_ns);
}

inline MaterialPair build_material_pair(const MaterialConfig& cfg) {
    MaterialPair out{};
    Vec3 field = cfg.external_field;
    if (cfg.preset == "free") {
        if (!cfg.field_explicit) field = {0, 0, 0};
        out.mu_plus = MuPlusMaterial{};
        out.mu_plus.external_field = field;
        out.mu_minus.reset();
    } else if (cfg.preset == "copper") {
        out.mu_plus = make_copper_mu_plus_example(field);
        out.mu_minus = make_copper_mu_minus_for_mode(
            cfg.copper_isotope, field, cfg.cu_map_B, cfg.m1_rate_ns.value_or(0.5));
    } else if (cfg.preset == "iron" || cfg.preset == "iron_pure") {
        out.mu_plus = make_iron_mu_plus_example(field, false);
        out.mu_minus = make_iron_mu_minus_conditional(field);
    } else if (cfg.preset == "iron_commercial") {
        out.mu_plus = make_iron_mu_plus_example(field, true);
        out.mu_minus = make_iron_mu_minus_conditional(field);
    } else if (cfg.preset == "sulfur") {
        out.mu_plus = make_sulfur_mu_plus_conditional(field);
        out.mu_minus = make_sulfur_mu_minus_conditional(field);
    } else if (cfg.preset == "quartz") {
        out.mu_plus = make_quartz_mu_plus_example(field);
        out.mu_minus.reset();
    } else if (cfg.preset == "glass") {
        out.mu_plus = MuPlusMaterial{};
        out.mu_plus.external_field = field;
        out.mu_plus.muonium_fraction = 0.30;
        out.mu_plus.ordinary_sigma_T = 0.5e-3 / gamma_muon;
        out.mu_plus.ordinary_hop_rate = 0.0;
        out.mu_plus.muonium_sigma_T = 0.5e-3 / gamma_muon;
        out.mu_plus.muonium_hop_rate = 0.0;
        out.mu_minus.reset();
    } else {
        throw std::invalid_argument(
            "未知材料预设: " + cfg.preset +
            "（可选: copper, iron, iron_pure, iron_commercial, sulfur, quartz, glass, free）");
    }

    if (cfg.mu_plus_sigma_T) out.mu_plus.ordinary_sigma_T = *cfg.mu_plus_sigma_T;
    if (cfg.mu_plus_hop_rate) out.mu_plus.ordinary_hop_rate = *cfg.mu_plus_hop_rate;
    if (cfg.muonium_fraction) out.mu_plus.muonium_fraction = *cfg.muonium_fraction;
    if (cfg.muonium_sigma_T) out.mu_plus.muonium_sigma_T = *cfg.muonium_sigma_T;
    if (cfg.muonium_hop_rate) out.mu_plus.muonium_hop_rate = *cfg.muonium_hop_rate;
    if (cfg.muonium_hyperfine) out.mu_plus.hyperfine = *cfg.muonium_hyperfine;
    if (cfg.muonium_dipolar_D) out.mu_plus.muonium_dipolar_D = *cfg.muonium_dipolar_D;
    if (cfg.muonium_state2_hyperfine) out.mu_plus.muonium_state2_hyperfine = *cfg.muonium_state2_hyperfine;
    if (cfg.muonium_state2_dipolar_D) out.mu_plus.muonium_state2_dipolar_D = *cfg.muonium_state2_dipolar_D;
    if (cfg.muonium_exchange_w12) out.mu_plus.muonium_exchange_w12 = *cfg.muonium_exchange_w12;
    if (cfg.muonium_exchange_w21) out.mu_plus.muonium_exchange_w21 = *cfg.muonium_exchange_w21;
    if (cfg.muonium_spin_exchange_rate) out.mu_plus.muonium_spin_exchange_rate = *cfg.muonium_spin_exchange_rate;
    if (cfg.muonium_reorientation_rate) out.mu_plus.muonium_reorientation_rate = *cfg.muonium_reorientation_rate;
    if (cfg.max_hops) out.mu_plus.max_hops = *cfg.max_hops;
    if (cfg.domain_mode) out.mu_plus.domain_field.mode = *cfg.domain_mode;
    if (cfg.domain_dir) out.mu_plus.domain_field.magnetization_dir = *cfg.domain_dir;
    if (cfg.mu_plus_internal_field_T) out.mu_plus.domain_field.internal_field_T = *cfg.mu_plus_internal_field_T;

    const bool has_cu_override = cfg.f1_pop.has_value() || cfg.f1_orient.has_value() ||
        cfg.f2_orient.has_value() || cfg.cu_hyperfine.has_value() ||
        cfg.m1_rate_ns.has_value() || cfg.dio_rate_ns.has_value() || cfg.capture_rate_ns.has_value();
    if (!out.mu_minus.has_value() && has_cu_override) {
        out.mu_minus = make_copper_mu_minus_for_mode(
            cfg.copper_isotope, field, cfg.cu_map_B, cfg.m1_rate_ns.value_or(0.5));
    }

    auto apply_cu_overrides = [&](MuMinusMaterial& cu) {
        cu.spin.field_T = field;
        cu.spin.gamma_nucleus_rad_ns_T =
            (cu.isotope == 63) ? 0.07100118490581576 : 0.0760556158609357;
        if (cfg.f1_pop) cu.preparation.population_F1 = *cfg.f1_pop;
        if (cfg.f1_orient) cu.preparation.orientation_F1 = *cfg.f1_orient;
        if (cfg.f2_orient) cu.preparation.orientation_F2 = *cfg.f2_orient;
        if (cfg.cu_hyperfine) cu.spin.hyperfine_rad_ns = *cfg.cu_hyperfine;
        if (cfg.m1_rate_ns) cu.conversion_rate_ns = *cfg.m1_rate_ns;
        if (cfg.dio_rate_ns) cu.rates.dio_ns = *cfg.dio_rate_ns;
        if (cfg.capture_rate_ns) cu.rates.capture_ns = *cfg.capture_rate_ns;
    };

    if (out.mu_minus.has_value()) {
        if (auto* cu = std::get_if<MuMinusMaterial>(&*out.mu_minus)) {
            apply_cu_overrides(*cu);
        } else if (auto* nat_cu = std::get_if<NaturalCopperMuMinusMaterial>(&*out.mu_minus)) {
            apply_cu_overrides(nat_cu->cu63);
            apply_cu_overrides(nat_cu->cu65);
        } else if (auto* s0 = std::get_if<MuMinusSpin0Material>(&*out.mu_minus)) {
            s0->external_field = field;
            if (cfg.domain_mode) s0->domain_field.mode = *cfg.domain_mode;
            if (cfg.domain_dir) s0->domain_field.magnetization_dir = *cfg.domain_dir;
            if (cfg.mu_minus_internal_field_T) s0->domain_field.internal_field_T = *cfg.mu_minus_internal_field_T;
            if (cfg.residual_pol) s0->residual_polarization = *cfg.residual_pol;
            if (cfg.mu_minus_sigma_T) s0->local_sigma_T = *cfg.mu_minus_sigma_T;
            if (cfg.mu_minus_hop_rate) s0->hop_rate = *cfg.mu_minus_hop_rate;
            if (cfg.dio_rate_ns) s0->rates.dio_ns = *cfg.dio_rate_ns;
            if (cfg.capture_rate_ns) s0->rates.capture_ns = *cfg.capture_rate_ns;
            if (cfg.max_hops) s0->max_hops = *cfg.max_hops;
        }
    }

    out.effective_field = field;
    return out;
}

} // namespace stopped_muon
