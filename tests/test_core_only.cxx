#include "mu_spin_dynamic/material_config.hxx"
#include "mu_spin_dynamic/spin_dynamic_decay.hxx"

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

// Verify that including MuSpinDynamic::Core does NOT pollute the global
// namespace with CLHEP/Geant4-conflicting unit names or vector aliases.
namespace {
constexpr double ns = 123.0;
constexpr double us = 456.0;
constexpr double T = 789.0;
constexpr double s = 321.0;
constexpr double pi = 3.0;
struct Vec3 { int sentinel = 42; };

void check(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}
} // namespace

int main() {
    // Silence unused warnings while proving the global/anonymous symbols compile cleanly.
    check(ns + us + T + s + pi > 0.0 && Vec3{}.sentinel == 42, "namespace isolation check");

    std::mt19937_64 engine(20261120ULL);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    auto rng = mu_spin_dynamic::make_uniform_rng([&]() { return dist(engine); });

    const mu_spin_dynamic::Vec3 B_ext{0.0, 100.0e-6, 0.0};
    const mu_spin_dynamic::StoppedMuonState mu_plus{
        1, {1.0, -2.0, 3.0}, {0.0, 0.0, 1.0}, +1, 10.0};
    const mu_spin_dynamic::StoppedMuonState mu_minus{
        2, {1.0, -2.0, 3.0}, {0.0, 0.0, 1.0}, -1, 10.0};

    std::size_t assertions = 0;

    // 1. mu+ presets: copper, quartz (muonium), iron (domain field), sulfur (Cox 2-state exchange)
    const mu_spin_dynamic::MuPlusMaterial mu_plus_presets[] = {
        mu_spin_dynamic::make_copper_mu_plus_example(B_ext),
        mu_spin_dynamic::make_quartz_mu_plus_example(B_ext),
        mu_spin_dynamic::make_iron_mu_plus_example(B_ext, false),
        mu_spin_dynamic::make_sulfur_mu_plus_conditional(B_ext),
    };
    for (const auto& mat : mu_plus_presets) {
        for (int i = 0; i < 25; ++i) {
            const auto r = mu_spin_dynamic::simulate_mu_plus_decay(mu_plus, mat, rng, i);
            check(r.has_decay_lepton, "mu+ must emit Michel positron");
            check(r.lepton_pdg == -11 && r.lepton_charge == +1, "mu+ lepton must be e+");
            check(r.elapsed_time_ns >= 0.0 && std::isfinite(r.decay_time_ns), "valid decay time");
            check(r.decay_polarization.R() <= 1.0 + 1e-6, "|P_decay| <= 1");
            check(r.total_energy_MeV > mu_spin_dynamic::electron_mass_MeV &&
                  r.total_energy_MeV <= mu_spin_dynamic::michel_endpoint_MeV + 1e-12,
                  "valid finite-me Michel energy");
            const double expected_p = std::sqrt(
                (r.total_energy_MeV - mu_spin_dynamic::electron_mass_MeV) *
                (r.total_energy_MeV + mu_spin_dynamic::electron_mass_MeV));
            check(std::abs(r.momentum_mag_MeV_c - expected_p) < 1e-11,
                  "finite-me Michel mass shell violated");
            check(std::abs(r.kinetic_energy_MeV - (r.total_energy_MeV - mu_spin_dynamic::electron_mass_MeV)) < 1e-12,
                  "finite-me Michel kinetic energy mismatch");
            check(std::abs(r.direction.R() - 1.0) < 1e-9, "unit lepton direction");
            assertions += 8;
        }
    }

    // 2. mu- presets: Cu-63, natural Cu (63/65), Fe-56 (I=0, ferromagnetic), S-32 (I=0)
    const mu_spin_dynamic::AnyMuMinusMaterial mu_minus_presets[] = {
        mu_spin_dynamic::make_copper_mu_minus_conditional(63, B_ext),
        mu_spin_dynamic::make_natural_copper_mu_minus_conditional(B_ext),
        mu_spin_dynamic::make_iron_mu_minus_conditional(B_ext),
        mu_spin_dynamic::make_sulfur_mu_minus_conditional(B_ext),
    };
    const auto dummy_plus = mu_spin_dynamic::make_copper_mu_plus_example(B_ext);
    for (const auto& mat : mu_minus_presets) {
        int dio_count = 0;
        int capture_count = 0;
        for (int i = 0; i < 200; ++i) {
            const auto r = mu_spin_dynamic::simulate_stopped_muon_decay(
                mu_minus, dummy_plus, mat, rng, i);
            check(r.charge == -1, "mu- charge preserved");
            check(r.elapsed_time_ns >= 0.0 && std::isfinite(r.decay_time_ns), "valid mu- time");
            check(r.decay_polarization.R() <= 1.0 + 1e-6, "valid mu- polarization");
            if (r.has_decay_lepton) {
                ++dio_count;
                check(r.lepton_pdg == 11 && r.lepton_charge == -1, "DIO lepton must be e-");
                check(r.total_energy_MeV > 0.511 && r.total_energy_MeV <= r.endpoint_energy_MeV + 1e-9,
                      "valid DIO energy");
                check(std::abs(r.direction.R() - 1.0) < 1e-9, "unit DIO direction");
                assertions += 3;
            } else {
                ++capture_count;
            }
            assertions += 3;
        }
        check(dio_count > 0 && capture_count > 0, "both DIO and nuclear_capture channels sampled");
        ++assertions;
    }

    // 3. Spin-evolution-only interface (NLO hook) exact match with full simulation up to decay time
    {
        const auto cu_plus = mu_spin_dynamic::make_copper_mu_plus_example(B_ext);
        const auto nat_cu_minus = mu_spin_dynamic::make_natural_copper_mu_minus_conditional(B_ext);
        for (int i = 0; i < 40; ++i) {
            auto rng_spin = mu_spin_dynamic::make_event_replica_rng(9001ULL, 500 + i, 2);
            auto rng_full = mu_spin_dynamic::make_event_replica_rng(9001ULL, 500 + i, 2);
            const auto term = mu_spin_dynamic::evolve_stopped_muon_to_termination(
                mu_plus, cu_plus, nat_cu_minus, rng_spin, i, 2, 0.25);
            const auto full = mu_spin_dynamic::simulate_stopped_muon_decay(
                mu_plus, cu_plus, nat_cu_minus, rng_full, i, 2, 0.25);
            check(term.will_emit_lepton && full.has_decay_lepton, "mu+ must emit lepton");
            check(term.replica_id == 2 && full.replica_id == 2, "replica_id propagation");
            check(std::abs(term.weight - 0.25) < 1e-15 && std::abs(full.weight - 0.25) < 1e-15, "weight scaling");
            check(term.decay_time_ns == full.decay_time_ns &&
                  term.decay_polarization.X() == full.decay_polarization.X() &&
                  term.decay_polarization.Y() == full.decay_polarization.Y() &&
                  term.decay_polarization.Z() == full.decay_polarization.Z(),
                  "spin-only interface must match full decay time and polarization bit-for-bit");
            assertions += 4;
        }
        for (int i = 0; i < 40; ++i) {
            auto rng_spin = mu_spin_dynamic::make_event_replica_rng(9002ULL, 600 + i, 0);
            auto rng_full = mu_spin_dynamic::make_event_replica_rng(9002ULL, 600 + i, 0);
            const auto term = mu_spin_dynamic::evolve_stopped_muon_to_termination(
                mu_minus, cu_plus, nat_cu_minus, rng_spin, i, 0, 1.0);
            const auto full = mu_spin_dynamic::simulate_stopped_muon_decay(
                mu_minus, cu_plus, nat_cu_minus, rng_full, i, 0, 1.0);
            check(term.will_emit_lepton == full.has_decay_lepton &&
                  term.decay_time_ns == full.decay_time_ns &&
                  std::string(term.spin_category) == std::string(full.spin_category) &&
                  std::string(term.termination) == std::string(full.termination),
                  "mu- spin-only interface must match full termination outcome");
            ++assertions;
        }
    }

    // 4. Shared MaterialConfig / build_material_pair covers all presets (default copper is natural)
    {
        for (const char* name : {"copper", "iron", "iron_pure", "iron_commercial", "sulfur", "quartz", "glass", "free"}) {
            mu_spin_dynamic::MaterialConfig cfg{};
            cfg.preset = name;
            const auto pair = mu_spin_dynamic::build_material_pair(cfg);
            if (std::string(name) == "copper") {
                check(pair.mu_minus.has_value() &&
                      std::holds_alternative<mu_spin_dynamic::NaturalCopperMuMinusMaterial>(*pair.mu_minus),
                      "default copper preset must use NaturalCopperMuMinusMaterial");
                ++assertions;
            }
            const auto r = mu_spin_dynamic::simulate_mu_plus_decay(mu_plus, pair.mu_plus, rng, 0);
            check(r.has_decay_lepton, "MaterialConfig preset mu+ simulation failed");
            ++assertions;
        }
    }

    // 5. UniformRngAdapter statistical sanity (Exp mean & Gaus moments)
    {
        constexpr int N = 10000;
        double sum_exp = 0.0;
        double sum_g = 0.0;
        double sum_g2 = 0.0;
        for (int i = 0; i < N; ++i) {
            sum_exp += rng.Exp(2.5);
            const double g = rng.Gaus(1.0, 2.0);
            sum_g += g;
            sum_g2 += (g - 1.0) * (g - 1.0);
        }
        const double mean_exp = sum_exp / N;
        const double mean_g = sum_g / N;
        const double var_g = sum_g2 / N;
        check(std::abs(mean_exp - 2.5) < 0.15, "Exp(tau) mean sanity");
        check(std::abs(mean_g - 1.0) < 0.10, "Gaus mean sanity");
        check(std::abs(var_g - 4.0) < 0.25, "Gaus variance sanity");
        assertions += 3;
    }

    std::cout << "[PASS] test_core_only: " << assertions
              << " assertions passed (zero TFile/TTree/TRandom3 dependency).\n";
    return 0;
}
