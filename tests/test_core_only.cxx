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
            check(r.total_energy_MeV > 0.0 && r.total_energy_MeV <= r.endpoint_energy_MeV + 1e-9,
                  "valid Michel energy");
            check(std::abs(r.direction.R() - 1.0) < 1e-9, "unit lepton direction");
            assertions += 6;
        }
    }

    // 2. mu- presets: Cu-63 (I=3/2), Fe-56 (I=0, ferromagnetic), S-32 (I=0)
    const mu_spin_dynamic::AnyMuMinusMaterial mu_minus_presets[] = {
        mu_spin_dynamic::make_copper_mu_minus_conditional(63, B_ext),
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

    // 3. UniformRngAdapter statistical sanity (Exp mean & Gaus moments)
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
