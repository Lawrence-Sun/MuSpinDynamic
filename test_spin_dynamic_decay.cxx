#include "spin_dynamic_decay.hxx"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

using namespace stopped_muon;

void check(bool cond, const char* msg) {
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

template <typename Fn>
void check_throws(Fn&& fn, const char* msg) {
    bool threw = false;
    try {
        fn();
    } catch (const std::exception&) {
        threw = true;
    }
    if (!threw) {
        throw std::runtime_error(msg);
    }
}

} // namespace

int test_spin_dynamic_decay() {
    std::size_t assertions = 0;

    // 1. Input validation checks
    {
        TRandom3 rng(101);
        const auto mat_plus = make_copper_mu_plus_example();
        const auto mat_minus = make_copper_mu_minus_conditional(63);

        check_throws([&]() {
            simulate_stopped_muon_decay({0, {NAN, 0, 0}, {0, 0, 1}, +1, 0}, mat_plus, mat_minus, rng);
        }, "non-finite stop_position should throw");
        ++assertions;

        check_throws([&]() {
            simulate_stopped_muon_decay({0, {0, 0, 0}, {0, 0, 1.01}, +1, 0}, mat_plus, mat_minus, rng);
        }, "|P0| > 1 should throw");
        ++assertions;

        check_throws([&]() {
            simulate_stopped_muon_decay({0, {0, 0, 0}, {0, 0, 1}, 0, 0}, mat_plus, mat_minus, rng);
        }, "charge == 0 should throw");
        ++assertions;

        check_throws([&]() {
            simulate_stopped_muon_decay({0, {0, 0, 0}, {0, 0, -1}, -1, 0}, mat_plus, std::nullopt, rng);
        }, "mu- without MuMinusMaterial should throw");
        ++assertions;
    }

    // 2. Deterministic replay & position/time propagation for mu+ (copper and quartz with muonium)
    {
        const StoppedMuonState st{42, {12.5, -3.25, 48.0}, {0.1, -0.2, 0.8}, +1, 150.0};
        const auto mat_cu = make_copper_mu_plus_example({0, 100e-6, 0});
        const auto mat_qz = make_quartz_mu_plus_example({0, 100e-6, 0});

        TRandom3 rng1(2026), rng2(2026);
        const auto r1 = simulate_stopped_muon_decay(st, mat_cu, std::nullopt, rng1, 7);
        const auto r2 = simulate_stopped_muon_decay(st, mat_cu, std::nullopt, rng2, 7);

        check(r1.muon_entry == 7 && r1.event_id == 42, "muon_entry or event_id mismatch");
        check(r1.stop_position.X() == 12.5 && r1.stop_position.Y() == -3.25 && r1.stop_position.Z() == 48.0,
              "stop_position must be preserved");
        check(r1.charge == +1 && std::string(r1.spin_category) == "mu_plus", "unexpected mu+ category");
        check(std::abs(r1.decay_time_ns - (150.0 + r1.elapsed_time_ns)) < 1e-12, "decay_time_ns mismatch");
        check(r1.has_decay_lepton && r1.lepton_pdg == -11 && r1.lepton_charge == +1,
              "mu+ must emit positron (PDG -11)");
        check(r1.total_energy_MeV > 0 && r1.total_energy_MeV <= 0.5 * muon_mass_MeV, "positron energy out of range");
        check(std::abs(r1.momentum_mag_MeV_c - r1.total_energy_MeV) < 1e-12, "massless Michel |p| == E violated");
        check(std::abs(r1.direction.R() - 1.0) < 1e-12, "direction must be unit vector");
        check(r1.total_energy_MeV == r2.total_energy_MeV &&
              r1.momentum_MeV_c.X() == r2.momentum_MeV_c.X() &&
              r1.decay_time_ns == r2.decay_time_ns,
              "fixed-seed replay mismatch for mu+");
        assertions += 9;

        bool saw_muonium = false, saw_free = false;
        for (int i = 0; i < 100; ++i) {
            const auto rq = simulate_stopped_muon_decay(st, mat_qz, std::nullopt, rng1, i);
            if (std::string(rq.spin_category) == "muonium") saw_muonium = true;
            if (std::string(rq.spin_category) == "mu_plus") saw_free = true;
            check(rq.has_decay_lepton && rq.lepton_pdg == -11, "quartz mu+ event missing positron");
            ++assertions;
        }
        check(saw_muonium && saw_free, "quartz preset should sample both muonium and free mu+");
        ++assertions;
    }

    // 3. Copper mu- (Cu63 and Cu65): competitive termination & massive DIO electron kinematics
    {
        for (int isotope : {63, 65}) {
            TRandom3 rng(3000 + isotope);
            const StoppedMuonState st{99, {-4.0, 7.5, -10.0}, {0, 0, -1}, -1, 25.0};
            const auto mat_minus = make_copper_mu_minus_conditional(isotope, {0, 100e-6, 0}, false, 0.5);
            const auto mat_plus = make_copper_mu_plus_example({0, 100e-6, 0});

            int n_dio = 0, n_cap = 0;
            for (int i = 0; i < 2000; ++i) {
                const auto r = simulate_stopped_muon_decay(st, mat_plus, mat_minus, rng, i);
                check(r.stop_position.X() == -4.0 && r.stop_position.Y() == 7.5 && r.stop_position.Z() == -10.0,
                      "mu- stop_position mismatch");
                check(r.charge == -1, "mu- charge mismatch");
                check(std::abs(r.decay_time_ns - (25.0 + r.elapsed_time_ns)) < 1e-12,
                      "mu- decay_time_ns mismatch");
                if (std::string(r.termination) == "DIO") {
                    ++n_dio;
                    check(r.has_decay_lepton && r.lepton_pdg == 11 && r.lepton_charge == -1,
                          "DIO must emit electron (PDG 11)");
                    check(r.total_energy_MeV > cu_dio_detail::electron_mass_MeV &&
                          r.total_energy_MeV < r.endpoint_energy_MeV,
                          "DIO energy out of [m_e, W]");
                    check(std::abs(r.kinetic_energy_MeV - (r.total_energy_MeV - cu_dio_detail::electron_mass_MeV)) < 1e-12,
                          "DIO kinetic energy mismatch");
                    const double expected_p = std::sqrt(
                        r.total_energy_MeV * r.total_energy_MeV -
                        cu_dio_detail::electron_mass_MeV * cu_dio_detail::electron_mass_MeV);
                    check(std::abs(r.momentum_mag_MeV_c - expected_p) < 1e-11,
                          "DIO relativistic energy-momentum relation violated");
                    check(std::abs(r.direction.R() - 1.0) < 1e-12, "DIO direction unit norm violated");
                } else {
                    ++n_cap;
                    check(std::string(r.termination) == "nuclear_capture", "unexpected mu- termination label");
                    check(!r.has_decay_lepton && r.lepton_pdg == 0 && std::isnan(r.total_energy_MeV),
                          "nuclear_capture must not emit decay electron");
                }
                assertions += 5;
            }
            check(n_dio > 80 && n_dio < 230 && n_cap + n_dio == 2000,
                  "Cu mu- DIO branching count outside expected statistical range");
            ++assertions;
        }
    }

    // 4. Mode B ROOT TTree round-trip:
    //    Input .root (1 entry = 1 muon) -> Output .root (1 entry = 1 Michel/DIO electron)
    {
        const std::filesystem::path tmp_dir = "spin_dynamic_decay/build/tmp_root_test";
        std::filesystem::create_directories(tmp_dir);
        const std::string in_root = (tmp_dir / "input_muons.root").string();
        const std::string out_root = (tmp_dir / "output_electrons.root").string();
        std::filesystem::remove(in_root);
        std::filesystem::remove(out_root);

        std::vector<StoppedMuonState> input_muons;
        constexpr int n_plus = 120;
        constexpr int n_minus = 280;
        for (int i = 0; i < n_plus + n_minus; ++i) {
            const bool is_plus = (i < n_plus);
            StoppedMuonState m{};
            m.event_id = 1000 + i;
            m.stop_position = {0.5 * i, -1.25 * i, 2.0 * i + 1.0};
            m.initial_polarization = is_plus ? Vec3{0.0, 0.0, 0.9} : Vec3{0.0, 0.0, -0.9};
            m.charge = is_plus ? +1 : -1;
            m.stop_time_ns = 10.0 * i;
            input_muons.push_back(m);
        }

        write_stopped_muons_root(in_root, input_muons, true);
        check_throws([&]() {
            write_stopped_muons_root(in_root, input_muons, false);
        }, "overwrite=false should reject existing input root file");
        ++assertions;

        const auto loaded_muons = read_stopped_muons_root(in_root);
        check(loaded_muons.size() == input_muons.size(), "loaded muon entry count mismatch");
        ++assertions;

        const auto mat_plus = make_copper_mu_plus_example({0, 100e-6, 0});
        const auto mat_minus = make_copper_mu_minus_conditional(63, {0, 100e-6, 0}, false, 0.5);
        TRandom3 rng(777);

        std::vector<DecaySimulationResult> sim_results;
        sim_results.reserve(loaded_muons.size());
        Long64_t expected_electrons = 0;
        for (std::size_t idx = 0; idx < loaded_muons.size(); ++idx) {
            check(loaded_muons[idx].event_id == input_muons[idx].event_id, "loaded event_id mismatch");
            const auto res = simulate_stopped_muon_decay(
                loaded_muons[idx], mat_plus, mat_minus, rng, static_cast<Long64_t>(idx));
            if (res.has_decay_lepton) ++expected_electrons;
            sim_results.push_back(res);
            ++assertions;
        }

        const Long64_t written_electrons = write_michel_electrons_root(
            out_root, sim_results, true, "michel_electrons", 777, "copper", {0, 100e-6, 0});
        check(written_electrons == expected_electrons, "written electron count mismatch");
        check(written_electrons > n_plus && written_electrons < n_plus + n_minus,
              "written electron count should include all mu+ plus DIO subset of mu-");
        assertions += 2;

        std::unique_ptr<TFile> f(TFile::Open(out_root.c_str(), "READ"));
        check(f && !f->IsZombie(), "failed to open output ROOT file");
        auto* tree = dynamic_cast<TTree*>(f->Get("michel_electrons"));
        check(tree != nullptr, "missing michel_electrons TTree in output ROOT file");
        check(tree->GetEntries() == written_electrons,
              "output TTree GetEntries() must equal number of emitted Michel/DIO electrons");

        auto* p_in = dynamic_cast<TParameter<Long64_t>*>(f->Get("n_input_muons"));
        auto* p_el = dynamic_cast<TParameter<Long64_t>*>(f->Get("n_michel_electrons"));
        auto* p_nc = dynamic_cast<TParameter<Long64_t>*>(f->Get("n_nuclear_captures"));
        check(p_in && p_el && p_nc, "missing summary TParameters in output ROOT file");
        check(p_in->GetVal() == static_cast<Long64_t>(input_muons.size()), "n_input_muons mismatch");
        check(p_el->GetVal() == written_electrons, "n_michel_electrons mismatch");
        check(p_in->GetVal() == p_el->GetVal() + p_nc->GetVal(), "branching closure failed in output ROOT file");
        assertions += 7;

        Long64_t electron_id = -1, muon_entry = -1, event_id = -1;
        Int_t pdg = 0, charge = 0, muon_charge = 0;
        Double_t x = 0, y = 0, z = 0;
        Double_t total_energy_MeV = 0, kinetic_energy_MeV = 0, momentum_MeV_c = 0;
        Double_t px = 0, py = 0, pz = 0;

        tree->SetBranchAddress("electron_id", &electron_id);
        tree->SetBranchAddress("muon_entry", &muon_entry);
        tree->SetBranchAddress("event_id", &event_id);
        tree->SetBranchAddress("pdg", &pdg);
        tree->SetBranchAddress("charge", &charge);
        tree->SetBranchAddress("muon_charge", &muon_charge);
        tree->SetBranchAddress("x", &x);
        tree->SetBranchAddress("y", &y);
        tree->SetBranchAddress("z", &z);
        tree->SetBranchAddress("total_energy_MeV", &total_energy_MeV);
        tree->SetBranchAddress("kinetic_energy_MeV", &kinetic_energy_MeV);
        tree->SetBranchAddress("momentum_MeV_c", &momentum_MeV_c);
        tree->SetBranchAddress("px_MeV_c", &px);
        tree->SetBranchAddress("py_MeV_c", &py);
        tree->SetBranchAddress("pz_MeV_c", &pz);

        for (Long64_t e = 0; e < tree->GetEntries(); ++e) {
            check(tree->GetEntry(e) > 0, "failed to read entry from michel_electrons");
            check(electron_id == e, "electron_id must match output entry index");
            check(muon_entry >= 0 && muon_entry < static_cast<Long64_t>(input_muons.size()),
                  "muon_entry out of bounds");
            const auto& parent = input_muons[static_cast<std::size_t>(muon_entry)];
            check(event_id == parent.event_id, "parent event_id mismatch");
            check(muon_charge == parent.charge, "parent muon_charge mismatch");
            check(x == parent.stop_position.X() &&
                  y == parent.stop_position.Y() &&
                  z == parent.stop_position.Z(),
                  "electron birth vertex (x,y,z) must match parent muon stop_position");
            check((pdg == -11 && charge == +1 && muon_charge == +1) ||
                  (pdg == +11 && charge == -1 && muon_charge == -1),
                  "inconsistent pdg/charge/muon_charge in output entry");
            check(std::isfinite(total_energy_MeV) && total_energy_MeV > 0 &&
                  std::isfinite(kinetic_energy_MeV) && kinetic_energy_MeV > 0 &&
                  std::isfinite(momentum_MeV_c) && momentum_MeV_c > 0 &&
                  std::isfinite(px) && std::isfinite(py) && std::isfinite(pz),
                  "non-finite kinematics in michel_electrons entry");
            assertions += 8;
        }
        f->Close();
        std::filesystem::remove_all(tmp_dir);
    }

    // 5. Iron (pure & commercial) and Sulfur mu+ / mu- end-to-end kinematics & ROOT Mode B round-trip
    {
        // 5a. Iron pure & commercial mu+ and Fe56 mu-
        for (bool commercial : {false, true}) {
            const auto fe_plus = make_iron_mu_plus_example({0, 0, 0}, commercial);
            const auto fe_minus = make_iron_mu_minus_conditional({0, 0, 0});
            TRandom3 rng(commercial ? 5002 : 5001);
            const StoppedMuonState st_plus{10, {1.0, 2.0, 3.0}, {0, 0, 1}, +1, 5.0};
            const StoppedMuonState st_minus{11, {-1.0, -2.0, -3.0}, {0, 0, -1}, -1, 5.0};

            for (int i = 0; i < 50; ++i) {
                const auto rp = simulate_stopped_muon_decay(st_plus, fe_plus, fe_minus, rng, i);
                check(rp.has_decay_lepton && rp.lepton_pdg == -11 && std::string(rp.spin_category) == "mu_plus",
                      "iron mu+ event invalid");
                check(std::abs(rp.direction.R() - 1.0) < 1e-12, "iron mu+ direction unit norm violated");
                assertions += 2;
            }
            int n_dio = 0, n_cap = 0;
            for (int i = 0; i < 2000; ++i) {
                const auto rm = simulate_stopped_muon_decay(st_minus, fe_plus, fe_minus, rng, i);
                check(std::string(rm.spin_category) == "fe56_mu_minus", "expected fe56_mu_minus spin_category");
                if (std::string(rm.termination) == "DIO") {
                    ++n_dio;
                    check(rm.has_decay_lepton && rm.lepton_pdg == 11 &&
                          rm.total_energy_MeV > bound_dio_detail::electron_mass_MeV &&
                          rm.total_energy_MeV < rm.endpoint_energy_MeV,
                          "Fe56 DIO kinematics out of range");
                } else {
                    ++n_cap;
                    check(!rm.has_decay_lepton && rm.lepton_pdg == 0, "Fe56 capture must not emit lepton");
                }
                assertions += 2;
            }
            // Expected Fe56 DIO fraction = 9.14% -> ~183 / 2000
            check(n_dio > 120 && n_dio < 260 && n_dio + n_cap == 2000,
                  "Fe56 mu- DIO branching count outside expected statistical range");
            ++assertions;
        }

        // 5b. Solid sulfur mu+ (Cox 2011 two-state chemical exchange, f_D = 0) and S32 mu- (Mamedov 2001)
        {
            const auto s_plus = make_sulfur_mu_plus_conditional({0, 0, 0});
            const auto s_minus = make_sulfur_mu_minus_conditional({0, 0, 0});
            check(s_plus.muonium_fraction == 1.0 && s_plus.muonium_exchange_w12 == 28.0 &&
                  s_plus.muonium_exchange_w21 == 15.0,
                  "sulfur mu+ preset must use Cox 2011 f_D=0 two-state exchange parameters");
            check(std::abs(s_minus.residual_polarization - 0.073) < 1e-12 &&
                  std::abs(s_minus.gyromagnetic_ratio - gamma_muon * (1.0 - 4.24e-3)) < 1e-12,
                  "sulfur mu- preset must use Mamedov 2001 f_res=0.073 and bound g-factor");
            assertions += 2;

            TRandom3 rng(6001);
            const StoppedMuonState st_plus{20, {0.5, 0.5, 0.5}, {0, 0, 1}, +1, 0.0};
            const StoppedMuonState st_minus{21, {0.5, 0.5, 0.5}, {0, 0, -1}, -1, 0.0};

            int n_rad = 0;
            for (int i = 0; i < 300; ++i) {
                const auto rp = simulate_stopped_muon_decay(st_plus, s_plus, s_minus, rng, i);
                check(rp.has_decay_lepton && rp.lepton_pdg == -11, "sulfur mu+ missing positron");
                check(std::string(rp.spin_category) == "anisotropic_muonium" && rp.hop_count == 0,
                      "sulfur mu+ Cox 2011 two-state exchange should have 100% anisotropic_muonium and 0 stochastic hops");
                ++n_rad;
                ++assertions;
            }
            check(n_rad == 300, "sulfur mu+ count mismatch");
            ++assertions;

            int n_dio = 0, n_cap = 0;
            for (int i = 0; i < 2000; ++i) {
                const auto rm = simulate_stopped_muon_decay(st_minus, s_plus, s_minus, rng, i);
                check(std::string(rm.spin_category) == "s32_mu_minus", "expected s32_mu_minus spin_category");
                check(std::abs(rm.decay_polarization.Z() - (-0.073)) < 1e-12,
                      "S32 zero-field mu- decay_polarization should equal -0.073");
                if (std::string(rm.termination) == "DIO") {
                    ++n_dio;
                    check(rm.has_decay_lepton && rm.lepton_pdg == 11 &&
                          rm.total_energy_MeV > bound_dio_detail::electron_mass_MeV &&
                          rm.total_energy_MeV < rm.endpoint_energy_MeV,
                          "S32 DIO kinematics out of range");
                } else {
                    ++n_cap;
                    check(!rm.has_decay_lepton && rm.lepton_pdg == 0, "S32 capture must not emit lepton");
                }
                assertions += 3;
            }
            // Expected S32 DIO fraction = 25.00% -> ~500 / 2000
            check(n_dio > 420 && n_dio < 580 && n_dio + n_cap == 2000,
                  "S32 mu- DIO branching count outside expected statistical range");
            ++assertions;
        }
    }

    std::cout << "PASS test_spin_dynamic_decay (" << assertions << " assertions)\n";
    return 0;
}

#ifndef __CLING__
int main() {
    try {
        return test_spin_dynamic_decay();
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
#endif
