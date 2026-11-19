#pragma once

#include "mu_minus_event.hxx"
#include "mu_minus_spin0.hxx"
#include "stopped_decay.hxx"

#include "TBranch.h"
#include "TFile.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

#include <cmath>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace stopped_muon {

struct StoppedMuonState {
    Long64_t event_id = 0;              // Muon event identifier (defaults to entry index).
    Vec3 stop_position{0, 0, 0};        // Stopped muon position in target (e.g. mm), passed through as decay vertex.
    Vec3 initial_polarization{0, 0, 1}; // Laboratory-frame spin polarization P0 at stop time (|P0| <= 1).
    int charge = +1;                    // +1 for mu+, -1 for mu-.
    double stop_time_ns = 0;            // Stopping time t_stop [ns].
};

struct DecaySimulationResult {
    Long64_t muon_entry = 0;
    Long64_t event_id = 0;
    Vec3 stop_position{0, 0, 0};
    int charge = +1;
    const char* spin_category = "mu_plus";
    const char* termination = "decay";
    double stop_time_ns = 0;
    double elapsed_time_ns = 0;
    double decay_time_ns = 0;
    Vec3 initial_polarization{0, 0, 0};
    Vec3 decay_polarization{0, 0, 0};
    std::size_t hop_count = 0;

    bool has_decay_lepton = false;
    int lepton_pdg = 0;
    int lepton_charge = 0;
    double total_energy_MeV = NAN;
    double kinetic_energy_MeV = NAN;
    double momentum_mag_MeV_c = NAN;
    Vec3 momentum_MeV_c{NAN, NAN, NAN};
    Vec3 direction{NAN, NAN, NAN};
    double endpoint_energy_MeV = NAN;
    double reduced_energy_x = NAN;
    const char* spectrum_model = "none";
};

inline void validate_stopped_muon_state(const StoppedMuonState& state) {
    if (!finite(state.stop_position)) {
        throw std::invalid_argument("stop_position must have finite (x, y, z) coordinates");
    }
    if (!finite(state.initial_polarization) || !std::isfinite(state.initial_polarization.R()) ||
        state.initial_polarization.R() > 1) {
        throw std::invalid_argument("initial_polarization must be finite with |P0| <= 1");
    }
    if (!std::isfinite(state.stop_time_ns)) {
        throw std::invalid_argument("stop_time_ns must be finite");
    }
    if (state.charge != +1 && state.charge != -1) {
        throw std::invalid_argument("charge must be +1 (mu+) or -1 (mu-)");
    }
}

inline DecaySimulationResult simulate_mu_plus_decay(
        const StoppedMuonState& muon, const MuPlusMaterial& material, TRandom3& rng,
        Long64_t muon_entry = 0) {
    validate_stopped_muon_state(muon);
    if (muon.charge != +1) {
        throw std::invalid_argument("simulate_mu_plus_decay requires charge == +1");
    }
    const auto ev = simulate_stopped_muon(material, +1, muon.initial_polarization, rng);
    const double decay_t = muon.stop_time_ns + ev.decay_time_ns;
    if (!std::isfinite(decay_t)) {
        throw std::overflow_error("decay_time_ns overflow");
    }
    const double endpoint = 0.5 * muon_mass_MeV;
    const double e_tot = ev.positron.energy_MeV;
    const double p_mag = ev.positron.momentum_MeV_c.R();
    const Vec3 dir = p_mag > 0 ? ev.positron.momentum_MeV_c / p_mag : Vec3{NAN, NAN, NAN};

    DecaySimulationResult out{};
    out.muon_entry = muon_entry;
    out.event_id = muon.event_id;
    out.stop_position = muon.stop_position;
    out.charge = +1;
    if (std::holds_alternative<AnisotropicMuoniumTrajectory>(ev.trajectory) ||
        std::holds_alternative<TwoStateExchangeTrajectory>(ev.trajectory)) {
        out.spin_category = "anisotropic_muonium";
    } else if (ev.spin_state == SpinState::Muonium) {
        out.spin_category = "muonium";
    } else {
        out.spin_category = "mu_plus";
    }
    out.termination = ev.termination;
    out.stop_time_ns = muon.stop_time_ns;
    out.elapsed_time_ns = ev.decay_time_ns;
    out.decay_time_ns = decay_t;
    out.initial_polarization = ev.initial_polarization;
    out.decay_polarization = ev.decay_polarization;
    out.hop_count = ev.hop_count();
    out.has_decay_lepton = true;
    out.lepton_pdg = ev.positron.pdg;
    out.lepton_charge = +1;
    out.total_energy_MeV = e_tot;
    out.kinetic_energy_MeV = e_tot; // Massless LO Michel baseline: E_kin == E_tot.
    out.momentum_mag_MeV_c = p_mag;
    out.momentum_MeV_c = ev.positron.momentum_MeV_c;
    out.direction = dir;
    out.endpoint_energy_MeV = endpoint;
    out.reduced_energy_x = e_tot / endpoint;
    out.spectrum_model = ev.positron.model;
    return out;
}

inline DecaySimulationResult simulate_copper_mu_minus_decay(
        const StoppedMuonState& muon, const MuMinusMaterial& material, TRandom3& rng,
        Long64_t muon_entry = 0) {
    validate_stopped_muon_state(muon);
    if (muon.charge != -1) {
        throw std::invalid_argument("simulate_copper_mu_minus_decay requires charge == -1");
    }
    const CopperDioSpectrum spectrum(material.isotope);
    const auto ev = simulate_stopped_mu_minus(material, muon.initial_polarization, rng);
    const double decay_t = muon.stop_time_ns + ev.termination.time_ns;
    if (!std::isfinite(decay_t)) {
        throw std::overflow_error("decay_time_ns overflow");
    }

    DecaySimulationResult out{};
    out.muon_entry = muon_entry;
    out.event_id = muon.event_id;
    out.stop_position = muon.stop_position;
    out.charge = -1;
    out.spin_category = (ev.isotope == 63) ? "cu63_mu_minus" : "cu65_mu_minus";
    out.termination = (ev.termination.channel == MuMinusTerminationChannel::DIO)
        ? "DIO"
        : "nuclear_capture";
    out.stop_time_ns = muon.stop_time_ns;
    out.elapsed_time_ns = ev.termination.time_ns;
    out.decay_time_ns = decay_t;
    out.initial_polarization = muon.initial_polarization;
    out.decay_polarization = ev.termination.polarization();
    out.hop_count = 0;

    if (ev.electron.has_value()) {
        const auto& el = *ev.electron;
        const double e_tot = el.energy_MeV;
        const double e_kin = e_tot - cu_dio_detail::electron_mass_MeV;
        const double p_mag = el.momentum_MeV_c.R();
        const Vec3 dir = p_mag > 0 ? el.momentum_MeV_c / p_mag : Vec3{NAN, NAN, NAN};
        const double endpoint = spectrum.endpoint_MeV();

        out.has_decay_lepton = true;
        out.lepton_pdg = el.pdg;
        out.lepton_charge = -1;
        out.total_energy_MeV = e_tot;
        out.kinetic_energy_MeV = e_kin;
        out.momentum_mag_MeV_c = p_mag;
        out.momentum_MeV_c = el.momentum_MeV_c;
        out.direction = dir;
        out.endpoint_energy_MeV = endpoint;
        out.reduced_energy_x = e_tot / endpoint;
        out.spectrum_model = el.model;
    }
    return out;
}

inline DecaySimulationResult simulate_spin0_mu_minus_decay(
        const StoppedMuonState& muon, const MuMinusSpin0Material& material, TRandom3& rng,
        Long64_t muon_entry = 0) {
    validate_stopped_muon_state(muon);
    if (muon.charge != -1) {
        throw std::invalid_argument("simulate_spin0_mu_minus_decay requires charge == -1");
    }
    const BoundDioSpectrum spectrum(material.nuclide);
    const auto ev = simulate_stopped_mu_minus_spin0(material, muon.initial_polarization, rng);
    const double decay_t = muon.stop_time_ns + ev.time_ns;
    if (!std::isfinite(decay_t)) {
        throw std::overflow_error("decay_time_ns overflow");
    }

    DecaySimulationResult out{};
    out.muon_entry = muon_entry;
    out.event_id = muon.event_id;
    out.stop_position = muon.stop_position;
    out.charge = -1;
    out.spin_category = dio_nuclide_spin_category(ev.nuclide);
    out.termination = (ev.channel == MuMinusTerminationChannel::DIO)
        ? "DIO"
        : "nuclear_capture";
    out.stop_time_ns = muon.stop_time_ns;
    out.elapsed_time_ns = ev.time_ns;
    out.decay_time_ns = decay_t;
    out.initial_polarization = muon.initial_polarization;
    out.decay_polarization = ev.decay_polarization;
    out.hop_count = ev.hop_count();

    if (ev.electron.has_value()) {
        const auto& el = *ev.electron;
        const double e_tot = el.energy_MeV;
        const double e_kin = e_tot - bound_dio_detail::electron_mass_MeV;
        const double p_mag = el.momentum_MeV_c.R();
        const Vec3 dir = p_mag > 0 ? el.momentum_MeV_c / p_mag : Vec3{NAN, NAN, NAN};
        const double endpoint = spectrum.endpoint_MeV();

        out.has_decay_lepton = true;
        out.lepton_pdg = el.pdg;
        out.lepton_charge = -1;
        out.total_energy_MeV = e_tot;
        out.kinetic_energy_MeV = e_kin;
        out.momentum_mag_MeV_c = p_mag;
        out.momentum_MeV_c = el.momentum_MeV_c;
        out.direction = dir;
        out.endpoint_energy_MeV = endpoint;
        out.reduced_energy_x = e_tot / endpoint;
        out.spectrum_model = el.model;
    }
    return out;
}

using AnyMuMinusMaterial = std::variant<MuMinusMaterial, MuMinusSpin0Material>;

inline DecaySimulationResult simulate_stopped_muon_decay(
        const StoppedMuonState& muon,
        const MuPlusMaterial& mu_plus_material,
        const std::optional<AnyMuMinusMaterial>& mu_minus_material,
        TRandom3& rng,
        Long64_t muon_entry = 0) {
    validate_stopped_muon_state(muon);
    if (muon.charge == +1) {
        return simulate_mu_plus_decay(muon, mu_plus_material, rng, muon_entry);
    }
    if (!mu_minus_material.has_value()) {
        throw std::invalid_argument("charge == -1 requires an explicit MuMinusMaterial or MuMinusSpin0Material configuration");
    }
    if (std::holds_alternative<MuMinusMaterial>(*mu_minus_material)) {
        return simulate_copper_mu_minus_decay(
            muon, std::get<MuMinusMaterial>(*mu_minus_material), rng, muon_entry);
    }
    return simulate_spin0_mu_minus_decay(
        muon, std::get<MuMinusSpin0Material>(*mu_minus_material), rng, muon_entry);
}

// --- ROOT TTree I/O for Mode B (1 input entry = 1 muon; 1 output entry = 1 Michel/DIO electron) ---

namespace root_io_detail {

inline TLeaf* find_leaf_any(TTree* tree, std::initializer_list<const char*> names) {
    for (const char* name : names) {
        if (TBranch* b = tree->GetBranch(name)) {
            if (TLeaf* l = b->GetLeaf(name)) return l;
            if (b->GetListOfLeaves() && b->GetListOfLeaves()->GetEntries() > 0) {
                return static_cast<TLeaf*>(b->GetListOfLeaves()->At(0));
            }
        }
        if (TLeaf* l = tree->GetLeaf(name)) return l;
    }
    return nullptr;
}

inline TTree* resolve_input_tree(TFile& file, const std::string& requested_name) {
    if (!requested_name.empty()) {
        auto* t = dynamic_cast<TTree*>(file.Get(requested_name.c_str()));
        if (!t) {
            throw std::runtime_error("输入 ROOT 文件中找不到指定的 TTree: " + requested_name);
        }
        return t;
    }
    if (auto* t = dynamic_cast<TTree*>(file.Get("stopped_muons"))) {
        return t;
    }
    TTree* single_tree = nullptr;
    int tree_count = 0;
    TIter next(file.GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (auto* obj = key->ReadObj()) {
            if (auto* t = dynamic_cast<TTree*>(obj)) {
                single_tree = t;
                ++tree_count;
            }
        }
    }
    if (tree_count == 1 && single_tree != nullptr) {
        return single_tree;
    }
    throw std::runtime_error(
        "输入 ROOT 文件中未找到默认 TTree \"stopped_muons\"（且文件中不是唯一 TTree，请用 --input-tree 指定）");
}

} // namespace root_io_detail

inline std::vector<StoppedMuonState> read_stopped_muons_root(
        const std::string& path,
        const std::string& tree_name = "",
        int default_charge = +1,
        double default_stop_time_ns = 0.0) {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("输入 ROOT 文件不存在: " + path);
    }
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        throw std::runtime_error("无法打开输入 ROOT 文件: " + path);
    }
    TTree* tree = root_io_detail::resolve_input_tree(*file, tree_name);
    const Long64_t n_entries = tree->GetEntries();
    if (n_entries <= 0) {
        throw std::invalid_argument("输入 ROOT TTree 为空（0 entries）: " + path);
    }

    TLeaf* l_x = root_io_detail::find_leaf_any(tree, {"stop_x", "x", "pos_x"});
    TLeaf* l_y = root_io_detail::find_leaf_any(tree, {"stop_y", "y", "pos_y"});
    TLeaf* l_z = root_io_detail::find_leaf_any(tree, {"stop_z", "z", "pos_z"});
    if (!l_x || !l_y || !l_z) {
        throw std::invalid_argument(
            "输入 TTree 缺少停止位置分支（需包含 stop_x/stop_y/stop_z 或 x/y/z）");
    }

    TLeaf* l_px = root_io_detail::find_leaf_any(tree, {"pol_x", "P0_x", "px0"});
    TLeaf* l_py = root_io_detail::find_leaf_any(tree, {"pol_y", "P0_y", "py0"});
    TLeaf* l_pz = root_io_detail::find_leaf_any(tree, {"pol_z", "P0_z", "pz0"});
    if (!l_px || !l_py || !l_pz) {
        throw std::invalid_argument(
            "输入 TTree 缺少初始极化分支（需包含 pol_x/pol_y/pol_z 或 P0_x/P0_y/P0_z）");
    }

    TLeaf* l_id = root_io_detail::find_leaf_any(tree, {"event_id", "muon_id", "id"});
    TLeaf* l_q = root_io_detail::find_leaf_any(tree, {"charge", "muon_charge", "q"});
    TLeaf* l_t = root_io_detail::find_leaf_any(tree, {"stop_time_ns", "t_stop", "time_ns"});

    std::vector<StoppedMuonState> muons;
    muons.reserve(static_cast<std::size_t>(n_entries));
    for (Long64_t i = 0; i < n_entries; ++i) {
        if (tree->GetEntry(i) <= 0) {
            throw std::runtime_error("读取输入 TTree 第 " + std::to_string(i) + " 个 entry 失败");
        }
        StoppedMuonState st{};
        st.event_id = l_id ? static_cast<Long64_t>(std::llround(l_id->GetValue())) : i;
        st.stop_position = {l_x->GetValue(), l_y->GetValue(), l_z->GetValue()};
        st.initial_polarization = {l_px->GetValue(), l_py->GetValue(), l_pz->GetValue()};
        st.charge = l_q ? static_cast<int>(std::llround(l_q->GetValue())) : default_charge;
        st.stop_time_ns = l_t ? l_t->GetValue() : default_stop_time_ns;
        validate_stopped_muon_state(st);
        muons.push_back(st);
    }
    return muons;
}

inline void write_stopped_muons_root(
        const std::string& path,
        const std::vector<StoppedMuonState>& muons,
        bool overwrite = false,
        const std::string& tree_name = "stopped_muons") {
    if (!overwrite && std::filesystem::exists(path)) {
        throw std::runtime_error("文件已存在（如需覆盖请启用 overwrite）: " + path);
    }
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), overwrite ? "RECREATE" : "NEW"));
    if (!file || file->IsZombie() || !file->IsOpen()) {
        throw std::runtime_error("无法创建输入缪子 ROOT 文件: " + path);
    }
    TTree tree(tree_name.c_str(), "Stopped muons (1 entry = 1 stopped muon)");

    Long64_t event_id = 0;
    Double_t stop_x = 0, stop_y = 0, stop_z = 0;
    Double_t pol_x = 0, pol_y = 0, pol_z = 0;
    Int_t charge = +1;
    Double_t stop_time_ns = 0;

    tree.Branch("event_id", &event_id, "event_id/L");
    tree.Branch("stop_x", &stop_x, "stop_x/D");
    tree.Branch("stop_y", &stop_y, "stop_y/D");
    tree.Branch("stop_z", &stop_z, "stop_z/D");
    tree.Branch("pol_x", &pol_x, "pol_x/D");
    tree.Branch("pol_y", &pol_y, "pol_y/D");
    tree.Branch("pol_z", &pol_z, "pol_z/D");
    tree.Branch("charge", &charge, "charge/I");
    tree.Branch("stop_time_ns", &stop_time_ns, "stop_time_ns/D");

    for (const auto& m : muons) {
        validate_stopped_muon_state(m);
        event_id = m.event_id;
        stop_x = m.stop_position.X();
        stop_y = m.stop_position.Y();
        stop_z = m.stop_position.Z();
        pol_x = m.initial_polarization.X();
        pol_y = m.initial_polarization.Y();
        pol_z = m.initial_polarization.Z();
        charge = m.charge;
        stop_time_ns = m.stop_time_ns;
        tree.Fill();
    }
    file->cd();
    tree.Write();
    file->Close();
}

// Writes emitted Michel positrons / DIO electrons to a ROOT file:
// 1 TTree entry strictly corresponds to 1 emitted Michel/DIO electron (nuclear_capture events are not added as entries).
inline Long64_t write_michel_electrons_root(
        const std::string& path,
        const std::vector<DecaySimulationResult>& results,
        bool overwrite = false,
        const std::string& tree_name = "michel_electrons",
        UInt_t seed = 0,
        const std::string& preset_name = "custom",
        const Vec3& external_field_T = {0, 0, 0}) {
    if (!overwrite && std::filesystem::exists(path)) {
        throw std::runtime_error("输出 ROOT 文件已存在（如需覆盖请显式添加 --overwrite）: " + path);
    }
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), overwrite ? "RECREATE" : "NEW"));
    if (!file || file->IsZombie() || !file->IsOpen()) {
        throw std::runtime_error("无法创建输出 ROOT 文件: " + path);
    }

    TTree tree(tree_name.c_str(), "Michel positrons and DIO electrons (1 entry = 1 emitted electron)");

    Long64_t electron_id = 0;
    Long64_t muon_entry = 0;
    Long64_t event_id = 0;
    Int_t pdg = 0;
    Int_t charge = 0;
    Int_t muon_charge = 0;
    Double_t x = 0, y = 0, z = 0;
    Double_t stop_time_ns = 0, elapsed_time_ns = 0, decay_time_ns = 0;
    Double_t total_energy_MeV = 0, kinetic_energy_MeV = 0, momentum_MeV_c = 0;
    Double_t px_MeV_c = 0, py_MeV_c = 0, pz_MeV_c = 0;
    Double_t dir_x = 0, dir_y = 0, dir_z = 0;
    Double_t endpoint_MeV = 0, reduced_energy_x = 0;
    Double_t P0_x = 0, P0_y = 0, P0_z = 0;
    Double_t Pdecay_x = 0, Pdecay_y = 0, Pdecay_z = 0;
    ULong64_t hop_count = 0;
    std::string spin_category;
    std::string termination;
    std::string spectrum_model;

    tree.Branch("electron_id", &electron_id, "electron_id/L");
    tree.Branch("muon_entry", &muon_entry, "muon_entry/L");
    tree.Branch("event_id", &event_id, "event_id/L");
    tree.Branch("pdg", &pdg, "pdg/I");
    tree.Branch("charge", &charge, "charge/I");
    tree.Branch("muon_charge", &muon_charge, "muon_charge/I");
    tree.Branch("x", &x, "x/D");
    tree.Branch("y", &y, "y/D");
    tree.Branch("z", &z, "z/D");
    tree.Branch("stop_time_ns", &stop_time_ns, "stop_time_ns/D");
    tree.Branch("elapsed_time_ns", &elapsed_time_ns, "elapsed_time_ns/D");
    tree.Branch("decay_time_ns", &decay_time_ns, "decay_time_ns/D");
    tree.Branch("total_energy_MeV", &total_energy_MeV, "total_energy_MeV/D");
    tree.Branch("kinetic_energy_MeV", &kinetic_energy_MeV, "kinetic_energy_MeV/D");
    tree.Branch("momentum_MeV_c", &momentum_MeV_c, "momentum_MeV_c/D");
    tree.Branch("px_MeV_c", &px_MeV_c, "px_MeV_c/D");
    tree.Branch("py_MeV_c", &py_MeV_c, "py_MeV_c/D");
    tree.Branch("pz_MeV_c", &pz_MeV_c, "pz_MeV_c/D");
    tree.Branch("dir_x", &dir_x, "dir_x/D");
    tree.Branch("dir_y", &dir_y, "dir_y/D");
    tree.Branch("dir_z", &dir_z, "dir_z/D");
    tree.Branch("endpoint_MeV", &endpoint_MeV, "endpoint_MeV/D");
    tree.Branch("reduced_energy_x", &reduced_energy_x, "reduced_energy_x/D");
    tree.Branch("P0_x", &P0_x, "P0_x/D");
    tree.Branch("P0_y", &P0_y, "P0_y/D");
    tree.Branch("P0_z", &P0_z, "P0_z/D");
    tree.Branch("Pdecay_x", &Pdecay_x, "Pdecay_x/D");
    tree.Branch("Pdecay_y", &Pdecay_y, "Pdecay_y/D");
    tree.Branch("Pdecay_z", &Pdecay_z, "Pdecay_z/D");
    tree.Branch("hop_count", &hop_count, "hop_count/l");
    tree.Branch("spin_category", &spin_category);
    tree.Branch("termination", &termination);
    tree.Branch("spectrum_model", &spectrum_model);

    Long64_t n_captures = 0;
    for (const auto& r : results) {
        if (!r.has_decay_lepton) {
            ++n_captures;
            continue;
        }
        muon_entry = r.muon_entry;
        event_id = r.event_id;
        pdg = r.lepton_pdg;
        charge = r.lepton_charge;
        muon_charge = r.charge;
        x = r.stop_position.X();
        y = r.stop_position.Y();
        z = r.stop_position.Z();
        stop_time_ns = r.stop_time_ns;
        elapsed_time_ns = r.elapsed_time_ns;
        decay_time_ns = r.decay_time_ns;
        total_energy_MeV = r.total_energy_MeV;
        kinetic_energy_MeV = r.kinetic_energy_MeV;
        momentum_MeV_c = r.momentum_mag_MeV_c;
        px_MeV_c = r.momentum_MeV_c.X();
        py_MeV_c = r.momentum_MeV_c.Y();
        pz_MeV_c = r.momentum_MeV_c.Z();
        dir_x = r.direction.X();
        dir_y = r.direction.Y();
        dir_z = r.direction.Z();
        endpoint_MeV = r.endpoint_energy_MeV;
        reduced_energy_x = r.reduced_energy_x;
        P0_x = r.initial_polarization.X();
        P0_y = r.initial_polarization.Y();
        P0_z = r.initial_polarization.Z();
        Pdecay_x = r.decay_polarization.X();
        Pdecay_y = r.decay_polarization.Y();
        Pdecay_z = r.decay_polarization.Z();
        hop_count = static_cast<ULong64_t>(r.hop_count);
        spin_category = r.spin_category ? r.spin_category : "";
        termination = r.termination ? r.termination : "";
        spectrum_model = r.spectrum_model ? r.spectrum_model : "";

        tree.Fill();
        ++electron_id;
    }

    file->cd();
    tree.Write();
    TParameter<Long64_t>("n_input_muons", static_cast<Long64_t>(results.size())).Write();
    TParameter<Long64_t>("n_michel_electrons", electron_id).Write();
    TParameter<Long64_t>("n_nuclear_captures", n_captures).Write();
    TParameter<Long64_t>("seed", static_cast<Long64_t>(seed)).Write();
    TNamed("material_preset", preset_name.c_str()).Write();
    TVector3(external_field_T.X(), external_field_T.Y(), external_field_T.Z()).Write("external_field_T");
    file->Close();
    return electron_id;
}

// --- Reference / example material presets (parameters remain caller-overridable) ---

inline MuPlusMaterial make_copper_mu_plus_example(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9}) {
    MuPlusMaterial m{};
    m.external_field = external_field_T;
    m.muonium_fraction = 0.0;
    m.ordinary_sigma_T = 0.252e-3 / gamma_muon; // delta = 0.252e-3 ns^-1 from batch_crmusr_examples.sh
    m.ordinary_hop_rate = 5.06e-3;              // ns^-1
    return m;
}

inline MuPlusMaterial make_quartz_mu_plus_example(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9}) {
    MuPlusMaterial m{};
    m.external_field = external_field_T;
    m.muonium_fraction = 0.65;
    m.ordinary_sigma_T = 0.04e-3 / gamma_muon;
    m.ordinary_hop_rate = 0.0;
    m.muonium_sigma_T = 0.04e-3 / gamma_muon;
    m.muonium_hop_rate = 0.0;
    return m;
}

inline MuMinusMaterial make_copper_mu_minus_conditional(
        int isotope = 63,
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9},
        bool map_B = false,
        double m1_conversion_rate_ns = 0.5) {
    if (isotope != 63 && isotope != 65) {
        throw std::invalid_argument("Copper mu- preset requires isotope 63 or 65");
    }
    const MuMinusIsotropicInitialMap prep = map_B
        ? MuMinusIsotropicInitialMap{3.0 / 8.0, 1.0 / 8.0, 1.0 / 20.0}
        : MuMinusIsotropicInitialMap{3.0 / 8.0, 0.0, 1.0 / 40.0};
    const double gamma_nuc = (isotope == 63) ? 0.07100118490581576 : 0.0760556158609357;
    const double hfs_rad_ns = 925.0 / (2.0 * 6.582119569e-7); // Conditional Winston scale (~702.66 rad/ns)
    const MuMinusNuclearParameters spin{external_field_T, hfs_rad_ns, -gamma_muon, gamma_nuc};
    const MuMinusTerminationRates rates{0.0004401494396105638, 0.005676058511459773}; // tau = 163.5 ns, capture = 0.92803
    return MuMinusMaterial{isotope, prep, spin, m1_conversion_rate_ns, rates};
}

// Ferromagnetic alpha-Fe (300 K) mu+ preset (Foy 1973 / Nishida 1977 / Moslang 1983).
// Net internal field B_mu is antiparallel to domain magnetization m_hat:
//   pure single-crystal Fe: |B_mu| = 0.37 T (~50.1 MHz), lambda_T ~ 0.80 us^-1
//   commercial iron / steel: |B_mu| = 0.35 T (~47.4 MHz), lambda_T ~ 7.0 us^-1
inline MuPlusMaterial make_iron_mu_plus_example(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9},
        bool commercial = false,
        DomainFieldMode domain_mode = DomainFieldMode::IsotropicDomains,
        const Vec3& saturation_dir = {0, 0, 1}) {
    MuPlusMaterial m{};
    m.external_field = external_field_T;
    m.muonium_fraction = 0.0;
    m.ordinary_sigma_T = (commercial ? 7.0e-3 : 0.80e-3) / gamma_muon;
    m.ordinary_hop_rate = 0.0;
    m.domain_field = {domain_mode, commercial ? -0.35 : -0.37, saturation_dir};
    return m;
}

// Ferromagnetic alpha-Fe (56Fe, I=0, 300 K) bound 1s mu- preset.
// Verified rates from Suzuki et al. PRC 35, 2212 (1987): tau = 206.0 ns, capture = 4.411e-3 ns^-1,
// Q = 0.975 => Lambda_DIO = 0.4438e-3 ns^-1 (DIO 9.14% / capture 90.86%).
// Verified spin parameters from Keller et al. PRB 35, 2008 (1987):
//   B_mu(300 K) = -18.8 T (Table III 322.9 K = -18.298 T extrapolated via Eq. (9) + Delta = -0.9%),
//   bound 1s relativistic g-factor shift (g - g0)/g0 = -0.85(5)% (Table I => 134.38 MHz/T),
//   f_res = 1/6 (Sec. III: 14.8(3.0)% vs I=0 cascade limit 1/6), lambda_mu(322.9 K) = 1.5(2.7) us^-1 (~0).
inline MuMinusSpin0Material make_iron_mu_minus_conditional(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9},
        DomainFieldMode domain_mode = DomainFieldMode::IsotropicDomains,
        const Vec3& saturation_dir = {0, 0, 1}) {
    MuMinusSpin0Material m{};
    m.nuclide = DioNuclide::Fe56;
    m.residual_polarization = 1.0 / 6.0;
    m.external_field = external_field_T;
    m.domain_field = {domain_mode, -18.8, saturation_dir};
    m.local_sigma_T = 0.0;
    m.hop_rate = 0.0;
    m.gyromagnetic_ratio = gamma_muon * (1.0 - 0.0085); // Keller 1987 Table I bound 1s g-factor shift
    m.rates = {0.4438e-3, 4.411e-3};
    return m;
}

// Solid elemental sulfur (alpha-S8, 300 K) mu+ preset.
// Verified primary literature: Cox et al. J. Phys.: Condens. Matter 23, 315801 (2011) &
// Eisenstein et al. Phys. Rev. 142, 217 (1966):
//   no prompt diamagnetic fraction at 300 K (f_D = 0 => muonium_fraction = 1.0);
//   two-state chemical exchange with electron spin memory between:
//     State 1 (closed-ring bond-centre S8Mu_BC, Eqs. (6)-(7)):
//       A1/(2*pi) = -184 MHz, D1_Cox/(2*pi) = +97 MHz => D1_code = (3/2)*D1_Cox = 2*pi*(+145.5 MHz)
//       (A_parallel = -87.0 MHz, A_perp = -232.5 MHz);
//     State 2 (open-ring / weakly coupled S8Mu, Sec. 6): A2 = D2 = 0;
//   300 K bidirectional exchange rates (Fig. 9(b)): W12 = 28.0 ns^-1, W21 = 15.0 ns^-1 (p1_eq = 0.349).
inline MuPlusMaterial make_sulfur_mu_plus_conditional(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9}) {
    MuPlusMaterial m{};
    m.external_field = external_field_T;
    m.muonium_fraction = 1.0;  // f_D = 0 (Eisenstein 1966 Table I; Cox 2011 Appendix A)
    m.ordinary_sigma_T = 0.0;  // 99.2% I=0 nuclei
    m.ordinary_hop_rate = 0.0;
    m.muonium_sigma_T = 0.0;
    m.muonium_hop_rate = 0.0;
    m.hyperfine = -2.0 * pi * 0.184;         // State 1 A1 = -184 MHz in rad/ns (Cox 2011 Eq. 6)
    m.muonium_dipolar_D = 2.0 * pi * 0.1455; // State 1 D1_code = (3/2)*97 = +145.5 MHz in rad/ns (Cox 2011 Eq. 7)
    m.muonium_state2_hyperfine = 0.0;        // State 2 A2 = 0 (Cox 2011 Sec. 6)
    m.muonium_state2_dipolar_D = 0.0;        // State 2 D2 = 0 (Cox 2011 Sec. 6)
    m.muonium_exchange_w12 = 28.0;           // 300 K exit rate W12 [ns^-1] (Cox 2011 Fig. 9(b))
    m.muonium_exchange_w21 = 15.0;           // 300 K return rate W21 [ns^-1] (Cox 2011 Fig. 9(b))
    m.muonium_spin_exchange_rate = 0.0;      // Electron spin retains coherence across 1 <-> 2 exchange
    m.muonium_reorientation_rate = 0.0;
    m.muonium_initial_p1 = std::nullopt;     // Equilibrium p1(0) = W21 / (W12 + W21)
    return m;
}

// Solid elemental sulfur (32S, I=0, 300 K) bound 1s mu- preset.
// Verified rates from Suzuki et al. PRC 35, 2212 (1987): tau = 554.7 ns, capture = 1.352e-3 ns^-1,
// Q = 0.990 => Lambda_DIO = 0.4506e-3 ns^-1 (DIO 25.00% / capture 75.00%).
// Verified 300 K spin parameters from Mamedov et al. JETP 93(5), 941 (2001):
//   a(mu-, S) / a(mu-, C) = 0.0213(2) / 0.0486(3) = 0.4383 => f_res(S) ~ 0.4383/6 = 0.073 (Sec. 4),
//   bound 1s relativistic g-factor shift (g_free - g_1s)/g_free = 42.4(2.1)e-4 (Tables 1 & 2).
inline MuMinusSpin0Material make_sulfur_mu_minus_conditional(
        const Vec3& external_field_T = {-1762.0e-9, 37638.2e-9, -26120.2e-9}) {
    MuMinusSpin0Material m{};
    m.nuclide = DioNuclide::S32;
    m.residual_polarization = 0.073;
    m.external_field = external_field_T;
    m.domain_field = {DomainFieldMode::None, 0.0, {0, 0, 1}};
    m.local_sigma_T = 0.0;
    m.hop_rate = 0.0;
    m.gyromagnetic_ratio = gamma_muon * (1.0 - 4.24e-3); // Mamedov 2001 Tables 1 & 2 bound 1s g-factor shift
    m.rates = {0.4506e-3, 1.352e-3};
    return m;
}

} // namespace stopped_muon
