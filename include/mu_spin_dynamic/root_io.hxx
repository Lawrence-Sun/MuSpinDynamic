#pragma once

#include "spin_dynamic_decay.hxx"

#include "TBranch.h"
#include "TFile.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace stopped_muon {

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
        st.event_id = l_id ? static_cast<std::int64_t>(std::llround(l_id->GetValue()))
                           : static_cast<std::int64_t>(i);
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
        event_id = static_cast<Long64_t>(m.event_id);
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
        muon_entry = static_cast<Long64_t>(r.muon_entry);
        event_id = static_cast<Long64_t>(r.event_id);
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

} // namespace stopped_muon
