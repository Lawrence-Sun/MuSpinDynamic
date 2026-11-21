#pragma once

#include "rng.hxx"
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
#include <unordered_set>
#include <vector>

namespace stopped_muon {

// --- ROOT TTree I/O for Mode 3 / Mode B batch processing ---

struct BatchOptions {
    std::string input_tree;
    std::string output_tree = "michel_electrons";
    bool overwrite = false;
    int default_charge = +1;
    double default_stop_time_ns = 0.0;
    double default_weight = 1.0;
    std::int32_t replicas = 1;
    bool keep_captures = true;
    bool event_seeded = true;
    bool check_unique_event_ids = true;
    std::uint64_t seed = 42;
    std::string preset_name = "custom";
    Vec3 external_field_T{0, 0, 0};
};

struct BatchSummary {
    Long64_t n_input_muons = 0;
    Long64_t n_simulated_events = 0;
    Long64_t n_written_entries = 0;
    Long64_t n_emitted_leptons = 0;
    Long64_t n_nuclear_captures = 0;
    double sum_input_weights = 0.0;
    double sum_output_weights = 0.0;
};

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

class StoppedMuonBranchReader {
public:
    StoppedMuonBranchReader(
            TTree* tree,
            int default_charge,
            double default_stop_time_ns,
            double default_weight)
        : tree_(tree),
          default_charge_(default_charge),
          default_stop_time_ns_(default_stop_time_ns),
          default_weight_(default_weight) {
        l_x_ = find_leaf_any(tree_, {"stop_x", "x", "pos_x"});
        l_y_ = find_leaf_any(tree_, {"stop_y", "y", "pos_y"});
        l_z_ = find_leaf_any(tree_, {"stop_z", "z", "pos_z"});
        if (!l_x_ || !l_y_ || !l_z_) {
            throw std::invalid_argument(
                "输入 TTree 缺少停止位置分支（需包含 stop_x/stop_y/stop_z 或 x/y/z）");
        }
        l_px_ = find_leaf_any(tree_, {"pol_x", "P0_x", "px0"});
        l_py_ = find_leaf_any(tree_, {"pol_y", "P0_y", "py0"});
        l_pz_ = find_leaf_any(tree_, {"pol_z", "P0_z", "pz0"});
        if (!l_px_ || !l_py_ || !l_pz_) {
            throw std::invalid_argument(
                "输入 TTree 缺少初始极化分支（需包含 pol_x/pol_y/pol_z 或 P0_x/P0_y/P0_z）");
        }
        l_id_ = find_leaf_any(tree_, {"event_id", "EventID", "muon_id", "id"});
        l_q_ = find_leaf_any(tree_, {"charge", "muon_charge", "q"});
        l_t_ = find_leaf_any(tree_, {"stop_time_ns", "t_stop", "time_ns"});
        l_w_ = find_leaf_any(tree_, {"weight", "Weight", "event_weight"});
    }

    StoppedMuonState read_entry(Long64_t i) const {
        if (tree_->GetEntry(i) <= 0) {
            throw std::runtime_error("读取输入 TTree 第 " + std::to_string(i) + " 个 entry 失败");
        }
        StoppedMuonState st{};
        st.event_id = l_id_ ? static_cast<std::int64_t>(std::llround(l_id_->GetValue()))
                            : static_cast<std::int64_t>(i);
        st.stop_position = {l_x_->GetValue(), l_y_->GetValue(), l_z_->GetValue()};
        st.initial_polarization = {l_px_->GetValue(), l_py_->GetValue(), l_pz_->GetValue()};
        st.charge = l_q_ ? static_cast<int>(std::llround(l_q_->GetValue())) : default_charge_;
        st.stop_time_ns = l_t_ ? l_t_->GetValue() : default_stop_time_ns_;
        st.weight = l_w_ ? l_w_->GetValue() : default_weight_;
        validate_stopped_muon_state(st);
        return st;
    }

private:
    TTree* tree_ = nullptr;
    int default_charge_ = +1;
    double default_stop_time_ns_ = 0.0;
    double default_weight_ = 1.0;
    TLeaf* l_x_ = nullptr;
    TLeaf* l_y_ = nullptr;
    TLeaf* l_z_ = nullptr;
    TLeaf* l_px_ = nullptr;
    TLeaf* l_py_ = nullptr;
    TLeaf* l_pz_ = nullptr;
    TLeaf* l_id_ = nullptr;
    TLeaf* l_q_ = nullptr;
    TLeaf* l_t_ = nullptr;
    TLeaf* l_w_ = nullptr;
};

// Disk-backed streaming TTree writer shared by vector and streaming ROOT output paths.
class DecayTreeWriter {
public:
    DecayTreeWriter(
            const std::string& path,
            bool overwrite,
            const std::string& tree_name,
            bool keep_captures)
        : keep_captures_(keep_captures) {
        if (!overwrite && std::filesystem::exists(path)) {
            throw std::runtime_error("输出 ROOT 文件已存在（如需覆盖请显式添加 --overwrite）: " + path);
        }
        file_.reset(TFile::Open(path.c_str(), overwrite ? "RECREATE" : "NEW"));
        if (!file_ || file_->IsZombie() || !file_->IsOpen()) {
            throw std::runtime_error("无法创建输出 ROOT 文件: " + path);
        }
        file_->cd();
        const char* title = keep_captures_
            ? "Stopped-muon spin-dynamics decay and nuclear-capture outcomes (check has_decay_lepton)"
            : "Michel positrons and DIO electrons (1 entry = 1 emitted electron)";
        tree_ = std::make_unique<TTree>(tree_name.c_str(), title);
        tree_->SetDirectory(file_.get());

        tree_->Branch("electron_id", &electron_id_, "electron_id/L");
        tree_->Branch("muon_entry", &muon_entry_, "muon_entry/L");
        tree_->Branch("event_id", &event_id_, "event_id/L");
        tree_->Branch("replica_id", &replica_id_, "replica_id/I");
        tree_->Branch("weight", &weight_, "weight/D");
        tree_->Branch("has_decay_lepton", &has_decay_lepton_, "has_decay_lepton/O");
        tree_->Branch("pdg", &pdg_, "pdg/I");
        tree_->Branch("charge", &charge_, "charge/I");
        tree_->Branch("muon_charge", &muon_charge_, "muon_charge/I");
        tree_->Branch("x", &x_, "x/D");
        tree_->Branch("y", &y_, "y/D");
        tree_->Branch("z", &z_, "z/D");
        tree_->Branch("stop_time_ns", &stop_time_ns_, "stop_time_ns/D");
        tree_->Branch("elapsed_time_ns", &elapsed_time_ns_, "elapsed_time_ns/D");
        tree_->Branch("decay_time_ns", &decay_time_ns_, "decay_time_ns/D");
        tree_->Branch("total_energy_MeV", &total_energy_MeV_, "total_energy_MeV/D");
        tree_->Branch("kinetic_energy_MeV", &kinetic_energy_MeV_, "kinetic_energy_MeV/D");
        tree_->Branch("momentum_MeV_c", &momentum_MeV_c_, "momentum_MeV_c/D");
        tree_->Branch("px_MeV_c", &px_MeV_c_, "px_MeV_c/D");
        tree_->Branch("py_MeV_c", &py_MeV_c_, "py_MeV_c/D");
        tree_->Branch("pz_MeV_c", &pz_MeV_c_, "pz_MeV_c/D");
        tree_->Branch("dir_x", &dir_x_, "dir_x/D");
        tree_->Branch("dir_y", &dir_y_, "dir_y/D");
        tree_->Branch("dir_z", &dir_z_, "dir_z/D");
        tree_->Branch("endpoint_MeV", &endpoint_MeV_, "endpoint_MeV/D");
        tree_->Branch("reduced_energy_x", &reduced_energy_x_, "reduced_energy_x/D");
        tree_->Branch("P0_x", &P0_x_, "P0_x/D");
        tree_->Branch("P0_y", &P0_y_, "P0_y/D");
        tree_->Branch("P0_z", &P0_z_, "P0_z/D");
        tree_->Branch("Pdecay_x", &Pdecay_x_, "Pdecay_x/D");
        tree_->Branch("Pdecay_y", &Pdecay_y_, "Pdecay_y/D");
        tree_->Branch("Pdecay_z", &Pdecay_z_, "Pdecay_z/D");
        tree_->Branch("hop_count", &hop_count_, "hop_count/l");
        tree_->Branch("spin_category", &spin_category_);
        tree_->Branch("termination", &termination_);
        tree_->Branch("spectrum_model", &spectrum_model_);
    }

    DecayTreeWriter(const DecayTreeWriter&) = delete;
    DecayTreeWriter& operator=(const DecayTreeWriter&) = delete;

    void push(const DecaySimulationResult& r) {
        if (!r.has_decay_lepton) {
            ++n_captures_;
            if (!keep_captures_) {
                return;
            }
            electron_id_ = -1;
            has_decay_lepton_ = kFALSE;
        } else {
            electron_id_ = n_emitted_++;
            has_decay_lepton_ = kTRUE;
        }

        muon_entry_ = static_cast<Long64_t>(r.muon_entry);
        event_id_ = static_cast<Long64_t>(r.event_id);
        replica_id_ = static_cast<Int_t>(r.replica_id);
        weight_ = r.weight;
        pdg_ = r.lepton_pdg;
        charge_ = r.lepton_charge;
        muon_charge_ = r.charge;
        x_ = r.stop_position.X();
        y_ = r.stop_position.Y();
        z_ = r.stop_position.Z();
        stop_time_ns_ = r.stop_time_ns;
        elapsed_time_ns_ = r.elapsed_time_ns;
        decay_time_ns_ = r.decay_time_ns;
        total_energy_MeV_ = r.total_energy_MeV;
        kinetic_energy_MeV_ = r.kinetic_energy_MeV;
        momentum_MeV_c_ = r.momentum_mag_MeV_c;
        px_MeV_c_ = r.momentum_MeV_c.X();
        py_MeV_c_ = r.momentum_MeV_c.Y();
        pz_MeV_c_ = r.momentum_MeV_c.Z();
        dir_x_ = r.direction.X();
        dir_y_ = r.direction.Y();
        dir_z_ = r.direction.Z();
        endpoint_MeV_ = r.endpoint_energy_MeV;
        reduced_energy_x_ = r.reduced_energy_x;
        P0_x_ = r.initial_polarization.X();
        P0_y_ = r.initial_polarization.Y();
        P0_z_ = r.initial_polarization.Z();
        Pdecay_x_ = r.decay_polarization.X();
        Pdecay_y_ = r.decay_polarization.Y();
        Pdecay_z_ = r.decay_polarization.Z();
        hop_count_ = static_cast<ULong64_t>(r.hop_count);
        spin_category_ = r.spin_category ? r.spin_category : "";
        termination_ = r.termination ? r.termination : "";
        spectrum_model_ = r.spectrum_model ? r.spectrum_model : "";

        tree_->Fill();
        ++n_written_;
        sum_output_weights_ += r.weight;
    }

    BatchSummary finalize(
            Long64_t n_input_muons,
            Long64_t n_simulated_events,
            double sum_input_weights,
            Int_t n_replicas,
            bool event_seeded,
            std::uint64_t seed,
            const std::string& preset_name,
            const Vec3& external_field_T) {
        file_->cd();
        tree_->Write();
        TParameter<Long64_t>("n_input_muons", n_input_muons).Write();
        TParameter<Long64_t>("n_simulated_events", n_simulated_events).Write();
        TParameter<Long64_t>("n_written_entries", n_written_).Write();
        TParameter<Long64_t>("n_michel_electrons", n_emitted_).Write();
        TParameter<Long64_t>("n_nuclear_captures", n_captures_).Write();
        TParameter<Int_t>("n_replicas", n_replicas).Write();
        TParameter<Int_t>("keep_captures", keep_captures_ ? 1 : 0).Write();
        TParameter<Int_t>("event_seeded", event_seeded ? 1 : 0).Write();
        TParameter<Long64_t>("seed", static_cast<Long64_t>(seed)).Write();
        TNamed("material_preset", preset_name.c_str()).Write();
        TVector3(external_field_T.X(), external_field_T.Y(), external_field_T.Z()).Write("external_field_T");

        // Detach tree before closing file so destructor does not double-free it.
        tree_->SetDirectory(nullptr);
        file_->Close();

        BatchSummary summary{};
        summary.n_input_muons = n_input_muons;
        summary.n_simulated_events = n_simulated_events;
        summary.n_written_entries = n_written_;
        summary.n_emitted_leptons = n_emitted_;
        summary.n_nuclear_captures = n_captures_;
        summary.sum_input_weights = sum_input_weights;
        summary.sum_output_weights = sum_output_weights_;
        return summary;
    }

private:
    bool keep_captures_ = true;
    std::unique_ptr<TFile> file_;
    std::unique_ptr<TTree> tree_;

    Long64_t n_emitted_ = 0;
    Long64_t n_captures_ = 0;
    Long64_t n_written_ = 0;
    double sum_output_weights_ = 0.0;

    Long64_t electron_id_ = -1;
    Long64_t muon_entry_ = 0;
    Long64_t event_id_ = 0;
    Int_t replica_id_ = 0;
    Double_t weight_ = 1.0;
    Bool_t has_decay_lepton_ = kFALSE;
    Int_t pdg_ = 0;
    Int_t charge_ = 0;
    Int_t muon_charge_ = 0;
    Double_t x_ = 0, y_ = 0, z_ = 0;
    Double_t stop_time_ns_ = 0, elapsed_time_ns_ = 0, decay_time_ns_ = 0;
    Double_t total_energy_MeV_ = 0, kinetic_energy_MeV_ = 0, momentum_MeV_c_ = 0;
    Double_t px_MeV_c_ = 0, py_MeV_c_ = 0, pz_MeV_c_ = 0;
    Double_t dir_x_ = 0, dir_y_ = 0, dir_z_ = 0;
    Double_t endpoint_MeV_ = 0, reduced_energy_x_ = 0;
    Double_t P0_x_ = 0, P0_y_ = 0, P0_z_ = 0;
    Double_t Pdecay_x_ = 0, Pdecay_y_ = 0, Pdecay_z_ = 0;
    ULong64_t hop_count_ = 0;
    std::string spin_category_;
    std::string termination_;
    std::string spectrum_model_;
};

} // namespace root_io_detail

inline std::vector<StoppedMuonState> read_stopped_muons_root(
        const std::string& path,
        const std::string& tree_name = "",
        int default_charge = +1,
        double default_stop_time_ns = 0.0,
        double default_weight = 1.0) {
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

    root_io_detail::StoppedMuonBranchReader reader(
        tree, default_charge, default_stop_time_ns, default_weight);
    std::vector<StoppedMuonState> muons;
    muons.reserve(static_cast<std::size_t>(n_entries));
    for (Long64_t i = 0; i < n_entries; ++i) {
        muons.push_back(reader.read_entry(i));
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
    file->cd();
    TTree tree(tree_name.c_str(), "Stopped muons (1 entry = 1 stopped muon)");

    Long64_t event_id = 0;
    Double_t stop_x = 0, stop_y = 0, stop_z = 0;
    Double_t pol_x = 0, pol_y = 0, pol_z = 0;
    Int_t charge = +1;
    Double_t stop_time_ns = 0;
    Double_t weight = 1.0;

    tree.Branch("event_id", &event_id, "event_id/L");
    tree.Branch("stop_x", &stop_x, "stop_x/D");
    tree.Branch("stop_y", &stop_y, "stop_y/D");
    tree.Branch("stop_z", &stop_z, "stop_z/D");
    tree.Branch("pol_x", &pol_x, "pol_x/D");
    tree.Branch("pol_y", &pol_y, "pol_y/D");
    tree.Branch("pol_z", &pol_z, "pol_z/D");
    tree.Branch("charge", &charge, "charge/I");
    tree.Branch("stop_time_ns", &stop_time_ns, "stop_time_ns/D");
    tree.Branch("weight", &weight, "weight/D");

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
        weight = m.weight;
        tree.Fill();
    }
    file->cd();
    tree.Write();
    file->Close();
}

// Writes decay simulation results to a ROOT file.
// When keep_captures == true (default), nuclear_capture outcomes are preserved as placeholder entries
// with has_decay_lepton == false and electron_id == -1 so Stage 1 / Stage 2 mergers keep 1:K alignment.
// Returns the total number of entries written to the output TTree.
inline Long64_t write_michel_electrons_root(
        const std::string& path,
        const std::vector<DecaySimulationResult>& results,
        bool overwrite = false,
        const std::string& tree_name = "michel_electrons",
        std::uint64_t seed = 0,
        const std::string& preset_name = "custom",
        const Vec3& external_field_T = {0, 0, 0},
        bool keep_captures = true,
        Int_t n_replicas = 1,
        bool event_seeded = false) {
    root_io_detail::DecayTreeWriter writer(path, overwrite, tree_name, keep_captures);
    double sum_weights = 0.0;
    for (const auto& r : results) {
        sum_weights += r.weight;
        writer.push(r);
    }
    const Long64_t n_sim = static_cast<Long64_t>(results.size());
    const Long64_t n_in = (n_replicas > 0 && n_sim % n_replicas == 0) ? (n_sim / n_replicas) : n_sim;
    const auto summary = writer.finalize(
        n_in, n_sim, sum_weights, n_replicas, event_seeded, seed, preset_name, external_field_T);
    return summary.n_written_entries;
}

// O(1)-memory streaming batch processor: reads stopped muons entry-by-entry from input_path,
// runs `opt.replicas` spin-dynamics simulations per stopped muon, and fills output_path directly.
template <typename SequentialRng>
inline BatchSummary process_stopped_muons_root_streaming(
        const std::string& input_path,
        const std::string& output_path,
        const MuPlusMaterial& mu_plus_material,
        const std::optional<AnyMuMinusMaterial>& mu_minus_material,
        const BatchOptions& opt,
        SequentialRng& sequential_rng) {
    if (opt.replicas <= 0) {
        throw std::invalid_argument("replicas 必须为正整数");
    }
    if (!std::filesystem::exists(input_path)) {
        throw std::runtime_error("输入 ROOT 文件不存在: " + input_path);
    }
    if (std::filesystem::exists(output_path) &&
        std::filesystem::equivalent(input_path, output_path)) {
        throw std::invalid_argument("输入与输出 ROOT 文件路径不能指向同一文件: " + output_path);
    }

    std::unique_ptr<TFile> in_file(TFile::Open(input_path.c_str(), "READ"));
    if (!in_file || in_file->IsZombie()) {
        throw std::runtime_error("无法打开输入 ROOT 文件: " + input_path);
    }
    TTree* in_tree = root_io_detail::resolve_input_tree(*in_file, opt.input_tree);
    const Long64_t n_entries = in_tree->GetEntries();
    if (n_entries <= 0) {
        throw std::invalid_argument("输入 ROOT TTree 为空（0 entries）: " + input_path);
    }

    root_io_detail::StoppedMuonBranchReader reader(
        in_tree, opt.default_charge, opt.default_stop_time_ns, opt.default_weight);
    root_io_detail::DecayTreeWriter writer(
        output_path, opt.overwrite, opt.output_tree, opt.keep_captures);

    std::unordered_set<std::int64_t> seen_event_ids;
    if (opt.event_seeded && opt.check_unique_event_ids) {
        seen_event_ids.reserve(static_cast<std::size_t>(n_entries));
    }

    const double weight_scale = 1.0 / static_cast<double>(opt.replicas);
    double sum_input_weights = 0.0;
    for (Long64_t i = 0; i < n_entries; ++i) {
        const StoppedMuonState st = reader.read_entry(i);
        if (st.charge == -1 && !mu_minus_material.has_value()) {
            throw std::invalid_argument(
                "输入 ROOT 文件包含 mu- 事例，但当前材料未配置 mu- 参数（请使用 copper|iron|sulfur 或显式提供参数）");
        }
        if (opt.event_seeded && opt.check_unique_event_ids) {
            if (!seen_event_ids.insert(st.event_id).second) {
                throw std::invalid_argument(
                    "检测到重复的 event_id = " + std::to_string(st.event_id) +
                    "（按事例派生随机流模式要求输入 event_id 唯一；请修复上游 event_id 或显式改用 --sequential-rng）");
            }
        }
        sum_input_weights += st.weight;
        for (std::int32_t rep = 0; rep < opt.replicas; ++rep) {
            if (opt.event_seeded) {
                auto ev_rng = make_event_replica_rng(opt.seed, st.event_id, rep);
                writer.push(simulate_stopped_muon_decay(
                    st, mu_plus_material, mu_minus_material, ev_rng, i, rep, weight_scale));
            } else {
                writer.push(simulate_stopped_muon_decay(
                    st, mu_plus_material, mu_minus_material, sequential_rng, i, rep, weight_scale));
            }
        }
    }

    const Long64_t n_simulated = n_entries * static_cast<Long64_t>(opt.replicas);
    return writer.finalize(
        n_entries, n_simulated, sum_input_weights,
        static_cast<Int_t>(opt.replicas), opt.event_seeded,
        opt.seed, opt.preset_name, opt.external_field_T);
}

inline BatchSummary process_stopped_muons_root_streaming(
        const std::string& input_path,
        const std::string& output_path,
        const MuPlusMaterial& mu_plus_material,
        const std::optional<AnyMuMinusMaterial>& mu_minus_material,
        const BatchOptions& opt) {
    auto fallback_rng = make_event_replica_rng(opt.seed, 0, 0);
    return process_stopped_muons_root_streaming(
        input_path, output_path, mu_plus_material, mu_minus_material, opt, fallback_rng);
}

} // namespace stopped_muon
