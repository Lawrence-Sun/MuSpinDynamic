#include "mu_spin_dynamic/material_config.hxx"
#include "mu_spin_dynamic/root_io.hxx"

#include "TRandom3.h"

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using namespace stopped_muon;

struct CliOptions {
    // Mode selection: empty input_file => Mode A (default direct input); non-empty => Mode B (.root -> .root streaming batch).
    std::string input_file;
    std::string input_tree;
    std::string output_file;
    std::string output_tree = "michel_electrons";
    bool overwrite = false;

    // Mode A muon stopping state
    Vec3 stop_position{0, 0, 0};
    Vec3 initial_polarization{0, 0, 1};
    bool pol_explicit = false;
    int charge = +1;
    double stop_time_ns = 0.0;
    double weight = 1.0;
    std::size_t events = 1;
    std::int32_t replicas = 1;
    bool keep_captures = true;
    bool event_seeded = true;
    UInt_t seed = 42;

    // Shared material baseline & overrides (same config used by Mode 1 Geant4 integration)
    MaterialConfig material;
};

bool has_root_extension(const std::string& path) {
    return path.size() > 5 && path.substr(path.size() - 5) == ".root";
}

double parse_double(const std::string& s, const char* flag) {
    std::size_t idx = 0;
    double val = 0;
    try {
        val = std::stod(s, &idx);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("无法解析数值参数 ") + flag + ": " + s);
    }
    while (idx < s.size() && std::isspace(static_cast<unsigned char>(s[idx]))) ++idx;
    if (idx != s.size() || !std::isfinite(val)) {
        throw std::invalid_argument(std::string("非法或非有限数值参数 ") + flag + ": " + s);
    }
    return val;
}

long long parse_ll(const std::string& s, const char* flag) {
    std::size_t idx = 0;
    long long val = 0;
    try {
        val = std::stoll(s, &idx, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("无法解析整数参数 ") + flag + ": " + s);
    }
    while (idx < s.size() && std::isspace(static_cast<unsigned char>(s[idx]))) ++idx;
    if (idx != s.size()) {
        throw std::invalid_argument(std::string("非法整数参数 ") + flag + ": " + s);
    }
    return val;
}

void print_usage(std::ostream& os, const char* prog) {
    os << "用法:\n"
       << "  方案 A（默认，命令行直接输入单个停止缪子与材料参数）:\n"
       << "    " << prog << " [--pos x y z] [--pol px py pz] [--charge +1|-1] [--material copper|iron|iron_pure|iron_commercial|sulfur|quartz|glass|free] [材料细调参数...]\n\n"
       << "  方案 B（ROOT 文件流式批处理，支持事件分裂与核俘获占位保留）:\n"
       << "    " << prog << " --input stopped_muons.root --output michel_electrons.root [--replicas K] [--material copper|iron|sulfur] [材料细调参数...]\n\n"
       << "主要参数说明:\n"
       << "  --pos <x> <y> <z>            缪子停止空间位置（如 mm），原样作为衰变电子出生顶点（默认: 0 0 0）\n"
       << "  --pol <px> <py> <pz>         停止时刻实验室系初始极化 P0，|P0|<=1（默认: mu+ 为 0 0 1，mu- 为 0 0 -1）\n"
       << "  --charge <+1|-1>             缪子电荷，+1 为正缪子，-1 为负缪子（默认: +1）\n"
       << "  --stop-time <t_ns>           缪子停止时刻 [ns]（默认: 0）\n"
       << "  --weight <w>                 默认事例权重（默认: 1.0；输入 TTree 含 weight 分支时优先读入）\n"
       << "  -n, --events <N>             方案 A 下连续独立模拟的事例数（默认: 1）\n"
       << "  --replicas <K>               每个停止缪子复用的分裂倍数 K，输出权重自动缩放为 weight/K（默认: 1）\n"
       << "  --seed <uint>                基准随机数种子（默认: 42）\n"
       << "  --event-seeded               按 (seed, event_id, replica_id) 派生独立随机流（默认开启，并检查 event_id 唯一性）\n"
       << "  --sequential-rng             回退到单流顺序 TRandom3(seed) 随机数模式\n\n"
       << "材料预设与磁场参数:\n"
       << "  --material <name>            基准材料预设: copper（默认）| iron / iron_pure | iron_commercial | sulfur | quartz | glass | free\n"
       << "  --B <Bx> <By> <Bz>           外加静态磁场矢量 [T]（覆盖预设外场）\n"
       << "  --domain-mode <mode>         铁磁磁畴模型: none | isotropic（未磁化多畴）| saturated（单畴饱和）\n"
       << "  --domain-dir <mx> <my> <mz>  饱和磁化方向矢量（仅 --domain-mode saturated 时生效）\n"
       << "  --internal-field <T>         覆盖 mu+ 磁畴有符号静态内场 [T]（如纯铁默认 -0.37）\n\n"
       << "正缪子（mu+ / Muonium / 各向异性自由基）材料细调参数:\n"
       << "  --delta <ns^-1>              普通 mu+ 核偶极局域场频率宽度 Delta [ns^-1]（等价设置 sigma_B = Delta/gamma_mu）\n"
       << "  --sigma-B <T>                普通 mu+ 每笛卡尔分量局域场标准差 [T]\n"
       << "  --hop-rate <ns^-1>           普通 mu+ 格点跳跃率 [ns^-1]\n"
       << "  --muonium-frac <0..1>        缪子素 / 顺磁自由基形成比例（1 - f_D）\n"
       << "  --muonium-sigma-B <T>        缪子素 / 自由基每分量局域场标准差 [T]\n"
       << "  --muonium-hop-rate <ns^-1>   缪子素 / 自由基局域场跳跃率 [ns^-1]\n"
       << "  --hyperfine <rad/ns>         有符号各向同性超精细耦合常数 A_iso（态 1）[rad/ns]\n"
       << "  --dipolar-D <rad/ns>         有符号轴对称偶极超精细常数 D = A_parallel - A_perp（态 1）[rad/ns]\n"
       << "  --state2-hyperfine <rad/ns>  两态化学交换态 2 各向同性超精细常数 A2 [rad/ns]（默认: 0）\n"
       << "  --state2-dipolar-D <rad/ns>  两态化学交换态 2 偶极超精细常数 D2 [rad/ns]（默认: 0）\n"
       << "  --exchange-w12 <ns^-1>       两态化学交换跃出率 W12（态 1 -> 态 2）[ns^-1]\n"
       << "  --exchange-w21 <ns^-1>       两态化学交换返回率 W21（态 2 -> 态 1）[ns^-1]\n"
       << "  --spin-exchange-rate <ns^-1> 自由基电子自旋随机化率 [ns^-1]\n"
       << "  --reorient-rate <ns^-1>      自由基分子对称轴重取向率 [ns^-1]\n"
       << "  --max-hops <N>               单事例最大允许跳跃步数预算（默认: 100000）\n\n"
       << "负缪子（mu-，支持自然铜/Cu63/Cu65 [I=3/2] 及 Fe56/S32 [I=0]）材料细调参数:\n"
       << "  --isotope <natural|63|65>    铜同位素模式（默认: natural，按 69.15% Cu63 / 30.85% Cu65 抽样）\n"
       << "  --cu-map <A|B>               铜靶级联后有效初态预设映射 A 或 B（默认: A）\n"
       << "  --f1-pop <p1>                铜靶自定义 F=1 布居数\n"
       << "  --f1-orient <c1>             铜靶自定义 F=1 取向系数\n"
       << "  --f2-orient <c2>             铜靶自定义 F=2 取向系数\n"
       << "  --cu-hyperfine <rad/ns>      铜缪原子 1s 超精细常数 [rad/ns]\n"
       << "  --m1-rate <ns^-1>            铜靶 F+ -> F- 单向 M1 跃迁率 [ns^-1]（默认: 0.5）\n"
       << "  --residual-pol <0..1>        I=0 核素（Fe56/S32）级联后剩余极化因子 f_res（铁默认 1/6，硫默认 0.073）\n"
       << "  --mu-minus-internal-field <T> I=0 铁磁靶 mu- 有符号超精细内场 [T]（铁靶 300 K 默认 -18.8）\n"
       << "  --mu-minus-sigma-B <T>       I=0 靶 mu- 局域场展宽标准差 [T]\n"
       << "  --mu-minus-hop-rate <ns^-1>  I=0 靶 mu- 局域场跳跃率 [ns^-1]\n"
       << "  --dio-rate <ns^-1>           束缚轨道衰变（DIO）率 [ns^-1]\n"
       << "  --capture-rate <ns^-1>       原子核俘获率 [ns^-1]\n\n"
       << "ROOT 文件批处理控制:\n"
       << "  --input <in.root>            启用方案 B：输入停止缪子 ROOT 文件（1 entry = 1 个停止缪子）\n"
       << "  --input-tree <name>          输入 TTree 名称（默认查找 stopped_muons，若文件仅有唯一 TTree 则自动选用）\n"
       << "  --output <out.root>          输出衰变末态 ROOT 文件（方案 B 必填，方案 A 可选）\n"
       << "  --output-tree <name>         输出 TTree 名称（默认: michel_electrons）\n"
       << "  --keep-captures              在输出 TTree 中保留 mu- 核俘获占位记录（默认开启，has_decay_lepton 标记）\n"
       << "  --electrons-only             仅写出产生 Michel/DIO 电子的事例（跳过核俘获记录）\n"
       << "  --overwrite                  允许覆盖已存在的 --output ROOT 文件\n"
       << "  -h, --help                   显示本帮助信息\n";
}

void write_pretty(std::ostream& os, Long64_t muon_entry, const DecaySimulationResult& r) {
    os << std::fixed << std::setprecision(6);
    os << "=== Muon Entry #" << muon_entry
       << " (event_id = " << r.event_id
       << ", replica_id = " << r.replica_id
       << ", weight = " << r.weight << ") ===\n"
       << "  停止位置 / 出生顶点 : (" << r.stop_position.X() << ", "
       << r.stop_position.Y() << ", " << r.stop_position.Z() << ")\n"
       << "  入射电荷 / 自旋态   : " << (r.charge > 0 ? "+1" : "-1")
       << " (" << r.spin_category << ")\n"
       << "  终止通道            : " << r.termination << '\n'
       << "  停止时刻 / 演化寿命 : t_stop = " << r.stop_time_ns
       << " ns, dt = " << r.elapsed_time_ns
       << " ns -> t_end = " << r.decay_time_ns << " ns\n"
       << "  初始极化 P0         : (" << r.initial_polarization.X() << ", "
       << r.initial_polarization.Y() << ", " << r.initial_polarization.Z() << ")\n"
       << "  终止极化 P(t_end)   : (" << r.decay_polarization.X() << ", "
       << r.decay_polarization.Y() << ", " << r.decay_polarization.Z()
       << ")  [hops = " << r.hop_count << "]\n";
    if (r.has_decay_lepton) {
        os << "  衰变出射轻子 (PDG)  : " << r.lepton_pdg
           << (r.lepton_pdg == -11 ? " (e+)" : " (e-)")
           << "  [model: " << r.spectrum_model << "]\n"
           << "  总能量 / 动能       : E_tot = " << r.total_energy_MeV
           << " MeV, E_kin = " << r.kinetic_energy_MeV
           << " MeV  (x = E/E_max = " << r.reduced_energy_x << ")\n"
           << "  三维动量 [MeV/c]    : (" << r.momentum_MeV_c.X() << ", "
           << r.momentum_MeV_c.Y() << ", " << r.momentum_MeV_c.Z()
           << "), |p| = " << r.momentum_mag_MeV_c << " MeV/c\n"
           << "  出射方向单位矢量    : (" << r.direction.X() << ", "
           << r.direction.Y() << ", " << r.direction.Z() << ")\n";
    } else {
        os << "  衰变出射轻子        : 无（发生原子核俘获 nuclear_capture）\n";
    }
}

CliOptions parse_cli(int argc, char** argv) {
    CliOptions opt;
    auto require_args = [&](int i, int count, const char* flag) {
        if (i + count >= argc) {
            throw std::invalid_argument(std::string("参数 ") + flag + " 缺少所需的后续值");
        }
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(std::cout, argv[0]);
            std::exit(0);
        } else if (arg == "--pos") {
            require_args(i, 3, "--pos");
            opt.stop_position = {
                parse_double(argv[i + 1], "--pos"),
                parse_double(argv[i + 2], "--pos"),
                parse_double(argv[i + 3], "--pos")
            };
            i += 3;
        } else if (arg == "--pol") {
            require_args(i, 3, "--pol");
            opt.initial_polarization = {
                parse_double(argv[i + 1], "--pol"),
                parse_double(argv[i + 2], "--pol"),
                parse_double(argv[i + 3], "--pol")
            };
            opt.pol_explicit = true;
            i += 3;
        } else if (arg == "--charge") {
            require_args(i, 1, "--charge");
            const auto q = parse_ll(argv[++i], "--charge");
            if (q != 1 && q != -1) {
                throw std::invalid_argument("--charge 仅支持 +1 或 -1");
            }
            opt.charge = static_cast<int>(q);
        } else if (arg == "--stop-time") {
            require_args(i, 1, "--stop-time");
            opt.stop_time_ns = parse_double(argv[++i], "--stop-time");
        } else if (arg == "--weight") {
            require_args(i, 1, "--weight");
            opt.weight = parse_double(argv[++i], "--weight");
            if (opt.weight < 0.0) throw std::invalid_argument("--weight 必须 >= 0");
        } else if (arg == "-n" || arg == "--events") {
            require_args(i, 1, arg.c_str());
            const auto n = parse_ll(argv[++i], arg.c_str());
            if (n <= 0) throw std::invalid_argument("--events 必须为正整数");
            opt.events = static_cast<std::size_t>(n);
        } else if (arg == "--replicas") {
            require_args(i, 1, "--replicas");
            const auto k = parse_ll(argv[++i], "--replicas");
            if (k <= 0 || k > std::numeric_limits<std::int32_t>::max()) {
                throw std::invalid_argument("--replicas 必须为正整数");
            }
            opt.replicas = static_cast<std::int32_t>(k);
        } else if (arg == "--seed") {
            require_args(i, 1, "--seed");
            const auto s = parse_ll(argv[++i], "--seed");
            if (s < 0 || s > 0xffffffffLL) {
                throw std::invalid_argument("--seed 必须在 [0, 2^32-1] 范围内");
            }
            opt.seed = static_cast<UInt_t>(s);
        } else if (arg == "--event-seeded") {
            opt.event_seeded = true;
        } else if (arg == "--sequential-rng") {
            opt.event_seeded = false;
        } else if (arg == "--material") {
            require_args(i, 1, "--material");
            opt.material.preset = argv[++i];
        } else if (arg == "--B") {
            require_args(i, 3, "--B");
            opt.material.external_field = {
                parse_double(argv[i + 1], "--B"),
                parse_double(argv[i + 2], "--B"),
                parse_double(argv[i + 3], "--B")
            };
            opt.material.field_explicit = true;
            i += 3;
        } else if (arg == "--delta") {
            require_args(i, 1, "--delta");
            const double d = parse_double(argv[++i], "--delta");
            require_nonnegative(d, "--delta must be >= 0");
            opt.material.mu_plus_sigma_T = d / gamma_muon;
        } else if (arg == "--sigma-B") {
            require_args(i, 1, "--sigma-B");
            opt.material.mu_plus_sigma_T = parse_double(argv[++i], "--sigma-B");
        } else if (arg == "--hop-rate") {
            require_args(i, 1, "--hop-rate");
            opt.material.mu_plus_hop_rate = parse_double(argv[++i], "--hop-rate");
        } else if (arg == "--muonium-frac") {
            require_args(i, 1, "--muonium-frac");
            opt.material.muonium_fraction = parse_double(argv[++i], "--muonium-frac");
        } else if (arg == "--muonium-sigma-B") {
            require_args(i, 1, "--muonium-sigma-B");
            opt.material.muonium_sigma_T = parse_double(argv[++i], "--muonium-sigma-B");
        } else if (arg == "--muonium-hop-rate") {
            require_args(i, 1, "--muonium-hop-rate");
            opt.material.muonium_hop_rate = parse_double(argv[++i], "--muonium-hop-rate");
        } else if (arg == "--hyperfine") {
            require_args(i, 1, "--hyperfine");
            opt.material.muonium_hyperfine = parse_double(argv[++i], "--hyperfine");
        } else if (arg == "--dipolar-D") {
            require_args(i, 1, "--dipolar-D");
            opt.material.muonium_dipolar_D = parse_double(argv[++i], "--dipolar-D");
        } else if (arg == "--state2-hyperfine") {
            require_args(i, 1, "--state2-hyperfine");
            opt.material.muonium_state2_hyperfine = parse_double(argv[++i], "--state2-hyperfine");
        } else if (arg == "--state2-dipolar-D") {
            require_args(i, 1, "--state2-dipolar-D");
            opt.material.muonium_state2_dipolar_D = parse_double(argv[++i], "--state2-dipolar-D");
        } else if (arg == "--exchange-w12") {
            require_args(i, 1, "--exchange-w12");
            opt.material.muonium_exchange_w12 = parse_double(argv[++i], "--exchange-w12");
        } else if (arg == "--exchange-w21") {
            require_args(i, 1, "--exchange-w21");
            opt.material.muonium_exchange_w21 = parse_double(argv[++i], "--exchange-w21");
        } else if (arg == "--spin-exchange-rate") {
            require_args(i, 1, "--spin-exchange-rate");
            opt.material.muonium_spin_exchange_rate = parse_double(argv[++i], "--spin-exchange-rate");
        } else if (arg == "--reorient-rate") {
            require_args(i, 1, "--reorient-rate");
            opt.material.muonium_reorientation_rate = parse_double(argv[++i], "--reorient-rate");
        } else if (arg == "--domain-mode") {
            require_args(i, 1, "--domain-mode");
            const std::string dm = argv[++i];
            if (dm == "none") opt.material.domain_mode = DomainFieldMode::None;
            else if (dm == "isotropic") opt.material.domain_mode = DomainFieldMode::IsotropicDomains;
            else if (dm == "saturated") opt.material.domain_mode = DomainFieldMode::SaturatedDomain;
            else throw std::invalid_argument("--domain-mode 仅支持 none, isotropic 或 saturated");
        } else if (arg == "--domain-dir") {
            require_args(i, 3, "--domain-dir");
            opt.material.domain_dir = {
                parse_double(argv[i + 1], "--domain-dir"),
                parse_double(argv[i + 2], "--domain-dir"),
                parse_double(argv[i + 3], "--domain-dir")
            };
            i += 3;
        } else if (arg == "--internal-field") {
            require_args(i, 1, "--internal-field");
            opt.material.mu_plus_internal_field_T = parse_double(argv[++i], "--internal-field");
        } else if (arg == "--residual-pol") {
            require_args(i, 1, "--residual-pol");
            opt.material.residual_pol = parse_double(argv[++i], "--residual-pol");
        } else if (arg == "--mu-minus-internal-field") {
            require_args(i, 1, "--mu-minus-internal-field");
            opt.material.mu_minus_internal_field_T = parse_double(argv[++i], "--mu-minus-internal-field");
        } else if (arg == "--mu-minus-sigma-B") {
            require_args(i, 1, "--mu-minus-sigma-B");
            opt.material.mu_minus_sigma_T = parse_double(argv[++i], "--mu-minus-sigma-B");
        } else if (arg == "--mu-minus-hop-rate") {
            require_args(i, 1, "--mu-minus-hop-rate");
            opt.material.mu_minus_hop_rate = parse_double(argv[++i], "--mu-minus-hop-rate");
        } else if (arg == "--max-hops") {
            require_args(i, 1, "--max-hops");
            const auto mh = parse_ll(argv[++i], "--max-hops");
            if (mh <= 0) throw std::invalid_argument("--max-hops 必须为正整数");
            opt.material.max_hops = static_cast<std::size_t>(mh);
        } else if (arg == "--isotope") {
            require_args(i, 1, "--isotope");
            opt.material.copper_isotope = parse_copper_isotope_mode(argv[++i]);
        } else if (arg == "--cu-map") {
            require_args(i, 1, "--cu-map");
            const std::string m = argv[++i];
            if (m == "A") opt.material.cu_map_B = false;
            else if (m == "B") opt.material.cu_map_B = true;
            else throw std::invalid_argument("--cu-map 仅支持 A 或 B");
        } else if (arg == "--f1-pop") {
            require_args(i, 1, "--f1-pop");
            opt.material.f1_pop = parse_double(argv[++i], "--f1-pop");
        } else if (arg == "--f1-orient") {
            require_args(i, 1, "--f1-orient");
            opt.material.f1_orient = parse_double(argv[++i], "--f1-orient");
        } else if (arg == "--f2-orient") {
            require_args(i, 1, "--f2-orient");
            opt.material.f2_orient = parse_double(argv[++i], "--f2-orient");
        } else if (arg == "--cu-hyperfine") {
            require_args(i, 1, "--cu-hyperfine");
            opt.material.cu_hyperfine = parse_double(argv[++i], "--cu-hyperfine");
        } else if (arg == "--m1-rate") {
            require_args(i, 1, "--m1-rate");
            opt.material.m1_rate_ns = parse_double(argv[++i], "--m1-rate");
        } else if (arg == "--dio-rate") {
            require_args(i, 1, "--dio-rate");
            opt.material.dio_rate_ns = parse_double(argv[++i], "--dio-rate");
        } else if (arg == "--capture-rate") {
            require_args(i, 1, "--capture-rate");
            opt.material.capture_rate_ns = parse_double(argv[++i], "--capture-rate");
        } else if (arg == "--input") {
            require_args(i, 1, "--input");
            opt.input_file = argv[++i];
        } else if (arg == "--input-tree") {
            require_args(i, 1, "--input-tree");
            opt.input_tree = argv[++i];
        } else if (arg == "--output") {
            require_args(i, 1, "--output");
            opt.output_file = argv[++i];
        } else if (arg == "--output-tree") {
            require_args(i, 1, "--output-tree");
            opt.output_tree = argv[++i];
        } else if (arg == "--keep-captures") {
            opt.keep_captures = true;
        } else if (arg == "--electrons-only") {
            opt.keep_captures = false;
        } else if (arg == "--overwrite") {
            opt.overwrite = true;
        } else {
            throw std::invalid_argument("未知命令行选项: " + arg + "（使用 --help 查看帮助）");
        }
    }

    if (opt.charge == -1 && !opt.pol_explicit) {
        opt.initial_polarization = {0, 0, -1};
    }
    if (!opt.input_file.empty()) {
        if (!has_root_extension(opt.input_file)) {
            throw std::invalid_argument("方案 B 的 --input 必须为 .root 文件: " + opt.input_file);
        }
        if (opt.output_file.empty()) {
            throw std::invalid_argument("方案 B（ROOT 批处理模式）必须同时指定 --output <out.root>");
        }
    }
    if (!opt.output_file.empty() && !has_root_extension(opt.output_file)) {
        throw std::invalid_argument("--output 必须为 .root 文件: " + opt.output_file);
    }
    return opt;
}

int run_cli(const CliOptions& opt) {
    const MaterialPair mats = build_material_pair(opt.material);
    const bool batch_mode = !opt.input_file.empty();

    if (batch_mode) {
        BatchOptions bopt{};
        bopt.input_tree = opt.input_tree;
        bopt.output_tree = opt.output_tree;
        bopt.overwrite = opt.overwrite;
        bopt.default_charge = opt.charge;
        bopt.default_stop_time_ns = opt.stop_time_ns;
        bopt.default_weight = opt.weight;
        bopt.replicas = opt.replicas;
        bopt.keep_captures = opt.keep_captures;
        bopt.event_seeded = opt.event_seeded;
        bopt.check_unique_event_ids = true;
        bopt.seed = opt.seed;
        bopt.preset_name = opt.material.preset;
        bopt.external_field_T = mats.effective_field;

        TRandom3 seq_rng(opt.seed);
        const BatchSummary summary = process_stopped_muons_root_streaming(
            opt.input_file, opt.output_file, mats.mu_plus, mats.mu_minus, bopt, seq_rng);

        std::cout << "[INFO] spin_dynamic_decay (streaming): 输入缪子 = " << summary.n_input_muons
                  << ", replicas = " << opt.replicas
                  << ", 总演化事例 = " << summary.n_simulated_events
                  << ", 写出 entries = " << summary.n_written_entries
                  << " (Michel/DIO 电子 = " << summary.n_emitted_leptons
                  << ", 核俘获 = " << summary.n_nuclear_captures
                  << ", keep_captures = " << (opt.keep_captures ? "true" : "false")
                  << ", event_seeded = " << (opt.event_seeded ? "true" : "false")
                  << ") -> " << opt.output_file
                  << " (TTree: " << opt.output_tree << ")\n";
        return 0;
    }

    StoppedMuonState single{
        0, opt.stop_position, opt.initial_polarization, opt.charge, opt.stop_time_ns, opt.weight
    };
    validate_stopped_muon_state(single);
    if (single.charge == -1 && !mats.mu_minus.has_value()) {
        throw std::invalid_argument(
            "charge == -1 需要支持 mu- 的材料配置（请使用 --material copper|iron|sulfur 或显式设置参数）");
    }

    TRandom3 seq_rng(opt.seed);
    const double weight_scale = 1.0 / static_cast<double>(opt.replicas);
    std::vector<DecaySimulationResult> results;
    results.reserve(opt.events * static_cast<std::size_t>(opt.replicas));
    for (std::size_t idx = 0; idx < opt.events; ++idx) {
        single.event_id = static_cast<std::int64_t>(idx);
        for (std::int32_t rep = 0; rep < opt.replicas; ++rep) {
            if (opt.event_seeded) {
                auto ev_rng = make_event_replica_rng(opt.seed, single.event_id, rep);
                results.push_back(simulate_stopped_muon_decay(
                    single, mats.mu_plus, mats.mu_minus, ev_rng,
                    static_cast<std::int64_t>(idx), rep, weight_scale));
            } else {
                results.push_back(simulate_stopped_muon_decay(
                    single, mats.mu_plus, mats.mu_minus, seq_rng,
                    static_cast<std::int64_t>(idx), rep, weight_scale));
            }
        }
    }

    if (!opt.output_file.empty()) {
        const Long64_t n_written = write_michel_electrons_root(
            opt.output_file, results, opt.overwrite, opt.output_tree,
            opt.seed, opt.material.preset, mats.effective_field,
            opt.keep_captures, static_cast<Int_t>(opt.replicas), opt.event_seeded);
        Long64_t n_emitted = 0;
        for (const auto& r : results) {
            if (r.has_decay_lepton) ++n_emitted;
        }
        const Long64_t n_sim = static_cast<Long64_t>(results.size());
        std::cout << "[INFO] spin_dynamic_decay: 模拟事例 = " << n_sim
                  << ", 写出 entries = " << n_written
                  << " (Michel/DIO 电子 = " << n_emitted
                  << ", 核俘获 = " << (n_sim - n_emitted)
                  << ") -> " << opt.output_file
                  << " (TTree: " << opt.output_tree << ")\n";
    } else {
        for (const auto& r : results) {
            write_pretty(std::cout, static_cast<Long64_t>(r.muon_entry), r);
        }
    }
    return 0;
}

} // namespace

// ROOT macro entry point — Mode A (default direct single-event call)
int spin_dynamic_decay(
        double stop_x = 0.0,
        double stop_y = 0.0,
        double stop_z = 0.0,
        double pol_x = 0.0,
        double pol_y = 0.0,
        double pol_z = 1.0,
        int charge = +1,
        const char* material_preset = "copper",
        UInt_t seed = 42,
        std::size_t events = 1,
        const char* output_root_file = "") {
    try {
        CliOptions opt;
        opt.stop_position = {stop_x, stop_y, stop_z};
        opt.initial_polarization = {pol_x, pol_y, pol_z};
        opt.pol_explicit = true;
        opt.charge = charge;
        opt.material.preset = material_preset ? material_preset : "copper";
        opt.seed = seed;
        opt.events = events;
        if (output_root_file && output_root_file[0] != '\0') {
            opt.output_file = output_root_file;
        }
        return run_cli(opt);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] spin_dynamic_decay: " << e.what() << '\n';
        return 1;
    }
}

// ROOT macro entry point — Mode B (.root stopped muons -> .root streaming batch)
int spin_dynamic_decay_batch(
        const char* input_root_file,
        const char* output_root_file,
        const char* material_preset = "copper",
        UInt_t seed = 42,
        bool overwrite = false,
        std::int32_t replicas = 1,
        bool keep_captures = true,
        bool event_seeded = true) {
    try {
        CliOptions opt;
        opt.input_file = input_root_file ? input_root_file : "";
        opt.output_file = output_root_file ? output_root_file : "";
        opt.material.preset = material_preset ? material_preset : "copper";
        opt.seed = seed;
        opt.overwrite = overwrite;
        opt.replicas = replicas;
        opt.keep_captures = keep_captures;
        opt.event_seeded = event_seeded;
        return run_cli(opt);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] spin_dynamic_decay_batch: " << e.what() << '\n';
        return 1;
    }
}

int spin_dynamic_decay(
        const char* input_root_file,
        const char* output_root_file,
        const char* material_preset = "copper",
        UInt_t seed = 42,
        bool overwrite = false,
        std::int32_t replicas = 1,
        bool keep_captures = true,
        bool event_seeded = true) {
    return spin_dynamic_decay_batch(
        input_root_file, output_root_file, material_preset, seed,
        overwrite, replicas, keep_captures, event_seeded);
}

#ifndef __CLING__
int main(int argc, char** argv) {
    try {
        const auto opt = parse_cli(argc, argv);
        return run_cli(opt);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] spin_dynamic_decay: " << e.what() << '\n';
        return 1;
    }
}
#endif
