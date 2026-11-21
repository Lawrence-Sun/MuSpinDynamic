# MuSpinDynamic

独立维护的**逐事例（Event-by-event）停止缪子自旋动力学演化与衰变/俘获抽样库及命令行工具**。给定停止缪子的空间位置、初始自旋极化矢量、电荷以及靶材料参数，即可模拟其在靶内的自旋演化，并抽样输出衰变产生的 Michel 正电子或束缚轨道衰变（DIO）电子（以及 $\mu^-$ 原子核俘获终止分支）。

本仓库从 `crmusr-toymc/spin_dynamic_decay` 拆分独立并保留完整提交历史，设计为**纯物理核心与 Geant4 解耦**：
- **`MuSpinDynamic::Core`**（纯头文件 `INTERFACE` target）：仅依赖 `ROOT::GenVector` 与 `ROOT::Matrix`，无 Geant4、无 Mustard、无 `TFile`/`TTree` I/O、无 `TRandom3` 硬绑定，也不会污染全局命名空间；可直接接入 `crmusr-lts` 等 Geant4 全模拟框架（通过 `make_uniform_rng([] { return G4UniformRand(); })` 接入线程/MPI 安全的随机数流）。
- **`MuSpinDynamic::RootIO`**（可选纯头文件 `INTERFACE` target）：封装方案 B 的 ROOT `TTree` 批处理读写接口。
- **`spin_dynamic_decay`**（独立 CLI / ROOT 宏）：支持单事例模拟与 `.root -> .root` 批处理，便于独立物理验证和参数扫描。

---

## 1. 目录结构

```text
MuSpinDynamic/
├── CMakeLists.txt
├── cmake/MuSpinDynamicConfig.cmake.in
├── include/mu_spin_dynamic/
│   ├── types.hxx                      # 单位、物理常量（含 muon_mass_MeV / electron_mass_MeV）与 Vec2/Vec3
│   ├── rng.hxx                        # UniformRngAdapter / make_uniform_rng() 与按 (seed, event_id, replica_id) 确定性派生的 EventReplicaRng
│   ├── spin_dynamic_decay.hxx         # 纯单事例统一 API（含仅自旋演化轻接口 evolve_* 与自然铜同位素封装）与参考材料预设
│   ├── material_config.hxx            # 共享纯头文件材料配置与工厂 build_material_pair()（供 Geant4 模式 1 与 CLI 模式 3 共用）
│   ├── root_io.hxx                    # 可选 ROOT TTree 流式批处理读写（支持事件分裂 --replicas 与核俘获占位）
│   ├── stopped_muon.hxx               # 普通 mu± 格点跳跃强碰撞轨迹 HoppingTrajectory 与铁磁磁畴内场模型 DomainFieldModel
│   ├── muonium_hopping.hxx            # 各向同性缪子素（mu+ ⊗ e-）4×4 密度矩阵相干跳跃演化
│   ├── anisotropic_muonium.hxx        # 轴对称偶极自由基 4×4 轨迹与 Cox (2011) 两态化学交换 30×30 刘维尔主方程（alpha-S8）
│   ├── stopped_decay.hxx              # mu+ 单事例寿命自旋演化与精确有限电子质量 LO Michel 正电子抽样
│   ├── bound_dio_data.hxx             # 14 段切比雪夫 + Bernstein 包络哈希认证多核素有质量 DIO 谱表（63/65Cu, 56Fe, 32S）
│   ├── bound_dio.hxx / copper_dio*.hxx
│   ├── mu_minus_spin0.hxx             # 零核自旋（I=0）束缚 1s mu- 演化、磁畴内场与 DIO/俘获竞争引擎（56Fe, 32S）
│   └── mu_minus_{nuclear,initial_state,m1_conversion,secular_trajectory,termination,event}.hxx
│                                      # mu-(1/2) ⊗ Cu(3/2) 8×8 级联初态、超精细+双自旋塞曼演化及 F+ -> F- M1 跃迁轨迹
├── apps/
│   └── spin_dynamic_decay.cxx         # 独立命令行程序与 ROOT 宏入口（方案 A / 方案 B）
├── notebooks/
│   └── spin_dynamic_decay_analysis.ipynb
│                                      # 上下游 Geant4 联调与铜/铁/硫多靶材物理分析与出图 Notebook
└── tests/
    ├── test_core_only.cxx             # 纯核心隔离测试（仅链接 MuSpinDynamic::Core，验证无 I/O 与全局污染及自定义 RNG）
    └── test_spin_dynamic_decay.cxx    # 完整物理断言与 .root -> .root 回读闭合测试
```

---

## 2. 编译与测试

要求：C++17（或更高标准，已验证兼容 C++23）、CMake ≥ 3.18、ROOT（核心库仅需 `GenVector` + `Matrix`；构建 CLI 和完整测试额外需要 `RIO`、`Tree`、`Physics`、`MathCore`）。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

---

## 3. 作为第三方库接入外部工程（例如 `crmusr-lts` / Geant4）

### CMake 集成（`FetchContent` / `add_subdirectory` / `find_package`）

当作为子工程引入时，`MUSPINDYNAMIC_BUILD_APPS` 与 `MUSPINDYNAMIC_BUILD_TESTING` 默认自动关闭，仅导出轻量纯头文件目标：

```cmake
include(FetchContent)
FetchContent_Declare(
    MuSpinDynamic
    GIT_REPOSITORY <your-repo-url-or-local-path>
    GIT_TAG main
)
FetchContent_MakeAvailable(MuSpinDynamic)

target_link_libraries(YourGeant4Target PRIVATE MuSpinDynamic::Core)
```

### C++ 调用示例（桥接 Geant4 `G4UniformRand()`）

```cpp
#include "mu_spin_dynamic/material_config.hxx"
#include "Randomize.hh" // Geant4 G4UniformRand()

auto rng = mu_spin_dynamic::make_uniform_rng([] { return G4UniformRand(); });

mu_spin_dynamic::StoppedMuonState muon{};
muon.event_id = event_id;
muon.stop_position = {stop_x_mm, stop_y_mm, stop_z_mm};
muon.initial_polarization = {pol_x, pol_y, pol_z};
muon.charge = muon_charge; // +1 或 -1
muon.stop_time_ns = stop_time_ns;

mu_spin_dynamic::MaterialConfig cfg{};
cfg.preset = "copper"; // 铜靶默认自然铜（69.15% 63Cu / 30.85% 65Cu）
cfg.external_field = B_ext_T;
cfg.field_explicit = true;
const auto mats = mu_spin_dynamic::build_material_pair(cfg);

// 方式 A：由 MuSpinDynamic 直接输出有限电子质量 Michel / DIO 末态
const auto res = mu_spin_dynamic::simulate_stopped_muon_decay(
    muon, mats.mu_plus, mats.mu_minus, rng);

if (res.has_decay_lepton) {
    // 出生时刻: res.decay_time_ns [ns]
    // 动量: res.momentum_MeV_c [MeV/c], 动能: res.kinetic_energy_MeV [MeV]
    // PDG: res.lepton_pdg (-11 为 e+, +11 为 e-)
}

// 方式 B（仅自旋演化，供 Geant4/Mustard NLO 辐射衰变等外部末态生成器接管）:
const auto term = mu_spin_dynamic::evolve_mu_plus_to_decay(muon, mats.mu_plus, rng);
// 使用 term.decay_time_ns 与 term.decay_polarization 交由外部衰变通道生成末态
```

---

## 4. 独立 CLI 使用方法

### 方案 A（默认）：命令行或 ROOT 直接输入单个缪子停止位置与材料参数

#### （1）模拟单个正缪子（$\mu^+$）在铜靶或自定义材料参数下的自旋演化与 Michel 衰变

```bash
./build/spin_dynamic_decay \
  --pos 1.5 -2.0 3.2 \
  --pol 0 0 1 \
  --charge +1 \
  --material copper \
  --B 0 100e-6 0 \
  --delta 0.252e-3 \
  --hop-rate 5.06e-3 \
  --seed 42
```

若在方案 A 下附带 `--events N --output out.root`，也会将产生的所有 Michel/DIO 电子写入 `out.root`（每个 entry 对应 1 个衰变电子）。

#### （2）通过 ROOT 宏直接运行方案 A

```bash
root -l -b -q -e 'gInterpreter->AddIncludePath("include")' \
  'apps/spin_dynamic_decay.cxx(1.5, -2.0, 3.2, 0.0, 0.0, 1.0, +1, "copper", 42, 1)'
```

#### （3）铁靶（$\alpha\text{-Fe}$）与固体硫靶（$\alpha\text{-S}_8$）示例

```bash
# 未磁化多畴纯铁靶 mu+（~50.1 MHz，2/3 振荡 + 1/3 纵向尾）
./build/spin_dynamic_decay \
  --material iron_pure --charge +1 --B 0 0 0 --events 5 --seed 42

# 商业铁 / 钢板（--material iron_commercial，~47.4 MHz，更快横向弛豫）与单畴饱和外场覆盖
./build/spin_dynamic_decay \
  --material iron_commercial --domain-mode saturated --domain-dir 0 0 1 --charge -1 --events 5 --seed 42

# 固体硫靶（300 K 一手文献参数：mu+ 为 Cox 2011 两态化学交换模型 f_D=0；mu- 为 32S DIO/俘获竞争与 Mamedov 2001 f_res=0.073）
./build/spin_dynamic_decay \
  --material sulfur --charge +1 --B 0 0 0.3 --events 5 --seed 42
./build/spin_dynamic_decay \
  --material sulfur --charge -1 --B 0 0 0 --events 5 --seed 42
```

> **物理参数声明**：
> - $^{56}\text{Fe}$ 与 $^{32}\text{S}$ 的核电荷半径（IAEA LiveChart）、$\mu^-$ 有效寿命与核俘获率（Suzuki et al., *PRC* **35**, 2212 (1987)）、静态有限核 DIO 谱表（Watanabe 1993 近似）及 $\alpha$-Fe $\mu^+$ 内场均来自已核实一手资料。
> - **铁靶 $^{56}\text{Fe}\ \mu^-$（300 K）**：已核对 Keller et al., *PRB* **35**, 2008 (1987) 全文，默认 $B_\mu = -18.8\text{ T}$（Table III + Eq. (9) 平均场与 $\Delta = -0.9\%$ 超精细反常外推至 300 K）、束缚 $1s$ 相对论 $g$ 因子修正 $(g-g_0)/g_0 = -0.85\%$（Table I）、$f_{\text{res}} = 1/6$（Sec. III）、$\lambda_\mu = 0$（Table II）。
> - **固体硫靶 $\alpha\text{-S}_8\ \mu^+$（300 K）**：已核对 Eisenstein et al., *PR* **142**, 217 (1966) 与 Cox et al., *JPCM* **23**, 315801 (2011) 全文，默认无瞬时抗磁分量（$f_D = 0$），采用带电子自旋记忆的两态化学交换主方程：态 1（闭环键心 $\text{S}_8\text{Mu}_{\text{BC}}$）$A_1/(2\pi) = -184\text{ MHz}$、$D_{1,\text{code}}/(2\pi) = \frac{3}{2} D_{1,\text{Cox}}/(2\pi) = +145.5\text{ MHz}$（$A_\parallel = -87.0\text{ MHz},\ A_\perp = -232.5\text{ MHz}$），态 2（开环/弱耦合）$A_2 = D_2 = 0$，300 K 交换速率 $W_{12} = 28.0\text{ ns}^{-1},\ W_{21} = 15.0\text{ ns}^{-1}$（Fig. 9(b)）。
> - **固体硫靶 $^{32}\text{S}\ \mu^-$（300 K）**：已核对 Mamedov et al., *JETP* **93**(5), 941 (2001) 全文，默认 $f_{\text{res}} = 0.073$（Sec. 4 与石墨碳相对不对称度比值 $0.4383/6$）及束缚 $1s$ 相对论 $g$ 因子修正 $(g^{\text{free}}-g^{1s})/g^{\text{free}} = 42.4\times 10^{-4}$（Tables 1 & 2）。
> - 以上所有参数均可通过 `--internal-field`、`--mu-minus-internal-field`、`--residual-pol`、`--muonium-frac`、`--hyperfine`、`--dipolar-D`、`--state2-hyperfine`、`--state2-dipolar-D`、`--exchange-w12`、`--exchange-w21`、`--spin-exchange-rate`、`--reorient-rate` 等选项自由覆盖（调整 $W_{12}, W_{21}$ 即可推广至任意温度）。

---

### 方案 B：ROOT 文件流式批处理模式（`.root` 输入 $\to$ `.root` 输出）

- **输入 `.root` 文件规范（1 个 `TTree` entry = 1 个停止缪子）**
  - 默认读取名为 `stopped_muons` 的 `TTree`（若文件中仅存在唯一一个 `TTree` 会自动识别，也可通过 `--input-tree <name>` 显式指定）。
  - 必需分支：
    - 停止位置：`stop_x`, `stop_y`, `stop_z`（兼容 `x`, `y`, `z` 或 `pos_x`, `pos_y`, `pos_z`）
    - 停止时实验室系初始极化：`pol_x`, `pol_y`, `pol_z`（兼容 `P0_x`, `P0_y`, `P0_z` 或 `px0`, `py0`, `pz0`）
  - 可选分支（缺省时使用 entry 下标或命令行 `--charge` / `--stop-time` / `--weight` 默认值）：
    - `event_id`（整数，母缪子事例编号；默认 `--event-seeded` 模式下要求唯一）
    - `charge`（整数，`+1` 为 $\mu^+$，`-1` 为 $\mu^-$）
    - `stop_time_ns`（浮点，停止时刻 [ns]）
    - `weight`（浮点，母事例权重，兼容 `Weight` / `event_weight`，默认 `1.0`）

- **输出 `.root` 文件规范（$O(1)$ 内存流式写出，支持事件分裂与核俘获占位）**
  - 默认输出到名为 `michel_electrons` 的 `TTree`（可通过 `--output-tree <name>` 修改）。
  - **默认保留 $\mu^-$ 核俘获占位（`--keep-captures` 默认开启）**：每条输入停缪对应严格 `K` 条输出（`--replicas K`，权重自动缩放为 `weight / K`）。通过 `has_decay_lepton`（`Bool_t`）区分出射电子与核俘获占位；核俘获行 `electron_id = -1`、电子能量/动量分支填 `NAN`。若仅需纯电子清单，可显式传 `--electrons-only`。
  - **默认按事例派生确定性随机流（`--event-seeded` 默认开启）**：由 `(seed, event_id, replica_id)` 经 SplitMix64 + Xoshiro256** 派生独立随机流，支持集群分片并行结果逐事例一致，并自动拦截重复 `event_id`；如需传统单流顺序 `TRandom3`，可传 `--sequential-rng`。
  - 文件顶层保存汇总元数据 `n_input_muons`、`n_simulated_events`、`n_written_entries`、`n_michel_electrons`、`n_nuclear_captures`、`n_replicas`、`keep_captures`、`event_seeded`、`seed`、`material_preset` 与 `external_field_T`。
  - 主要分支：
    - 溯源与权重：`electron_id`、`muon_entry`、`event_id`、`replica_id`、`weight`、`has_decay_lepton`、`muon_charge`
    - 电子属性与时空顶点：`pdg`（$e^+$ 为 `-11`，$e^-$ 为 `11`，核俘获为 `0`）、`charge`（$\pm 1$ 或 `0`）、`x`, `y`, `z`、`stop_time_ns`、`elapsed_time_ns`、`decay_time_ns`
    - 运动学参数：`total_energy_MeV`、`kinetic_energy_MeV`、`momentum_MeV_c`、`px_MeV_c`, `py_MeV_c`, `pz_MeV_c`、`dir_x`, `dir_y`, `dir_z`、`endpoint_MeV`、`reduced_energy_x`
    - 自旋演化信息：`P0_x/y/z`、`Pdecay_x/y/z`、`hop_count`、`spin_category`、`termination`、`spectrum_model`

- **命令行运行方案 B（含事件分裂示例）**：

  ```bash
  ./build/spin_dynamic_decay \
    --input stopped_muons.root \
    --output michel_electrons.root \
    --material copper \
    --replicas 10 \
    --B 0 100e-6 0 \
    --seed 123
  ```

- **ROOT 宏运行方案 B**：

  ```bash
  root -l -b -q -e 'gSystem->AddIncludePath("-Iinclude")' \
    'apps/spin_dynamic_decay.cxx+("stopped_muons.root", "michel_electrons.root", "copper", 123)'
  ```
