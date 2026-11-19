# 逐事例停止缪子自旋动力学演化与衰变抽样（`spin_dynamic_decay`）

本模块在不修改根目录原有 `crmusr_simulation.cxx` 主程序的前提下，提供一套完全独立的**逐事例（Event-by-event）停止缪子自旋动力学演化与衰变/俘获抽样**工具。给定停止缪子的空间位置、初始自旋极化矢量、电荷以及靶材料参数，即可模拟其在靶内的随机跳跃自旋演化，并输出衰变产生的 Michel 正电子 / 束缚轨道衰变（DIO）电子参数。

---

## 1. 目录结构

- `spin_dynamic_decay.hxx`：统一的 C++ 单事例调用接口与 ROOT `TTree` 读写封装（`StoppedMuonState`、`DecaySimulationResult`、`simulate_stopped_muon_decay()`、`read_stopped_muons_root()`、`write_michel_electrons_root()`）。
- `spin_dynamic_decay.cxx`：独立可执行程序与 ROOT 宏入口（默认**方案 A** 命令行/宏直接输入单缪子位置与材料参数；同时支持**方案 B** `.root` 文件批处理输入输出）。
- `test_spin_dynamic_decay.cxx`：端到端运动学与 `.root -> .root` 回读闭合测试。
- 正缪子（$\mu^+$）、缪子素（Muonium）与各向异性顺磁自由基核心头文件：
  - `stopped_muon.hxx`：普通 $\mu^\pm$ 格点间随机跳跃的动态强碰撞轨迹 `HoppingTrajectory` 与铁磁磁畴内场模型 `DomainFieldModel`。
  - `muonium_hopping.hxx`：各向同性缪子素（$\mu^+\otimes e^-$）$4\times 4$ 联合密度矩阵随机跳跃相干演化 `MuoniumHoppingTrajectory`。
  - `anisotropic_muonium.hxx`：轴对称偶极超精细张量（$A_{\text{iso}}, D$）、粉末取向、强碰撞电子自旋交换（$\rho \mapsto \mathrm{Tr}_e(\rho)\otimes\frac{1}{2}\mathbb{I}_e$）、构型重取向 $4\times 4$ 轨迹 `AnisotropicMuoniumTrajectory`，以及 Cox (2011) 两态化学交换主方程（$30\times 30$ 实无迹 Pauli 基耦合刘维尔算符，带电子自旋记忆与细致平衡 Padé(13,13) 矩阵指数）`TwoStateExchangeEvolution` / `TwoStateExchangeTrajectory`（用于固体硫 $\alpha\text{-S}_8$ 顺磁自由基）。
  - `stopped_decay.hxx`：$\mu^+$ 单事例寿命抽样与无质量最低阶自由 Michel 正电子抽样。
- 负缪子（$\mu^-$，支持铜靶 $^{63}\text{Cu} / {}^{65}\text{Cu}$ [$I=3/2$] 与铁靶 $^{56}\text{Fe}$ / 固体硫靶 $^{32}\text{S}$ [$I=0$]）核心头文件：
  - `bound_dio_data.hxx`、`bound_dio.hxx`（及兼容封装 `copper_dio.hxx`）：基于 14 段切比雪夫 + Bernstein 严格包络的哈希认证多核素有质量 DIO 电子谱表（$^{63}\text{Cu}$、$^{65}\text{Cu}$、$^{56}\text{Fe}$、$^{32}\text{S}$）。
  - `mu_minus_spin0.hxx`：通用零核自旋（$I=0$）束缚 $1s$ $\mu^-$ 2 态演化、磁畴内场与竞争终止引擎（供 $^{56}\text{Fe}$ 与 $^{32}\text{S}$ 共用）。
  - `mu_minus_nuclear.hxx`、`mu_minus_initial_state.hxx`、`mu_minus_m1_conversion.hxx`、`mu_minus_secular_trajectory.hxx`、`mu_minus_termination.hxx`、`mu_minus_event.hxx`：$\mu^-(1/2)\otimes\text{Cu}(3/2)$ 的 $8\times 8$ 级联后有效初态、超精细 + 双自旋塞曼演化及 $F_+\to F_-$ 单向 M1 跃迁随机轨迹。

---

## 2. 编译与测试

在项目根目录下执行：

```bash
mkdir -p spin_dynamic_decay/build

# 编译主程序
$(root-config --cxx) $(root-config --cflags) -O2 -Wall -Wextra -Werror \
  spin_dynamic_decay/spin_dynamic_decay.cxx \
  -o spin_dynamic_decay/build/spin_dynamic_decay $(root-config --libs)

# 编译并运行自检测试
$(root-config --cxx) $(root-config --cflags) -O2 -Wall -Wextra -Werror \
  spin_dynamic_decay/test_spin_dynamic_decay.cxx \
  -o spin_dynamic_decay/build/test_spin_dynamic_decay $(root-config --libs)
./spin_dynamic_decay/build/test_spin_dynamic_decay
```

---

## 3. 使用方法

### 方案 A（默认）：命令行或 ROOT 直接输入单个缪子停止位置与材料参数

#### （1）模拟单个正缪子（$\mu^+$）在铜靶或自定义材料参数下的自旋演化与 Michel 衰变
```bash
./spin_dynamic_decay/build/spin_dynamic_decay \
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
root -l -b -q 'spin_dynamic_decay/spin_dynamic_decay.cxx(1.5, -2.0, 3.2, 0.0, 0.0, 1.0, +1, "copper", 42, 1)'
```

#### （3）铁靶（$\alpha\text{-Fe}$）与固体硫靶（$\alpha\text{-S}_8$）示例
```bash
# 未磁化多畴纯铁靶 mu+（~50.1 MHz，2/3 振荡 + 1/3 纵向尾）
./spin_dynamic_decay/build/spin_dynamic_decay \
  --material iron_pure --charge +1 --B 0 0 0 --events 5 --seed 42

# 商业铁 / 钢板（--material iron_commercial，~47.4 MHz，更快横向弛豫）与单畴饱和外场覆盖
./spin_dynamic_decay/build/spin_dynamic_decay \
  --material iron_commercial --domain-mode saturated --domain-dir 0 0 1 --charge -1 --events 5 --seed 42

# 固体硫靶（300 K 一手文献参数：mu+ 为 Cox 2011 两态化学交换模型 f_D=0；mu- 为 32S DIO/俘获竞争与 Mamedov 2001 f_res=0.073）
./spin_dynamic_decay/build/spin_dynamic_decay \
  --material sulfur --charge +1 --B 0 0 0.3 --events 5 --seed 42
./spin_dynamic_decay/build/spin_dynamic_decay \
  --material sulfur --charge -1 --B 0 0 0 --events 5 --seed 42
```

> **物理参数声明**：
> - $^{56}\text{Fe}$ 与 $^{32}\text{S}$ 的核电荷半径（IAEA LiveChart）、$\mu^-$ 有效寿命与核俘获率（Suzuki et al., *PRC* **35**, 2212 (1987)）、静态有限核 DIO 谱表（Watanabe 1993 近似）及 $\alpha$-Fe $\mu^+$ 内场均来自已核实一手资料。
> - **铁靶 $^{56}\text{Fe}\ \mu^-$（300 K）**：已核对 Keller et al., *PRB* **35**, 2008 (1987) 全文，默认 $B_\mu = -18.8\text{ T}$（Table III + Eq. (9) 平均场与 $\Delta = -0.9\%$ 超精细反常外推至 300 K）、束缚 $1s$ 相对论 $g$ 因子修正 $(g-g_0)/g_0 = -0.85\%$（Table I）、$f_{\text{res}} = 1/6$（Sec. III）、$\lambda_\mu = 0$（Table II）。
> - **固体硫靶 $\alpha\text{-S}_8\ \mu^+$（300 K）**：已核对 Eisenstein et al., *PR* **142**, 217 (1966) 与 Cox et al., *JPCM* **23**, 315801 (2011) 全文，默认无瞬时抗磁分量（$f_D = 0$），采用带电子自旋记忆的两态化学交换主方程：态 1（闭环键心 $\text{S}_8\text{Mu}_{\text{BC}}$）$A_1/(2\pi) = -184\text{ MHz}$、$D_{1,\text{code}}/(2\pi) = \frac{3}{2} D_{1,\text{Cox}}/(2\pi) = +145.5\text{ MHz}$（$A_\parallel = -87.0\text{ MHz},\ A_\perp = -232.5\text{ MHz}$），态 2（开环/弱耦合）$A_2 = D_2 = 0$，300 K 交换速率 $W_{12} = 28.0\text{ ns}^{-1},\ W_{21} = 15.0\text{ ns}^{-1}$（Fig. 9(b)）。
> - **固体硫靶 $^{32}\text{S}\ \mu^-$（300 K）**：已核对 Mamedov et al., *JETP* **93**(5), 941 (2001) 全文，默认 $f_{\text{res}} = 0.073$（Sec. 4 与石墨碳相对不对称度比值 $0.4383/6$）及束缚 $1s$ 相对论 $g$ 因子修正 $(g^{\text{free}}-g^{1s})/g^{\text{free}} = 42.4\times 10^{-4}$（Tables 1 & 2）。
> - 以上所有参数均可通过 `--internal-field`、`--mu-minus-internal-field`、`--residual-pol`、`--muonium-frac`、`--hyperfine`、`--dipolar-D`、`--state2-hyperfine`、`--state2-dipolar-D`、`--exchange-w12`、`--exchange-w21`、`--spin-exchange-rate`、`--reorient-rate` 等选项自由覆盖（调整 $W_{12}, W_{21}$ 即可推广至任意温度）。

---

### 方案 B：ROOT 文件批处理模式（`.root` 输入 $\to$ `.root` 输出）

- **输入 `.root` 文件规范（1 个 `TTree` entry = 1 个停止缪子）**
  - 默认读取名为 `stopped_muons` 的 `TTree`（若文件中仅存在唯一一个 `TTree` 会自动识别，也可通过 `--input-tree <name>` 显式指定）。
  - 必需分支：
    - 停止位置：`stop_x`, `stop_y`, `stop_z`（兼容 `x`, `y`, `z` 或 `pos_x`, `pos_y`, `pos_z`）
    - 停止时实验室系初始极化：`pol_x`, `pol_y`, `pol_z`（兼容 `P0_x`, `P0_y`, `P0_z` 或 `px0`, `py0`, `pz0`）
  - 可选分支（缺省时使用 entry 下标或命令行 `--charge` / `--stop-time` 默认值）：
    - `event_id`（整数，母缪子事例编号）
    - `charge`（整数，`+1` 为 $\mu^+$，`-1` 为 $\mu^-$）
    - `stop_time_ns`（浮点，停止时刻 [ns]）

- **输出 `.root` 文件规范（1 个 `TTree` entry = 1 个 Michel 正电子 / DIO 电子）**
  - 默认输出到名为 `michel_electrons` 的 `TTree`（可通过 `--output-tree <name>` 修改）。
  - **严格保证每个 entry 对应一个真实出射的衰变正/负电子**：当负缪子发生原子核俘获（`nuclear_capture`，无出射衰变电子）时，不向 `michel_electrons` 写入空占位行；文件顶层同时保存 `TParameter<Long64_t>` 统计量 `n_input_muons`、`n_michel_electrons` 和 `n_nuclear_captures` 以保证计数闭合。
  - 每个电子 entry 包含的主要分支：
    - 溯源索引：`electron_id`（当前电子序号）、`muon_entry`（对应输入缪子树的 entry 下标）、`event_id`（母缪子 `event_id`）、`muon_charge`
    - 电子属性与时空顶点：`pdg`（$e^+$ 为 `-11`，$e^-$ 为 `11`）、`charge`（$\pm 1$）、`x`, `y`, `z`（出生顶点位置）、`stop_time_ns`、`elapsed_time_ns`、`decay_time_ns`
    - 运动学参数：`total_energy_MeV`、`kinetic_energy_MeV`、`momentum_MeV_c`、`px_MeV_c`, `py_MeV_c`, `pz_MeV_c`、`dir_x`, `dir_y`, `dir_z`、`endpoint_MeV`、`reduced_energy_x`
    - 自旋演化信息：`P0_x/y/z`、`Pdecay_x/y/z`、`hop_count`、`spin_category`、`termination`、`spectrum_model`

- **命令行运行方案 B**：
  ```bash
  ./spin_dynamic_decay/build/spin_dynamic_decay \
    --input stopped_muons.root \
    --output michel_electrons.root \
    --material copper \
    --B 0 100e-6 0 \
    --seed 123
  ```

- **ROOT 宏运行方案 B**：
  ```bash
  root -l -b -q 'spin_dynamic_decay/spin_dynamic_decay.cxx+("stopped_muons.root", "michel_electrons.root", "copper", 123)'
  ```
  （或调用 `spin_dynamic_decay_batch("stopped_muons.root", "michel_electrons.root", "copper", 123)`。）
