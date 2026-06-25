# 农业生态系统过程模型深度补充报告：可运行代码、严格率定、机器学习混合建模与临界点动力学

## TL;DR
- 本报告把先前的理论拆解为**可运行、可调试、可率定、可扩展**的工程实践：提供 RothC 五库模型的完整 Python 实现并复现 Rothamsted 官方算例、可插拔环境修正函数库、最小作物模型、DNDC "厌氧气球" 氧化还原核心的伪代码到代码，以及五类数值陷阱（质量不守恒、时间步不收敛、显式 Euler 刚性发散、负库、spin-up 未平衡）的检测与修复。
- 在率定层面，系统讲解 Morris/Sobol 敏感性分析（SALib）、PEST/PEST++ 的 Gauss-Marquardt-Levenberg 算法与 Tikhonov 正则化、贝叶斯 MCMC（Metropolis-Hastings、DREAM(ZS)、emcee、NUTS/PyMC）、等效性（equifinality）与 GLUE，以及 RMSE/NSE/KGE 等目标函数与多目标 Pareto 率定。
- 在前沿与理论层面，覆盖物理信息/知识引导机器学习（PINN、KGML、可微分建模）、过程模型的 ML 仿真器与残差混合建模（含 KGML-ag 的 N₂O 案例），以及把迭代模拟与突变临界点联系起来的动力系统理论（不动点稳定性、鞍结/叉式/Hopf 分岔、迟滞、临界慢化与早期预警信号），并明确指出 DNDC 的 Eh 阈值、土壤碳的多稳态、作物物候 GDD 阈值如何产生真实的 regime shift 而非代码 bug。

---

## Key Findings（核心结论）

1. **RothC 是入门复现的最佳起点**：它只有 5 个库、月时间步、解析的一阶衰减解 `C(t+1)=C(t)·exp(-abck·t)`，且 Rothamsted 官方文档给出了可逐位复现的算例（DPM 0.1533→0.1140 等）。掌握它就掌握了所有过程模型共有的 `f(T)·f(W)·f(substrate)` 乘性骨架。

2. **数值陷阱大多来自显式 Euler 与多尺度刚性**：tipping-bucket 之所以稳定，是因为它用"溢出"逻辑而非扩散方程，避免了 `Δt ≤ Δz²/(2D)` 的 von Neumann 稳定性约束；DNDC 必须用小时步 + 算子分裂（operator splitting）正是因为硝化/反硝化/发酵是刚性快过程。

3. **率定的黄金顺序是"先筛选、再优化、最后量化不确定性"**：用 Morris 廉价筛选出 top-5~10 参数，用 Sobol 量化方差贡献与交互，用 PEST（频率派梯度法）或 DREAM/emcee/PyMC（贝叶斯）得到参数后验分布。单一"最优值"会掩盖 equifinality，分布远胜点估计。

4. **混合建模的三条路线各有取舍**：仿真器（emulator/surrogate）最易上手且回报最高；残差/参数化混合（用 NN 替换某个 `f(W)`）保留 ODE 骨架与可解释性；可微分建模（differentiable modeling）是终极形态但需要把整个模型用 PyTorch/JAX 重写为可自动微分。

5. **农业生态模型中的"突变"是真实涌现现象**：DNDC 在 Eh<500 mV 触发反硝化、Eh<−150 mV 触发产甲烷，本质是鞍结分岔式的阈值切换；土壤碳在管理/土地利用变化下存在多稳态；乘性限制因子与饱和动力学（Michaelis-Menten）是非线性的来源。临界慢化（自相关上升、方差上升）可用来在时间序列中提前预警。

---

## Details（详细内容）

# 区域一：分阶段可运行代码 + 测试基线 + 数值调试

## STAGE 1 — RothC 的完整 Python 复现

### 1.1 模型结构与方程

RothC-26.3（Coleman & Jenkinson 1996；模型描述与用户指南更新于 2014）把土壤有机碳（SOC）分为 5 个库：DPM（Decomposable Plant Material 易分解植物物料）、RPM（Resistant Plant Material 抗性植物物料）、BIO（Microbial Biomass 微生物量）、HUM（Humified Organic Matter 腐殖化有机质）、IOM（Inert Organic Matter 惰性有机质）。前 4 个库各以一阶过程分解，年分解速率常数为：

- k_DPM = 10.0 yr⁻¹
- k_RPM = 0.3 yr⁻¹
- k_BIO = 0.66 yr⁻¹
- k_HUM = 0.02 yr⁻¹

（这些常数由 Rothamsted 长期试验调出，使用时一般不改动。）IOM 不参与周转。

每个活性库在一个月内的解析衰减为：

`C(t+Δt) = C(t) · exp(−a·b·c·k·Δt)`，其中 Δt = 1/12 年。

三个速率修正因子（rate-modifying factors）：

**(a) 温度修正因子 a（Coleman & Jenkinson 形式）：**

```
a = 47.91 / (1 + exp(106.06 / (T + 18.27)))
```
T 为月平均气温（°C）。注意该函数在低温处不为零、在高温处饱和。

**(b) 土壤水分修正因子 b（基于累积 TSMD）：**
先算最大土壤水分亏缺 Max_TSMD（0–23 cm 层、%clay 为黏粒百分比）：
```
Max_TSMD = −(20 + 1.3·%clay − 0.01·%clay²) · (depth/23)
```
（对 Rothamsted %clay=23，Max_TSMD≈−44 mm。）累积 TSMD 从开放水面蒸发首次超过降雨的月份起累积，封顶于 Max_TSMD。然后：
```
若 Acc.TSMD < 0.444·Max_TSMD:  b = 0.2 + (1−0.2)·(Max_TSMD − Acc.TSMD)/(Max_TSMD − 0.444·Max_TSMD)
否则:                          b = 1.0
```
最小值为 0.2。（该分段形式见于 arXiv:2108.00077 对 RothC 的形式化推导，与官方用户指南一致。）

**(c) 土壤覆盖修正因子 c：** 有植被覆盖时 c=0.6，裸地 c=1.0。

**分解产物的去向（黏粒控制）：** 每个库分解出的碳按比例 x 分给 CO₂ 与 (BIO+HUM)：
```
x = 1.67·(1.85 + 1.60·exp(−0.0786·%clay))
CO2/(BIO+HUM) = x
```
即每分解 1 单位碳，去 (BIO+HUM) 的比例为 `1/(x+1)`，去 CO₂ 的比例为 `x/(x+1)`。(BIO+HUM) 再按 **BIO:HUM = 0.46:0.54** 拆分。

**IOM 初始化（Falloon et al. 1998, *Soil Biology & Biochemistry* 30:1207）：**
```
IOM = 0.049 · SOC^1.139   (t C/ha)
```

**进入碳的 DPM/RPM 划分：** 多数农作物默认 DPM/RPM = 1.44（即 59% DPM、41% RPM）；非改良草地/灌丛 0.67；落叶/热带林 0.25。

### 1.2 可运行 Python 代码

```python
import numpy as np

def rate_temp(T):
    "温度修正因子 a (Coleman & Jenkinson)"
    return 47.91 / (1.0 + np.exp(106.06 / (T + 18.27)))

def rate_moist(acc_tsmd, max_tsmd):
    "水分修正因子 b, 最小 0.2"
    thr = 0.444 * max_tsmd          # max_tsmd 为负值
    if acc_tsmd > thr:              # 亏缺较小(更接近0)
        return 1.0
    return 0.2 + (1.0 - 0.2) * (max_tsmd - acc_tsmd) / (max_tsmd - thr)

def clay_x(clay):
    "CO2/(BIO+HUM) 比值"
    return 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay))

def iom_falloon(soc):
    return 0.049 * soc**1.139

K = np.array([10.0, 0.3, 0.66, 0.02])   # DPM,RPM,BIO,HUM (yr^-1)

def rothc_step(pools, a, b, c, clay, c_input=(0.0,0.0), fym=0.0, dt=1/12):
    """单月步进。pools=[DPM,RPM,BIO,HUM,IOM]。返回新库与本月CO2。"""
    DPM, RPM, BIO, HUM, IOM = pools
    rm = a * b * c                        # 综合速率修正
    active = np.array([DPM, RPM, BIO, HUM])
    # 解析衰减（每库分解掉的量）
    decayed = active * (1.0 - np.exp(-rm * K * dt))
    total_dec = decayed.sum()
    x = clay_x(clay)
    to_bh = total_dec / (x + 1.0)         # 去 BIO+HUM
    to_co2 = total_dec - to_bh            # 去 CO2
    dBIO = 0.46 * to_bh
    dHUM = 0.54 * to_bh
    # 衰减后的剩余
    DPM2, RPM2, BIO2, HUM2 = active - decayed
    BIO2 += dBIO
    HUM2 += dHUM
    # 外部碳输入（残体按 DPM/RPM=1.44; FYM 另有划分）
    cin = sum(c_input)
    DPM2 += c_input[0]; RPM2 += c_input[1]
    return np.array([DPM2, RPM2, BIO2, HUM2, IOM]), to_co2
```

### 1.3 单元测试：复现 Rothamsted 官方算例

官方文档给出某月综合修正因子 abc = 0.3561，时间步 Δt = 1/12，各库变化为：
- DPM: 0.1533 · exp(−10·0.3561/12) = 0.1140
- RPM: 4.4852 · exp(−0.3·0.3561/12) = 4.4455
- BIO: 0.6671 · exp(−0.66·0.3561/12) = 0.6542
- HUM: 25.8576 · exp(−0.02·0.3561/12) = 25.8423

```python
def test_rothc_official():
    abc = 0.3561; dt = 1/12
    pools = np.array([0.1533, 4.4852, 0.6671, 25.8576])
    expected = np.array([0.1140, 4.4455, 0.6542, 25.8423])
    got = pools * np.exp(-K * abc * dt)
    assert np.allclose(got, expected, atol=1e-4), got
    print("PASS", np.round(got, 4))

test_rothc_official()
```
此测试只验证衰减项（官方算例就是逐库的纯衰减），通过即说明速率常数、abc 因子、解析解的实现正确。完整模型还需把 (BIO+HUM) 回流与碳输入加上。

> **Fortran/C# 提示**：官方 RothC 用 Fortran（Shell.for），库以 `REAL` 数组存储；C#/.NET 移植时注意 `Math.Exp` 与按月循环。R 生态中 **SoilR** 包（Sierra, Müller & Trumbore 2012, "Models of soil organic matter decomposition: the SoilR package, version 1.0," *Geoscientific Model Development* 5(4):1045–1060, DOI:10.5194/gmd-5-1045-2012）把 RothC 写成线性库模型 `dC/dt = ξ(t)·A·C + I`，矩阵 A 即转移矩阵，是验证你 Python 实现的极佳交叉参照。

---

## STAGE 2 — 可插拔环境修正函数库

过程模型的核心是 `f(T)·f(W)·f(substrate)` 乘性框架。把每个因子写成可互换函数：

```python
import numpy as np

def f_Q10(T, Tref=20.0, Q10=2.0):
    return Q10 ** ((T - Tref) / 10.0)

def f_arrhenius(T_celsius, Ea=60000.0):
    "Ea J/mol; R=8.314"
    Tk = T_celsius + 273.15
    R = 8.314
    return np.exp(-Ea / (R * Tk))   # 常需对参考温度归一化

def f_beta_cardinal(T, Tmin=5, Topt=25, Tmax=40):
    "三基点 beta 函数, 0..1"
    T = np.asarray(T, float)
    out = np.zeros_like(T)
    mask = (T > Tmin) & (T < Tmax)
    a = np.log(2.0) / np.log((Tmax - Tmin) / (Topt - Tmin))
    num = 2*((T-Tmin)**a)*((Topt-Tmin)**a) - ((T-Tmin)**(2*a))
    den = (Topt - Tmin)**(2*a)
    out[mask] = (num/den)[mask]
    return np.clip(out, 0, 1)

def f_wfps_bell(wfps, opt=0.6, width=0.25):
    "WFPS 钟形水分曲线 (water-filled pore space)"
    return np.exp(-((wfps - opt)/width)**2)

# 演示：交换 f(T) 对输出的影响
def decomp_rate(k, T, wfps, fT=f_Q10):
    return k * fT(T) * f_wfps_bell(wfps)

for fT in (f_Q10, f_beta_cardinal):
    print(fT.__name__, decomp_rate(0.5, 30.0, 0.6, fT))
```
把 `fT` 作为参数传入即可在 Q10、Arrhenius、beta 之间切换，对比同一温度下分解速率的差异——这正是模型结构不确定性的来源之一。

---

## STAGE 3 — 最小作物模型（GDD + Beer-Lambert + RUE + 水分胁迫 + tipping-bucket）

理论骨架（Monteith 辐射利用效率范式）：
- **物候**：积温 GDD = Σ max(0, (Tmax+Tmin)/2 − Tbase)，达到阶段阈值即切换发育期。
- **光截获（Beer-Lambert）**：截获比例 = 1 − exp(−k·LAI)，k 为消光系数（extinction coefficient，约 0.4–0.7）。
- **生物量（RUE）**：dDM = RUE · (1−exp(−k·LAI)) · PAR · SWFAC。
- **水分胁迫 SWFAC**：由供需比驱动，0（全胁迫）到 1（无胁迫）。
- **tipping-bucket 土壤水**：分层 DUL（drained upper limit 田间持水）、LL（lower limit 凋萎点）、SAT（饱和）、SWCON（排水系数），超过 DUL 的水按 SWCON 比例下渗到下一层。

```python
import numpy as np

def crop_step(state, wx, par, p):
    gdd = max(0.0, (wx['tmax']+wx['tmin'])/2 - p['tbase'])
    state['gdd'] += gdd
    # LAI 简化为随 GDD 线性增长到上限后衰减
    state['lai'] = min(p['laimax'], p['laimax']*state['gdd']/p['gdd_mat'])
    fint = 1.0 - np.exp(-p['k'] * state['lai'])
    swfac = water_balance(state, wx, p, fint)
    dDM = p['rue'] * fint * par * swfac
    state['dm'] += dDM
    return state

def water_balance(state, wx, p, fint):
    sw = state['sw']                 # 各层含水(mm)
    dul, ll, sat = p['dul'], p['ll'], p['sat']
    # 入渗
    infil = wx['rain']
    for L in range(len(sw)):
        sw[L] += infil
        excess = max(0.0, sw[L] - dul[L])
        drain = excess * p['swcon']
        sw[L] -= drain
        infil = drain
    # 蒸腾需求与供给
    pet = wx['pet'] * fint
    avail = np.maximum(0.0, sw - ll).sum()
    uptake = min(pet, avail)
    swfac = 1.0 if pet <= 0 else min(1.0, uptake/pet)
    # 按比例从各层取水
    if avail > 0:
        for L in range(len(sw)):
            sw[L] -= uptake * max(0.0, sw[L]-ll[L]) / avail
    state['sw'] = sw
    return swfac
```
这与 DSSAT-CSM 的 tipping-bucket（SOILDYN/SPAM 模块）和 APSIM 的 SoilWat 思路一致。tipping-bucket 的关键优势在 STAGE 4 后的数值讨论中说明。

---

## STAGE 4 — DNDC "厌氧气球"（anaerobic balloon）氧化还原核心

### 4.1 关键方程（来自 DNDC v9.5 *Scientific Basis and Processes*, UNH 2017）

**Nernst 方程（官方 v9.5 印刷形式，公式 1）：**
```
Eh = E0 + (R·T)/(n·F) · ln([oxidant]/[reductant])
```
其中 Eh 为氧化还原电位（V）；E0 为标准电动势（V，随当前主导氧化还原电对取值）；R = 8.314 J/mol/K；T = 273 + t（绝对温度，t 为 °C）；n = 转移电子数；F = 96485 C/mol；[oxidant]、[reductant] 为主导氧化剂/还原剂浓度（mol/L）。

> **重要符号约定提醒**：先前报告与许多教材写作 `Eh = E0 − (RT/nF)·ln([red]/[ox])`，与 v9.5 印刷形式符号相反、比值倒置，但二者代数等价（因 ln(ox/red) = −ln(red/ox)）。实现时务必与你所参照的版本一致。

**顺序还原序列（氧化剂按 Gibbs 自由能依次被消耗，依据 Stumm & Morgan 1981；Li et al. 2004）：**
```
O2 → NO3⁻ → Mn⁴⁺ → Fe³⁺ → SO4²⁻ → CO2(产甲烷终端受体)
```

**厌氧气球体积分数（anvf，来自 DNDC 代码级方程，归于 Li, Aber, Stange, Butterbach-Bahl & Papen 2000）：**
```
anvf_L = a · (1 − (b − pO2[L]/pO2_air))
```
anvf 为第 L 层厌氧微位点体积分数；pO2[L] 为该层氧分压（由一维氧扩散-消耗平衡算出）；pO2_air 为大气氧分压；a、b 为经验系数。氧扩散系数 `D_s[L] = (afps[L]^3.33 / afps_max[L]^2.0)·D_air`，D_air = 0.07236 m²/h（Millington & Quirk 1961）。当 pO2→0，anvf 达最大值——气球"爆裂"，下一个氧化剂（NO3⁻）接管，新气球诞生。

**底物分配**：DOC、NO3⁻、NH4⁺ 按 anvf 分配——anvf 内的底物进入还原反应（反硝化、产甲烷），(1−anvf) 部分进入氧化反应（硝化、甲烷氧化）。

**Michaelis-Menten / dual-Monod（公式 2，氧化剂还原分数，Paul & Clark 1989）：**
```
F[oxidant] = a · [DOC/(b+DOC)] · [oxidant/(c+oxidant)]
```
Nernst 与 Michaelis-Menten **共享 [oxidant] 这一公因子**，anvf 是连接二者的动力学桥梁（这正是 Li 2007 与 Gilhespy et al. 2014 所强调的"两方程合并构成 DNDC 核心"的数学实质）。

**反硝化菌 dual-Monod 生长（Table 3, 公式 33–35）：**
```
u_NOx = u_NOx,max · [DOC/(Kc+DOC)] · [NOx/(Kn+NOx)]
```
代码级参数值：Kc(可溶性碳半饱和)=0.017 kg C/m³；Kn(N氧化物)=0.083；最大生长率 GR_NO3,max=GR_NO2,max=0.67 h⁻¹，GR_NO,max=GR_N2O,max=0.34 h⁻¹；反硝化菌 C/N=3.45（基于 Paracoccus denitrificans，Verseveld & Stouthamer 1978）。

**Eh 阈值（来自 v9.5 正文）：**
- 反硝化启动：Eh ≤ **500 mV**（氧耗尽后）。
- 产甲烷启动：Eh < **−150 mV**（在 NO3⁻、Mn⁴⁺、Fe³⁺、SO4²⁻ 依次耗尽之后）。
- 好氧区 Eh 范围 100–650 mV（有利分解、硝化、甲烷氧化）；厌氧区 Eh −300~0 mV。
- 反硝化内部序列：NO3⁻→NO2⁻→NO→N2O→N2（公式 32）。
- 甲烷氧化的 Eh 依赖：`CH4ox = CH4[L]·exp(8.6711·Eh[L]/1000)`（Schütz et al. 1989）。

### 4.2 时间步与算子分裂结构

DNDC 分两大组件、六个子模型：第一组件（土壤气候、作物生长、分解）在**日步**预测土壤温度、湿度、pH、Eh 与底物剖面；第二组件（硝化、反硝化、发酵）预测 C/N 气体通量。Eh 用 Nernst 在**日步**计算，但内部微生物动力学子模型以**小时步**积分（反硝化菌动力学系数单位均为 per hour，每日积分 24 个子步）。各过程按算子分裂顺序求解，作为"连续反应（consecutive reactions）"，每个反应在特定 Eh 条件下发生。土层厚约 2 cm；CH₄ 产生与氧化可在同一土层内同时发生，但分处由 Eh 决定的好氧/厌氧子体积（Li 2007）。

### 4.3 伪代码到代码

```python
R, F = 8.314, 96485.0
def nernst_Eh(E0, T_c, n, ox, red):
    Tk = 273.15 + T_c
    return E0 + (R*Tk)/(n*F) * np.log(max(ox,1e-12)/max(red,1e-12))

def anaerobic_fraction(pO2, pO2_air=0.209, a=1.0, b=1.0):
    return np.clip(a*(1.0 - (b - pO2/pO2_air)), 0.0, 1.0)

def dndc_day(layer, T, dt_hours=24):
    """日步内对 redox 序列做算子分裂 + 小时子步"""
    # 1) 分解 -> 产生 DOC, NH4
    decompose(layer, T)
    # 2) 计算 Eh 与厌氧分数
    Eh = nernst_Eh(layer['E0'], T, layer['n'], layer['ox'], layer['red'])
    anvf = anaerobic_fraction(layer['pO2'])
    # 3) 底物分配
    aer = {k: layer[k]*(1-anvf) for k in ('DOC','NH4','NO3')}
    ana = {k: layer[k]*anvf      for k in ('DOC','NH4','NO3')}
    for h in range(dt_hours):              # 小时子步
        nitrify(aer, T)                    # 仅好氧区
        if Eh <= 0.5:                      # 500 mV
            denitrify(ana, T)              # NO3->NO2->NO->N2O->N2
        if Eh < -0.15:                     # -150 mV
            methanogenesis(ana, T)
    recombine(layer, aer, ana)
```
为什么是小时步 + 算子分裂：硝化/反硝化是相对分解快几个数量级的"快过程"，若与慢的分解共用日步显式积分会发散；算子分裂让每个过程在自己合适的步长上顺序求解。

---

## 数值陷阱调试（NUMERICAL PITFALL DEBUGGING）

### (a) 质量不守恒（mass-balance non-closure）
在每步加入守恒检查：所有库之和 + 累积 CO₂ + 累积淋失，应恒等于初始库 + 累积输入。
```python
def check_mass_balance(pools, cum_co2, cum_leach, init_total, cum_input, tol=1e-9):
    lhs = pools.sum() + cum_co2 + cum_leach
    rhs = init_total + cum_input
    err = abs(lhs - rhs)
    assert err < tol, f"质量不守恒: 误差={err:.3e}"
```
误差随步数线性累积通常说明某个通量被算了两次或漏算；突然跳变通常是输入事件（施肥、残体）未计入 cum_input。

### (b) 时间步不收敛（step-size convergence test）
把 Δt 减半，若解显著变化，说明当前步长太大。
```python
def converged(run, dt, atol=1e-3):
    a = run(dt); b = run(dt/2)
    return np.max(np.abs(a-b)) < atol
```
显式 Euler 的局部截断误差 O(Δt²)、全局 O(Δt)；步长减半误差应约减半，否则可能进入不稳定区。

### (c) 刚性/快过程发散与振荡（CFL / von Neumann 稳定性）
对扩散型方程（如 Richards 方程数值解、热扩散），显式格式的稳定性约束为：
```
Δt ≤ Δz² / (2·D)
```
D 为扩散系数，Δz 为空间步长。违反则解出现振荡放大并发散。**tipping-bucket 之所以稳定**：它不解扩散方程，而用"库满即溢"的代数级联，没有 D，因而不受 CFL 约束——这是它在日步下稳健的根本原因。对策：(1) 减小到小时步；(2) 用算子分裂把快慢过程分开；(3) 对刚性项改用隐式方法（implicit/backward Euler 无条件稳定，但每步需解方程）。

### (d) 负库/负浓度（negative pools）
显式步过大时 `C − k·C·Δt` 可能变负。两种防护：(1) 用解析解 `C·exp(−k·Δt)`（恒正，RothC 即如此）；(2) 限制 `Δt < 1/k_max`（CFL 式约束），或对通量做 `min(flux, available)` 截断。切忌简单地 `max(C,0)` 截断——那会悄悄破坏质量守恒（见 (a)）。

### (e) 初始化/spin-up 未达平衡（慢库漂移）
慢库（HUM、IOM、CENTURY 的 passive 库）周转以百年计。检测：监控慢库年际变化率，若 `|ΔC_slow/C_slow| > ε`（如 0.1%/yr）则尚未平衡。对策：(1) 循环驱动气象数据跑数千年到稳态；(2) 用解析稳态解直接初始化 `C* = I/(ξ·k)`（线性库模型）或 Falloon 式经验公式；(3) 用谱/特征值法估计最慢模态的弛豫时间 τ = 1/(min eigenvalue)，跑 3–5τ。

> **数值方法参考**：Press et al. *Numerical Recipes*（显式/隐式、刚性）、Hairer & Wanner *Solving ODEs II: Stiff Problems*、LeVeque *Finite Difference Methods*（von Neumann 稳定性）。模型描述论文：Coleman & Jenkinson（RothC）、Parton et al. 1987/1988（CENTURY）、Li et al. 1992 + Li 2007 + Gilhespy et al. 2014（DNDC）、Jones et al. 2003 + Hoogenboom et al.（DSSAT）、Keating et al. 2003 + Holzworth et al. 2014（APSIM）。

---

# 区域二：参数率定实践（敏感性分析 + PEST + 贝叶斯/MCMC）

## 2.1 敏感性分析（Sensitivity Analysis）

### Morris 法（elementary effects，初等效应）
对每个参数沿轨迹做一次扰动 Δ，计算初等效应 `EE_i = [f(x+Δe_i) − f(x)]/Δ`。在 r 条轨迹上汇总：
- **μ\*** = mean(|EE_i|)：参数总体重要性（含交互与非线性）。
- **σ** = std(EE_i)：交互/非线性强度。
Morris 极廉价（成本约 r·(k+1) 次模型运行），用于**筛选**——快速识别哪些参数可忽略。

### Sobol 方差分解（variance-based）
把输出方差分解为各参数及其交互的贡献：
- **一阶指数 S_i** = V[E(Y|X_i)]/V(Y)：参数 i 单独的贡献。
- **总效应 STi** = 1 − V[E(Y|X_~i)]/V(Y)：参数 i 及其所有交互的贡献。
若 STi ≫ Si，说明存在强交互。Sobol 精确但昂贵（Saltelli 采样需 N·(2k+2) 次运行）。

**何时用哪个**：先 Morris 筛掉不重要参数（几百次运行），再对保留的 5–10 个参数做 Sobol（数千~数万次运行）量化。

### SALib 实操（以 RothC 或作物模型为例）
SALib（Herman & Usher 2017, *Journal of Open Source Software* 2(9):97）的典型流程是 sample → 运行模型 → analyze 三步：
```python
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
import numpy as np

problem = {
    'num_vars': 4,
    'names': ['k_DPM','k_RPM','k_HUM','clay'],
    'bounds': [[8,12],[0.2,0.4],[0.01,0.03],[10,40]]
}
X = sobol_sample.sample(problem, 1024)      # N*(2D+2) 行
Y = np.array([run_rothc_to_SOC(x) for x in X])   # 你的模型封装
Si = sobol_analyze.analyze(problem, Y, print_to_console=True)
# Si 含 'S1','S1_conf','ST','ST_conf'(95% bootstrap 置信区间)
```
Morris 版把 `from SALib.sample import morris` / `from SALib.analyze import morris`，返回 `mu_star`、`sigma`。解读：μ\* 排序取前几名作为率定参数集；若 ST 显著大于 S1，说明参数交互主导，需保留这些参数共同率定。

## 2.2 频率派/梯度率定：PEST 与 PEST++

PEST（Doherty）是**模型无关（model-independent）**率定器，核心是 **Gauss-Marquardt-Levenberg (GML)** 算法：
```
(JᵀQJ + λ·I)·Δp = JᵀQ·(observed − simulated)
```
J 为 Jacobian（用有限差分逐参数扰动数值估计），Q 为观测权重矩阵，λ 为 Marquardt lambda（在最速下降与 Gauss-Newton 之间插值）。

**模型无关的实现机制**：
- **模板文件 (.tpl)**：模型输入文件的"挖空"版，PEST 把参数值写入。
- **指令文件 (.ins)**：告诉 PEST 如何从模型输出文件中读取模拟值。
- **控制文件 (.pst)**：参数、观测、权重、算法设置。
这样 PEST 无需改动模型源码即可率定任何模型（DSSAT、APSIM、DNDC 皆可）。

**Tikhonov 正则化**：对欠定问题（参数多于可辨识信息）加"软知识"约束（如倾向于先验值或平滑），保持数值稳定。其哲学被 Doherty 总结为"如有疑问，就把它包括进来（if in doubt, include it）"。PEST++（Welter, White, Hunt & Doherty 2015, USGS Techniques and Methods 7-C12）加入 SVD-Assist（奇异值分解降维）、并行运行管理、线性（FOSM）不确定性分析，PESTPP-GLM/PESTPP-IES 是现代实现，实现 subspace Gauss-Levenberg-Marquardt 算法。

**应用到 DSSAT/APSIM**：DSSAT 内置 **GLUE**（Generalized Likelihood Uncertainty Estimation，伪贝叶斯，率定品种遗传系数 .CUL/.ECO，通常 phenology 与 growth 各需约 6000 次运行）与 **GenCalc/GENCALC**（Genotype Coefficient Calculator，Hunt et al. 1993，逐参数迭代细调）；新版 **GLUEP**（Ferreira et al. 2024, *Computers and Electronics in Agriculture*）支持并行计算并可率定 .ECO 生态型参数。R 包 **CroptimizR** 提供 Nelder-Mead、DREAM 等多种算法接口给 DSSAT/APSIM/STICS。

## 2.3 贝叶斯率定 + MCMC

**贝叶斯定理**：`p(θ|D) ∝ p(D|θ)·p(θ)`。后验 ∝ 似然 × 先验。
- **似然（高斯误差模型）**：`p(D|θ) = Π N(y_obs,i | model_i(θ), σ²)`，对数似然 `ℓ = −0.5·Σ[(y_i−m_i)²/σ² + ln(2πσ²)]`。
- **先验**：由生物学约束设定（如 RUE>0、k 在文献范围）。
- **后验**：通常高维、无解析解，需 MCMC 采样。

**Metropolis-Hastings**：从提议分布抽 θ'，以 `α = min(1, p(θ'|D)/p(θ|D))` 接受。简单但在相关/多峰后验中混合慢。

**现代采样器**：
- **DREAM / DREAM(ZS)**（Vrugt et al. 2009；Laloy & Vrugt 2012, "High-dimensional posterior exploration of hydrologic models using multiple-try DREAM(ZS) and high-performance computing," *Water Resources Research* 48:W01526）：differential evolution adaptive Metropolis，多链并行，从历史状态档案 + snooker 更新中生成提议，维持 detailed balance，对多峰、高维问题表现优异，是水文/生物地球化学率定的主力。Python 实现 PyDREAM（Shockley et al. 2018, *Bioinformatics* 34:695）。
- **emcee**（Foreman-Mackey et al.；Goodman & Weare 仿射不变 ensemble sampler）：用一群 walker 互相生成提议，对参数尺度不敏感，易用。
- **NUTS/HMC**（PyMC、Stan）：用梯度（Hamiltonian Monte Carlo + No-U-Turn Sampler）高效探索连续后验，但需模型可微或用数值梯度。

**把模型封装为黑盒似然（emcee 示例）**：
```python
import emcee, numpy as np

def log_prior(theta):
    k_dpm, rue, sigma = theta
    if 8<k_dpm<12 and 1.0<rue<4.0 and sigma>0:
        return 0.0
    return -np.inf

def log_like(theta, y_obs):
    *params, sigma = theta
    y_sim = run_model(params)
    return -0.5*np.sum((y_obs - y_sim)**2/sigma**2 + np.log(2*np.pi*sigma**2))

def log_post(theta, y_obs):
    lp = log_prior(theta)
    return lp + log_like(theta, y_obs) if np.isfinite(lp) else -np.inf

ndim, nwalkers = 3, 32
p0 = initial + 1e-3*np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_post, args=(y_obs,))
sampler.run_mcmc(p0, 5000, progress=True)
chain = sampler.get_chain(discard=1000, flat=True)
# 后验中位数与 95% 可信区间
print(np.percentile(chain, [2.5,50,97.5], axis=0))
```
**后验预测检验（posterior predictive check）**：从后验抽样跑模型，看观测是否落在预测分布内。

**等效性（equifinality）与 GLUE**：Beven & Binley 的 GLUE 方法承认**多组参数可同样好地拟合数据**——存在"行为参数集（behavioral parameter sets）"的集合而非唯一最优。因此后验分布（而非单点）才是诚实的表达。参数**可辨识性/相关性**可通过后验的协方差/相关矩阵诊断：强相关的参数对（如 RUE 与 k）无法被数据单独约束。Vrugt et al. 2009 还专门讨论了形式（DREAM）与非形式（GLUE）贝叶斯方法在 equifinality 上的关系。

## 2.4 目标函数与指标
- **RMSE** = √(Σ(sim−obs)²/n)：量纲同变量，对大值敏感。
- **NSE**（Nash-Sutcliffe Efficiency, Nash & Sutcliffe 1970）= 1 − Σ(obs−sim)²/Σ(obs−mean_obs)²。NSE=1 完美，NSE=0 等于用均值预测，<0 比均值还差。
- **KGE**（Kling-Gupta Efficiency, Gupta et al. 2009）= 1 − √[(r−1)² + (β−1)² + (γ−1)²]，其中 r 为相关系数，β=μ_sim/μ_obs（偏差比），γ=CV_sim/CV_obs（变率比，CV=标准差/均值）。KGE 把相关、偏差、变率三者解耦，比 NSE 更均衡（NSE 对高值不成比例敏感）。注意：用均值预测时 KGE = 1−√2 ≈ −0.41（而非 0），所以 KGE>−0.41 才算优于均值基准（Knoben et al. 2019, *HESS* 23:4323）；NSE 与 KGE 不可直接比较，其关系依赖观测的变异系数。
- **R²**、**一致性指数 d**（index of agreement, Willmott）。
- **多目标率定**：同时拟合产量 + N₂O + 土壤水 + SOC。各目标常冲突，解为 **Pareto 前沿（Pareto front）**——无法在不牺牲一个目标的情况下改善另一个。**NSGA-II**（非支配排序遗传算法）是求 Pareto 前沿的标准进化算法。

---

# 区域三：机理 × 机器学习混合建模（现代前沿）

## 3.1 物理信息 / 知识引导机器学习（PINN, PGML/KGML）

**PINN（Raissi, Perdikaris & Karniadakis 2019, "Physics-informed neural networks," *Journal of Computational Physics* 378:686–707, DOI:10.1016/j.jcp.2018.10.045）**：把微分方程残差放进神经网络损失函数。对 PDE/ODE `N[u]=0`，损失 = 数据拟合损失 + λ·物理残差损失（在配点上用自动微分算 ∂u/∂t、∂u/∂x 并代入方程）。既能"求解"正问题，也能"反演"未知参数（inverse problem）。综述见 Karniadakis et al. 2021, "Physics-informed machine learning," *Nature Reviews Physics* 3:422–440。

**知识引导机器学习（KGML / Theory-Guided Data Science，Karpatne, Atluri, Faghmous, Steinbach, Banerjee, Ganguly, Shekhar, Samatova & Kumar 2017, "Theory-guided data science," *IEEE Transactions on Knowledge and Data Engineering* 29(10):2318–2331；专著 Karpatne, Kannan & Kumar 2022 *Knowledge Guided Machine Learning: Accelerating Discovery Using Scientific Knowledge and Data*, CRC Press）**：更广义地把科学知识注入 ML——可通过损失函数惩罚违反物理（如质量/能量守恒）的预测、用物理引导网络结构/初始化、或物理引导的预训练。Karpatne, Watkins, Read & Kumar 2017 的湖泊温度建模 PGNN（Physics-Guided Neural Network, arXiv:1710.11431）是经典案例。

## 3.2 ML 仿真器 / 代理模型（emulator / surrogate）

用快速 NN 或高斯过程（Gaussian Process, GP）仿真慢过程模型，使昂贵的敏感性分析、率定、大尺度集合运行成为可能。**GP 仿真器**天然给出预测的均值 + 方差，特别适合不确定性量化（UQ）。工作流：(1) 用过程模型生成训练样本（Latin Hypercube 采样参数空间）；(2) 训练 GP/NN 仿真器；(3) 用仿真器做数千万次 Sobol/MCMC。这是把 DNDC/APSIM 嵌入大尺度优化的主流手段。

## 3.3 混合/残差建模与可微分建模

- **残差建模（residual modeling）**：ML 学习过程模型的误差 `error = obs − model`，预测时 `corrected = model + ML(features)`。保留机理可解释性，ML 只补结构缺陷。
- **参数化替换**：用 NN 替换某个不确定子过程（如经验 `f(W)` 水分响应函数），保留 ODE 骨架。
- **可微分建模（differentiable modeling，Shen et al. 2023, *Nature Reviews Earth & Environment*）**：把整个过程模型用支持自动微分的框架（PyTorch/JAX）重写，使输出对参数的梯度可由 AD/伴随法（adjoint）高效计算，从而端到端用梯度下降训练嵌入的 NN。Feng et al. 2022（*Water Resources Research* 58(10):e2022WR032404；arXiv:2203.14827）用可微分 HBV 水文模型（δ models）作为骨架、嵌入神经网络参数化，原文报告"δ models can obtain a median Nash-Sutcliffe efficiency of 0.732 for 671 basins across the USA for the Daymet forcing data set, compared to 0.748 from a state-of-the-art LSTM model with the same setup"——即在 671 个美国流域中位 NSE 达 0.732，逼近纯 LSTM 的 0.748（NLDAS 驱动下为 0.715 vs 0.722），同时保留物理可解释性。**Neural ODE**（把 ODE 右端用 NN 表示）是其重要形式（Höge et al. 2022 在流量预测中达到 DL 性能同时保留概念模型可解释性）。

## 3.4 农业/生物地球化学应用与原型路线

**旗舰案例 KGML-ag（Liu, Xu, Tang, Guan, Griffis, Erickson, Frie, Jia, Kim, Miller, Peng, Wu, Yang, Zhou, Kumar & Jin 2022, "KGML-ag: a modeling framework of knowledge-guided machine learning to simulate agroecosystems: a case study of estimating N₂O emission using data from mesocosm experiments," *Geoscientific Model Development* 15(7):2839–2858, DOI:10.5194/gmd-15-2839-2022）**：用知识引导 ML 估计农田 N₂O 排放。以 GRU（门控循环单元）为骨架，关键创新：(1) 用中间变量（IMV：CO₂通量、土壤 NO3⁻、NH4⁺、体积含水量）的初值而非时间序列作输入以降低数据需求；(2) 分层结构显式估计 IMV；(3) 多任务学习；(4) 用先进过程模型 **ecosys**（Grant et al. 2001）生成的数百万合成数据预训练，再用 mesocosm 观测微调。原文报告："KGML-ag did an excellent job in reproducing the mesocosm N2O fluxes (overall r2=0.81, and RMSE=3.6 mgNm-2d-1 from cross validation)"，并稳定优于纯过程模型和纯 ML 模型，尤其在复杂时间动态与排放峰值上。

**农学研究生的原型路线（推荐由易到难）**：
1. 先建**仿真器**（GP/NN 仿真你的 RothC/作物模型）——回报最高、风险最低。
2. 再做**残差建模**——拿过程模型预测 + 现场观测，训练 ML 学残差。
3. 最后尝试**参数化替换/可微分**——把一个 `f(T)` 或 `f(W)` 换成小 NN，用 PyTorch 重写该模块端到端训练。

**权衡（tradeoffs）**：
- **可解释性**：残差/参数化混合 > 仿真器 > 纯 ML。
- **外推（extrapolation）**：保留机理骨架的混合模型在分布外更稳健；纯 ML 易"out-of-sample failure"。
- **数据需求**：纯 ML 最贪婪；KGML 用过程模型合成数据预训练可大幅降低对观测的需求。

---

# 区域四：临界点 / 状态转换 / 复杂系统动力学（深度）

## 4.1 动力系统基础

考虑 ODE 系统 `dx/dt = f(x; p)`。
- **不动点/平衡（fixed point/equilibrium）**：f(x\*)=0。
- **线性稳定性（linear stability）**：在 x\* 处线性化，Jacobian J = ∂f/∂x。**所有特征值实部 < 0 → 稳定（吸引子）**；任一实部 > 0 → 不稳定。一维情形 `dx/dt=f(x)`，稳定 ⟺ f'(x\*)<0。
- **吸引子（attractor）、吸引域（basin of attraction）、相图（phase portrait）**：多稳态系统中，初始条件落在哪个吸引域决定终态——这正是 equifinality 与 spin-up 初值敏感性的几何根源。

应用到土壤碳：线性库模型 `dC/dt = I − ξkC` 有唯一稳定平衡 `C\*=I/(ξk)`；但**微生物显式模型**（含 Michaelis-Menten 的 SOC-DOC-MBC-酶模型）可出现多个平衡点甚至不稳定平衡与振荡。Georgiou 等与 *Biogeosciences* 2024（21:3441，"When and why microbial-explicit soil organic carbon models can be unstable"）系统分析了此类模型何时出现"存在解析解但不稳定、瞬态模拟无法到达"的平衡点，并指出参数相关性可缓解不稳定平衡的发生。

## 4.2 分岔理论（Bifurcation Theory）

随参数 p 变化，平衡的数量/稳定性发生定性改变：
- **鞍结（saddle-node / fold）分岔**——**临界点的典范**。正规形 `dx/dt = r − x²`（或 `r + x²`）。r 越过 0 时，一对稳定+不稳定平衡碰撞湮灭，系统被迫跳到远处的另一吸引子。这是"灾变式（catastrophic）"突变的核心机制（codimension-1，只需调一个参数）。判据：Jacobian 有一个简单零特征值，且 ∂²f/∂x²≠0。
- **跨临界（transcritical）**：`dx/dt = rx − x²`，两平衡交换稳定性（如植被 0 态与非 0 态的交换；Amazon 森林死亡模型即被刻画为接近跨临界分岔）。
- **叉式（pitchfork）**：`dx/dt = rx − x³`（超临界）或 `rx + x³`（亚临界），对称系统中一个平衡分裂为三。
- **Hopf 分岔**：一对共轭复特征值实部越过 0，平衡失稳并诞生极限环（limit cycle）——**振荡的起源**（微生物显式 SOC 模型的振荡即可由此解释）。

**迟滞（hysteresis）与多稳态（alternative stable states）**：折叠型（fold）响应曲线意味着系统状态依赖历史。**Scheffer 浅湖富营养化模型（Scheffer 1990, "Alternative stable states in eutrophic, shallow freshwater systems: a minimal model," *Hydrobiologia*；Scheffer et al. 1993, *Trends in Ecology & Evolution* 8:275）**是经典：清水（沉水植物主导）与浊水（浮游植物主导）两个稳态在一定营养盐范围内共存；营养盐升高使清水态消失（系统突跳到浊水），但降低营养盐到同一水平**不能**恢复——必须降到更低的阈值，形成迟滞回线。这就是**折叠灾变（fold catastrophe）**。同类还有草地-荒漠化、森林-草原模型。

## 4.3 临界慢化与早期预警信号（Critical Slowing Down & EWS）

**理论**：接近 fold 分岔时，主导特征值（dominant eigenvalue）趋近 0，系统受扰后恢复变慢（return rate→0），即**临界慢化（critical slowing down）**。Wissel (1984, *Oecologia* 65:101) 给出特征恢复时间在阈值附近发散的普适律。鞍结正规形带噪声的随机微分方程展开可解析推出这些指标的发散行为。

**统计指标（可从时间序列计算）**：
- **滞后-1 自相关（lag-1 autocorrelation）上升**：系统"记忆"变长。
- **方差/标准差上升**（Carpenter & Brock 2006, *Ecology Letters* 9:311："rising variance" 是先兆）。
- **偏度（skewness）变化**（Guttal & Jayaprakash 2008, *Ecology Letters* 11:450）。
- **闪烁（flickering）**：在双稳态间来回跳。
- 空间指标：空间方差、空间自相关上升（Dakos et al. 2010）。

**权威文献**：Scheffer, Bascompte, Brock, Brovkin, Carpenter, Dakos, Held, van Nes, Rietkerk & Sugihara 2009, "Early-warning signals for critical transitions," *Nature* 461(7260):53–59, DOI:10.1038/nature08227（综述）；Dakos et al. 2008, *PNAS* 105:14308（8 个古气候突变前自相关上升）；Dakos et al. 2012, *PLoS ONE* 7:e41010（方法学，配套 R 包 **earlywarnings** 的 `generic_ews` 函数）。理论参考教材 Strogatz *Nonlinear Dynamics and Chaos*（该书亦被 Scheffer 2009 在临界慢化推导中引用）。

**计算示例（概念）**：对去趋势后的状态变量时间序列，用滑动窗口计算滚动方差与 lag-1 AR(1) 系数，若两者随时间显著上升（Kendall τ 检验），即为接近临界点的预警。

```python
import numpy as np
from statsmodels.tsa.stattools import acf

def rolling_ews(x, win=50):
    var, ar1 = [], []
    for i in range(len(x)-win):
        seg = x[i:i+win]
        seg = seg - np.polyval(np.polyfit(np.arange(win), seg, 1), np.arange(win)) # 去趋势
        var.append(np.var(seg))
        ar1.append(acf(seg, nlags=1, fft=False)[1])
    return np.array(var), np.array(ar1)
```

## 4.4 与农业生态模型的关联

**DNDC 中的阈值/开关非线性**：Eh<500 mV 触发反硝化、Eh<−150 mV 触发产甲烷——这是分段切换，数学上接近鞍结/不连续分岔。"气球"的诞生-膨胀-爆裂是离散事件序列：当一个氧化剂耗尽，主导电对切换，Eh 阶跃，下一组反应突然启动。淹水后 N₂O/CH₄ 通量的脉冲式爆发正是这种阈值动力学的涌现，**不是数值 bug**。判别方法：(1) 检查是否对应物理阈值（Eh、底物耗尽）；(2) 做步长收敛测试（区域一(b)）——若突变在 Δt 减半后位置/幅度稳定，则是真实阈值而非数值振荡；(3) 检查质量守恒（区域一(a)）。

**土壤碳模型中的多稳态**：土地利用/管理变化下，SOC 可能存在多个稳定平衡（耕作 vs 免耕、草地 vs 农田）；微生物显式模型有碳饱和（C saturation）与替代稳态。实证上，Clayton, Lemanski & Bonkowski 2021（"Shifts in soil microbial stoichiometry and metabolic quotient provide evidence for a critical tipping point at 1% soil organic carbon in an agricultural post-mining chronosequence," *Biology and Fertility of Soils* 57:435–446）在一条 52 年褐煤矿区复垦（初始 0.2% SOC）的"空间换时间"序列中报告："Our data revealed sudden shifts in microbial stoichiometry and metabolic quotient with increasing SOC at a critical value of 1% SOC"——即在 **1% SOC** 处微生物化学计量与代谢商发生突变（临界点）。线性库模型（RothC/CENTURY）本身只有唯一平衡，**不会**自发产生多稳态——若你的线性模型出现多稳态，那确实是 bug；但若引入非线性微生物动力学或 CUE（碳利用效率）反馈，多稳态是真实的。

**作物物候的 GDD 阈值切换**：发育阶段在积温越过阈值时切换（如开花、灌浆），是离散事件，会使产量对播期/温度的响应出现"悬崖"——这是阈值非线性，可用早期预警视角理解但通常是确定性切换。

**与 equifinality 和乘性限制因子的连接**：`f(T)·f(W)·f(substrate)` 的**乘性**结构与 Michaelis-Menten 的**饱和**动力学是非线性的根本来源——乘积使任一因子趋零即整体趋零（开关效应），饱和使响应在高底物处变平（产生折叠曲线所需的非单调性）。多稳态意味着多个吸引域，不同初值（spin-up）或不同参数集（equifinality）可能收敛到不同终态——研究者必须区分"模型有多个真实吸引子"与"率定问题欠定"这两种不同的非唯一性。

---

## Recommendations（建议）

**阶段 0（1–2 周）——打地基**：先把 STAGE 1 的 RothC Python 代码跑通并通过官方单元测试。基准：复现 0.1533→0.1140 等四个数到 1e-4。这验证你理解了解析衰减、abc 因子、黏粒分配。

**阶段 1（2–4 周）——建可调试的管线**：实现 STAGE 2 函数库 + STAGE 3 最小作物模型，并把区域一的五个数值检查（质量守恒、步长收敛、负库守护、CFL、spin-up 漂移）做成每次运行自动跑的断言。基准：质量守恒误差 <1e-9，步长减半解变化 <0.1%。

**阶段 2（1–2 月）——率定**：先用 SALib Morris 筛选（几百次运行）→ 取 μ\* 前 5–10 个参数 → Sobol 量化（数千次运行）→ 用 emcee 或 PyMC 做贝叶斯率定得到后验。基准：用 NSE 与 KGE 双指标评估，报告参数后验的 95% 可信区间而非单点；若两参数后验相关系数 |ρ|>0.9，说明不可辨识，需固定其一或增加观测类型。切换阈值：若 KGE<−0.41，模型还不如用均值，重查结构。

**阶段 3（2–3 月）——混合 ML**：从 GP 仿真器起步（最高回报）。当仿真器 R²>0.95 时，再尝试残差建模；只有当你需要端到端学习某子过程时才投入可微分重写（PyTorch/JAX）。参考 KGML-ag 用过程模型合成数据预训练以降低观测需求。

**阶段 4（持续）——临界点分析**：对 DNDC 的 Eh/通量时间序列、SOC 长期轨迹用 `earlywarnings`/`generic_ews` 计算滚动方差与 AR(1)。判别突变真伪的三步法：物理阈值对应 → 步长收敛稳定 → 质量守恒——三者皆满足才是真实 regime shift。改变建议的基准：若突变随 Δt 变化而漂移，停止解释为临界点，先修数值格式（改隐式或减步长）。

---

## Caveats（注意事项与不确定性）

1. **DNDC Nernst 符号约定**：v9.5 官方印刷形式为 `Eh = E0 + (RT/nF)·ln([ox]/[red])`，与许多教材/先前报告的 `−(RT/nF)·ln([red]/[ox])` 符号相反、比值倒置，二者代数等价但实现时必须与所参照版本一致。

2. **DNDC 厌氧气球的显式 anvf 公式**（`anvf = a·(1−(b−pO2/pO2_air))`）及氧扩散方程**未印在 v9.5 Scientific Basis PDF**（该文档在此点为概念性描述），来自代码级复现（归于 Li et al. 2000），系数 a、b 无公开数值。生产用途须对照实际 DNDC 源码核实。Li (2007) 与 Gilhespy et al. (2014) 全文获取受限，部分引文来自检索片段；v9.5 PDF 与代码级复现中的半饱和常数、产率、阈值 mV 大体一致，但 DNDC 8.5/9.3/9.5 版本间可能有细微差异。

3. **RothC 算例测试**只验证纯衰减项；完整模型还需加 (BIO+HUM) 回流、碳输入、IOM。官方 Excel/Fortran 版（Rothamsted Research 发布，Zenodo v2.1.1）是最终校验基准。

4. **CFL 稳定性界 `Δt ≤ Δz²/(2D)`** 是线性扩散显式格式的经典结果；tipping-bucket 因不解扩散方程而不受此约束，但其代数级联在极端入渗下也可能产生非物理的瞬时跳变，需检查。

5. **早期预警信号并非万无一失**：Scheffer et al. 自己指出，在快速/非线性强迫下（如亚马逊森林死亡的 HadCM3 集合实验）自相关趋势可能不出现；方差上升也可能仅源于强迫本身的变率增加，而非临界慢化。EWS 是概率性提示，非确定性预测。

6. **线性 vs 非线性多稳态**：RothC/CENTURY 的线性库结构只有唯一平衡，不会自发多稳态；观察到的"多稳态"若来自线性模型则是 bug，若来自非线性微生物/CUE 反馈或外部驱动（管理、土地利用）切换则可能真实。务必先判断模型类别。

7. **KGE 基准值**：用均值预测得 KGE≈−0.41 而非 0（Knoben et al. 2019, *HESS* 23:4323），NSE 与 KGE 不可直接比较（其关系依赖观测变异系数）。

8. **DSSAT GLUE 的局限**：GLUE 不估计 PHINT、不支持固定部分系数（GenCalc 可），且单季只能用一个观测优化目标变量；多次季内观测的时间序列率定需借助 CroptimizR 等外部工具。
