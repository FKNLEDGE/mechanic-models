# Round 5｜动力系统 + 临界转变（中级）：用动力系统语言把"平衡、反馈、阈值、临界点"讲透

> 框架⑤（动力系统 / 稳定性 / 临界点）为本轮唯一核心。三个经典模型：Lotka–Volterra 捕食–被捕食、SIR 传染病、Scheffer 浅湖多稳态。

## TL;DR（三句话回答核心问题）
- **本轮把前四轮反复出现的"平衡、反馈、阈值、翻转"用严格数学语言重写**：固定点 f(x*)=0、用 Jacobian 矩阵特征值判稳定性、相图、分岔（saddle-node/fold、transcritical、pitchfork、Hopf）、滞后（hysteresis）、多稳态、临界慢化与早期预警信号（EWS）。
- **三个经典模型分别示范三种核心现象**：Lotka–Volterra 给出"中心/极限环"的振荡（boom-bust 循环），SIR 给出"阈值/transcritical 分岔"（R0=1 翻转），Scheffer 浅湖给出"fold 分岔/滞后/多稳态"（清水↔浊水不可逆翻转）。
- **农业生态系统是多稳态的、可能突然且不可逆地翻转**——害虫周期、病害暴发、湖泊富营养化、草地荒漠化都是同一套数学；预防远比恢复便宜，EWS 可提前预警但并非万无一失。

---

## Key Findings（核心要点）
1. **稳定性的唯一判据**：一维看 f'(x*) 的正负；多维看 Jacobian 矩阵 J 的特征值实部——全负则稳定，有正则不稳定，复数则振荡。这是"小球在山谷里"的严格版本，特征值就是"小球滚回谷底有多快"。
2. **fold（saddle-node）分岔是"临界点/tipping"的数学心脏**：两个平衡点（一稳一不稳）相撞并消失，系统被迫跳到远处的另一个状态。S 形（折叠）平衡曲线 → 两个 fold → 滞后环。
3. **R0 = β/γ 是传染病阈值**：R0>1 暴发、R0<1 熄灭，本质是无病平衡点在 R0=1 处失稳的 transcritical 型分岔；群体免疫阈值 = 1 − 1/R0。
4. **临界慢化（critical slowing down）**：逼近 fold 时主导特征值趋于 0，系统从扰动恢复越来越慢 → lag-1 自相关上升、方差上升、偏度变化、闪烁（flickering）。这是 EWS 的理论根基（Scheffer et al. 2009）。
5. **统一视角**：同一套工具箱（固定点、Jacobian/特征值、分岔、滞后、EWS）解释害虫循环、病害暴发、生态系统崩溃，并支撑第 2–4 轮过程模型中的阈值行为（DNDC 氧化还原开关、微生物显式土壤碳多稳态、regime shift）。

---

## Section 0：从第 1–4 轮架桥——给老朋友起个正式的名字

亲爱的小白同学，如果你跟着走完了前四轮，你其实已经无数次撞见今天的主角了，只是当时还没人告诉你它的学名。

- **第 1 轮**：Logistic 增长 dN/dt = rN(1 − N/K)。当种群长到承载力 K 时就稳定下来，不再变化——这个"稳定下来的 K"就是一个**稳定固定点**。我们当时用"小球滚到山谷底"的比喻，今天要把这个比喻变成可以算的数学。
- **第 2 轮**：DNDC 里土壤淹水后氧化还原电位（redox）跨过某个阈值，反硝化突然启动、N2O 出现脉冲式释放——这是一个**阈值切换**。系统不是平滑过渡，而是"啪"地翻到另一种行为。
- **第 4 轮**：tipping-bucket（倾倒水桶）土壤水分模型、land-surface 里的 regime 概念——都是"系统突然换挡"。

这些现象（平衡、反馈、阈值、翻转）有一个共同的数学母语，叫做**动力系统理论（dynamical systems theory）**。本轮就是给这门母语开蒙。

**核心工具箱预览**（后面每一项都会从第一性原理讲）：
1. 一个系统写成 dx/dt = f(x)（x 可以是多维向量）；
2. **固定点 / 平衡点**：让 f(x*) = 0 的那些 x*；
3. **稳定性**靠**线性化**判断——把 f 在 x* 附近近似成直线/线性映射，看 **Jacobian 矩阵**特征值实部的符号；
4. **相图（phase portrait）**：在状态平面上画轨迹；
5. **分岔（bifurcation）**：参数跨过某个临界值时，系统行为发生质变；
6. **滞后（hysteresis）与多稳态（alternative stable states）**；
7. **临界慢化与早期预警信号（EWS）**：在系统真正翻车之前从数据里嗅到危险。

**统一信息（请记住这一句）**：很多农业生态系统的行为不是渐变的——它们可能突然、甚至不可逆地翻转；动力系统理论就是我们理解并提前预判这种翻转的工具。

---

## Section 1：动力系统工具箱（第一性原理）

### 1(a) 状态空间与"流"：dx/dt = f(x)

把 x 想成"系统现在的状态"（比如某种群密度），f(x) 告诉你"在这个状态下，x 下一刻往哪个方向变、变多快"。

- 若 f(x) > 0，x 在增大（箭头朝右）；
- 若 f(x) < 0，x 在减小（箭头朝左）；
- 若 f(x) = 0，x 不动——这就是固定点。

**一维"线上的流"（flow on a line）图像**：在一条数轴上，每个点画一个小箭头（f>0 朝右、f<0 朝左）。箭头从两侧"指向"的点是稳定的（像谷底），箭头从两侧"背离"的点是不稳定的（像山顶）。这是 Strogatz《Nonlinear Dynamics and Chaos》一书第 2 章的标志性教学图像，也是整个动力系统直觉的起点。

```python
import numpy as np
import matplotlib.pyplot as plt

# 例子：Logistic 增长 f(N) = r N (1 - N/K)
r, K = 1.0, 10.0
def f(N): return r * N * (1 - N/K)

N = np.linspace(-2, 14, 400)
plt.figure(figsize=(8,3))
plt.axhline(0, color='gray', lw=0.8)
plt.plot(N, f(N), 'b', label="f(N)=rN(1-N/K)")

# 在数轴上画"流"的箭头：f>0朝右, f<0朝左
for N0 in np.linspace(-1, 13, 25):
    v = f(N0)
    plt.arrow(N0, 0, 0.0008*np.sign(v), 0, head_width=0.25,
              head_length=0.4, fc='r', ec='r')
# 固定点
for Nstar in [0, K]:
    stable = (f(Nstar+1e-3) < 0)  # 右侧朝左→稳定
    plt.plot(Nstar, 0, 'o', ms=11,
             mfc=('k' if stable else 'white'), mec='k')
plt.xlabel("N（种群密度）"); plt.ylabel("dN/dt")
plt.title("Flow on a line：实心=稳定固定点(K)，空心=不稳定(0)")
plt.legend(); plt.tight_layout(); plt.show()
```

### 1(b) 固定点 / 平衡点：f(x*) = 0

固定点就是"系统不再变化"的状态。**小球–山谷比喻**：把 f(x) 想成地形的坡度，固定点是坡度为 0 的地方——可能是谷底（稳定），也可能是山顶（不稳定）。

### 1(c) 稳定性 = 线性化（这是本轮最重要的一节）

**一维情形**：在固定点 x* 附近，给一个小扰动 η = x − x*。泰勒展开得 dη/dt ≈ f'(x*)·η。这是一条直线，解是 η(t) = η(0)·e^{f'(x*)·t}。于是：
- **f'(x*) < 0**：扰动指数衰减 → **稳定**（谷底）；
- **f'(x*) > 0**：扰动指数放大 → **不稳定**（山顶）；
- **f'(x*) = 0**：线性判据失效（往往正是分岔点）。

**已解算的一维例子（Logistic）**：f(N) = rN(1 − N/K)，f'(N) = r − 2rN/K。
- 在 N* = 0：f'(0) = r > 0 → **不稳定**（空地总会被占领）；
- 在 N* = K：f'(K) = r − 2r = −r < 0 → **稳定**（种群稳定在承载力）。

这正是第 1 轮"K 是稳定的"的严格证明！

**多维情形——Jacobian 矩阵**：当状态是向量 x = (x₁, x₂, …)，f 也是向量函数。线性化得到的不是一个数 f'，而是一个偏导数矩阵——**Jacobian 矩阵 J**：

J = [ ∂f₁/∂x₁  ∂f₁/∂x₂ ; ∂f₂/∂x₁  ∂f₂/∂x₂ ]（二维为例）

稳定性由 **J 的特征值（eigenvalues）λ** 决定：
- **所有特征值实部 < 0** → 稳定（任何小扰动都会回落）；
- **任一特征值实部 > 0** → 不稳定；
- **特征值为复数 a ± bi** → 扰动会**振荡**（虚部 b 是振荡频率，实部 a 决定振荡是衰减 a<0 还是放大 a>0）。

**给小白的特征值直觉**：特征值就是"小球被推一下之后回到谷底（或冲出谷顶）的速率与方式"。实部是"快慢与方向"（负=回得来、回多快；正=跑得掉、跑多快）；虚部是"回来的路上要不要绕圈/打转（振荡）"。**当实部趋近 0，小球回谷底越来越慢——这就是后面要讲的"临界慢化"。**

二维系统还有两个超好用的快捷判据（对 2×2 Jacobian）：迹 tr(J) = λ₁+λ₂，行列式 det(J) = λ₁λ₂。**稳定的充要条件是 tr(J) < 0 且 det(J) > 0。**

### 1(d) 相图（phase portrait，二维系统）

二维系统在 (x, y) 平面上画轨迹。关键元素：
- **零线（nullclines）**：dx/dt=0 的曲线和 dy/dt=0 的曲线；它们的交点就是固定点。
- **固定点类型**（由特征值决定）：
  - **节点（node）**：特征值都是实数、同号（都负=稳定结点，都正=不稳定结点）；
  - **螺旋/焦点（spiral/focus）**：特征值是复数（实部负=稳定螺旋=阻尼振荡；实部正=不稳定螺旋）；
  - **鞍点（saddle）**：特征值实数异号（一进一出，总是不稳定）；
  - **中心（center）**：特征值是纯虚数（±bi），轨迹是闭合圈——"中性稳定"，是个刀锋上的特例（Lotka–Volterra 经典模型就是这种）。

### 1(e) 分岔（bifurcation）：参数跨过阈值时行为发生质变

分岔是本轮的灵魂。下面四个"正规型（normal form）"是 Strogatz 教材里的标准范式，每个配一句直觉（μ 是分岔参数，临界值 μc = 0）：

| 分岔 | 正规型 | 一句话直觉 |
|---|---|---|
| **saddle-node / fold** | dx/dt = μ − x² | μ>0 时有两个平衡点（一稳一不稳），μ 减到 0 时两者相撞并**同时消失**——平衡点"凭空出现/凭空消失"。**这就是 tipping 的数学心脏。** |
| **transcritical** | dx/dt = μx − x² | 两个平衡点穿过彼此并**交换稳定性**——一个本来稳定的变不稳定，另一个反之。SIR 的 R0=1 阈值就是这一类。 |
| **pitchfork** | dx/dt = μx − x³ | μ 跨过 0，一个平衡点分裂成三个（一不稳+两稳，超临界）——"一变三"，对称破缺。 |
| **Hopf** | dr/dt = μr − r³, dθ/dt = 1 | μ 跨过 0，一个稳定焦点失稳并**诞生一个极限环（limit cycle）**——系统开始持续振荡。捕食–被捕食的"富营养化悖论"就是它。 |

**特别强调 fold/saddle-node**：在 saddle-node、transcritical、pitchfork 中，是一个**实特征值穿过 0**；在 Hopf 中，是一对**共轭复特征值的实部穿过 0**（越过虚轴）。fold 是最"稳健"的分岔（加扰动不会消失），也是真实世界"突然翻车"最常见的机制。

```python
# fold（saddle-node）一维分岔图：dx/dt = mu - x^2
import numpy as np, matplotlib.pyplot as plt
mu = np.linspace(-3, 3, 600)
xs = np.sqrt(mu[mu>=0])     # 稳定分支 x*=+sqrt(mu)
xu = -np.sqrt(mu[mu>=0])    # 不稳定分支 x*=-sqrt(mu)
plt.figure(figsize=(6,4))
plt.plot(mu[mu>=0], xs, 'b', label='稳定分支 x*=+√μ')
plt.plot(mu[mu>=0], xu, 'r--', label='不稳定分支 x*=-√μ')
plt.plot(0,0,'ko'); plt.axvline(0, color='gray', lw=0.6)
plt.xlabel('μ（参数）'); plt.ylabel('固定点 x*')
plt.title('Saddle-node/fold 分岔：μ<0 无平衡点，μ>0 出现两个')
plt.legend(); plt.tight_layout(); plt.show()
```

### 1(f) 滞后（hysteresis）与多稳态

当平衡曲线呈 **S 形（折叠）** 时，会出现两个稳定分支（上支、下支）被中间一段不稳定分支隔开。随参数缓慢增大，系统沿下支走，到达**上 fold 点**才被迫跳到上支；要想退回去，参数必须减小到**另一个更低的下 fold 点**才会跳回——**前向阈值与后向阈值不相等，形成滞后环（hysteresis loop）**。

**滞后比喻**：翻船容易、扶正难；把湖搅浑容易、再让它变清难。这是本轮反复出现的核心直觉，也是 Scheffer 浅湖模型的精髓。

### 1(g) 临界慢化与早期预警信号（EWS）

逼近 fold 时，主导特征值的实部趋近 0 → 系统从扰动中恢复的速度越来越慢（**临界慢化 critical slowing down**）。可观测的后果（Scheffer et al. 2009, *Nature* 461:53–59 "Early-warning signals for critical transitions"；Dakos et al. 2008, *PNAS* 105:14308–14312 "Slowing down as an early warning signal for abrupt climate change"）：
- **lag-1 自相关上升**（系统"记性变长"，这一刻越来越像上一刻）；
- **方差上升**（扰动后回不去，波动累积变大）；
- **偏度（skewness）变化**；
- **闪烁（flickering）**：在两个状态之间反复短暂跳动。

Dakos 等人的标准做法是：去趋势后在滑动窗口内拟合一阶自回归模型 AR(1)（x_{t+1} = α·x_t + ε），用系数 α 的上升趋势（Kendall's τ 检验）作为慢化指标。

**比喻**：一把快要散架的椅子，你每推它一下，它晃悠回正的时间越来越长——这就是临界慢化。**重要警告**：EWS 不是万无一失的，可能误报或漏报（噪声、突变型 non-bifurcation 翻转、数据长度都会影响），应作为"提示灯"而非"判决书"。

---

## Section 2：Lotka–Volterra 捕食–被捕食模型

### 2.1 第一性原理与历史
Alfred J. Lotka（1925，《Elements of Physical Biology》）与 Vito Volterra（1926）各自独立提出。Volterra 的动机来自生物学家 Umberto D'Ancona 的观察：第一次世界大战期间（1914–1918）亚得里亚海捕鱼强度下降，捕食性鱼类（掠食者）所占比例反而上升——为什么少捕鱼反而对掠食者更有利？

### 2.2 方程与每个符号
dx/dt = αx − βxy
dy/dt = δxy − γy

- **x**：被捕食者（猎物，如野兔）密度；
- **y**：捕食者（如猞猁）密度；
- **α**（alpha）：猎物在无捕食者时的内禀增长率（食物无限 → 指数增长 αx）；
- **β**（beta）：捕食率系数（相遇项 xy 表示"猎物遇上捕食者"的频率，−βxy 是被吃掉的损失）；
- **δ**（delta）：捕食者因吃到猎物而增殖的转化率（+δxy）；
- **γ**（gamma）：捕食者在无猎物时的自然死亡率（−γy）。

**boom-bust 比喻**：猎物多 → 捕食者吃得饱、繁殖快 → 捕食者变多 → 猎物被吃光 → 捕食者饿死 → 猎物又恢复……一场永不停歇的追逐/盛衰循环。

### 2.3 固定点、中心与守恒量
令右边为 0，得两个固定点：
- (0, 0)：双双灭绝（鞍点，不稳定）；
- **共存点 (x*, y*) = (γ/δ, α/β)**。

对共存点求 Jacobian，其特征值为**纯虚数 ±i√(αγ)**——这是一个**中心（center）**：轨迹是闭合圈，振荡既不放大也不衰减，**中性稳定（neutrally stable）**。经典模型还有一个守恒量 H = δx − γ ln x + βy − α ln y，沿轨迹不变（像能量守恒），轨迹周期约 T ≈ 2π/√(αγ)。

**警告**：中心是"刀锋上"的非典型情形——只要模型稍作现实修改，它就会变成螺旋或极限环。这恰恰说明经典 LV 模型是一个理想化起点。

### 2.4 现实修改①：Logistic 猎物 → 稳定螺旋（阻尼振荡）
给猎物加上承载力 K：dx/dt = αx(1 − x/K) − βxy。共存点变成**稳定焦点（spiral）**：扰动后做**衰减振荡**，最终回到平衡。Jacobian 的特征值变成"实部为负的复数"。

### 2.5 现实修改②：Rosenzweig–MacArthur + Holling II → Hopf 分岔与富营养化悖论
Rosenzweig & MacArthur（1963，*American Naturalist* 97:209–223，"Graphical representation and stability conditions of predator-prey interactions"）引入**Holling type II 功能反应** g(x) = ax/(1+ahx)（捕食者会被"处理时间"h 限制，吃饱了就吃不动了）：

dx/dt = rx(1 − x/K) − [a x/(1+ahx)] y
dy/dt = e [a x/(1+ahx)] y − m y

**富营养化悖论（paradox of enrichment）**——Rosenzweig（1971，*Science* 171:385–387，"Paradox of Enrichment: Destabilization of Exploitation Ecosystems in Ecological Time"）：提高猎物承载力 K（相当于给系统"加营养"），本应让系统更稳，结果反而**通过 Hopf 分岔使稳定共存点失稳，生出一个稳定极限环（持续大幅振荡）**，振幅随富营养化加剧而增大，极端时把种群压到极低值、招致灭绝。Rosenzweig 在原文警告：人类在试图"富集"生态系统以提高产量时必须非常小心，否则可能反而导致想要的物种被消灭。

```python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 经典 Lotka–Volterra
alpha, beta, delta, gamma = 1.0, 0.1, 0.075, 1.5
def LV(t, z):
    x, y = z
    return [alpha*x - beta*x*y, delta*x*y - gamma*y]

t = np.linspace(0, 60, 4000)
sol = solve_ivp(LV, [0,60], [10,5], t_eval=t, rtol=1e-9, atol=1e-9)
x, y = sol.y

fig, ax = plt.subplots(1, 2, figsize=(11,4))
ax[0].plot(t, x, 'g', label='猎物 x'); ax[0].plot(t, y, 'r', label='捕食者 y')
ax[0].set_xlabel('时间'); ax[0].set_ylabel('密度')
ax[0].set_title('时间序列：永不衰减的振荡（中心）'); ax[0].legend()
ax[1].plot(x, y, 'b'); ax[1].plot(gamma/delta, alpha/beta, 'ko')
ax[1].set_xlabel('猎物 x'); ax[1].set_ylabel('捕食者 y')
ax[1].set_title('相平面：闭合轨道（neutral cycle）')
plt.tight_layout(); plt.show()
```

### 2.6 农业案例（重头戏）
- **吹绵蚧 × 澳洲瓢虫（vedalia beetle）——经典生物防治**：吹绵蚧（cottony cushion scale, *Icerya purchasi*）在 1880 年代成为加州柑橘业的毁灭性害虫，几乎摧毁整个产业。USDA 的 C. V. Riley 派 Albert Koebele 于 1888 年赴澳，运回约 524 头澳洲瓢虫（*Rodolia cardinalis*）以及数千头寄生蝇 *Cryptochaetum iceryae*（据 Paul DeBach 史料）。1888–89 年冬释放后，到 1889 年底吹绵蚧即被基本控制，洛杉矶柑橘外运车皮数从 1888 年的约 400 车皮回升至约 2000 车皮——这是美国典型生物防治的开端（"捕食者追上猎物"的活教材）。
- **落叶松芽蛾（larch budmoth, *Zeiraphera diniana/griseana*）周期**：据瑞士联邦森林、雪与景观研究所（WSL）对上恩加丁（Engadine）山谷的长期研究，"暴发平均约每 8.5 年一次；在 4–5 个世代内，种群密度可波动高达约 30000 倍"（Wermelinger, Forster & Nievergelt 2018, WSL Fact Sheet 61）。树轮重建更显示这种规律周期由来已久——Esper, Büntgen, Frank, Nievergelt & Liebhold（2007, *Proc. R. Soc. B* 274:671–679, "1200 years of regular outbreaks in alpine insects"）用 47,513 个最大晚材密度测量重建出"过去 1173 年里持续存在的高度规律的 LBM 波动，种群峰值平均每 9.3 年一次"，规律周期持续到 1981 年后中断。Turchin et al.（2003, *Ecology* 84:1207–1214, "Dynamical Effects of Plant Quality and Parasitism on Population Cycles of Larch Budmoth"）用非线性时间序列分析判定**寄主–寄生蜂互作**（而非食物质量）为主导机制——"食物质量对芽蛾密度的影响很弱"，简单的芽蛾–寄生蜂模型可解释约 90% 的种群增长率方差。这是教科书级的"种群周期"（极限环）实例。
- **IPM（综合虫害治理）中的捕食–被捕食思维**：理解振荡与极限环，有助于避免"喷药把天敌也杀光 → 害虫报复性反弹"的悖论；引入/保育天敌相当于把系统拉回稳定共存。近期还有研究指出化肥可能通过"富营养化悖论"机制使作物–食草动物动态从稳定转为周期，加大歉收风险。

**连接框架⑤**：中心（经典 LV）、稳定螺旋（Logistic 猎物）、极限环（Rosenzweig–MacArthur 经 Hopf 分岔）——一条从"中性振荡"到"持续振荡"的完整谱系。

---

## Section 3：SIR 传染病模型

### 3.1 第一性原理与历史
Kermack & McKendrick（1927，*Proceedings of the Royal Society A* 115:700–721，"A contribution to the mathematical theory of epidemics"）。这是传染病数学建模的奠基之作。

### 3.2 方程与每个符号
dS/dt = −βSI/N
dI/dt = βSI/N − γI
dR/dt = γI

- **S**：易感者（susceptible）人数；
- **I**：感染者（infectious）人数；
- **R**：康复/移出者（recovered/removed）人数；
- **N = S + I + R**：总人口（封闭、守恒）；
- **β**（beta）：传播率（一个感染者单位时间有效接触并传染的强度）；
- **γ**（gamma）：康复/移出率，1/γ 是平均传染期。

**比喻**：β 是"火苗蹿到隔壁的速度"，γ 是"火堆烧完熄灭的速度"，R 区像"已经烧过的隔离带（firebreak）"。

### 3.3 基本再生数 R0 与阈值（transcritical 型分岔）
在疫情初期 S≈N，dI/dt ≈ (β − γ)I。于是 I 增长当且仅当 β>γ，即定义

**R0 = β/γ**

- **R0 > 1**：一个病人平均传染 >1 人 → 疫情暴发、I 增长；
- **R0 < 1**：平均传染 <1 人 → 疫情熄灭。

**比喻**：R0 = "一个病人平均会传染给几个人"。

**与分岔的联系**：无病平衡点（DFE：S=N, I=0, R=0）的稳定性在 **R0=1 处翻转**——R0<1 时 DFE 稳定（病自己消失），R0>1 时 DFE 失稳、出现"地方病/暴发"分支。这与 transcritical 分岔（dx/dt=μx−x²，μ 跨 0 时两平衡点交换稳定性）在结构上同型，R0−1 扮演 μ 的角色。

### 3.4 群体免疫阈值与最终规模
- **群体免疫阈值（herd immunity threshold, HIT）= 1 − 1/R0**：推导——当易感比例 S/N 降到 1/R0 时，有效再生数 Rₑ = R0·(S/N) = 1，疫情见顶；故需免疫比例达到 1 − 1/R0 才能从一开始阻断传播。
- **峰值条件**：当 S/N = 1/R0 时 dI/dt = 0，I 达到峰值。
- **最终规模（final size）**：即便疫情结束仍会有一部分 S 从未被感染；最终染病比例由超越方程（final size equation）决定，并随 R0 增大而增大。

```python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

N = 1000.0; gamma = 0.1   # 平均传染期 10 天
def sir(t, z, beta):
    S, I, R = z
    return [-beta*S*I/N, beta*S*I/N - gamma*I, gamma*I]

plt.figure(figsize=(7,4.5))
for R0 in [1.5, 2.5, 4.0]:
    beta = R0*gamma
    sol = solve_ivp(sir, [0,160], [N-1,1,0],
                    t_eval=np.linspace(0,160,1600), args=(beta,))
    plt.plot(sol.t, sol.y[1], label=f'R0={R0}, HIT={1-1/R0:.0%}')
plt.xlabel('时间（天）'); plt.ylabel('感染者 I')
plt.title('R0 越大，疫情峰值越高、越早'); plt.legend()
plt.tight_layout(); plt.show()
```

**已解算数值例子**：设 β=0.3/天、γ=0.1/天 → R0 = β/γ = 3。则群体免疫阈值 = 1 − 1/3 ≈ 66.7%；意味着需约 2/3 人口获得免疫才能阻断传播。

### 3.5 农业案例（重头戏）
- **植物病害流行学——SIR 的 HLIR 改编**：植物不会"走动接触"，但病害靠孢子/介体在田间一轮轮（polycyclic）传播。Jeger & van den Bosch（1994）把 SIR 改成 **HLIR（Healthy–Latent–Infectious–Removed，健康–潜伏–传染–移出）**，并给出阈值判据，分两部分发表：Jeger, M.J. & van den Bosch, F. 1994. "Threshold criteria for model plant disease epidemics. I. Asymptotic results." *Phytopathology* 84:24–27；以及 "II. Persistence and endemicity." *Phytopathology* 84:28–30。
- **植物病害的 R0 定义**（van den Bosch, McRoberts, van den Berg & Madden 2008, *Phytopathology* 98(2):239, "The Basic Reproduction Number of Plant Pathogens: Matrix Approaches to Complex Dynamics"）：原文定义为"The basic reproduction number, R0, is defined as the total number of infections arising from one newly infected individual introduced into a healthy (disease-free) host population."——和人类流行病学完全同构，阈值仍是 R0>1 才暴发。系统化处理见 Madden, Hughes & van den Bosch（2007，*The Study of Plant Disease Epidemics*, APS Press, St. Paul, MN, 421 pp.）。
- **具体作物/家畜疫情**：马铃薯晚疫病（*Phytophthora infestans*，多循环流行，Van der Plank 1963《Plant Diseases: Epidemics and Control》奠基）、小麦锈病、柑橘黄龙病（HLB，介体亚洲柑橘木虱，大量模型用 next-generation matrix 推 R0 阈值并证明 R0<1 时无病平衡全局稳定，如 Zhang et al. 2020, *Mathematical Biosciences and Engineering* 17(3):2048–2069）。家畜方面，2001 年英国口蹄疫（FMD）：Ferguson, Donnelly & Anderson（2001, *Science* 292:1155–1160；及 *Nature* 413:542–548）估计——"the basic reproduction number of FMD was estimated to be in the ranges of 3.5–4.5 among livestock farms in Great Britain"，即农场间 R0 约 **3.5–4.5**（疫情早期、全国移动禁令前），随扑杀/移动限制使有效再生数降到 1 以下。
- **管理含义**：选用抗病品种（降 β）、铲除病株 roguing（升 γ）、田间卫生/隔离带、调整种植密度——本质都是"把 R0 压到 1 以下"。

**连接框架⑤**：R0=1 阈值 = transcritical 型分岔；无病平衡点稳定性翻转。

---

## Section 4：Scheffer 浅湖 / 多稳态模型

### 4.1 第一性原理与历史
Scheffer（1990, *Hydrobiologia* 200/201:475–486，"Multiplicity of stable states in freshwater systems"）；Scheffer, Hosper, Meijer, Moss & Jeppesen（1993, *Trends in Ecology & Evolution* 8:275–279，"Alternative equilibria in shallow lakes"）；Scheffer, Carpenter, Foley, Folke & Walker（2001, *Nature* 413:591–596，"Catastrophic shifts in ecosystems"）。

### 4.2 两个交替稳态
浅湖在一定营养盐范围内存在**两个交替稳态**：
- **清水态（clear）**：沉水植物（macrophytes）主导——植物固定底泥、抑制悬浮、与藻类竞争营养 → 水清（正反馈维持清水）；
- **浊水态（turbid）**：浮游植物（phytoplankton）主导——藻类遮光使沉水植物消失、鱼搅动底泥 → 水浊（正反馈维持浑浊）。

随农业面源等导致的营养盐**缓慢升高**，湖会沿清水分支走，到达临界浊度时沉水植物崩溃，**跳到浊水态**；而要恢复，光把营养盐降回原阈值远远不够——必须降到一个**低得多的阈值**藻类才会被营养限制住、清水态才回来。这就是**滞后**。

### 4.3 最小模型与折叠平衡曲线
一个常用的通用最小模型（Scheffer et al. 1993/2001 风格，带 Hill 型正反馈）：

dx/dt = a − b·x + r · xᵖ/(xᵖ + hᵖ)

- **x**：状态变量（如藻类/浊度）；
- **a**：外部营养输入（控制参数，缓慢变化）；
- **b**：自然衰减/移除率；
- **r**：正反馈最大强度；
- **h**：半饱和常数（Hill 函数中点）；
- **p**：Hill 指数（p 越大反馈越"开关化"，p 足够大时平衡曲线呈 S 形折叠）。

当 p 较大时，x 对 a 的平衡响应呈 **S 形（折叠）曲线** → 出现**两个 fold 分岔点** → 上跳阈值 ≠ 下跳阈值 → 滞后环。

```python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

b, r, h, p = 1.0, 1.0, 1.0, 8        # p 大 → 折叠
def f(x, a): return a - b*x + r*x**p/(x**p + h**p)

def equilibrate(a, x0):
    sol = solve_ivp(lambda t,x: f(x[0],a), [0,500], [x0],
                    t_eval=[500], rtol=1e-8, atol=1e-8)
    return sol.y[0,-1]

a_up = np.linspace(0, 1.2, 120)      # 营养盐缓慢升
a_dn = a_up[::-1]                     # 再缓慢降
x = 0.05; up=[]
for a in a_up: x = equilibrate(a, x); up.append(x)   # 沿清水分支→跳浊
dn=[]
for a in a_dn: x = equilibrate(a, x); dn.append(x)   # 沿浊水分支→跳清

plt.figure(figsize=(6.5,4.5))
plt.plot(a_up, up, 'b-o', ms=3, label='营养盐上升（清→浊跳点高）')
plt.plot(a_dn, dn, 'r-o', ms=3, label='营养盐下降（浊→清跳点低）')
plt.xlabel('a（营养盐负荷）'); plt.ylabel('x（浊度/藻类）')
plt.title('折叠平衡曲线与滞后环：两个跳点不重合')
plt.legend(); plt.tight_layout(); plt.show()
```

### 4.4 管理含义
- **预防远比恢复便宜**：一旦越过上 fold，要退回必须把营养盐压到低得多的下 fold——成本极高甚至不可行。
- **biomanipulation（生物操纵）**：仅降营养往往不够；通过移除浮游食性鱼、放养食浮游动物的鱼等手段，给系统一个"大扰动"把它**推过盆地边界（basin boundary）**回到清水吸引盆——Scheffer 等（1993）及荷兰多湖（如 Veluwemeer、Wolderwijd）的实践印证了这一点。

### 4.5 农业/环境案例
- **农业面源营养（N、P）径流导致浅湖富营养化 regime shift**——最经典的 tipping 案例。
- **类比迁移**：草地荒漠化 / 灌丛侵入（**grass↔shrub 两稳态**）——半干旱牧场近几十年大量从草地翻转为灌丛态，正反馈（灌丛下"肥岛 islands of fertility"、火与放牧反馈）使其稳定难逆；Ratajczak 等（2014, *Journal of Ecology*）等强调"预防不期望的转变、局地管理优先"。土壤盐渍化、土壤退化亦可视为农业 tipping。
- **回连第 2/4 轮**：第 2 轮 MIMICS 等微生物显式模型的土壤碳多稳态、DNDC 氧化还原开关，与第 4 轮的 regime 思想，都是同一套 fold/多稳态数学。

**连接框架⑤**：fold 分岔、滞后、多稳态、吸引盆（basin of attraction）。

---

## Section 5：三者的共性与统一的动力系统图景

| 维度 | Lotka–Volterra | SIR | Scheffer 浅湖 |
|---|---|---|---|
| 状态变量 | 猎物 x、捕食者 y | S、I、R | 浊度/藻类 x（vs 沉水植物） |
| 关键参数 | α, β, δ, γ（及 K, Holling h） | β, γ（→ R0=β/γ） | 营养负荷 a（及 r, h, p） |
| 平衡/行为类型 | 中心（经典）/稳定螺旋/极限环 | 无病平衡稳定性在 R0=1 翻转 | 双稳态 + 折叠平衡曲线 |
| 标志性动力学现象 | 振荡（Hopf 生极限环） | 流行阈值（transcritical 型） | 多稳态 + 滞后（两个 fold） |
| 主要分岔 | Hopf | transcritical 型 | saddle-node / fold ×2 |
| 农业领域 | 害虫–天敌、生物防治、IPM | 植物病害/家畜疫情流行学 | 富营养化、荒漠化、盐渍化 |

**统一信息**：动力系统理论是一面**统一的透镜**——固定点、Jacobian/特征值、分岔、滞后、EWS 这一套工具，既解释害虫周期、疾病暴发、生态系统崩溃，也支撑第 2–4 轮过程模型里的阈值/翻转行为（DNDC 氧化还原开关、MIMICS 微生物显式土壤碳多稳态、各类 regime shift）。

**大结论**：农业生态系统可以是**多稳态**的，可以**突然且不可逆地翻转**。因此——
1. 用 EWS（临界慢化、自相关/方差上升）**提前预判**逼近的临界点；
2. 以**预防**为先（越过 fold 后恢复代价极高）；
3. 必要时用类似 biomanipulation 的"大扰动"把系统**推回**期望吸引盆。

---

## Section 6：Round 5 练习与里程碑

### 练习（动手做）
1. **固定点与稳定性**：对 Logistic f(N)=rN(1−N/K) 手算 f'(0)、f'(K)，判定稳定性；再画 flow-on-a-line 验证。
2. **LV 中性循环**：数值积分经典 Lotka–Volterra，画时间序列与相平面闭合轨道；验证周期 T≈2π/√(αγ)。
3. **加 Logistic 猎物**：在 LV 中给猎物加 (1−x/K)，观察振荡如何变成**阻尼振荡 / 稳定螺旋**；对共存点数值求 Jacobian 特征值，确认实部由 0 变负。
4. **R0 与群体免疫**：给定 β、γ 计算 R0 与 HIT=1−1/R0；用 SIR 代码改变 R0，观察峰值高度/时间变化，并验证 S/N=1/R0 时 I 达峰。
5. **浅湖滞后**：把营养盐 a 先升后降跑一遍 Scheffer 最小模型，测量上跳阈值与下跳阈值之差（滞后宽度）；改变 Hill 指数 p，观察 p 太小时滞后消失。
6. **EWS（进阶）**：构造一条缓慢逼近 fold 的时间序列（让参数缓慢漂移并加小噪声），用滑动窗口计算 **lag-1 自相关与方差的上升趋势**，复现临界慢化。

```python
# 练习6骨架：fold前的临界慢化（lag-1自相关上升）
import numpy as np, matplotlib.pyplot as plt
np.random.seed(0)
n=4000; x=np.zeros(n); x[0]=1.0
mu=np.linspace(1.0, 0.0, n)          # 参数缓慢逼近 fold(μ=0)
dt=0.01
for t in range(1,n):
    drift=(mu[t]-x[t-1]**2)*dt       # dx/dt=μ-x^2
    x[t]=x[t-1]+drift+0.02*np.sqrt(dt)*np.random.randn()
def ar1(seg): seg=seg-seg.mean(); return np.corrcoef(seg[:-1],seg[1:])[0,1]
w=300; ac=[ar1(x[i-w:i]) for i in range(w,n)]
fig,ax=plt.subplots(2,1,figsize=(7,5),sharex=True)
ax[0].plot(x,'b'); ax[0].set_ylabel('x'); ax[0].set_title('逼近fold：状态序列')
ax[1].plot(range(w,n),ac,'r'); ax[1].set_ylabel('lag-1 自相关')
ax[1].set_xlabel('时间步'); ax[1].set_title('临界慢化：自相关上升=预警')
plt.tight_layout(); plt.show()
```

### "你已掌握 Round 5"自查清单
- [ ] 能写出 dx/dt=f(x)，求固定点并用 f'(x*)（或 Jacobian 特征值）判稳定性；
- [ ] 能用"小球–山谷 + 特征值=回弹速率"解释稳定性，并说清复特征值=振荡；
- [ ] 能说出四种分岔的正规型与直觉，并指出 **fold 是 tipping 的核心**；
- [ ] 能解释滞后/多稳态为何让翻转"难以逆转"，并画出折叠平衡曲线；
- [ ] 能从 SIR 推出 R0=β/γ、HIT=1−1/R0，并把 R0=1 阈值理解为 transcritical 型分岔；
- [ ] 能解释临界慢化与 EWS（自相关/方差上升）及其局限；
- [ ] 能把三个模型映射到害虫周期、病害暴发、富营养化/荒漠化，并连回第 2–4 轮的阈值行为。

---

## Recommendations（建议的学习/应用路径）
1. **先跑代码再推公式**：本章每个 Python 块都可直接运行。建议顺序：flow-on-a-line → LV 时间序列/相图 → SIR 峰值 → 浅湖滞后 → EWS。看到现象后再回头啃 Jacobian/特征值，吸收率最高。
2. **把"特征值实部符号"当成贯穿全章的金线**：稳定=全负、振荡=复数、失稳/分岔=穿过 0。每遇一个新模型，先问"Jacobian 特征值在哪"。
3. **判断你面对的是哪类分岔**：要振荡找 Hopf；要阈值找 transcritical（如 R0）；要不可逆翻转找 fold/滞后。这一步决定你该用什么管理对策。
4. **管理上以阈值为锚**：传染病/病害把 R0 压到 <1（抗病品种↓β、铲除↑γ）；富营养化在越过上 fold 之前控源，越过后才考虑代价高昂的 biomanipulation。
5. **改变建议的触发阈值（benchmarks）**：若监测到 lag-1 自相关与方差持续同向上升（Kendall's τ 显著），视为逼近临界点的预警，应立即从"渐进式管理"切换为"主动干预/预防"；若 EWS 反复假阳性或数据过短，则回退到机理模型情景模拟而非单纯依赖统计指标。

## Caveats（重要提醒与不确定性）
- **EWS 并非万无一失**：临界慢化主要适用于 fold 型（分岔型）转变；噪声驱动的突变、rate-induced tipping、数据过短或去趋势不当都会导致误报/漏报。把 EWS 当提示灯，不当判决书。
- **经典 LV 的"中心"是理想化特例**：现实中几乎总会被修正为螺旋或极限环；不要把"中性循环"当作普遍规律。
- **R0=β/γ、HIT=1−1/R0 依赖均匀混合假设**：异质接触、空间结构、免疫衰减都会改变阈值；植物/家畜疫情中 R0 常需用 next-generation matrix 等更一般方法定义。本章给出的 FMD R0≈3.5–4.5 是 2001 年英国疫情早期、管控前的农场间估计，不同研究/地区差异很大。
- **吹绵蚧入侵加州的确切年份史料略有出入**（多在 1868–1880 年代之间），本章只保守陈述其在 1880 年代成为毁灭性害虫这一公认事实。
- **最小浅湖模型是定性示意**：真实湖泊涉及 N/P、鱼群、水深、风浪、气候等多因子；折叠/滞后的存在与宽度依赖参数（如 Hill 指数 p），p 太小则无多稳态。
- **跨系统类比要谨慎**：草地↔灌丛、土壤碳多稳态等"多稳态"证据强弱不一，部分仍有争议；务必区分"理论上可能多稳态"与"该系统已被实证为多稳态"。
