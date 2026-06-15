# Round 6（终章）：数据驱动 + 混合 + 不确定性——从机理到经验，再回到混合前沿

> 本章是整个 6 轮农业生态系统过程建模课程的"capstone（收官）"章节。面向农学专业的"小白"，强调 first-principles（第一性原理）、大量 case study（案例）和可运行的带注释 Python 代码。英文术语、模型名与文献按惯例保留英文。

---

## SECTION 0：从 Round 1–5 出发 + 机理↔经验光谱的总框架

### 0.1 我们走到哪里了

前五轮我们一直在搭建 **process-based（过程/机理）模型**，它们都是从 first-principles（守恒、限制因子、速率—状态—更新）推导出来的：

- **Round 1**：守恒定律 + explicit Euler（显式欧拉），用 Logistic、RothC、AquaCrop 把"库（pool）"和"通量（flux）"讲清楚。
- **Round 2**：土壤生物地球化学——CENTURY、DNDC、DAYCENT、MIMICS、Michaelis–Menten。
- **Round 3**：作物模型——WOFOST/SUCROS、DSSAT、APSIM、EPIC、STICS，光→生物量→产量，de Wit 生产层次。
- **Round 4**：水文 + 生态系统/陆面——SWAT、TOPMODEL、HYDRUS、Farquhar、LPJ/Biome-BGC、CLM，尺度层次、tipping-bucket vs Richards、能量平衡。
- **Round 5**：动力系统 + 临界转变——Lotka–Volterra、SIR、Scheffer 浅湖模型，不动点、Jacobian/特征值、分岔、滞后（hysteresis）、早期预警信号。

这些模型有一个共同基因：**先写下"自然界应该遵守的方程"，再求解**。它们可解释（interpretable）、能外推（extrapolate）到新条件，但需要大量过程知识，而且当某个子过程我们其实"不懂"时，模型就会系统性偏差。

### 0.2 Round 6 要去的两个新地方

Round 6 把指南针转向光谱的另一端和中间地带：

```
纯机理 (mechanistic)  ←——————————————→  纯经验 (empirical / ML)
  可解释、能外推             混合 (hybrid)            高精度(分布内)、弱外推
  需要过程知识            兼顾两者之长              数据饥渴、可解释性弱
  数据少也能跑          KGML / 可微建模 / 残差        需要大量数据
       Rounds 1–5              Round 6                  Round 6
```

- **Framework ⑥（标定—验证—不确定性 / 反问题 / equifinality / Bayesian / GP emulator）**：怎么把模型和数据对接，怎么诚实地说出"我有多不确定"。
- **Framework ⑦（mechanistic vs empirical vs hybrid；可解释性—精度权衡）**：什么时候用哪一端，怎么把两端拼起来。

三类主角：
1. **Random Forest / XGBoost**——纯数据驱动的产量预测（经验端）。
2. **Gaussian Process（GP）emulator**——既给预测又给不确定性的"桥梁工具"。
3. **KGML（Knowledge-Guided Machine Learning，知识引导机器学习）**——把机理"骨架"和 ML 融合的混合前沿。

**Capstone 核心信息（请记住这一句）**：成熟建模者要回答的从来不是"过程模型 vs ML 谁更好"，而是**针对你的问题、数据、尺度，选择或组合正确的工具**。

---

## SECTION 1：Framework ⑥——标定、验证、不确定性、equifinality、Bayesian 思维（第一性原理）

### 1.1 正问题 vs 反问题（forward vs inverse）

- **Forward modeling（正问题）**：参数 → 输出。给定参数 θ（如分解速率 k），跑模型得到输出 y。Round 1–5 我们一直在做正问题。
- **Inverse modeling（反问题）**：输出 → 参数。我们有观测数据 y_obs，**反推**最可能的参数 θ。

**关键认知转变**：参数不是"上帝给的常数"，而是**要从数据中估计的量，并且自带不确定性**。

> 🎚️ **类比（调收音机）**：calibration（标定）就像拧收音机旋钮去对准一个电台。旋钮 = 参数，听到的清晰度 = 模型与数据的拟合度。

### 1.2 Calibration vs Validation vs Evaluation（别在训练数据上考试）

- **Calibration（标定）**：用一部分数据调参数。
- **Validation（验证）**：用**没参与调参**的另一部分数据检验。
- **Evaluation（评估）**：在完全独立的场景报告最终性能。

**铁律：不要在训练数据上测试。** 否则就像"考试前老师把考卷给了你"，分数虚高。常用做法是 **split-sample（分样本）**：例如把多年数据分成"标定期"和"验证期"。

### 1.3 Equifinality（殊途同归）——没有唯一"最佳参数"

**Equifinality**：很多**不同的参数组合（甚至不同模型结构）对数据的拟合一样好**。因此并不存在唯一的"最优参数集"。

> 🎚️ **类比**：很多组旋钮位置听起来一样清楚。你没法只凭声音判断哪一组才是"真"的。

这正是 **Beven & Binley (1992)** 提出的 **GLUE（Generalized Likelihood Uncertainty Estimation，广义似然不确定性估计）** 的核心。该论文发表于 *Hydrological Processes* 6: 279–298。它是水文学领域影响极大的方法论文：作者本人在 2013 年的回顾文章《GLUE: 20 years on》中写到 "The paper has now received over 1200 citations (as of December 2012)"，此后被引持续攀升至 3000 次以上，Keith Beven 也被同行誉为"the world's most cited hydrologist"。GLUE 主张：**放弃"唯一全局最优参数集"的概念**，转而接受一组都"行为合格（behavioral）"的参数集，并用它们一起表达预测的不确定性。

### 1.4 Bayesian 思维——用证据更新信念

Bayesian（贝叶斯）公式：

```
posterior  ∝  likelihood  ×  prior
后验          似然          先验
P(θ|data)  ∝  P(data|θ)  ×  P(θ)
```

- **prior（先验）P(θ)**：看到数据前对参数的信念（如"k 大概在 0.01–0.1 之间"）。
- **likelihood（似然）P(data|θ)**：给定参数，数据出现的可能性——即"这组参数解释数据有多好"。
- **posterior（后验）P(θ|data)**：看到数据后更新的信念。

> 🔄 **类比（更新信念）**：先验 = 你的初始猜测；似然 = 新证据有多支持某个猜测；后验 = 综合后的新看法。看到越多证据，后验越被数据"拉"向真相。

**Bayesian calibration 给出的是参数的分布（distribution）而非单一数值**——这正是不确定性的来源。计算上常用 **MCMC（Markov Chain Monte Carlo，马尔可夫链蒙特卡洛）**：可以理解为一个"聪明的随机游走器"，在参数空间里多停留在后验概率高的地方，最终采样勾勒出后验分布的形状。

### 1.5 不确定性量化（UQ）——四种不确定性

1. **Parameter uncertainty（参数不确定性）**：θ 估不准。
2. **Structural uncertainty（结构不确定性）**：模型方程本身就是近似（"所有模型都是错的"）。
3. **Input uncertainty（输入不确定性）**：气象、土壤等驱动数据有误差。
4. **Observation uncertainty（观测不确定性）**：用来标定/验证的实测值也有误差。

**传播（propagation）**：把参数不确定性"灌"进模型，跑很多遍 → 得到一组输出（**ensemble，集合**）→ 输出的散布就是预测不确定性。

### 1.6 可运行代码：grid-search 标定 + equifinality 演示

```python
import numpy as np

# ---- 玩具"过程模型":一阶分解   C(t) = C0 * exp(-k * t) ----
def model(k, C0, t):
    return C0 * np.exp(-k * t)

# ---- 制造"观测数据"(真值 k=0.05, C0=100, 加噪声) ----
rng = np.random.default_rng(0)
t = np.linspace(0, 20, 12)
C_obs = model(0.05, 100, t) + rng.normal(0, 3, size=t.size)

# ---- grid-search:扫描 k 和 C0,计算每组的 RMSE ----
ks  = np.linspace(0.02, 0.09, 60)
C0s = np.linspace(80, 120, 60)
best = []
for k in ks:
    for C0 in C0s:
        rmse = np.sqrt(np.mean((model(k, C0, t) - C_obs)**2))
        best.append((rmse, k, C0))
best.sort()

print("最优:", round(best[0][0],3), "k=", round(best[0][1],4), "C0=", round(best[0][2],1))
# 打印"几乎一样好"的前若干组 —— 这就是 equifinality
print("\n前 5 组近乎等优的参数(注意 k 与 C0 互相补偿):")
for rmse, k, C0 in best[:5]:
    print(f"  RMSE={rmse:.3f}  k={k:.4f}  C0={C0:.1f}")
```

**预期现象**：前几组 RMSE 几乎相同，但 (k, C0) 各不相同——**k 偏大可被 C0 偏大补偿**。这就是 equifinality 的可视化：数据无法唯一确定参数。这也解释了为何要报告参数**分布**而非单点。

### 1.7 可运行代码：一次极简 Bayesian update

```python
import numpy as np
# 估计一个分解速率 k 的后验(网格贝叶斯,够小白直观)
k_grid = np.linspace(0.0, 0.15, 300)
prior  = np.exp(-0.5*((k_grid-0.04)/0.03)**2)   # 先验:大概 0.04 附近
# 似然:假设有 3 个独立观测点,残差服从高斯
t_obs = np.array([5,10,15]); y_obs = np.array([78, 61, 47]); C0=100; sigma=4
def loglik(k):
    pred = C0*np.exp(-k*t_obs)
    return -0.5*np.sum(((pred-y_obs)/sigma)**2)
like = np.exp(np.array([loglik(k) for k in k_grid]))
post = prior*like; post/=np.trapz(post,k_grid)   # 归一化
mean = np.trapz(k_grid*post,k_grid)
print("后验均值 k ≈", round(mean,4))              # 后验把先验"拉"向数据
```

---

## SECTION 2：Random Forest / Gradient Boosting（XGBoost）——数据驱动之极

### 2.1 第一性原理：从决策树到集成

**Decision tree（决策树）**：一系列"if-then"问题（如"降雨 > 500 mm？"→"氮肥 > 150 kg/ha？"），在叶子节点给出预测。单棵深树容易**过拟合（overfitting）**：把训练数据的噪声也背下来。

两条"集成（ensemble）"路线，对应 **bias–variance tradeoff（偏差—方差权衡）**：

**(A) Random Forest（Breiman 2001, *Machine Learning* 45:5–32）——降方差**
- **Bagging（bootstrap aggregating，自助聚合）**：从训练集有放回抽样，造出很多个略不同的数据集，各训一棵树。
- **随机特征子集**：每次分裂只看随机一部分特征 → 让树之间**去相关（de-correlated）**。
- 把许多去相关的树**平均**起来 → 方差大幅下降，预测更稳。
- **OOB（out-of-bag）误差**：每棵树没抽到的样本天然就是它的验证集，可免费估计泛化误差。
- **Feature importance（特征重要性）**：哪些输入对降低不纯度贡献最大 → 给出一点点可解释性。

> 🗳️ **类比**：Random Forest = 一群略有不同的专家投票。单个专家可能偏激，但众人取平均后稳健。

**(B) Gradient Boosting / XGBoost（Friedman 2001, *Annals of Statistics* 29:1189–1232；Chen & Guestrin 2016, KDD '16, pp. 785–794）——降偏差**
- **顺序**添加树，每棵新树去拟合当前集成的**残差/梯度**（还没解释好的部分）。
- 逐步纠错 → 偏差下降，精度通常很高。XGBoost 加了正则化、稀疏处理、并行，成为表格数据竞赛"事实标准"。

> 👥 **类比**：Boosting = 一个团队，每个新队员专门修补前面队伍犯的错。

### 2.2 关键限制：ML 不会外推

树模型的预测本质是"训练样本的加权平均"。**超出训练数据范围（out-of-distribution）时，它只能输出见过的边界值，无法外推。** 这是经验模型的根本软肋，也是后面 hybrid 要解决的痛点。

### 2.3 可运行代码：合成产量数据上的 RandomForest + 外推失败演示

```python
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

rng = np.random.default_rng(42)
n = 600
temp = rng.uniform(15, 30, n)        # 温度 ℃
rain = rng.uniform(200, 800, n)      # 降雨 mm
Nrate= rng.uniform(0, 250, n)        # 施氮 kg/ha
# "真实"产量(带饱和与噪声) —— 仅用于造数据
yield_ = (3 + 0.25*(temp-15) - 0.004*(temp-22)**2
          + 0.006*rain - 3e-6*rain**2
          + 0.02*Nrate - 4e-5*Nrate**2 + rng.normal(0,0.4,n))
X = np.column_stack([temp, rain, Nrate]); y = yield_

Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=1)
rf = RandomForestRegressor(n_estimators=300, oob_score=True, random_state=1)
rf.fit(Xtr, ytr)
print("Test RMSE :", round(mean_squared_error(yte, rf.predict(Xte))**0.5, 3))
print("OOB score :", round(rf.oob_score_, 3))
for name, imp in zip(["temp","rain","Nrate"], rf.feature_importances_):
    print(f"  importance {name}: {imp:.3f}")

# ---- 外推失败演示:把氮肥推到训练范围外(400 远超 250) ----
X_extrap = np.array([[22, 500, 400]])
print("外推预测(N=400, 训练最大仅250):", round(rf.predict(X_extrap)[0], 3),
      "→ 与 N=250 的预测几乎一样,模型拒绝外推")
print("对照 N=250:", round(rf.predict(np.array([[22,500,250]]))[0], 3))
```

### 2.4 农业案例

- **作物产量预测（旗舰案例：Jeong et al. 2016, *PLoS ONE* 11(6):e0156571）**：在全球小麦、美国县级玉米、美国东北部马铃薯/青贮玉米上比较 Random Forest 与 multiple linear regression（MLR）。摘要原文：**RF 的 RMSE 约为平均产量的 6–14%，而 MLR 为 14–49%**（"the root mean square errors (RMSE) ranged between 6 and 14% of the average observed yield with RF models in all test cases whereas these values ranged from 14% to 49% for MLR models"）。
  - **重要勘误**：该文**并未对玉米（任何案例）报告"R²"**。它用三项统计量：RMSE、Nash–Sutcliffe model efficiency（EF）、Willmott's d。对美国县级 30 年玉米：RF 的 **EF = 0.76、d = 0.92、RMSE = 1.13 t/ha（占均值 16.7%）**，观测—预测 **Pearson r = 0.88**（r² ≈ 0.77）；同案例 MLR 仅 EF = 0.30、RMSE = 1.94 t/ha（28.7%）。**唯一明确的"方差解释"是对全球小麦的 96%**（"The RF model explained 96% of yield variance"），不是玉米。引用时切勿把玉米的 EF=0.76 误标为"R²"。
  - 另需注意：逐案例精确百分比为 RF 5.8%–16.7%、MLR 13.8%–49.2%，摘要中的"6–14% vs 14–49%"是四舍五入后的近似概括，方向性结论（RF 大幅优于 MLR）完全成立。
- **遥感 + ML 大区域估产**：以 Jeong et al. 2016 的全球小麦案例为具名实例——用气候与生物物理变量作特征，Random Forest 解释了约 96% 的产量方差，是大尺度统计估产的典型范式。
- **物种分布 / 土地利用分类**：RF 在遥感分类里是经典基线。

### 2.5 连接两大框架

- **Framework ⑦（经验端）**：高精度、弱外推、可解释性有限，但 feature importance 给一点洞察。
- **Framework ⑥（不确定性）**：用 cross-validation 估泛化；用 **quantile regression forests（Meinshausen 2006, *JMLR* 7:983–999）** 或集成散布给出预测区间。

---

## SECTION 3：Gaussian Process（GP）emulator——桥梁工具

### 3.1 第一性原理：函数上的概率分布

**Gaussian Process（高斯过程）= 一个"函数的分布"**。它不只给一条预测曲线，而是在**每个点同时给出均值预测和校准过的不确定性（方差）**。这一点让它天生适合做"仿真器/不确定性"工作（Rasmussen & Williams 2006, *Gaussian Processes for Machine Learning*, MIT Press）。

工作流程（小白版）：
1. **prior（先验）**：先假设函数是"平滑"的，平滑程度由 **kernel（核 / covariance function，协方差函数）** 决定，如 **RBF（径向基 / squared-exponential，平方指数）核**：距离近的点函数值相关性高。
2. **conditioning（条件化）**：用观测点"钉住"这些函数。
3. **posterior（后验）**：得到后验均值 + credible interval（可信区间）。**离数据越远，不确定性越大**——这是 GP 最美的性质。

> 🟦 **类比**：GP = 一张柔软的"橡皮膜（rubber sheet）"，既穿过数据点，又会告诉你"在没数据的地方我没把握"（膜在那里晃得厉害）。

### 3.2 Emulation / Surrogate modeling（仿真器/代理模型）

一个昂贵的过程模型（如 Round 4 的 CLM、或 DSSAT）跑一次要很久。**emulator（仿真器）** 的思路：
1. 在精心设计的参数点（如 Latin hypercube 采样）上跑少数几次昂贵模型；
2. 训练一个**快速 GP** 去模仿"输入→输出"映射；
3. 之后用 GP 廉价地评估成千上万次，用于 sensitivity analysis（敏感性分析）、calibration（标定）、UQ。

> 🎬 **类比**：emulator = 慢模型的"快替身/特技演员（stunt double）"。危险又昂贵的镜头让替身上，省时省钱。

### 3.3 可运行代码：稀疏 1-D 数据上的 GP + 置信带变宽

```python
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# 真函数(假装是昂贵模型) f(x)=x*sin(x)
f = lambda x: x*np.sin(x)
X_train = np.array([1., 3., 5., 6., 7.]).reshape(-1,1)   # 稀疏样本
y_train = f(X_train).ravel()

kernel = C(1.0)*RBF(length_scale=1.0)
gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-6, normalize_y=True)
gp.fit(X_train, y_train)

X_test = np.linspace(0, 10, 200).reshape(-1,1)
mu, sd = gp.predict(X_test, return_std=True)
lo, hi = mu - 1.96*sd, mu + 1.96*sd      # 95% 置信带

# 检查:远离训练点处 sd 更大
print("x=6.0 附近 sd:", round(float(gp.predict([[6.0]], return_std=True)[1]),3))
print("x=9.5 (外推区) sd:", round(float(gp.predict([[9.5]], return_std=True)[1]),3))
# 画图(可选):plt.fill_between(X_test.ravel(), lo, hi, alpha=.3)
```

**预期**：在 x=9.5（没有训练点）处的标准差远大于 x=6.0 附近——置信带向无数据区张开，这正是 GP 自带的诚实不确定性。

### 3.4 案例

- **陆面模型参数敏感性（旗舰案例：Gao et al. 2021, *J. Hydrometeorology* 22(2):259–278）**：对 CLM5 的五个土壤相关参数（saturated hydraulic conductivity、porosity、saturated matric potential、shape parameter、organic matter fraction）用 **GP emulator** 在五维参数空间做 variance-based 敏感性分析。结论：**土壤水分方差主要由 porosity（孔隙度）和 shape parameter（持水曲线形状参数）主导**；其中表层土壤水分对 shape parameter 更敏感，根区土壤水分则 porosity 更重要。emulator 只用少量 CLM5 模拟（Maximin Latin hypercube 采样）就能高精度复现 CLM5。
- **标定昂贵作物/生态模型**：先建 emulator，再在其上做 Bayesian calibration，把原本不可行的上万次评估变得可行。
- **管理的 Bayesian optimization**：用 GP 指导"下一个最该试的施肥/灌溉方案"。

### 3.5 连接框架

- **Framework ⑥**：GP 是 UQ 和 emulator 加速标定/敏感性分析的核心工具。
- **Framework ⑦**：GP 是连接机理世界和数据世界的桥梁——它既能拟合数据，又把昂贵机理模型"压缩"成快替身。

---

## SECTION 4：KGML——Knowledge-Guided Machine Learning（混合前沿）

### 4.1 第一性原理：把过程知识"焊进"ML

**KGML / physics-informed / knowledge-guided ML** = 把过程知识与机器学习融合。四种主流混合策略（小白版，保持自洽）：

1. **Emulation / surrogate（仿真，见 Section 3）**：ML 模仿过程模型以提速。
2. **Residual modeling（残差建模）**：ML 学习过程模型的**误差**，`corrected = process + ML(features)`。机理骨架保留 → 保留可解释性。
   > 🔧 **类比**：给过程模型外挂一个 ML"纠错器（error-corrector）"。
3. **Parameter replacement / differentiable modeling（参数替换 / 可微建模）**：把某个薄弱子过程（如胁迫函数）换成 ML 组件，端到端训练。代表作 **Shen et al. 2023, *Nature Reviews Earth & Environment* 4(8):552–567**（可微建模统一 ML 与物理模型）。
4. **Physics-informed loss / constraints（物理约束损失）**：在 ML 损失里惩罚违反守恒律的行为。代表作 **Raissi et al. 2019, *J. Computational Physics* 378:686–707（PINN）** 和 **Karpatne et al. 2017, *IEEE TKDE* 29(10):2318–2331（theory-guided data science）**。
   > 🚓 **类比**：物理约束损失 = "学生答题若违反物理定律就扣分"，逼它给出物理上说得通的答案。

### 4.2 为什么混合更好

- **比纯 ML 更能外推**：保留机理骨架，离开训练分布时仍受物理约束。
- **比纯过程更准**：在某个子过程"我们不懂"的地方让数据说话。
- **数据需求更低**：过程模型可生成合成数据用于 pretraining（预训练）。

### 4.3 旗舰案例

- **KGML-ag-Carbon（Liu et al. 2024, *Nature Communications* 15:357，DOI 10.1038/s41467-023-43860-5）**：开发分三步——(1) 用农业过程模型导出的因果关系搭 ML 架构；(2) 用过程模型（**ecosys**）生成的合成数据 **pre-train**；(3) 用稀疏分布的涡度相关站点观测与低分辨率产量数据 **fine-tune**；预训练和微调阶段都加入源自过程模型的 knowledge-guided losses 来约束响应。结论：以美国 Corn Belt 为试验场，KGML **在精度上超过传统过程模型和纯黑箱 ML 模型（尤其在数据有限时）**。据 University of Minnesota（College of Science and Engineering）2024-02-21 官方新闻稿，该方法"is 10,000 times faster than current systems"（比传统过程模型快 1 万倍以上）；论文原文称其高分辨率方法 "quantitatively reveals 86% more spatial detail of soil organic carbon changes than conventional coarse-resolution approaches"（多揭示 86% 的土壤有机碳变化空间细节）。
- **APSIM + KGML 渍水—产量（Yangtze River Basin 案例，2024, ScienceDirect S2666916124000380）**：用 transfer learning 把 waterlogging-enabled APSIM 的渍水过程迁移到八个网格作物模型。原文："KGML could accurately replicate the behavior of the improved APSIM model under waterlogging conditions, achieving an R2 of 0.83 and an RMSE of 272.3 kg/ha for yield loss simulations"，并指出 "Soil properties were identified as the primary factors influencing yield losses under waterlogging"。
- **可微水文（Feng et al. 2022, *Water Resources Research* 58(10):e2022WR032404）**：以 HBV 为骨架嵌入神经网络做参数化/增强。原文："δ models can obtain a median Nash-Sutcliffe efficiency of 0.732 for 671 basins across the USA for the Daymet forcing data set, compared to 0.748 from a state-of-the-art LSTM model with the same setup"（另一驱动数据集 NLDAS 下为 0.715 vs LSTM 0.722），逼近 state-of-the-art LSTM，同时保留过程清晰度，还能输出未训练的物理变量（蒸散、基流）。

### 4.4 可运行代码：残差校正混合（概念示意）

```python
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

rng = np.random.default_rng(7)
n = 500
temp = rng.uniform(15,30,n); rain = rng.uniform(200,800,n); N = rng.uniform(0,250,n)
# "真值"(含一个过程模型没有的交互项 0.0008*rain*N/100)
truth = (3 + 0.006*rain - 3e-6*rain**2 + 0.02*N - 4e-5*N**2
         + 0.0008*rain*N/100 + 0.25*(temp-15) + rng.normal(0,0.3,n))

# ---- 简化"过程模型":漏掉了 rain×N 交互项 ----
def process_model(temp, rain, N):
    return 3 + 0.006*rain - 3e-6*rain**2 + 0.02*N - 4e-5*N**2 + 0.25*(temp-15)

X = np.column_stack([temp,rain,N])
proc = process_model(temp,rain,N)
resid = truth - proc                       # 过程模型的"误差"

Xtr,Xte,rtr,rte,ptr,pte,ttr,tte = train_test_split(
    X, resid, proc, truth, test_size=0.3, random_state=0)

# 混合:用 RF 学残差
rf = RandomForestRegressor(n_estimators=300, random_state=0).fit(Xtr, rtr)
hybrid_pred = pte + rf.predict(Xte)

rmse_proc   = mean_squared_error(tte, pte)**0.5
rmse_hybrid = mean_squared_error(tte, hybrid_pred)**0.5
print("纯过程模型 RMSE :", round(rmse_proc,3))
print("残差混合   RMSE :", round(rmse_hybrid,3), "← 通常显著更低")
```

**预期**：混合模型 RMSE 明显低于纯过程模型，因为 RF 学到了过程模型遗漏的 rain×N 交互——同时保留了机理骨架的可解释性与外推性。

### 4.5 回到整门课

Round 1–5 的过程模型，**正是引导 ML 的那份"知识"**。没有 Round 1–5，你就没有 KGML 里的"K"。这就是 capstone 的闭环。

---

## SECTION 5：选择 mechanistic vs empirical vs hybrid（决策框架）

### 5.1 决策表

| 情况 | 首选 | 原因 |
|---|---|---|
| 需要外推到新气候/新管理 | **纯过程** | 守恒律在分布外仍成立 |
| 要机理理解、情景分析 | **纯过程** | 可解释、可问"为什么" |
| 数据稀少但过程知识好 | **纯过程** | 不需要大数据 |
| 分布内数据充足、只要预测 | **纯 ML** | 高精度、易上手 |
| 过程机理不清楚 | **纯 ML** | 让数据说话 |
| 需要极快评估（实时/大批量） | **纯 ML 或 GP emulator** | 速度 |
| 有已知骨架 + 某子过程薄弱 | **hybrid（残差/可微）** | 兼顾精度与可解释 |
| 既要可解释又要高精度 | **hybrid / KGML** | 两全 |
| 想用过程模型降低数据需求 | **hybrid（合成预训练）** | 数据高效 |

### 5.2 文字流程图

```
你需要外推到训练数据没覆盖的条件吗?
 ├─ 是 → 过程知识够好吗?
 │        ├─ 够 → 纯过程模型(必要时 GP emulator 加速标定/UQ)
 │        └─ 不够(某子过程薄弱) → hybrid / KGML
 └─ 否(分布内预测) → 数据多吗?
          ├─ 多 → 纯 ML(RF/XGBoost);要预测区间 → quantile RF / GP
          └─ 少 → hybrid(用过程模型合成数据预训练)
```

核心张力是 **interpretability–accuracy tradeoff（可解释性—精度权衡）**，再叠加 **数据量、外推需求、计算预算** 三个维度。

---

## SECTION 6：三者异同 + 课程 capstone

### 6.1 三类方法对照表

| 维度 | Random Forest / XGBoost | GP emulator | KGML（混合） |
|---|---|---|---|
| 是什么 | 树的集成 | 函数的概率分布 | 机理 + ML 融合 |
| 核心原理 | bagging 降方差 / boosting 降偏差 | 核函数定义的先验 + 条件化 | 用知识约束/引导 ML |
| 可解释性 | 中（feature importance） | 中（核、长度尺度） | 高（保留机理骨架） |
| 外推能力 | 弱 | 弱—中（但诚实报告不确定） | 较强 |
| 数据需求 | 大（分布内） | 小—中 | 中（可用合成数据降需求） |
| 不确定性量化 | 间接（quantile RF/集成） | **内生、校准良好** | 可继承机理+ML 两者 |
| 主要农业用途 | 估产、分类 | 敏感性/标定/UQ、代理 | 碳/氮通量、产量、渍水 |
| 上手难度 | 低 | 中 | 高 |

### 6.2 整门课的 capstone 综合

**7 个 thinking frameworks 回顾**：
1. **守恒定律 + rate–state–update**（Round 1）——库与通量、Euler 积分。
2. **限制因子 / 多过程耦合**（Round 2–3）——Michaelis–Menten、Liebig、de Wit 层次。
3. **光→生物量→产量**（Round 3）——作物生长引擎。
4. **尺度层次 + 能量/水量平衡**（Round 4）——点→田→流域→全球；tipping-bucket vs Richards。
5. **动力系统思维**（Round 5）——不动点、Jacobian/特征值、分岔、hysteresis、早期预警。
6. **标定—验证—不确定性 / 反问题 / equifinality / Bayesian**（Round 6 Framework ⑥）。
7. **mechanistic vs empirical vs hybrid + 可解释性—精度权衡**（Round 6 Framework ⑦）。

**20 个模型在光谱上的映射**：
- **机理端**：Logistic、RothC、AquaCrop、CENTURY、DNDC、DAYCENT、MIMICS、WOFOST/SUCROS、DSSAT、APSIM、EPIC、STICS、SWAT、TOPMODEL、HYDRUS、Farquhar、LPJ/Biome-BGC、CLM、Lotka–Volterra、SIR、Scheffer 浅湖（全部 Round 1–5）。
- **混合端**：KGML、可微 HBV、残差校正、PINN——把上面的机理当骨架。
- **经验端**：Random Forest、XGBoost；**GP** 横跨混合—经验，既是 ML 又是 emulator 桥梁。

**统一信息**：守恒律 + 限制因子 + rate–state–update + 尺度层次 + 动力系统思维 + 标定/不确定性 + 机理/经验/混合思维——这七件套合起来，就是一套完整的农业生态系统建模"心智工具箱（mental toolkit）"。

### 6.3 下一步建议

1. 用真实数据集（如 eddy-covariance 通量、县级产量）复现一个 Section 案例。
2. 学 PyTorch/JAX 的 automatic differentiation，动手做一个可微模型（参考 Feng et al. 2022 的 δHBV）。
3. 深入一种 UQ 方法（MCMC 用 PyMC，或 GP 用 GPyTorch）。
4. 读旗舰论文原文（本章引用），关注"它如何标定、如何报告不确定性"。

---

## SECTION 7：Round 6 练习 / 里程碑 + 全课程总清单

### 7.1 动手练习（4–6 个）

1. **RF 估产 + 外推失败**：在合成产量数据上训 `RandomForestRegressor`，打印 feature importance，然后把某特征推到训练范围外，验证预测"卡住"不外推。
2. **grid-search 显 equifinality**：扫描 1–2 参数玩具模型，列出 RMSE 几乎相同但参数不同的多组解，解释 k 与 C0 的补偿。
3. **GP 稀疏拟合**：用 `GaussianProcessRegressor` + RBF 核拟合 5 个点，画 95% 置信带，确认远离数据处带宽变大。
4. **残差校正混合**：玩具过程模型（故意漏交互项）+ RF 学残差，比较纯过程、纯 RF、混合三者 RMSE。
5. **极简 Bayesian update**：网格贝叶斯估一个参数，改变先验宽度，观察后验如何被数据"拉动"。
6. **设计一个 KGML 工作流（纸面）**：为你关心的问题（如 N₂O 排放）写出"过程模型生成合成数据 → 预训练 → 观测微调 → 物理约束损失"的流程。

### 7.2 "你掌握了 Round 6，如果你能……"

- [ ] 用自己的话解释正问题 vs 反问题，并说出参数为何带不确定性。
- [ ] 解释 calibration/validation 的区别和"别在训练集上考试"。
- [ ] 用调收音机的类比讲清 equifinality，并说出 GLUE 的主张。
- [ ] 写出 `posterior ∝ likelihood × prior`，解释每一项。
- [ ] 区分 bagging（降方差）与 boosting（降偏差）。
- [ ] 演示并解释 RF 为何不会外推。
- [ ] 说清 GP 为何能同时给均值和不确定性，以及置信带为何远离数据变宽。
- [ ] 解释 emulator 是什么、为什么能加速标定与敏感性分析。
- [ ] 列出 KGML 四种策略并各举一例。
- [ ] 用决策表为给定问题选出 mechanistic / empirical / hybrid。

### 7.3 "你完成了整门 6 轮课程，如果你能……"（capstone 总清单）

- [ ] **F①** 守恒律 + Euler：写出库—通量方程并显式积分。
- [ ] **F②** 限制因子：用 Michaelis–Menten/Liebig 解释速率受限。
- [ ] **F③** 光→生物量→产量：讲清作物模型主干与 de Wit 层次。
- [ ] **F④** 尺度层次 + 能水平衡：区分 tipping-bucket 与 Richards，点到全球的尺度跃迁。
- [ ] **F⑤** 动力系统：求不动点、用 Jacobian/特征值判稳定性、解释分岔/hysteresis/早期预警。
- [ ] **F⑥** 标定—验证—不确定性：做 split-sample、识别 equifinality、用 Bayesian 思维量化不确定性。
- [ ] **F⑦** 机理/经验/混合：在光谱上定位任一模型，并按问题—数据—尺度选择或组合工具。
- [ ] 把 20 个模型在 mechanistic↔hybrid↔empirical 光谱上各归其位。
- [ ] 用一句话总结全课：**先用守恒律和限制因子搭骨架，用动力系统理解行为，用标定/不确定性诚实评估，用机理—经验—混合思维选对工具。**

🎓 **恭喜——你已经走完从纯机理到纯数据、再回到混合前沿的完整旅程。**

---

### 附：本章核心文献（均已核验）

- Breiman, L. (2001). Random Forests. *Machine Learning* 45:5–32. DOI 10.1023/A:1010933404324.
- Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. *Annals of Statistics* 29(5):1189–1232.
- Chen, T. & Guestrin, C. (2016). XGBoost: A Scalable Tree Boosting System. *KDD '16*, pp. 785–794. DOI 10.1145/2939672.2939785.
- Meinshausen, N. (2006). Quantile Regression Forests. *JMLR* 7:983–999.
- Rasmussen, C.E. & Williams, C.K.I. (2006). *Gaussian Processes for Machine Learning*. MIT Press.
- Beven, K. & Binley, A. (1992). The future of distributed models: model calibration and uncertainty prediction. *Hydrological Processes* 6:279–298.
- Jeong, J.H. et al. (2016). Random Forests for Global and Regional Crop Yield Predictions. *PLoS ONE* 11(6):e0156571.
- Gao, X. et al. (2021). Emulation of CLM5 to Quantify Sensitivity of Soil Moisture to Uncertain Parameters. *J. Hydrometeorology* 22(2):259–278. DOI 10.1175/JHM-D-20-0043.1.
- Raissi, M., Perdikaris, P. & Karniadakis, G.E. (2019). Physics-informed neural networks. *J. Computational Physics* 378:686–707.
- Karpatne, A. et al. (2017). Theory-Guided Data Science. *IEEE TKDE* 29(10):2318–2331.
- Shen, C. et al. (2023). Differentiable modelling to unify machine learning and physical models for geosciences. *Nature Reviews Earth & Environment* 4(8):552–567.
- Liu, L. et al. (2024). Knowledge-guided machine learning can improve carbon cycle quantification in agroecosystems. *Nature Communications* 15:357. DOI 10.1038/s41467-023-43860-5.
- Feng, D. et al. (2022). Differentiable, Learnable, Regionalized Process-Based Models… *Water Resources Research* 58(10):e2022WR032404.

**两处需向读者透明声明的准确性说明**：(1) Jeong et al. (2016) 原文并未报告玉米的"R²"，其报告的是 EF=0.76、d=0.92、Pearson r=0.88（r²≈0.77）；唯一明确的"方差解释 96%"针对全球小麦。课程项目中若曾写"corn-belt R²≈0.72–0.80"应据此修正。(2) Jeong 摘要的"RMSE 6–14% vs 14–49%"为四舍五入概括，逐案例实测为 RF 5.8–16.7%、MLR 13.8–49.2%。


---

## 🧪 动手练习（配套脚本）

读完这一轮，去跑配套练习，把核心机理亲手验证一遍——随机森林 + GP 仿真器（带不确定性）：

```bash
cd exercises && python3 round6_exercise.py
```

→ 脚本：[`exercises/round6_exercise.py`](exercises/round6_exercise.py) ｜ 全部练习说明见 [`exercises/README.md`](exercises/README.md)
