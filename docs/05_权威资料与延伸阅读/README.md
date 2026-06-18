# 05 · 权威资料与延伸阅读（避免重复造轮子）

> 本指南把概念用大白话讲了一遍、配了能跑的代码。但"**最权威的原始说法、官方能跑的完整版、
> 系统的公开课和教材**"在外面已经有人做得很好了——这一章就是把它们**精选、核实、按小白友好度分级**
> 收进来，给你指路，**不重复造轮子**。
>
> 这里**只收第一梯队**：各模型官方主页/文档、领域公认教材、知名学者公开课、官方代码包、权威机构数据集。

本章四个文件：

| 文件 | 内容 |
|---|---|
| **本页（总览）** | 怎么用 · 🥇 每个方向"从这里开始" · **从零到一的外部学习路线** |
| [01 · 分主题资料库](01_分主题资料库.md) | 按 6 大主题(数学地基/土壤碳/作物/水文光合/动力系统/ML-不确定性)，每个模型指向官方与经典资料 |
| [02 · 公开数据集与可复现工具](02_公开数据集与可复现工具.md) | 喂给模型的气象/土壤/通量/产量数据 + 让结果可复现的工具链 |
| [03 · 权威开源项目](03_权威开源项目.md) | 能 clone 即跑的模型实现 + 数据工具栈 + **"怎么判断项目靠不靠谱"小白清单** |

---

## 怎么用这一章（重要）

1. **先用本指南建立直觉、跑通对应 Lab**，再来这里点权威材料——顺序反了容易被官方文档的细节劝退。
2. 每条资料都标了 **适合（入门/进阶/高级）** 和 **🥇（该方向小白首选）**。**没标🥇的、标"高级"的，第一遍可以跳过。**
3. **💰 = 付费**（多为纸质教材，大学图书馆通常可借；公开课一般可免费旁听）。其余默认免费。
4. **打不开的链接**：部分官方站(FAO、出版社、Rothamsted、NREL、UNH、政府数据站)会拦截自动访问(403)，
   **这是反爬、不是失效**。深层 PDF 点不开，就回该模型**官方主页**站内搜标题。

---

## 🥇 每个方向"从这里开始"（选择困难时看这张）

把全章每个方向的小白首选挑出来，集中一张表。**只看这一列也能起步：**

| 方向 | 从这里开始 | 为什么是它 |
|---|---|---|
| 💻 先会跑代码 | [Software Carpentry · Python](https://swcarpentry.github.io/python-novice-inflammation/) | 半天上手，研究人员学 Python 的黄金起点 |
| ➗ 微分方程直觉 | [3Blue1Brown · 微分方程](https://www.3blue1brown.com/topics/differential-equations) | 动画建立 `dC/dt` 画面感，零基础可看 |
| 📘 通用建模教材 | [Soetaert《生态建模实用指南》+ deSolve](https://link.springer.com/book/10.1007/978-1-4020-8624-3) | 和本指南最对口，从守恒律推到 ODE 全程 R 实战 |
| 🟤 土壤碳 | [RothC 官方](https://www.rothamsted.ac.uk/rothamsted-carbon-model-rothc) / [SoilR First steps](https://www.bgc-jena.mpg.de/TEE/basics/2015/09/25/First-steps/) | 模型诞生地 + 一个 R 包跑通 RothC |
| 🌾 作物 | [FAO AquaCrop](https://www.fao.org/aquacrop/en) | 作物模型里参数最少、最适合上手，免费 GUI+视频 |
| 💧 水文/光合 | [SWAT+ 官网](https://swat.tamu.edu/) / [plantecophys(R)](https://cran.r-project.org/package=plantecophys) | 流域水文官方门户 / 动手玩 Farquhar 摩擦最低 |
| 🔄 动力系统/突变 | [Strogatz 公开课](https://www.youtube.com/playlist?list=PLbN57C5Zdl6j_qJA-pARJnKsmROzPnO9V) | 平衡/振荡/分岔直觉的最佳来源，数学友好 |
| 🤖 ML/不确定性 | [scikit-learn 指南](https://scikit-learn.org/stable/user_guide.html) / [VanderPlas《PDSH》](https://jakevdp.github.io/PythonDataScienceHandbook/) | 经典 ML 权威参考 + 最平易近人的免费教材 |
| 📊 数据 | [NASA POWER](https://power.larc.nasa.gov/data-access-viewer/) | 免注册一键下全球气象，直接喂作物模型 |
| 🧰 数据处理 | [xarray](https://github.com/pydata/xarray) / [Project Pythia](https://foundations.projectpythia.org/) | 把 NetCDF/遥感数据读进来、画出来的核心栈（详见 03） |
| 🧩 开源项目 | [03 · 权威开源项目](03_权威开源项目.md) | clone 即跑的模型实现 + "怎么挑靠谱项目"清单 |
| 🔁 可复现 | [The Turing Way](https://book.the-turing-way.org/) | 让你的模型实验别人能重跑的方法论底座 |

---

## 从零到一：一条全程免费的外部学习路线

如果你想"**有充分扎实的基础**"，按下面顺序走——**每一步都先用外部资源建立直觉，再回本指南做对应的
Round 练习 / Lab，形成"看懂 → 动手"闭环**。全程用上表的免费资源即可。

```
第 0 步 · 把工具备好（半天）
  Software Carpentry「Programming with Python」→ 会跑 .py、会用 Jupyter
  ↳ 配套：本指南《安装与常见报错》

第 1 步 · 微分方程直觉（1–2 天）
  3Blue1Brown「微分方程」前 3 集（看动画，不用做题）
  ↳ 配套：本指南《数学补给站》§1–§4（含手把手 Euler 例子）

第 2 步 · 会用工具求解（1 天）
  SciPy solve_ivp 教程 + Scientific Python Lectures 的 SciPy 章
  ↳ 配套：Round 1 + round1_exercise.py（亲眼看"步长减半、误差减半"）

第 3 步 · 第一个机理模型 = 土壤碳（2–3 天）
  RothC 官方说明 + SoilR「First steps」跑通一个 RothC 模拟
  ↳ 配套：Round 1–2 + Lab 1（复现 Rothamsted 官方四个数）

第 4 步 · 第一个作物模型（2–3 天）
  FAO AquaCrop（GUI，参数少）跑一季玉米；或 PCSE/WOFOST 文档
  ↳ 配套：Round 3 + Lab 2（NASA POWER 真实气象，潜在 vs 水分限制）

第 4.5 步 · 学会"喂数据"（可选但强烈推荐，2–3 天）
  Project Pythia Foundations 学 xarray → 把 ERA5/MODIS 的 NetCDF 读进来、画出来
  ↳ 配套：02 公开数据集 + 03 数据工具栈

第 5 步 · 动力系统与突变（2–3 天）
  Strogatz 公开课前 5 讲 + 玩 MIT Mathlets 分岔交互
  ↳ 配套：Round 5 + round5_exercise.py（SIR 阈值、Scheffer 迟滞、早期预警）

第 6 步 · 数据驱动与不确定性（3–5 天）
  VanderPlas《PDSH》第5章 + scikit-learn 入门 → SALib 做敏感性
  → PyMC / Statistical Rethinking 前几讲 看贝叶斯
  ↳ 配套：Round 6 + Lab 3（残差混合）+ Lab 4（标定/Sobol/GLUE/GP）

第 7 步 · 让工作可复现（1 天，贯穿始终）
  The Turing Way 通读关键章 + 用 Zenodo 给代码存档拿 DOI
  ↳ 配套：小白教程第4部分（从想法到可复现论文）
```

> **不必一次走完。** 只有一个周末？走第 0–3 步，能跑通 Lab 1 就已经入门了。
> 想深入某一块，再去 [01 · 分主题资料库](01_分主题资料库.md) 里那一节按"入门→进阶→高级"往下挖。

---

## 这一章和本指南其它部分的关系

- **不重复**：本指南正文负责"用大白话讲懂 + 给最小可跑代码"；这一章负责"指向**官方/权威的完整版**"。
- **要数据**？→ [02 · 公开数据集](02_公开数据集与可复现工具.md)。
- **想要能 clone 即跑的开源项目 / 代码**？→ [03 · 权威开源项目](03_权威开源项目.md)（含"怎么判断一个开源项目靠不靠谱"小白清单）。
- **某个词不懂**？→ [术语速查表](../00_新手上路/术语表.md)；**某个符号不懂**？→ [数学补给站](../00_新手上路/数学补给站.md)。
- **诚实提醒**：外部链接会随官网改版变动；若某条失效，用资源名称在官方主页或搜索引擎重搜即可——
  名称都是稳定的，且都选了各领域**公认第一梯队**的资料。
