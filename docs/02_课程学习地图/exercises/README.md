# 六轮课程 · 配套练习脚本

每一轮课程配一个**自包含、可直接运行**的 Python 脚本，把那一轮的核心 first-principles
亲手跑一遍。和 `03_动手实操包` 的 4 个 lab 一样，每个脚本都带**预期输出**和**"你的练习"**接口。

> 这些是"读完一轮 → 立刻动手 5 分钟"的小练习；想要更完整的科研工作流（标定 / 不确定性 /
> WOFOST 真实气象），去跑 [`03_动手实操包`](../../03_动手实操包/README.md) 的四个 lab。

## 环境

```bash
pip install numpy scipy scikit-learn      # round6 才用到 scikit-learn
```

全部纯本地、离线可跑（不像 Lab 2 需要联网）。

## 清单

| 脚本 | 对标章节 | 你会亲手看到 |
|---|---|---|
| [`round1_exercise.py`](round1_exercise.py) | Round 1 地基 | 显式 Euler 步长减半→误差减半（收敛）；一阶衰减精确步进 + 质量守恒查账 |
| [`round2_exercise.py`](round2_exercise.py) | Round 2 土壤生化 | Michaelis–Menten 的一阶/零阶两个极限；输入翻倍时**线性库翻倍、微生物显式不翻倍**（非线性） |
| [`round3_exercise.py`](round3_exercise.py) | Round 3 作物 | 30 行玩具作物模型跑出 de Wit 潜在 vs 水分限制产量差（Lab 2 的离线版） |
| [`round4_exercise.py`](round4_exercise.py) | Round 4 水文/生态 | SCS-CN 产流随 CN 变化；Farquhar 光合 Ac/Aj 交叉点（Rubisco 限制 ↔ 光限制） |
| [`round5_exercise.py`](round5_exercise.py) | Round 5 动力系统 | SIR 的 R0 阈值与群体免疫；Scheffer 浅湖迟滞（同一营养、两种命运）；早期预警信号上升 |
| [`round6_exercise.py`](round6_exercise.py) | Round 6 数据驱动 | 随机森林产量预测（过拟合信号 + 特征重要性）；GP 仿真器自带不确定性误差带 |

## 跑法

```bash
cd docs/02_课程学习地图/exercises
python3 round1_exercise.py        # 依次 round1 … round6
```

每个脚本都在本仓库环境里逐个运行验证过。哪一步报错或数字对不上，把输出贴出来即可定位。
