#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 3 · 残差混合 + 外推 + equifinality：把 Round 6 的能力落地
============================================================
目的
  1) 残差混合（hybrid）：机理"骨架" + ML 学其误差，对比 纯机理 / 纯ML / 混合 的 RMSE；
  2) 外推失败：亲手看到纯 ML 走出训练分布就"卡住"，而混合保留机理外推性；
  3) equifinality：grid-search 标定，看到多组不同参数拟合一样好。

依赖：numpy, scikit-learn（matplotlib 可选）
运行：python3 lab3_residual_hybrid.py
对标：第 6 轮教程 SECTION 1（equifinality/标定）+ SECTION 2（RF）+ SECTION 4（KGML 残差）

★ 换成你自己的数据：见文件底部 `load_your_csv()`，把一行注释打开即可。
  自学阶段先用本文件内置的、可复现的"拟真"产量数据练手。
"""

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error


# ======================================================================
# 0) 造一份"拟真"产量数据（可复现）。真值里藏了一个 rain×N 交互项，
#    我们故意让"机理模型"漏掉它——这正是现实中机理模型的典型缺陷。
# ======================================================================
def make_dataset(n=800, seed=42):
    rng = np.random.default_rng(seed)
    temp = rng.uniform(15, 30, n)          # 温度 ℃
    rain = rng.uniform(200, 800, n)        # 季降雨 mm
    Nrate = rng.uniform(0, 250, n)         # 施氮 kg/ha
    # "真实"产量 (t/ha)：含温度二次、降雨饱和、氮饱和，以及一个交互项
    truth = (3.0
             + 0.25 * (temp - 15) - 0.004 * (temp - 22) ** 2
             + 0.006 * rain - 3e-6 * rain ** 2
             + 0.02 * Nrate - 4e-5 * Nrate ** 2
             + 0.0009 * rain * Nrate / 100.0          # ← 交互项（机理模型不知道）
             + rng.normal(0, 0.35, n))                # 观测噪声
    X = np.column_stack([temp, rain, Nrate])
    return X, truth


# ======================================================================
# 1) "机理"玩具模型：基于已知农学关系，但缺了 rain×N 交互项
# ======================================================================
def process_model(X):
    temp, rain, Nrate = X[:, 0], X[:, 1], X[:, 2]
    return (3.0
            + 0.25 * (temp - 15) - 0.004 * (temp - 22) ** 2
            + 0.006 * rain - 3e-6 * rain ** 2
            + 0.02 * Nrate - 4e-5 * Nrate ** 2)        # 注意：没有交互项


def rmse(a, b):
    return float(mean_squared_error(a, b) ** 0.5)


# ======================================================================
# 2) 三方对比：纯机理 vs 纯ML vs 残差混合
# ======================================================================
def compare_three():
    X, y = make_dataset()
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=0)

    # (a) 纯机理
    proc_te = process_model(Xte)

    # (b) 纯 ML（直接学产量）
    ml = RandomForestRegressor(n_estimators=400, random_state=0).fit(Xtr, ytr)
    ml_te = ml.predict(Xte)

    # (c) 残差混合：ML 只学"机理模型的误差"，最终 = 机理 + ML(残差)
    resid_tr = ytr - process_model(Xtr)
    rf_resid = RandomForestRegressor(n_estimators=400, random_state=0).fit(Xtr, resid_tr)
    hybrid_te = process_model(Xte) + rf_resid.predict(Xte)

    print("=" * 60)
    print("三方对比（测试集 RMSE，越低越好）")
    print("-" * 60)
    print(f"  纯机理 (process-only) : {rmse(yte, proc_te):.3f} t/ha")
    print(f"  纯 ML  (RF on yield)  : {rmse(yte, ml_te):.3f} t/ha")
    print(f"  残差混合 (hybrid)     : {rmse(yte, hybrid_te):.3f} t/ha  ← 通常最低")
    print("解读：混合用 RF 补上了机理模型漏掉的 rain×N 交互，")
    print("      同时保留机理骨架的可解释性与外推性。\n")
    return ml, rf_resid


# ======================================================================
# 3) 外推失败演示：把氮肥推到训练范围外（>250），看谁还靠谱
# ======================================================================
def extrapolation_demo(ml, rf_resid):
    print("=" * 60)
    print("外推失败演示：N=400 远超训练上限 250")
    print("-" * 60)
    x250 = np.array([[22.0, 500.0, 250.0]])
    x400 = np.array([[22.0, 500.0, 400.0]])
    # 纯 ML：训练范围外预测会"卡住"（≈边界值），无法反映 N 过量减产
    print(f"  纯 ML   : N=250 → {ml.predict(x250)[0]:.2f};  "
          f"N=400 → {ml.predict(x400)[0]:.2f}  （几乎不变=拒绝外推）")
    # 混合：机理项 0.02*N - 4e-5*N^2 在 N=400 处会给出"过量减产"，外推有物理依据
    h250 = process_model(x250) + rf_resid.predict(x250)
    h400 = process_model(x400) + rf_resid.predict(x400)
    print(f"  残差混合: N=250 → {h250[0]:.2f};  "
          f"N=400 → {h400[0]:.2f}  （机理项给出过量减产，外推可信）\n")


# ======================================================================
# 4) equifinality：grid-search 标定一个 2 参数玩具模型，看多组解一样好
#    模型 C(t)=C0*exp(-k t)；真值 k=0.05, C0=100；故意只给少量含噪观测。
# ======================================================================
def equifinality_demo():
    rng = np.random.default_rng(0)
    t = np.linspace(0, 20, 12)
    C_obs = 100 * np.exp(-0.05 * t) + rng.normal(0, 3, t.size)
    ks = np.linspace(0.02, 0.09, 80)
    C0s = np.linspace(80, 120, 80)
    results = []
    for k in ks:
        for C0 in C0s:
            pred = C0 * np.exp(-k * t)
            results.append((rmse(C_obs, pred), k, C0))
    results.sort()
    print("=" * 60)
    print("equifinality：前 6 组'近乎等优'的参数（注意 k 与 C0 互相补偿）")
    print("-" * 60)
    print(f"{'RMSE':>8}{'k':>10}{'C0':>10}")
    for r, k, C0 in results[:6]:
        print(f"{r:>8.3f}{k:>10.4f}{C0:>10.1f}")
    print("解读：RMSE 几乎相同但 (k,C0) 各异——数据无法唯一确定参数，")
    print("      所以要报告参数'分布'而非单点（呼应 GLUE / Bayesian）。\n")


# ======================================================================
# ★ 换成你自己的数据（自学进阶时打开）
# ======================================================================
def load_your_csv(path):
    """期望 CSV 含列：temp, rain, Nrate, yield。返回 (X, y)。
    用法：把下面 main() 里的 make_dataset() 换成 load_your_csv('your.csv')。
    """
    import pandas as pd
    df = pd.read_csv(path)
    X = df[["temp", "rain", "Nrate"]].to_numpy()
    y = df["yield"].to_numpy()
    return X, y


def main():
    ml, rf_resid = compare_three()
    extrapolation_demo(ml, rf_resid)
    equifinality_demo()
    print("✅ Lab 3 跑通。下一步：把 make_dataset() 换成 load_your_csv('your.csv')，")
    print("   或用 Lab 2 的 WOFOST 输出当机理预测，真正接你自己的研究数据。")


if __name__ == "__main__":
    main()
