#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 6 练习 · 数据驱动与混合：随机森林产量预测 + GP 仿真器（带不确定性）
==========================================================================
对标章节：02_课程学习地图/Round6_数据驱动_混合_不确定性.md
框架：⑥ 校准-验证-不确定性、⑦ 机理 vs 经验 vs 混合

两件事（与 Lab 3 / Lab 4 互补，聚焦"纯数据驱动"这一极）：
  A) 随机森林预测产量：训练/测试 RMSE + 特征重要性（哪个因子最说了算）。
  B) 高斯过程(GP)仿真器：用很少的"昂贵"样本学一个慢函数的"快替身"，
     而且每次预测都自带"我有多确定"（误差带）——这正是 GP 比纯神经网络
     更适合做不确定性分析的原因。

依赖：numpy, scikit-learn
运行：python3 round6_exercise.py
"""
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as Cst


# ----------------------------------------------------------------------
# A) 随机森林：用 温度/降雨/施氮 预测产量
# ----------------------------------------------------------------------
def part_A(seed=0):
    rng = np.random.default_rng(seed)
    n = 1000
    temp = rng.uniform(15, 30, n)
    rain = rng.uniform(200, 800, n)
    Nrate = rng.uniform(0, 250, n)
    yield_ = (3 + 0.25 * (temp - 15) - 0.004 * (temp - 22) ** 2
              + 0.006 * rain - 3e-6 * rain ** 2
              + 0.02 * Nrate - 4e-5 * Nrate ** 2
              + rng.normal(0, 0.35, n))
    X = np.column_stack([temp, rain, Nrate])
    Xtr, Xte, ytr, yte = train_test_split(X, yield_, test_size=0.3, random_state=0)
    rf = RandomForestRegressor(n_estimators=400, random_state=0).fit(Xtr, ytr)
    rmse_tr = mean_squared_error(ytr, rf.predict(Xtr)) ** 0.5
    rmse_te = mean_squared_error(yte, rf.predict(Xte)) ** 0.5
    print("=" * 60)
    print("A) 随机森林产量预测")
    print("-" * 60)
    print(f"  训练集 RMSE = {rmse_tr:.3f} t/ha   测试集 RMSE = {rmse_te:.3f} t/ha")
    print("  特征重要性：", {n: round(v, 2) for n, v in
                       zip(["temp", "rain", "Nrate"], rf.feature_importances_)})
    print("  → 训练误差远小于测试误差是过拟合的信号；重要性告诉你哪个因子最说了算。")
    print('    但注意：RF 只在训练分布内可信，外推会"卡住"（见 Lab 3）。\n')


# ----------------------------------------------------------------------
# B) GP 仿真器：学一个"昂贵慢函数"的快替身，自带误差带
# ----------------------------------------------------------------------
def expensive(x):
    """假装这是一个跑一次很慢的过程模型。"""
    return np.sin(1.5 * x) + 0.3 * x


def part_B():
    rng = np.random.default_rng(1)
    # 只用 8 个"昂贵"样本训练
    Xd = np.sort(rng.uniform(0, 6, 8)).reshape(-1, 1)
    yd = expensive(Xd.ravel())
    gp = GaussianProcessRegressor(
        kernel=Cst(1.0) * RBF(1.0), normalize_y=True, alpha=1e-6,
        n_restarts_optimizer=2, random_state=0).fit(Xd, yd)
    print("=" * 60)
    print("B) GP 仿真器：8 个样本学慢函数，预测自带不确定性")
    print("-" * 60)
    print(f"{'x':>6}{'GP 预测':>12}{'±标准差':>10}{'真值':>10}{'落在带内?':>10}")
    for x in (0.5, 1.7, 3.0, 4.5, 5.8):
        mu, sd = gp.predict([[x]], return_std=True)
        true = expensive(x)
        inside = "是" if abs(true - mu[0]) <= 2 * sd[0] + 1e-9 else "否(外插)"
        print(f"{x:>6.1f}{mu[0]:>12.3f}{sd[0]:>10.3f}{true:>10.3f}{inside:>10}")
    print('  → 离训练点近 → 误差带窄、预测准；离得远 → 误差带变宽（GP 诚实地说"我不确定"）。')
    print('    把它当过程模型的"快替身"，就能毫秒级做上万次敏感性/标定。\n')
    print("★ 你的练习：把训练样本数从 8 调到 4 和 20，看误差带怎么收窄；")
    print("   把 part_A 的产量真值里加一个 rain×Nrate 交互项，看 RF 能否抓住它。")


if __name__ == "__main__":
    part_A()
    part_B()
