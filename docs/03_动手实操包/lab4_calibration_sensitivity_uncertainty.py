#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 4 · 标定 + 敏感性 + 不确定性：从"会跑模型"到"会做科研"
============================================================
目的（你目的里"用农学做灵魂"最吃功夫的一关）
  把 Lab 1 的 RothC 当被研究对象，完整走一遍科研工作流：
    A) 敏感性分析（Sobol，SALib）：先搞清"哪些参数才值得调"；
    B) 标定（calibration）：用观测把关键参数拟合出来（反问题）；
    C) GLUE 不确定性：看到 equifinality，并把参数不确定性传播成预测的"误差带"；
    D) GP 仿真器（emulator）：用快替身加速昂贵模型的评估。

依赖：numpy, scipy, scikit-learn, SALib（matplotlib 可选）
  pip install SALib scikit-learn scipy
运行：python3 lab4_calibration_sensitivity_uncertainty.py
对标：第 6 轮教程 SECTION 1（标定/equifinality/Bayesian）+ SECTION 3（GP emulator）

★ 这是把课程框架⑥（标定—验证—不确定性）真正落到可运行代码的一关。
"""

import warnings; warnings.filterwarnings("ignore")
import numpy as np

# ----------------------------------------------------------------------
# RothC 稳态 SOC（向量化版，衰减因子恒定→预计算一次，快 ~10x）
# 参数：annual_input 年碳输入(t C/ha/yr)，clay 黏粒%，abc 速率修正，ratio=DPM/RPM
# ----------------------------------------------------------------------
K = np.array([10.0, 0.3, 0.66, 0.02])      # DPM, RPM, BIO, HUM 年分解速率 (yr^-1)


def _co2_fraction(clay_pct):
    x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay_pct))
    return x / (x + 1.0)


def rothc_soc(annual_input, clay_pct, abc, ratio, years=600, return_traj=False):
    """跑 spin-up，返回稳态总 SOC（或逐年总 SOC 轨迹）。"""
    f_soil = 1.0 - _co2_fraction(clay_pct)
    f_dpm = ratio / (1.0 + ratio)
    dt = 1.0 / 12
    decay = np.exp(-abc * K * dt)                       # 恒定月衰减因子
    add = np.array([f_dpm, 1 - f_dpm, 0.0, 0.0]) * (annual_input / 12.0)
    p = np.zeros(4)
    traj = np.empty(years)
    for yr in range(years):
        for _m in range(12):
            a = p * decay
            dec = (p - a).sum() * f_soil
            a[2] += 0.46 * dec
            a[3] += 0.54 * dec
            a = a + add
            p = a
        traj[yr] = p.sum()
    return (traj if return_traj else p.sum())


# ======================================================================
# A) 敏感性分析（Sobol）：哪些参数主导稳态 SOC？
# ======================================================================
def part_A_sensitivity():
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol as sobol_analyze
    problem = {
        "num_vars": 4,
        "names": ["input", "clay", "abc", "ratio"],
        "bounds": [[1.0, 4.0], [5, 50], [0.2, 0.5], [0.25, 2.0]],
    }
    X = sobol_sample.sample(problem, 128, seed=42)   # 固定种子→可复现；128*(2*4+2)=1280 次运行
    Y = np.array([rothc_soc(*x) for x in X])
    Si = sobol_analyze.analyze(problem, Y, seed=42, print_to_console=False)
    print("=" * 60)
    print("A) Sobol 敏感性：稳态 SOC 对各参数的敏感度")
    print("-" * 60)
    print(f"{'参数':<8}{'一阶 S1':>10}{'总效应 ST':>12}")
    for n, s1, st in zip(problem["names"], Si["S1"], Si["ST"]):
        print(f"{n:<8}{s1:>10.3f}{st:>12.3f}")
    print("解读：S1/ST 越大越该重点标定。这里碳输入与 abc 主导，")
    print("      clay 次要、DPM/RPM 比几乎无影响→标定时可固定 ratio。\n")


# ======================================================================
# B) 标定：用"观测"反推关键参数（反问题）
#    造一条 synthetic 观测（真值 input=2.5, abc=0.28，固定 clay=24, ratio=1.44），
#    用 scipy least_squares 拟合 (input, abc)。
# ======================================================================
TRUE = dict(annual_input=2.5, clay=24.0, abc=0.28, ratio=1.44)
OBS_YEARS = np.array([20, 40, 60, 80, 100, 150])      # 有观测的年份


def _make_observations(seed=1):
    traj = rothc_soc(TRUE["annual_input"], TRUE["clay"], TRUE["abc"],
                     TRUE["ratio"], years=200, return_traj=True)
    rng = np.random.default_rng(seed)
    obs = traj[OBS_YEARS - 1] + rng.normal(0, 1.5, OBS_YEARS.size)   # 加观测噪声
    return obs


def part_B_calibration(obs):
    from scipy.optimize import least_squares

    def residuals(theta):
        ai, abc = theta
        traj = rothc_soc(ai, TRUE["clay"], abc, TRUE["ratio"],
                         years=200, return_traj=True)
        return traj[OBS_YEARS - 1] - obs

    res = least_squares(residuals, x0=[3.0, 0.35],
                        bounds=([1.0, 0.2], [4.0, 0.5]))
    ai_hat, abc_hat = res.x
    print("=" * 60)
    print("B) 标定（反问题）：用观测拟合 (annual_input, abc)")
    print("-" * 60)
    print(f"  真值   : input={TRUE['annual_input']:.2f}, abc={TRUE['abc']:.3f}")
    print(f"  标定得 : input={ai_hat:.2f}, abc={abc_hat:.3f}")
    print(f"  最终 RMSE = {np.sqrt(np.mean(res.fun**2)):.2f} t C/ha")
    print("  注意：标定值≈真值，但单点拟合不告诉你'有多确定'——看 C 部分。\n")
    return ai_hat, abc_hat


# ======================================================================
# C) GLUE：Monte Carlo 采样→保留'行为合格'参数集→看 equifinality + 误差带
# ======================================================================
def part_C_glue(obs, n_samples=1500, seed=2):
    rng = np.random.default_rng(seed)
    ais = rng.uniform(1.0, 4.0, n_samples)
    abcs = rng.uniform(0.2, 0.5, n_samples)
    sigma = 1.5
    rmses = np.empty(n_samples)
    for i in range(n_samples):
        traj = rothc_soc(ais[i], TRUE["clay"], abcs[i], TRUE["ratio"],
                         years=200, return_traj=True)
        rmses[i] = np.sqrt(np.mean((traj[OBS_YEARS - 1] - obs) ** 2))
    # 似然（高斯）→ 取最优 10% 为'行为合格(behavioral)'集
    thr = np.percentile(rmses, 10)
    beh = rmses <= thr
    print("=" * 60)
    print(f"C) GLUE：{n_samples} 组随机参数，保留最优 10% 为'行为合格'集")
    print("-" * 60)
    print(f"  行为合格集数量: {beh.sum()}")
    print(f"  input 行为区间 : [{ais[beh].min():.2f}, {ais[beh].max():.2f}]  "
          f"(真值 {TRUE['annual_input']:.2f})")
    print(f"  abc   行为区间 : [{abcs[beh].min():.3f}, {abcs[beh].max():.3f}]  "
          f"(真值 {TRUE['abc']:.3f})")
    print("  → 多组不同 (input,abc) 拟合都'合格' = equifinality（殊途同归）。\n")

    # 把参数不确定性传播成"未来 100 年情景"的预测误差带
    # 情景：碳输入提高 20%（如增施有机肥），看行为合格集给出的 SOC 增量分布
    fut = []
    for ai, abc in zip(ais[beh], abcs[beh]):
        base = rothc_soc(ai, TRUE["clay"], abc, TRUE["ratio"], years=200)
        scen = rothc_soc(ai * 1.2, TRUE["clay"], abc, TRUE["ratio"], years=200)
        fut.append(scen - base)
    fut = np.array(fut)
    print("  预测不确定性传播（情景：碳输入+20%，100→200 年稳态 SOC 增量）：")
    print(f"    中位数 = {np.median(fut):.2f} t C/ha")
    print(f"    5–95% 区间 = [{np.percentile(fut,5):.2f}, {np.percentile(fut,95):.2f}] t C/ha")
    print("  → 这才是诚实的结论：给区间，不是给单点。\n")
    return ais, abcs, rmses, beh


# ======================================================================
# D) GP 仿真器：训练快替身模仿 RothC(input,abc)->SOC，加速评估并自带不确定性
# ======================================================================
def part_D_emulator():
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
    rng = np.random.default_rng(0)
    # 设计点：30 组 (input, abc) 上跑真模型（"昂贵"的几十次）
    Xd = np.column_stack([rng.uniform(1, 4, 30), rng.uniform(0.2, 0.5, 30)])
    yd = np.array([rothc_soc(x[0], TRUE["clay"], x[1], TRUE["ratio"]) for x in Xd])
    gp = GaussianProcessRegressor(
        kernel=C(1.0) * RBF([1.0, 0.1]), normalize_y=True, alpha=1e-6
    ).fit(Xd, yd)

    print("=" * 60)
    print("D) GP 仿真器：用 30 次真模型训练'快替身'，之后预测近乎瞬时")
    print("-" * 60)
    test = np.array([[2.5, 0.28], [3.5, 0.45]])
    mu, sd = gp.predict(test, return_std=True)
    for (ai, abc), m, s in zip(test, mu, sd):
        true_v = rothc_soc(ai, TRUE["clay"], abc, TRUE["ratio"])
        print(f"  (input={ai}, abc={abc}): 仿真器 {m:.1f}±{s:.1f}  真模型 {true_v:.1f}")
    print("  → 仿真器可在毫秒内评估上万次，用于敏感性/标定/不确定性，自带误差估计。\n")


def main():
    part_A_sensitivity()
    obs = _make_observations()
    part_B_calibration(obs)
    part_C_glue(obs)
    part_D_emulator()
    print("✅ Lab 4 跑通。你已把'敏感性→标定→不确定性→仿真器'整套科研工作流跑了一遍。")
    print("   下一步：把 rothc_soc 换成你的模型（或 Lab 2 的 WOFOST），把 obs 换成你的实测。")


if __name__ == "__main__":
    main()
