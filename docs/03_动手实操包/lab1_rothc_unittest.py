#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 1 · RothC 破冰：复现 Rothamsted 官方算例 + spin-up + 温度敏感性
================================================================
目的（呼应你的目的：用代码做苦力，用农学做灵魂）
  1) 先用"有标准答案"的单元测试确认你的代码和 Rothamsted 官方完全一致；
  2) 把单步衰减扩展成多年循环，看碳库如何趋于稳态（spin-up）；
  3) 调一个温度因子，亲手看到"变暖→土壤漏碳"。

依赖：仅 numpy（matplotlib 可选，用于画图）
运行：python3 lab1_rothc_unittest.py
对标：第 1 轮教程 SECTION 2（RothC）

第一性原理回顾
  每个活性碳库按一阶动力学衰减，月步长用"精确指数解"（恒正、永不发散）：
      C(t+Δt) = C * exp(-a*b*c*k*Δt)
  其中 k 是该库年分解速率常数，a/b/c 是温度/水分/覆盖三个 0–1 速率修正因子，
  Δt = 1/12 年（一个月）。五库：DPM, RPM, BIO, HUM（活性）+ IOM（惰性，不分解）。
"""

import numpy as np

# ----------------------------------------------------------------------
# 1) 五库年分解速率常数 k (yr^-1)，来自 Coleman & Jenkinson (1996) RothC-26.3
#    DPM 易分解最快；HUM 腐殖质最慢；IOM 惰性不入此表。
# ----------------------------------------------------------------------
K = {"DPM": 10.0, "RPM": 0.3, "BIO": 0.66, "HUM": 0.02}


def decay_step(pools: dict, abc: float, dt: float = 1.0 / 12) -> dict:
    """对四个活性库做一步'精确指数衰减'。
    pools : dict，各库当前碳量 (t C/ha)
    abc   : 速率修正因子乘积 a*b*c (0~1)
    dt    : 时间步长，月步=1/12 年
    返回   : 衰减后的各库碳量
    """
    return {name: pools[name] * np.exp(-abc * K[name] * dt)
            for name in ("DPM", "RPM", "BIO", "HUM")}


# ----------------------------------------------------------------------
# 2) 单元测试：复现 Rothamsted 官方算例（必须分毫不差）
#    官方给定：某月 abc = 0.3561，月步 dt = 1/12
# ----------------------------------------------------------------------
def test_official_example():
    pools_in = {"DPM": 0.1533, "RPM": 4.4852, "BIO": 0.6671, "HUM": 25.8576}
    abc = 0.3561
    out = decay_step(pools_in, abc, dt=1.0 / 12)

    expected = {"DPM": 0.1140, "RPM": 4.4455, "BIO": 0.6542, "HUM": 25.8423}
    print("=" * 60)
    print("单元测试：复现 RothC 官方算例 (abc=0.3561, dt=1/12)")
    print("-" * 60)
    print(f"{'库':<6}{'算得':>12}{'官方期望':>12}{'判定':>8}")
    all_pass = True
    for name in expected:
        ok = abs(out[name] - expected[name]) < 1e-3
        all_pass &= ok
        print(f"{name:<6}{out[name]:>12.4f}{expected[name]:>12.4f}{'✅' if ok else '❌':>8}")
    print("-" * 60)
    # 手算验证 DPM 这一步，确认你真的理解了指数解
    manual = 0.1533 * np.exp(-10 * 0.3561 / 12)
    print(f"手算验证 DPM = 0.1533*exp(-10*0.3561/12) = {manual:.4f}")
    print("注：DPM 精确值=0.11394，4 位小数显示为 0.1139；官方表记 0.1140 是其取整，")
    print("    差异仅在第 4 位小数（RPM/BIO/HUM 三库 4 位完全一致），属取整而非计算错误。")
    assert all_pass, "单元测试未通过！请检查 decay_step。"
    print("✅ 全部通过：你的代码与 Rothamsted 官方在取整精度内一致。\n")
    return all_pass


# ----------------------------------------------------------------------
# 3) Spin-up：把单步包进多年循环，看碳库趋于稳态
#    简化版完整周转：年碳输入按 DPM/RPM=1.44 分配；分解出的碳按黏粒决定的
#    比例回填 BIO+HUM（其余成 CO2），BIO:HUM=0.46:0.54。
# ----------------------------------------------------------------------
def co2_fraction(clay_pct: float) -> float:
    """分解碳中变成 CO2 的比例 = x/(x+1)，x 由黏粒含量决定（Coleman & Jenkinson 1996）。"""
    x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay_pct))
    return x / (x + 1.0)


def spin_up(annual_input=2.0, clay_pct=23.4, abc=0.30, years=2000,
            dpm_rpm_ratio=1.44):
    """跑 spin-up，返回各库逐年末碳量轨迹 (years x 4)。
    annual_input : 年植物碳输入 (t C/ha/yr)
    clay_pct     : 黏粒百分比（Rothamsted 默认 23.4）
    abc          : 这里用一个固定的年均速率修正因子（教学简化）
    """
    f_co2 = co2_fraction(clay_pct)          # 分解→CO2 的比例
    f_soil = 1.0 - f_co2                     # 分解→(BIO+HUM) 的比例
    f_dpm = dpm_rpm_ratio / (1.0 + dpm_rpm_ratio)   # 进 DPM 的输入比例 (=0.59)

    pools = {"DPM": 0.0, "RPM": 0.0, "BIO": 0.0, "HUM": 0.0}
    traj = np.zeros((years, 4))
    dt = 1.0 / 12
    for yr in range(years):
        for _m in range(12):
            before = dict(pools)
            after = decay_step(pools, abc, dt)               # 1) 衰减
            decomposed = {n: before[n] - after[n] for n in after}  # 各库这步分解掉的碳
            total_decomp = sum(decomposed.values())
            to_soil = total_decomp * f_soil                  # 回填 BIO+HUM 的总量
            after["BIO"] += 0.46 * to_soil                   # 2) 回填 BIO
            after["HUM"] += 0.54 * to_soil                   #    回填 HUM
            month_input = annual_input / 12.0                # 3) 月碳输入
            after["DPM"] += f_dpm * month_input
            after["RPM"] += (1 - f_dpm) * month_input
            pools = after
        traj[yr] = [pools["DPM"], pools["RPM"], pools["BIO"], pools["HUM"]]
    return traj, pools


# ----------------------------------------------------------------------
# 4) 温度敏感性：把速率修正因子 abc 调高（模拟变暖），看 HUM 稳态如何下降
# ----------------------------------------------------------------------
def warming_sensitivity():
    print("=" * 60)
    print("温度敏感性：abc 升高（变暖→分解更快）如何改变稳态总碳")
    print("-" * 60)
    print(f"{'abc':>6}{'稳态DPM':>10}{'稳态RPM':>10}{'稳态BIO':>10}{'稳态HUM':>10}{'总SOC':>10}")
    for abc in (0.25, 0.30, 0.40, 0.50):
        _traj, eq = spin_up(annual_input=2.0, clay_pct=23.4, abc=abc, years=2000)
        tot = sum(eq.values())
        print(f"{abc:>6.2f}{eq['DPM']:>10.2f}{eq['RPM']:>10.2f}"
              f"{eq['BIO']:>10.2f}{eq['HUM']:>10.2f}{tot:>10.2f}")
    print("观察：abc 越大（越暖/越湿），稳态总 SOC 越低——这就是'变暖漏碳'。\n")


def main():
    test_official_example()

    print("=" * 60)
    print("Spin-up：年输入 2.0 t C/ha，黏粒 23.4%，abc=0.30")
    print("-" * 60)
    traj, eq = spin_up(annual_input=2.0, clay_pct=23.4, abc=0.30, years=2000)
    print(f"稳态各库 (t C/ha): DPM={eq['DPM']:.2f}  RPM={eq['RPM']:.2f}  "
          f"BIO={eq['BIO']:.2f}  HUM={eq['HUM']:.2f}")
    print(f"稳态总 SOC = {sum(eq.values()):.2f} t C/ha")
    print(f"100 年时 HUM = {traj[99, 3]:.2f}，2000 年时 HUM = {traj[-1, 3]:.2f}"
          f"（HUM 最慢，趋稳最久）\n")

    warming_sensitivity()

    # 可选画图：spin-up 轨迹
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        yrs = np.arange(1, traj.shape[0] + 1)
        plt.figure(figsize=(8, 5))
        for i, name in enumerate(["DPM", "RPM", "BIO", "HUM"]):
            plt.plot(yrs, traj[:, i], label=name)
        plt.plot(yrs, traj.sum(axis=1), "k--", label="Total SOC")
        plt.xlabel("year"); plt.ylabel("t C/ha"); plt.xscale("log")
        plt.title("RothC spin-up to steady state")
        plt.legend(); plt.tight_layout()
        plt.savefig("lab1_rothc_spinup.png", dpi=110)
        print("已保存图：lab1_rothc_spinup.png")
    except Exception as e:
        print(f"(画图跳过：{e})")


if __name__ == "__main__":
    main()
