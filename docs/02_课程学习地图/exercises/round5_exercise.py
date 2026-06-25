#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 5 练习 · 动力系统与临界转变：SIR 阈值 + Scheffer 浅湖迟滞 + 早期预警
==========================================================================
对标章节：02_课程学习地图/Round5_动力系统与临界转变.md
框架：⑤ 动力系统稳定性（阈值、分岔、迟滞、临界慢化）

三件事：
  A) SIR：基本再生数 R0=β/γ 决定疫情是否暴发；群体免疫阈值 = 1−1/R0。
  B) Scheffer 浅湖：缓慢加营养盐 → 清水突然变浊；再降回去，要降得更低才变清
     （迟滞 / 替代稳态）。打印"去程"和"回程"两个临界阈值。
  C) 早期预警：逼近临界点时，扰动恢复变慢 → 滞后-1 自相关与方差上升。

依赖：仅 numpy
运行：python3 round5_exercise.py
"""
import numpy as np


# ----------------------------------------------------------------------
# A) SIR（Euler 积分）
# ----------------------------------------------------------------------
def sir(beta, gamma, S0=0.999, I0=0.001, dt=0.05, T=160):
    S, I, R = S0, I0, 0.0
    peakI = I
    for _ in range(int(T / dt)):
        dS = -beta * S * I
        dI = beta * S * I - gamma * I
        S += dt * dS; I += dt * dI; R += dt * gamma * I
        peakI = max(peakI, I)
    return peakI, R


def part_A():
    print("=" * 64)
    print("A) SIR：R0=β/γ 决定暴发，群体免疫阈值=1−1/R0")
    print("-" * 64)
    gamma = 0.1
    print(f"{'β':>6}{'R0':>6}{'峰值感染比例':>14}{'最终感染过':>12}{'群体免疫阈值':>14}")
    for beta in (0.08, 0.15, 0.25, 0.4):
        R0 = beta / gamma
        peakI, Rfinal = sir(beta, gamma)
        hit = max(0.0, 1 - 1 / R0)
        print(f"{beta:>6.2f}{R0:>6.1f}{peakI:>14.3f}{Rfinal:>12.3f}{hit:>13.0%}")
    print("→ R0<1 不暴发；R0 越大，峰值越高、需接种到 1−1/R0 才能挡住。\n")


# ----------------------------------------------------------------------
# B) Scheffer 浅湖：dx/dt = a − b·x + r·x^q/(x^q+m^q)（x=藻类浊度/营养）
#    缓慢上调 / 下调营养输入 a，稳态走两条不同分支 → 迟滞回线
# ----------------------------------------------------------------------
def lake_equilibrium(a, x0, b=0.6, r=1.0, m=1.0, q=4, dt=0.05, steps=6000):
    x = x0
    for _ in range(steps):
        x += dt * (a - b * x + r * x**q / (x**q + m**q))
        x = max(x, 0.0)
    return x


def part_B():
    print("=" * 64)
    print("B) Scheffer 浅湖迟滞：去程(加营养)与回程(减营养)阈值不同")
    print("-" * 64)
    b, r, m, q = 0.6, 1.0, 1.0, 4
    # 折点 = 平衡驱动曲线 g(x)=b·x − r·x^q/(x^q+m^q) 的局部极大(去程)与极小(回程)
    P = np.linspace(0.001, 4.0, 200000)
    g = b * P - r * P**q / (P**q + m**q)
    turns = np.where(np.diff(np.sign(np.diff(g))) != 0)[0] + 1
    a_up, a_down = g[turns].max(), g[turns].min()
    # 在两阈值之间取一个 a，验证"从清水出发"与"从浊水出发"落到不同稳态(双稳)
    a_mid = 0.5 * (a_up + a_down)
    x_clear = lake_equilibrium(a_mid, 0.0)     # 从清水(低浊度)出发
    x_turbid = lake_equilibrium(a_mid, 3.0)    # 从浊水(高浊度)出发
    print(f"  去程：营养 a 升过 ≈ {a_up:.3f} → 清水突然崩成浊水")
    print(f"  回程：营养 a 必须降到 ≈ {a_down:.3f} 才由浊变清（低得多）")
    print(f"  迟滞宽度 ≈ {a_up - a_down:.3f}")
    print(f"  同一 a={a_mid:.3f}：从清水出发→浊度 {x_clear:.2f}（清），"
          f"从浊水出发→浊度 {x_turbid:.2f}（浊）")
    print('  → 同一营养水平、两种命运，全看历史 = 替代稳态 / 迟滞。\n')


# ----------------------------------------------------------------------
# C) 早期预警：逼近 fold 时恢复变慢 → AR1 与方差上升
# ----------------------------------------------------------------------
def part_C():
    print("=" * 64)
    print("C) 早期预警信号：逼近临界点，滞后-1 自相关与方差上升")
    print("-" * 64)
    rng = np.random.default_rng(0)
    n = 1200
    # 让"恢复速率" lam 随时间从 0.6 慢慢趋近 0（越来越接近临界慢化）
    lam = np.linspace(0.6, 0.03, n)
    x = np.zeros(n)
    for t in range(1, n):
        x[t] = (1 - lam[t]) * x[t - 1] + rng.normal(0, 0.1)   # AR(1)，系数=1−lam
    win = 200
    print(f"{'时段':>10}{'方差':>12}{'滞后-1自相关':>14}")
    for start in (0, 400, 800, n - win):
        seg = x[start:start + win]
        seg = seg - seg.mean()
        ar1 = np.corrcoef(seg[:-1], seg[1:])[0, 1]
        tag = "早" if start == 0 else ("晚(临界前)" if start >= n - win else "中")
        print(f"{tag:>10}{np.var(seg):>12.4f}{ar1:>14.3f}")
    print('→ 越接近临界点，方差和自相关越往上爬——这就是可观测的"预警"。\n')
    print("★ 你的练习：把 part_B 的 r 调小到 0.6（反馈变弱），看迟滞回线变窄甚至消失；")
    print("   把 part_C 的噪声 0.1 调大，看预警信号是否被噪声淹没。")


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
