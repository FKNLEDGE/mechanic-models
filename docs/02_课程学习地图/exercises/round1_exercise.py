#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 1 练习 · 地基：守恒律 + 显式 Euler（Logistic 增长 + 一阶衰减）
====================================================================
对标章节：02_课程学习地图/Round1_地基_守恒律与最简模型.md
框架：① 状态变量+通量（守恒"浴缸"）、③ 速率-状态-更新（Euler 离散化）

这个练习让你亲手验证两件"地基级"的事：
  A) 计算机不是"解方程"，是"一步步往前挪"——但步子太大就会算错；
     把步长 dt 减半，数值解会越来越靠近精确解（收敛）。
  B) 一阶衰减的精确指数解恒为正、永不发散，且严格守恒。

依赖：仅 numpy
运行：python3 round1_exercise.py
"""
import numpy as np


# ----------------------------------------------------------------------
# A) Logistic 增长：dN/dt = r·N·(1 − N/K)
#    有解析解，正好用来当"标准答案"检验 Euler 的误差随步长怎么变。
# ----------------------------------------------------------------------
def logistic_exact(t, N0, r, K):
    return K / (1.0 + ((K - N0) / N0) * np.exp(-r * t))


def logistic_euler(N0, r, K, dt, T):
    n = int(round(T / dt))
    N = N0
    for _ in range(n):
        N = N + dt * (r * N * (1.0 - N / K))      # 算速率 → 积分 → 更新
    return N


def part_A():
    N0, r, K, T = 1.0, 1.0, 100.0, 10.0
    truth = logistic_exact(T, N0, r, K)
    print("=" * 62)
    print("A) 显式 Euler 收敛性：步长减半，误差应随之缩小")
    print("-" * 62)
    print(f"{'步长 dt':>10}{'Euler 解':>14}{'精确解':>12}{'误差':>12}")
    prev_err = None
    for dt in (2.0, 1.0, 0.5, 0.25, 0.125):
        approx = logistic_euler(N0, r, K, dt, T)
        err = abs(approx - truth)
        ratio = "" if prev_err is None else f"  (误差≈上行的 {err/prev_err:.2f} 倍)"
        print(f"{dt:>10.3f}{approx:>14.4f}{truth:>12.4f}{err:>12.4f}{ratio}")
        prev_err = err
    print("→ 步长越小越准；Euler 是一阶方法，dt 减半误差大致也减半。\n")


# ----------------------------------------------------------------------
# B) 一阶衰减 + 守恒检查：dC/dt = −k·C，精确步进 C(t+dt)=C·exp(−k·dt)
#    "跑掉的"碳应等于"减少的"碳——账要平。
# ----------------------------------------------------------------------
def part_B():
    C, k, dt, steps = 100.0, 0.1, 1.0, 30
    init = C
    cum_out = 0.0
    print("=" * 62)
    print("B) 一阶衰减的精确指数步进 + 质量守恒查账")
    print("-" * 62)
    for _ in range(steps):
        after = C * np.exp(-k * dt)     # 恒正，永不出负数
        cum_out += (C - after)          # 这步"分解跑掉"的碳
        C = after
    left = C + cum_out                  # 现在账上：剩下的 + 跑掉的
    print(f"  初始碳            = {init:.4f}")
    print(f"  {steps} 步后剩余    = {C:.4f}")
    print(f"  累计分解跑掉      = {cum_out:.4f}")
    print(f"  剩余 + 跑掉        = {left:.4f}  (应 = 初始 {init:.4f})")
    assert abs(left - init) < 1e-9, "质量不守恒！"
    print("  ✅ 守恒成立（差 < 1e-9）。指数步进天然不会出负数。\n")


def your_turn():
    print("=" * 62)
    print("★ 你的练习")
    print("-" * 62)
    print("""  1) 把 part_A 的精确解换成你心算的预期，验证 dt=0.125 已足够准。
  2) 在 part_B 里把 exp 精确步进换成朴素 Euler  C = C - k*C*dt，
     再把 dt 调到 > 1/k（例如 dt=12, k=0.1→不会，试 k=1,dt=1.5），
     看 C 何时变成负数——这就是第二部分讲的"负库"坑。""")


if __name__ == "__main__":
    part_A()
    part_B()
    your_turn()
