#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 4 练习 · 水文与生态系统：SCS-CN 产流 + Farquhar 光合 A=min(Ac,Aj)
======================================================================
对标章节：02_课程学习地图/Round4_水文与生态系统.md
框架：④ 尺度层级 / 能量与水量平衡；"概念水桶 vs 物理机理"

两个经典 worked example：
  A) SCS 曲线数法（SWAT 等用它把降雨分成"产流 vs 下渗"）——一个查表式经验公式。
  B) Farquhar 叶片光合：受 Rubisco 限制(Ac) 与受光/电子传递限制(Aj) 的交叉，
     净光合 A = min(Ac, Aj) − Rd。找出两条曲线的交叉点 Ci。

依赖：仅 numpy
运行：python3 round4_exercise.py
"""
import numpy as np


# ----------------------------------------------------------------------
# A) SCS-CN 产流：S = 25400/CN − 254 (mm)，Q = (P−0.2S)²/(P+0.8S)
# ----------------------------------------------------------------------
def scs_runoff(P, CN):
    S = 25400.0 / CN - 254.0           # 最大滞蓄量 (mm)
    Ia = 0.2 * S                       # 初损
    return np.where(P > Ia, (P - Ia) ** 2 / (P - Ia + S), 0.0), S


def part_A():
    print("=" * 62)
    print("A) SCS-CN 产流：降雨 50 mm 落在不同地表(CN)上产多少流")
    print("-" * 62)
    print(f"{'CN':>6}{'滞蓄S(mm)':>12}{'P=50→产流Q(mm)':>18}")
    for CN in (60, 75, 85, 95):
        Q, S = scs_runoff(50.0, CN)
        print(f"{CN:>6}{S:>12.1f}{float(Q):>18.1f}")
    print("→ CN 越大（地越硬/越湿/不透水），同样的雨产流越多、下渗越少。\n")


# ----------------------------------------------------------------------
# B) Farquhar FvCB：Ac(Rubisco限制) 与 Aj(光限制) 随 Ci 变化，A=min−Rd
# ----------------------------------------------------------------------
def farquhar(Ci, Vcmax=60.0, J=120.0, Gstar=42.0, Kc=404.9, Ko=278.4, O=210.0, Rd=1.0):
    Ac = Vcmax * (Ci - Gstar) / (Ci + Kc * (1 + O / Ko))   # 羧化(Rubisco)限制
    Aj = J * (Ci - Gstar) / (4 * Ci + 8 * Gstar)           # 电子传递(光)限制
    A = np.minimum(Ac, Aj) - Rd
    return Ac, Aj, A


def part_B():
    print("=" * 62)
    print("B) Farquhar 叶片光合：Ac vs Aj 谁限制（µmol m⁻² s⁻¹）")
    print("-" * 62)
    print(f"{'Ci(µmol/mol)':>14}{'Ac':>9}{'Aj':>9}{'净A=min−Rd':>13}{'限制者':>10}")
    Cis = np.array([100, 200, 300, 400, 600, 800], float)
    for Ci in Cis:
        Ac, Aj, A = farquhar(Ci)
        who = "Rubisco" if Ac < Aj else "光/电子"
        print(f"{Ci:>14.0f}{Ac:>9.2f}{Aj:>9.2f}{A:>13.2f}{who:>10}")
    # 细扫找交叉点（Ac≈Aj）
    grid = np.linspace(50, 1000, 20000)
    Ac, Aj, _ = farquhar(grid)
    ci_cross = grid[np.argmin(np.abs(Ac - Aj))]
    print("-" * 62)
    print(f"  交叉点 Ci ≈ {ci_cross:.0f} µmol/mol：低于它 Rubisco 限制，高于它 光限制。")
    print('  这就是 A–Ci 曲线"先陡升后压平"的机理来源。\n')


def your_turn():
    print("=" * 62)
    print("★ 你的练习")
    print("-" * 62)
    print("""  1) part_A 里把 P 从 50 改成 10 和 120，看小雨几乎不产流、大雨产流陡增。
  2) part_B 里把 J（光照强度的代理）从 120 降到 60（阴天），
     看交叉点 Ci 左移、光限制更早接管。""")


if __name__ == "__main__":
    part_A()
    part_B()
    your_turn()
