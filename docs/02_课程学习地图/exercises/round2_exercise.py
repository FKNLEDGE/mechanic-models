#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 2 练习 · 土壤生物地球化学：Michaelis–Menten + 线性 vs 非线性稳态
======================================================================
对标章节：02_课程学习地图/Round2_土壤生物地球化学.md
框架：② 多重限制因子、⑤ 动力系统稳定性（线性唯一稳态 vs 微生物显式非线性）

两件事：
  A) Michaelis–Menten 速率 V = Vmax·S/(Km+S) 的两个极限——
     底物少时像"一阶"(∝S)，底物多时像"零阶"(饱和=Vmax)。
  B) "输入翻倍，稳态碳是否翻倍?"——
     线性一阶库模型：严格翻倍（叠加性）；
     微生物显式（M–M 反馈）模型：不翻倍（非线性）。这正是 Round 2 的核心跃迁。

依赖：仅 numpy
运行：python3 round2_exercise.py
"""
import numpy as np


# ----------------------------------------------------------------------
# A) Michaelis–Menten 的两个极限
# ----------------------------------------------------------------------
def mm(S, Vmax=10.0, Km=2.0):
    return Vmax * S / (Km + S)


def part_A():
    print("=" * 64)
    print("A) Michaelis–Menten：V = Vmax·S/(Km+S)，Vmax=10, Km=2")
    print("-" * 64)
    print(f"{'底物 S':>8}{'M-M 速率':>12}{'一阶近似 Vmax/Km·S':>20}{'占 Vmax':>10}")
    for S in (0.1, 0.5, 2.0, 10.0, 100.0):
        v = mm(S)
        lin = (10.0 / 2.0) * S
        print(f"{S:>8.1f}{v:>12.3f}{lin:>20.3f}{v/10.0*100:>9.0f}%")
    print("→ S≪Km 时 ≈ 一阶(∝S)；S≫Km 时饱和到 Vmax(零阶)。\n")


# ----------------------------------------------------------------------
# B) 线性一阶库 vs 微生物显式(M–M) —— 输入翻倍，稳态翻倍吗?
# ----------------------------------------------------------------------
def linear_pool_steady(I, k=0.1):
    """dC/dt = I − k·C  → 稳态 C* = I/k（与 I 成严格正比）。"""
    return I / k


def microbial_explicit_steady(I, years=4000):
    """极简微生物显式：
        dSOC/dt = I − Vmax·MIC·SOC/(Km+SOC)
        dMIC/dt = CUE·Vmax·MIC·SOC/(Km+SOC) − m·MIC
    用 Euler 跑到稳态，返回 (SOC*, MIC*)。非线性 → 稳态不随 I 线性。"""
    Vmax, Km, CUE, m = 0.8, 30.0, 0.4, 0.3
    SOC, MIC, dt = 50.0, 1.0, 0.05
    for _ in range(int(years / dt)):
        uptake = Vmax * MIC * SOC / (Km + SOC)
        SOC += dt * (I - uptake)
        MIC += dt * (CUE * uptake - m * MIC)
        SOC = max(SOC, 1e-9); MIC = max(MIC, 1e-9)
    return SOC, MIC


def part_B():
    print("=" * 64)
    print("B) 输入翻倍，稳态 SOC 是否翻倍?（线性叠加性 vs 非线性反馈）")
    print("-" * 64)
    I = 2.0
    c1 = linear_pool_steady(I)
    c2 = linear_pool_steady(2 * I)
    print(f"  线性一阶库   : I→稳态 {c1:.2f};  2I→稳态 {c2:.2f};  比值 = {c2/c1:.2f}  (应=2.00)")
    s1, _ = microbial_explicit_steady(I)
    s2, _ = microbial_explicit_steady(2 * I)
    print(f"  微生物显式   : I→稳态 {s1:.2f};  2I→稳态 {s2:.2f};  比值 = {s2/s1:.2f}  (≠2 → 非线性)")
    print("→ 线性模型严格翻倍（这也意味着它只有唯一稳态、不会自发多稳态）；")
    print("  微生物显式因 M–M 反馈，稳态对输入是非线性响应——这是 Round 2 的关键跃迁。\n")


def your_turn():
    print("=" * 64)
    print("★ 你的练习")
    print("-" * 64)
    print("""  1) 改 microbial_explicit_steady 里的 CUE（碳利用效率）从 0.4→0.6，
     看稳态 SOC/MIC 怎么变（CUE 越高，微生物越能"囤碳"）。
  2) 把 part_A 的 Km 调大到 20，重看哪一段还像"一阶"——
     Km 越大，越晚饱和。""")


if __name__ == "__main__":
    part_A()
    part_B()
    your_turn()
