#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Round 3 练习 · 作物生长：30 行玩具作物模型（de Wit 潜在 vs 水分限制）
======================================================================
对标章节：02_课程学习地图/Round3_作物生长与产量.md
框架：② 潜在×胁迫（de Wit 生产层级）、③ 速率-状态-更新

这是 Lab 2(WOFOST) 的**纯 numpy 离线版**：不依赖 pcse、不联网，
但把同一套 first-principles 跑出来——
  GDD 物候 → Beer–Lambert 截光 → RUE 转干物质 → 水分胁迫打折 → ×HI → 产量。
跑两遍：潜在（不缺水）vs 水分限制（tipping-bucket 土壤水），看 de Wit 产量差。

依赖：仅 numpy
运行：python3 round3_exercise.py
"""
import numpy as np


def make_weather(days=160, seed=0):
    """造一季逐日天气（可复现）：温度、辐射、降雨、潜在蒸散。"""
    rng = np.random.default_rng(seed)
    t = np.arange(days)
    tmean = 18 + 6 * np.sin(2 * np.pi * (t - 20) / 200) + rng.normal(0, 1.5, days)
    rad = np.clip(16 + 6 * np.sin(2 * np.pi * (t - 20) / 200) + rng.normal(0, 2, days), 4, None)  # MJ/m2/d
    pet = np.clip(0.18 * rad, 1.0, None)                 # 粗略潜在蒸散 mm/d
    rain = rng.choice([0, 0, 0, 0, 6, 12, 25], size=days) * (rng.random(days) > 0.55)  # mm/d
    return dict(tmean=tmean, rad=rad, pet=pet, rain=rain)


def run_crop(wx, water_limited, p):
    """逐日跑一季，返回最终产量 (kg/ha)。water_limited=False 即潜在生产。"""
    gdd = dm = 0.0
    soil_water = p["WAV"]                                 # 可用土壤水 (mm)
    for i in range(len(wx["tmean"])):
        gdd += max(0.0, wx["tmean"][i] - p["tbase"])      # 攒积温
        lai = min(p["laimax"], p["laimax"] * gdd / p["gdd_mat"])
        fint = 1.0 - np.exp(-p["k"] * lai)                # Beer–Lambert 截光比例
        par = 0.5 * wx["rad"][i]                          # 一半总辐射是 PAR

        if water_limited:
            soil_water += wx["rain"][i]                   # 下雨进水
            demand = wx["pet"][i] * fint                  # 作物想蒸腾的水
            uptake = min(demand, soil_water)              # 实际吸到的
            soil_water -= uptake
            soil_water = min(soil_water, p["WHC"])        # 超过持水量的溢流走掉
            swfac = uptake / demand if demand > 1e-9 else 1.0
        else:
            swfac = 1.0                                   # 潜在：永不缺水

        dm += p["rue"] * fint * par * swfac               # 今天长的干物质
    return dm * p["HI"] * 10.0                             # ×收获指数, t/ha→kg/ha


def main():
    wx = make_weather()
    p = dict(tbase=0.0, gdd_mat=1600.0, laimax=5.0, k=0.6,
             rue=1.4, HI=0.45, WAV=40.0, WHC=120.0)
    pp = run_crop(wx, water_limited=False, p=p)
    wlp = run_crop(wx, water_limited=True, p=p)
    gap = pp - wlp
    print("=" * 60)
    print("Round 3 · 玩具作物模型：de Wit 潜在 vs 水分限制")
    print("-" * 60)
    print(f"  潜在产量 PP        ≈ {pp:7.0f} kg/ha")
    print(f"  水分限制产量 WLP   ≈ {wlp:7.0f} kg/ha")
    print(f"  产量差（可灌溉弥补）≈ {gap:7.0f} kg/ha = {100*gap/pp:.0f}%")
    print("-" * 60)
    print("→ 同一作物、同一天气，水分把潜在天花板压低的那一截 = de Wit 层级。")
    print("  （这是 Lab 2 WOFOST 的离线最小版：同一 first-principles，无需联网。）\n")
    print("★ 你的练习：把 p['WAV'] 初始土壤水从 40 调到 10，或给 make_weather")
    print("   换个干旱的 seed，看产量差怎么变大；再把 rue/HI 换成你作物的值。")


if __name__ == "__main__":
    main()
