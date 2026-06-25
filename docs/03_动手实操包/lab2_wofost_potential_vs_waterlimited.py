#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 2 · WOFOST 跑作物：用真实公开气象，演示 de Wit 潜在 vs 水分限制产量差
========================================================================
目的
  1) 在你自己选的地点（经纬度）上，用真正的公开气象数据 NASA POWER 跑 WOFOST；
  2) 同一地块跑两个"生产水平"：潜在产量(PP) 与 水分限制产量(WLP)，
     亲眼看到 de Wit 生产层级的产量差（= 灌溉理论上可弥补的部分）；
  3) 学会改地点/作物/播期——这就是把模型"指向你自己田块"的第一步。

依赖：pip install pcse  （首次运行会自动下载作物参数并联网取 NASA POWER 气象）
运行：python3 lab2_wofost_potential_vs_waterlimited.py
对标：第 3 轮教程（WOFOST/SUCROS、de Wit 生产层级）

说明：NASA POWER 是 NASA 的全球免费气象再分析数据，覆盖 1984→今、任意经纬度，
      是名副其实的"公开数据"。换经纬度即可换成你关心的地点。
"""

import warnings
warnings.filterwarnings("ignore")

import yaml
from pcse.input import (YAMLCropDataProvider, WOFOST72SiteDataProvider,
                        NASAPowerWeatherDataProvider, DummySoilDataProvider)
from pcse.base import ParameterProvider
from pcse.models import Wofost72_PP, Wofost72_WLP_FD


# ----------------------------------------------------------------------
# ★ 你只需改这几行，就能把模型指向你自己的地点/作物/季节
# ----------------------------------------------------------------------
LATITUDE = 37.4          # 纬度（示例：西班牙 Sevilla 附近，地中海气候、雨养小麦受旱）
LONGITUDE = -5.9         # 经度
CROP = "wheat"           # 作物（可选见运行输出的 crop list）
VARIETY = "Winter_wheat_101"
SOW_DATE = "2006-11-01"  # 播种日
END_DATE = "2007-07-01"  # 收获日
WAV = 20                 # 初始土壤有效含水量 (cm)，偏小→更易缺水→产量差更明显


def build_agromanagement():
    text = f"""
- 2006-01-01:
    CropCalendar:
        crop_name: {CROP}
        variety_name: {VARIETY}
        crop_start_date: {SOW_DATE}
        crop_start_type: sowing
        crop_end_date: {END_DATE}
        crop_end_type: harvest
        max_duration: 300
    TimedEvents: null
    StateEvents: null
"""
    return yaml.safe_load(text)


def main():
    print("=" * 64)
    print("Lab 2 · WOFOST：潜在 vs 水分限制产量（de Wit 生产层级）")
    print("-" * 64)

    # 1) 作物参数（首次联网下载并缓存）
    cropd = YAMLCropDataProvider()
    print("可选作物：", ", ".join(sorted(cropd.get_crops_varieties().keys())))
    cropd.set_active_crop(CROP, VARIETY)

    # 2) 土壤 + 站点 + 参数打包
    soild = DummySoilDataProvider()
    sited = WOFOST72SiteDataProvider(WAV=WAV)
    params = ParameterProvider(cropdata=cropd, soildata=soild, sitedata=sited)

    # 3) 真实公开气象：NASA POWER（按经纬度取，全球可用）
    print(f"\n取 NASA POWER 公开气象：lat={LATITUDE}, lon={LONGITUDE} ...")
    wdp = NASAPowerWeatherDataProvider(latitude=LATITUDE, longitude=LONGITUDE)
    print(f"气象可用区间：{wdp.first_date} → {wdp.last_date}")

    # 4) 农事管理
    agro = build_agromanagement()

    # 5) 跑两个生产水平
    out = {}
    for label, Model in [("潜在 PP", Wofost72_PP),
                         ("水分限制 WLP", Wofost72_WLP_FD)]:
        wof = Model(params, wdp, agro)
        wof.run_till_terminate()
        r = wof.get_summary_output()[0]
        out[label] = r
        print(f"\n{label}:")
        print(f"  产量 TWSO   = {r['TWSO']:.0f} kg/ha")
        print(f"  地上生物量 TAGP = {r['TAGP']:.0f} kg/ha")
        print(f"  最大叶面积 LAImax = {r['LAIMAX']:.2f}")

    # 6) 产量差 = de Wit 层级里"灌溉理论上可弥补"的部分
    gap = out["潜在 PP"]["TWSO"] - out["水分限制 WLP"]["TWSO"]
    pct = 100 * gap / out["潜在 PP"]["TWSO"] if out["潜在 PP"]["TWSO"] else 0
    print("\n" + "-" * 64)
    print(f"产量差（潜在 − 水分限制）= {gap:.0f} kg/ha = {pct:.0f}%")
    print("这就是 de Wit 生产层级：水分把'潜在天花板'压低的那一截。")
    print("把 LATITUDE/LONGITUDE 改成你的田块经纬度，就跑你自己的地点。")
    print("=" * 64)


if __name__ == "__main__":
    main()
