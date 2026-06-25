#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成教程配图（PNG）到本目录 docs/img/。
可复现：python3 make_figures.py    依赖：numpy, matplotlib
图里用中文标签（需要 WenQuanYi Zen Hei 或任一 CJK 字体；缺字体会退化为方框，但不影响数据）。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 中文字体（环境里有 WenQuanYi Zen Hei）
for cand in ("WenQuanYi Zen Hei", "Noto Sans CJK SC", "Source Han Sans SC", "SimHei"):
    try:
        plt.rcParams["font.sans-serif"] = [cand]; break
    except Exception:
        pass
plt.rcParams["axes.unicode_minus"] = False
HERE = os.path.dirname(os.path.abspath(__file__))


def save(fig, name):
    p = os.path.join(HERE, name)
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# 1) RothC spin-up：五库随时间趋于稳态 -----------------------------------
def fig_rothc_spinup():
    K = np.array([10.0, 0.3, 0.66, 0.02])       # DPM,RPM,BIO,HUM
    abc, clay, ann_in, ratio = 0.30, 23.4, 2.0, 1.44
    x = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * clay)); f_soil = 1 - x/(x+1)
    f_dpm = ratio/(1+ratio); dt = 1/12
    decay = np.exp(-abc*K*dt); add = np.array([f_dpm, 1-f_dpm, 0, 0])*(ann_in/12)
    p = np.zeros(4); years = 2000; traj = np.empty((years, 4))
    for yr in range(years):
        for _m in range(12):
            a = p*decay; dec = (p-a).sum()*f_soil
            a[2] += 0.46*dec; a[3] += 0.54*dec; p = a+add
        traj[yr] = p
    yrs = np.arange(1, years+1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, n in enumerate(["DPM 易分解", "RPM 难分解", "BIO 微生物", "HUM 腐殖质"]):
        ax.plot(yrs, traj[:, i], label=n)
    ax.plot(yrs, traj.sum(1), "k--", label="总 SOC")
    ax.set_xscale("log"); ax.set_xlabel("年（对数轴）"); ax.set_ylabel("碳量 (t C/ha)")
    ax.set_title('RothC spin-up：五个"浴缸"如何各自趋于稳态')
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    save(fig, "rothc_spinup.png")


# 2) 精确指数解 vs 显式 Euler（步长太大会出错甚至变负）-------------------
def fig_exp_decay_euler():
    k = 1.0; C0 = 100; T = 6
    t = np.linspace(0, T, 400)
    exact = C0*np.exp(-k*t)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t, exact, "k", lw=2, label="精确解 C0·e^(-kt)")
    for dt, style in [(1.5, "o--"), (0.5, "s--")]:
        n = int(T/dt); tt = np.arange(n+1)*dt; C = C0; ys = [C]
        for _ in range(n):
            C = C - k*C*dt; ys.append(C)      # 朴素 Euler（故意展示其缺陷）
        ax.plot(tt, ys, style, ms=4, label=f"Euler dt={dt}")
    ax.axhline(0, color="r", lw=.8, alpha=.6)
    ax.set_xlabel("时间"); ax.set_ylabel("剩余量 C")
    ax.set_title("一步步往前挪：步子太大（dt=1.5）会冲过头甚至算出负数")
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "exp_decay_euler.png")


# 3) de Wit 产量层级：潜在 vs 水分限制 -----------------------------------
def fig_dewit_gap():
    rng = np.random.default_rng(0); days = 160
    t = np.arange(days)
    rad = np.clip(16+6*np.sin(2*np.pi*(t-20)/200)+rng.normal(0, 2, days), 4, None)
    tmean = 18+6*np.sin(2*np.pi*(t-20)/200)+rng.normal(0, 1.5, days)
    pet = np.clip(0.18*rad, 1, None)
    rain = rng.choice([0, 0, 0, 0, 6, 12, 25], days)*(rng.random(days) > 0.55)
    p = dict(tbase=0, gdd_mat=1600, laimax=5, k=0.6, rue=1.4)

    def run(wl):
        gdd = dm = 0; sw = 40.0; out = []
        for i in range(days):
            gdd += max(0, tmean[i]-p["tbase"])
            lai = min(p["laimax"], p["laimax"]*gdd/p["gdd_mat"])
            fint = 1-np.exp(-p["k"]*lai); par = 0.5*rad[i]
            if wl:
                sw += rain[i]; demand = pet[i]*fint; up = min(demand, sw)
                sw = min(sw-up, 120); swf = up/demand if demand > 1e-9 else 1
            else:
                swf = 1.0
            dm += p["rue"]*fint*par*swf; out.append(dm)
        return np.array(out)
    pp, wlp = run(False), run(True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t, pp, label="潜在生产 PP（不缺水）")
    ax.plot(t, wlp, label="水分限制 WLP")
    ax.fill_between(t, wlp, pp, alpha=.2, color="tab:red", label="产量差（灌溉可弥补）")
    ax.set_xlabel("播种后天数"); ax.set_ylabel("累积干物质 (g/m²)")
    ax.set_title('de Wit 生产层级：水分把"潜在天花板"压低的那一截')
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "dewit_yield_gap.png")


# 4) Michaelis–Menten：一阶/零阶两个极限 --------------------------------
def fig_michaelis_menten():
    S = np.linspace(0, 60, 400); Vmax, Km = 10, 8
    V = Vmax*S/(Km+S)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(S, V, "b", lw=2, label="M–M：V=Vmax·S/(Km+S)")
    ax.plot(S, (Vmax/Km)*S, "g--", label="底物少→一阶 (∝S)")
    ax.axhline(Vmax, color="r", ls=":", label="底物多→零阶 (饱和=Vmax)")
    ax.axvline(Km, color="gray", ls=":"); ax.text(Km+1, 1, "S=Km 时\nV=½Vmax", fontsize=8)
    ax.set_ylim(0, Vmax*1.15); ax.set_xlabel("底物浓度 S"); ax.set_ylabel("反应速率 V")
    ax.set_title('Michaelis–Menten：从"越多越快"到"吃饱封顶"')
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "michaelis_menten.png")


# 5) Scheffer 浅湖迟滞回线（折叠分岔）-----------------------------------
def fig_hysteresis():
    b, r, m, q = 0.6, 1.0, 1.0, 4
    x = np.linspace(0.001, 3.5, 4000)
    a = b*x - r*x**q/(x**q+m**q)          # 平衡关系 a=g(x)
    gp = np.gradient(a, x)                  # g'(x)：>0 稳定，<0 不稳定
    stable = gp > 0
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(a[stable], x[stable], ".", ms=2, color="tab:blue", label="稳定态")
    ax.plot(a[~stable], x[~stable], ".", ms=2, color="tab:red", label="不稳定（门槛）")
    a_up, a_down = a[gp > 0].max() if False else a.max(), a.min()
    a_up = a[np.where(np.diff(np.sign(np.diff(a))) != 0)[0]+1].max()
    a_dn = a[np.where(np.diff(np.sign(np.diff(a))) != 0)[0]+1].min()
    ax.axvline(a_up, color="gray", ls="--", lw=.8); ax.axvline(a_dn, color="gray", ls="--", lw=.8)
    ax.annotate("加营养到此\n突然变浊→", xy=(a_up, 0.5), fontsize=8, ha="right")
    ax.annotate("←降到此\n才变清", xy=(a_dn, 1.7), fontsize=8, ha="left")
    ax.set_xlabel("营养输入 a"); ax.set_ylabel("浊度 x")
    ax.set_title('迟滞：去程阈值 ≠ 回程阈值，"翻过去容易翻回来难"')
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "hysteresis_loop.png")


# 6) 鞍结分岔图 dx/dt=r−x² → x*=±√r ------------------------------------
def fig_saddle_node():
    r = np.linspace(0, 4, 300); xs = np.sqrt(r)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(r, xs, "b", lw=2, label="稳定平衡 (+√r)")
    ax.plot(r, -xs, "r--", lw=2, label="不稳定平衡 (-√r)")
    ax.plot(0, 0, "ko"); ax.annotate("临界点 r=0\n两平衡相撞湮灭", xy=(0, 0),
                                     xytext=(1.2, -1.2), fontsize=8,
                                     arrowprops=dict(arrowstyle="->"))
    ax.axvline(0, color="gray", lw=.6)
    ax.set_xlabel("控制参数 r"); ax.set_ylabel("平衡状态 x*")
    ax.set_title('鞍结分岔：参数越过 r=0，稳定态突然消失（系统"翻车"）')
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "saddle_node.png")


# 7) Farquhar A–Ci：Ac/Aj 双限制与交叉 ----------------------------------
def fig_a_ci():
    Ci = np.linspace(50, 1000, 400)
    Vcmax, J, Gstar, Kc, Ko, O, Rd = 60, 120, 42, 404.9, 278.4, 210, 1.0
    Ac = Vcmax*(Ci-Gstar)/(Ci+Kc*(1+O/Ko))
    Aj = J*(Ci-Gstar)/(4*Ci+8*Gstar)
    A = np.minimum(Ac, Aj)-Rd
    cross = Ci[np.argmin(np.abs(Ac-Aj))]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(Ci, Ac, "--", color="tab:orange", label="Ac：Rubisco 限制")
    ax.plot(Ci, Aj, "--", color="tab:green", label="Aj：光/电子限制")
    ax.plot(Ci, A, "b", lw=2, label="净光合 A=min(Ac,Aj)-Rd")
    ax.axvline(cross, color="gray", ls=":"); ax.text(cross+10, 2, f"交叉 Ci≈{cross:.0f}", fontsize=8)
    ax.set_xlabel("胞间 CO₂ 浓度 Ci (µmol/mol)"); ax.set_ylabel("光合速率 (µmol/m²/s)")
    ax.set_title("Farquhar A–Ci：低 Ci 受 Rubisco 限制，高 Ci 受光限制")
    ax.legend(); ax.grid(alpha=.3)
    save(fig, "a_ci_farquhar.png")


if __name__ == "__main__":
    fig_rothc_spinup()
    fig_exp_decay_euler()
    fig_dewit_gap()
    fig_michaelis_menten()
    fig_hysteresis()
    fig_saddle_node()
    fig_a_ci()
    print("全部图已生成到", HERE)
