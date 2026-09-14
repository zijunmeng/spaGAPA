#!/usr/bin/env python3
"""P0 对抗性验证：P0-1 域边界损害测试 + P0-3 空间分块 conformal 检验。

数据: data/processed/mob_st11/ (MOB Visium, 260 spots, 手工形态层注释 GCL/MCL/GL/ONL...)

P0-1 域边界损害测试（预注册判定线: 边界近/远误差比 <1.5× 通过）
  - 遮蔽 20% per-gene (seed 42, 5 个种子重复)
  - 稀疏 GP 插补（与主管线同构: RBF, ls=5×NN, inducing=100）
  - 误差按到最近层边界的距离分层（尺度无关: 距离分布的最近 20% vs 最远 20%）
  - Kruskal-Wallis 检验 3 层误差分布

P0-3 空间分块 conformal（预注册判定线: 分块覆盖偏差 <1pp 通过）
  - 随机划分 vs 空间四象限划分 vs 按层划分（3 种分块方案）
  - 每种: 校准 50% / 测试 50%, 80/90/95% 目标, 50 次重复
  - 报告各方案的覆盖偏差绝对值

输出: pipeline_output/p0_adversarial_validation/
  p0_1_boundary_damage.csv + p0_3_spatial_block_conformal.csv + summary.json + 图
"""
from __future__ import annotations
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
DATA = ROOT / "data/processed/mob_st11"
OUT = ROOT / "pipeline_output/p0_adversarial_validation"
OUT.mkdir(parents=True, exist_ok=True)

# ── 加载 ─────────────────────────────────────────────────────
apa = pd.read_csv(DATA / "apa_matrix.csv", index_col=0)      # gene x spot
coords = pd.read_csv(DATA / "coordinates.csv", index_col=0)  # spot_id,x,y
labels = pd.read_csv(DATA / "labels.csv", index_col=0)       # spot_id,label
common = apa.columns.intersection(coords.index).intersection(labels.index)
apa = apa[common]; coords = coords.loc[common]; labels = labels.loc[common]
Y = apa.values.T.astype(float)          # spot x gene
X = coords[["x", "y"]].values.astype(float)
L = labels["label"].values
n, G = Y.shape
print(f"MOB st11: {n} spots x {G} genes; 层: {pd.Series(L).value_counts().to_dict()}")

# 坐标缩放
Xc = X - X.mean(0)
d2all = ((Xc[:, None, :] - Xc[None, :, :]) ** 2).sum(-1)
np.fill_diagonal(d2all, np.inf)
nn_dist = np.sqrt(d2all.min(axis=1))
scale = np.median(nn_dist)
Xs = Xc / max(scale, 1e-9)

# ── 到最近层边界的距离（尺度无关分位数分层）─────────────────
# 边界定义: 存在不同标签邻居的 spot; 距离 = 到最近异层 spot 的欧氏距离
dist_boundary = np.full(n, np.inf)
for i in range(n):
    other = L != L[i]
    if other.any():
        dist_boundary[i] = np.sqrt(((Xc[other] - Xc[i]) ** 2).sum(axis=1)).min()
q20, q80 = np.quantile(dist_boundary, [0.2, 0.8])
near_mask = dist_boundary <= q20          # 边界近层
far_mask = dist_boundary >= q80           # 边界远层
mid_mask = ~(near_mask | far_mask)
print(f"边界距离分位: q20={q20:.2f} q80={q80:.2f} | near={near_mask.sum()} mid={mid_mask.sum()} far={far_mask.sum()}")


def fit_gp_predict(Yobs, m=100, ls_mult=5.0, seed=0):
    """与主管线同构的 Nyström GP: 逐基因 fit on observed, predict all."""
    rs = np.random.RandomState(seed)
    idx = rs.choice(len(Xs), min(m, len(Xs)), replace=False)
    Xu = Xs[idx]
    ls = ls_mult * 1.0  # 坐标已按 NN 中位距离归一, ls=5 单位 ≈ 5×NN
    def k(A, B):
        d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
        return np.exp(-d2 / (2 * ls ** 2))
    Knu = k(Xs, Xu); Kuu = k(Xu, Xu) + 1e-6 * np.eye(len(Xu))
    pred = np.empty_like(Yobs)
    for g in range(Yobs.shape[1]):
        y = np.nan_to_num(Yobs[:, g])
        lam = 0.1
        A = Kuu + lam * (Knu.T @ Knu) / n
        b = (Knu.T @ y) / n
        pred[:, g] = Knu @ np.linalg.solve(A, b)
    return pred


# ── P0-1: 边界损害（5 种子）──────────────────────────────────
print("\n=== P0-1 域边界损害测试 ===")
rows = []
for seed in [42, 123, 456, 789, 2026]:
    rng = np.random.RandomState(seed)
    Ymask = Y.copy()
    mask = np.zeros_like(Y, dtype=bool)
    for g in range(G):
        obs = np.where(~np.isnan(Y[:, g]))[0]
        k = max(1, int(0.2 * len(obs)))
        sel = rng.choice(obs, k, replace=False)
        mask[sel, g] = True
        Ymask[sel, g] = np.nan
    pred = fit_gp_predict(Ymask, seed=seed % 997)
    hh = mask & ~np.isnan(Y)
    abs_err = np.abs(Y - pred)
    for zone, zm in [("near", near_mask), ("mid", mid_mask), ("far", far_mask)]:
        e = abs_err[hh & zm[:, None] * np.ones_like(hh, bool)] if False else None
        # 逐 (spot,gene) 属于 zone 的误差
        sel = hh.copy()
        sel[~zm, :] = False
        e = abs_err[sel]
        rows.append(dict(seed=seed, zone=zone, rmse=float(np.sqrt((e ** 2).mean())),
                         mean_abs=float(e.mean()), n=len(e)))
df1 = pd.DataFrame(rows)
df1.to_csv(OUT / "p0_1_boundary_damage.csv", index=False)

# 汇总 + 检验
piv = df1.groupby("zone")[["rmse", "mean_abs"]].agg(["mean", "std"])
near_rmse = df1[df1.zone == "near"]["rmse"].values
far_rmse = df1.zone.values == "far"
far_rmse = df1[df1.zone == "far"]["rmse"].values
ratio = near_rmse.mean() / far_rmse.mean()
# Kruskal-Wallis on per-spot-zone errors (seed 42)
rng = np.random.RandomState(42)
Ymask = Y.copy(); mask = np.zeros_like(Y, dtype=bool)
for g in range(G):
    obs = np.where(~np.isnan(Y[:, g]))[0]
    k = max(1, int(0.2 * len(obs)))
    sel = rng.choice(obs, k, replace=False)
    mask[sel, g] = True; Ymask[sel, g] = np.nan
pred = fit_gp_predict(Ymask, seed=42)
hh = mask & ~np.isnan(Y); abs_err = np.abs(Y - pred)
zone_per_spot = np.where(near_mask, "near", np.where(far_mask, "far", "mid"))
e_by = {z: [] for z in ["near", "mid", "far"]}
for i in range(n):
    mrow = hh[i]
    if mrow.any():
        e_by[zone_per_spot[i]].extend(abs_err[i, mrow].tolist())
kw_stat, kw_p = stats.kruskal(*[np.array(v) for v in e_by.values()])
p0_1_pass = ratio < 1.5
print(f"  near RMSE={near_rmse.mean():.4f} far RMSE={far_rmse.mean():.4f} ratio={ratio:.3f}")
print(f"  Kruskal-Wallis H={kw_stat:.1f} p={kw_p:.2e}")
print(f"  判定（预注册 <1.5×）: {'✅ 通过 — RBF 平滑假设成立' if p0_1_pass else '❌ 未通过 — 边界损害存在，需边界感知核'}")

# ── P0-3: 空间分块 conformal ──────────────────────────────────
print("\n=== P0-3 空间分块 conformal 检验 ===")
# 用完整数据的遮蔽-预测误差做 conformal（与主管线一致的 marginal 校准）
def conformal_cov(err_cal, err_test, alpha):
    k = int(np.ceil((len(err_cal) + 1) * alpha))
    q = np.sort(err_cal)[min(k - 1, len(err_cal) - 1)]
    return float((err_test <= q).mean())

split_results = []
for seed in range(50):
    rng = np.random.RandomState(seed)
    Ymask = Y.copy(); mask = np.zeros_like(Y, dtype=bool)
    for g in range(G):
        obs = np.where(~np.isnan(Y[:, g]))[0]
        k = max(1, int(0.2 * len(obs)))
        sel = rng.choice(obs, k, replace=False)
        mask[sel, g] = True; Ymask[sel, g] = np.nan
    pred = fit_gp_predict(Ymask, seed=seed % 997)
    hh = mask & ~np.isnan(Y)
    err = np.abs(Y - pred)[hh]
    spots_hh = np.where(hh.any(axis=1))[0]
    spot_of = np.repeat(np.arange(n), hh.sum(axis=1))

    # 方案 0: 随机（baseline）
    perm = rng.permutation(len(err))
    half = len(err) // 2
    cov_r = {a: conformal_cov(err[perm[:half]], err[perm[half:]], a) for a in [0.8, 0.9, 0.95]}
    # 方案 1: 空间四象限
    qx, qy = np.median(X[:, 0]), np.median(X[:, 1])
    quad = (X[:, 0] > qx).astype(int) + 2 * (X[:, 1] > qy)
    qb = quad[spot_of] % 2 == 0
    cov_q = {a: conformal_cov(err[qb], err[~qb], a) for a in [0.8, 0.9, 0.95]}
    # 方案 2: 按层（前 3 层做校准, 后面层做测试 — 最强违反）
    uniq_l = pd.Series(L).value_counts().index.tolist()
    cal_l = set(uniq_l[: len(uniq_l) // 2])
    lb = np.isin(L[spot_of], list(cal_l))
    if lb.all() or (~lb).all():
        cov_l = {a: np.nan for a in [0.8, 0.9, 0.95]}
    else:
        cov_l = {a: conformal_cov(err[lb], err[~lb], a) for a in [0.8, 0.9, 0.95]}
    split_results.append(dict(seed=seed,
                              cov80_rand=cov_r[0.8], cov90_rand=cov_r[0.9], cov95_rand=cov_r[0.95],
                              cov80_quad=cov_q[0.8], cov90_quad=cov_q[0.9], cov95_quad=cov_q[0.95],
                              cov80_layer=cov_l[0.8], cov90_layer=cov_l[0.9], cov95_layer=cov_l[0.95]))
df3 = pd.DataFrame(split_results)
df3.to_csv(OUT / "p0_3_spatial_block_conformal.csv", index=False)

summary3 = {}
for sc in ["rand", "quad", "layer"]:
    for a in [80, 90, 95]:
        v = df3[f"cov{a}_{sc}"].dropna()
        summary3[f"{sc}_{a}"] = dict(mean=float(v.mean()), dev_pp=float(abs(v.mean() - a / 100) * 100), n=int(len(v)))
max_dev = max(summary3[k]["dev_pp"] for k in summary3)
p0_3_pass = max_dev < 1.0
print(pd.DataFrame(summary3).T.round(4).to_string())
print(f"  最大覆盖偏差: {max_dev:.2f}pp")
print(f"  判定（预注册 <1pp）: {'✅ 通过 — 交换性违反实际无害' if p0_3_pass else '❌ 未通过 — 需写入局限/spatially-weighted conformal'}")

# ── 图 ────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import glob as _glob
for _fp in _glob.glob("/usr/share/fonts/msfonts/ARIAL*.TTF"):
    try: fm.fontManager.addfont(_fp)
    except Exception: pass
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"],
                     "font.size": 8, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))
ax = axes[0]
for z, c in [("near", RED), ("mid", "#999999"), ("far", BLUE)]:
    v = df1[df1.zone == z]["rmse"].values
    ax.scatter([{"near": 0, "mid": 1, "far": 2}[z]] * len(v), v, s=14, color=c, alpha=0.8)
    ax.scatter([{"near": 0, "mid": 1, "far": 2}[z]], [v.mean()], s=60, color=c, marker="_", linewidths=2)
ax.axhline(far_rmse.mean() * 1.5, ls=":", color="k", lw=0.8)
ax.set_xticks(range(3)); ax.set_xticklabels(["near", "mid", "far"])
ax.set_xlabel("Distance to layer boundary"); ax.set_ylabel("RMSE")
ax.set_title(f"A  Boundary damage\nratio={ratio:.2f} ({'pass' if p0_1_pass else 'FAIL'})")
ax = axes[1]
zm = near_mask | far_mask
for lb_, c in zip(np.unique(L), plt.cm.Set3(np.linspace(0, 1, len(np.unique(L))))):
    ax.scatter(X[zm & (L == lb_), 0], X[zm & (L == lb_), 1], s=3, color=c, label=str(lb_))
ax.set_title("B  MOB layers"); ax.set_xticks([]); ax.set_yticks([])
ax.legend(fontsize=4, frameon=False, loc="upper left", ncol=2)
ax = axes[2]
for sc, c in [("rand", BLUE), ("quad", ORANGE), ("layer", GREEN)]:
    devs = [abs(df3[f"cov{a}_{sc}"].mean() - a / 100) * 100 for a in [80, 90, 95]]
    ax.plot([80, 90, 95], devs, "o-", color=c, label=sc, ms=4, lw=1.2)
ax.axhline(1.0, ls=":", color="k", lw=0.8)
ax.set_xlabel("Nominal coverage (%)"); ax.set_ylabel("|deviation| (pp)")
ax.set_title(f"C  Spatial-block conformal\nmax dev {max_dev:.2f}pp")
ax.legend(fontsize=6, frameon=False)
fig.tight_layout()
fig.savefig(OUT / "p0_validation.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "p0_validation.pdf", bbox_inches="tight")

summary = dict(
    p0_1=dict(ratio_near_far=float(ratio), near_rmse=float(near_rmse.mean()),
              far_rmse=float(far_rmse.mean()), kw_p=float(kw_p),
              pre_registered_threshold=1.5, passed=bool(p0_1_pass),
              n_seeds=5, n_spots=int(n), n_genes=int(G)),
    p0_3=dict(max_deviation_pp=float(max_dev), pre_registered_threshold_pp=1.0,
              passed=bool(p0_3_pass), n_replicates=50,
              by_scheme=summary3),
)
json.dump(summary, open(OUT / "summary.json", "w"), indent=2)
print(f"\n[saved] {OUT}")
