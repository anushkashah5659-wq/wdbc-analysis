#!/usr/bin/env python3
"""Which tumour-nucleus measurements differ between malignant and benign
breast masses? Wisconsin Diagnostic Breast Cancer (WDBC) dataset, 569 samples.

For each of 30 features: Mann-Whitney U test (malignant vs benign), effect size,
then multiple-testing control (Benjamini-Hochberg FDR and Bonferroni).
A permuted-label control shows what the correction is protecting against.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.datasets import load_breast_cancer
from pathlib import Path

ALPHA = 0.05
N_PERMUTATIONS = 200
SEED = 1
OUT = Path("results")
OUT.mkdir(exist_ok=True)
Path("data").mkdir(exist_ok=True)


def benjamini_hochberg(p):
    """BH-adjusted p-values (q-values), same order as input."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1.0)
    return q


def test_all(X, is_malignant):
    """Mann-Whitney U for every column; returns p-values and effect sizes."""
    pvals, effects = [], []
    for j in range(X.shape[1]):
        a, b = X[is_malignant, j], X[~is_malignant, j]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        pvals.append(p)
        effects.append(2 * u / (len(a) * len(b)) - 1)  # rank-biserial: +1 = malignant always higher
    return np.array(pvals), np.array(effects)


def main():
    data = load_breast_cancer(as_frame=True)
    df = data.frame.copy()
    # sklearn codes target 0 = malignant, 1 = benign
    df["diagnosis"] = np.where(df["target"] == 0, "malignant", "benign")
    df.drop(columns="target").to_csv("data/wdbc.csv", index=False)

    features = list(data.feature_names)
    X = df[features].to_numpy()
    is_malignant = (df["diagnosis"] == "malignant").to_numpy()
    n_mal, n_ben = int(is_malignant.sum()), int((~is_malignant).sum())

    p, effect = test_all(X, is_malignant)
    q = benjamini_hochberg(p)
    bonf = np.minimum(p * len(p), 1.0)

    res = pd.DataFrame({
        "feature": features,
        "median_benign": np.median(X[~is_malignant], axis=0),
        "median_malignant": np.median(X[is_malignant], axis=0),
        "effect_size_rank_biserial": effect,
        "p_value": p,
        "q_value_BH": q,
        "p_bonferroni": bonf,
        "significant_uncorrected": p < ALPHA,
        "significant_BH_FDR5": q < ALPHA,
        "significant_bonferroni": bonf < ALPHA,
    }).sort_values("p_value")
    res.to_csv(OUT / "feature_tests.tsv", sep="\t", index=False, float_format="%.4g")

    # Permuted-label control: shuffle diagnosis, repeat everything.
    rng = np.random.default_rng(SEED)
    uncorr, bh_hits = [], []
    for _ in range(N_PERMUTATIONS):
        perm = rng.permutation(is_malignant)
        pp, _ = test_all(X, perm)
        uncorr.append(int((pp < ALPHA).sum()))
        bh_hits.append(int((benjamini_hochberg(pp) < ALPHA).sum()))
    perm_summary = {
        "n_permutations": N_PERMUTATIONS,
        "mean_uncorrected_false_positives_per_run": float(np.mean(uncorr)),
        "runs_with_at_least_one_uncorrected_false_positive": int(sum(u > 0 for u in uncorr)),
        "mean_BH_false_positives_per_run": float(np.mean(bh_hits)),
        "runs_with_at_least_one_BH_false_positive": int(sum(b > 0 for b in bh_hits)),
    }

    # Plot 1: effect size vs significance
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(res["effect_size_rank_biserial"], -np.log10(res["q_value_BH"].clip(lower=1e-300)),
               c=np.where(res["significant_BH_FDR5"], "tab:red", "grey"))
    ax.axhline(-np.log10(ALPHA), ls="--", c="k", lw=0.8)
    ax.set_xlabel("Effect size (rank-biserial; + = higher in malignant)")
    ax.set_ylabel("-log10 BH q-value")
    ax.set_title("Malignant vs benign: 30 nucleus features")
    fig.tight_layout(); fig.savefig(OUT / "effect_vs_significance.png", dpi=150); plt.close(fig)

    # Plot 2: feature correlation (shows the tests are not independent)
    corr = pd.DataFrame(X, columns=features).corr(method="spearman")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(features))); ax.set_xticklabels(features, rotation=90, fontsize=6)
    ax.set_yticks(range(len(features))); ax.set_yticklabels(features, fontsize=6)
    fig.colorbar(im, label="Spearman correlation"); ax.set_title("Feature correlations")
    fig.tight_layout(); fig.savefig(OUT / "feature_correlations.png", dpi=150); plt.close(fig)

    upper = corr.where(np.triu(np.ones(corr.shape, dtype=bool), k=1)).stack().dropna()
    summary = {
        "dataset": "Wisconsin Diagnostic Breast Cancer (UCI), via scikit-learn",
        "n_malignant": n_mal, "n_benign": n_ben, "n_features_tested": len(features),
        "alpha": ALPHA, "test": "two-sided Mann-Whitney U",
        "significant_uncorrected": int(res["significant_uncorrected"].sum()),
        "significant_BH_FDR5": int(res["significant_BH_FDR5"].sum()),
        "significant_bonferroni": int(res["significant_bonferroni"].sum()),
        "feature_pairs_with_abs_spearman_ge_0.9": int((upper.abs() >= 0.9).sum()),
        "feature_pairs_total": int(len(upper)),
        "permutation_control": perm_summary,
        "random_seed": SEED,
    }
    import json
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(res[["feature", "effect_size_rank_biserial", "q_value_BH"]].head(8).to_string(index=False))
    print("...")
    print(res[["feature", "effect_size_rank_biserial", "q_value_BH", "significant_BH_FDR5"]].tail(8).to_string(index=False))


if __name__ == "__main__":
    main()
