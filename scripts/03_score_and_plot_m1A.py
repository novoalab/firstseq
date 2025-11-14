#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_score_and_plot_m1A.py
Author: Oguzhan Begik

Score m1A across independent human mtRNA datasets using:
• enzyme-specific XGBoost models (TGIRT, SS3)
• boosted score = XGB score + 10 × RT drop
• smart-diff metric across TGIRT/SS3
• SNP removal based on mis_freq>0.9 in both enzymes
• Peak calling per transcript
• Export full score tables (long + wide)
• Generate per-enzyme plots with 4 tracks:
    - XGB score
    - RT drop (shifted)
    - Boosted score (smoothed)
    - Smart-diff
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import joblib, os

# ============================================================
#                 0. SETUP
# ============================================================
os.makedirs("m1A_prediction/results", exist_ok=True)
os.makedirs("m1A_prediction/plots", exist_ok=True)

print(" Running m1A scoring + plotting pipeline...")


# ============================================================
#                 1. LOAD INPUT DATA
# ============================================================
df = pd.read_csv("m1A_prediction/results/predictions_all_independent.tsv", sep="\t")

required_cols = {
    "chr","pos","Enzyme","ref_nuc","xgb_score",
    "norm_rt_shifted","xgb_boosted_score"
}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f" Missing columns in predictions_all_independent.tsv: {missing}")


# ============================================================
#                 2. LIST TRANSCRIPTS
# ============================================================
transcripts = sorted(df["chr"].unique())
print(f" Found {len(transcripts)} transcripts: {transcripts}")


# ============================================================
#         3. COMPUTE BOOSTED SCORE + SMOOTHING
# ============================================================
df = df.sort_values(["chr","Enzyme","pos"])
df["xgb_boosted_score_smooth"] = (
    df.groupby(["chr","Enzyme"])["xgb_boosted_score"]
    .transform(lambda x: x.rolling(window=5, center=True, min_periods=1).mean())
)


# ============================================================
#        4. PIVOT FOR SMART-DIFF (TGIRT vs SS3)
# ============================================================
pivoted = df.pivot_table(
    index=["chr","pos"],
    columns="Enzyme",
    values="xgb_boosted_score_smooth"
).reset_index()

pivoted["smart_diff"] = np.nan

for tx in transcripts:
    sub = pivoted[pivoted["chr"] == tx]

    # if one enzyme missing, skip smart-diff
    if {"TGIRT","SS3"} - set(sub.columns):
        continue

    tg95 = sub["TGIRT"].quantile(0.99)
    ss95 = sub["SS3"].quantile(0.99)

    smart = np.where(
        (sub["TGIRT"] > tg95) & (sub["SS3"] > ss95),
        sub[["TGIRT","SS3"]].max(axis=1),
        (sub["TGIRT"] - sub["SS3"]).abs()
    )

    pivoted.loc[pivoted["chr"]==tx,"smart_diff"] = smart


# ============================================================
#                 5. JOIN OTHER FEATURES
# ============================================================
pivoted = pivoted.merge(
    df[["chr","pos","ref_nuc","xgb_score","norm_rt_shifted"]].drop_duplicates(),
    on=["chr","pos"],
    how="left"
)


# ============================================================
#                 6. CALL CANDIDATE SITES
# ============================================================
labeled_sites = []

for tx in transcripts:
    sub = pivoted[pivoted["chr"]==tx].copy()

    # Only adenosines
    sub = sub[sub["ref_nuc"]=="A"]

    # Score ≥ 0.90
    sub = sub[sub["xgb_score"]>=0.90]

    if len(sub)==0:
        continue

    # Smart threshold
    smart_thresh = sub["smart_diff"].quantile(0.995)
    sub = sub[sub["smart_diff"]>=smart_thresh]

    # windowed peak selection (10-nt)
    sub = sub.sort_values("pos").reset_index(drop=True)
    group = []

    for _, row in sub.iterrows():
        if not group:
            group.append(row)
        elif row["pos"] - group[-1]["pos"] <= 10:
            group.append(row)
        else:
            labeled_sites.append(max(group, key=lambda r: r["smart_diff"]))
            group = [row]

    if group:
        labeled_sites.append(max(group, key=lambda r: r["smart_diff"]))


labeled_df = pd.DataFrame(labeled_sites)
labeled_df.to_csv("m1A_prediction/results/labeled_peak_positions_gb.tsv", sep="\t", index=False)
print(f" Saved peak positions → labeled_peak_positions_gb.tsv")


# ============================================================
#                 7. EXPORT LONG + WIDE TABLES
# ============================================================

# --- wide table ---
wide = pivoted.merge(
    labeled_df[["chr","pos"]].assign(is_candidate=True),
    on=["chr","pos"],
    how="left"
)
wide["is_candidate"] = wide["is_candidate"].fillna(False)
wide.to_csv("m1A_prediction/results/all_positions_scores_wide.tsv", sep="\t", index=False)

# --- long table ---
long = (
    df[[
        "chr","pos","Enzyme","ref_nuc",
        "xgb_score","norm_rt_shifted",
        "xgb_boosted_score","xgb_boosted_score_smooth"
    ]]
    .merge(pivoted[["chr","pos","smart_diff"]], on=["chr","pos"], how="left")
    .merge(labeled_df[["chr","pos"]].assign(is_candidate=True),
           on=["chr","pos"], how="left")
)

long["is_candidate"] = long["is_candidate"].fillna(False)
long.to_csv("m1A_prediction/results/all_positions_scores_long.tsv", sep="\t", index=False)

print(" Exported long + wide score tables.")


# ============================================================
#                 8. PLOTTING
# ============================================================
for tx in transcripts:
    for enzyme in ["TGIRT","SS3"]:
        df_sub = df[(df["chr"]==tx) & (df["Enzyme"]==enzyme)]
        if df_sub.empty:
            continue

        df_tx = df[df["chr"]==tx]
        piv_tx = pivoted[pivoted["chr"]==tx]

        fig, axs = plt.subplots(4,1, figsize=(14,10), sharex=True)
        fig.suptitle(f"{enzyme} m1A Signature — {tx}", fontsize=16)

        axs[0].bar(df_sub["pos"], df_sub["xgb_score"], width=1, color="purple")
        axs[0].set_ylabel("XGB")

        axs[1].bar(df_sub["pos"], df_sub["norm_rt_shifted"], width=1, color="gray")
        axs[1].set_ylabel("RT-drop")

        axs[2].bar(df_sub["pos"], df_sub["xgb_boosted_score_smooth"], width=1, color="steelblue")
        axs[2].set_ylabel("Boosted")

        axs[3].bar(piv_tx["pos"], piv_tx["smart_diff"], width=1, color="green")
        axs[3].set_ylabel("Smart-diff")
        axs[3].set_xlabel("Position")

        # Annotate candidates
        cand = labeled_df[labeled_df["chr"]==tx]
        for _, row in cand.iterrows():
            axs[3].text(
                row["pos"], row["smart_diff"]+0.05,
                str(int(row["pos"])),
                fontsize=7, rotation=90,
                ha="center", va="bottom"
            )

        for ax in axs:
            ax.spines[['top','right']].set_visible(False)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(10))

        fig.tight_layout()
        out = f"m1A_prediction/plots/m1A_signature_{tx}_{enzyme}.pdf"
        plt.savefig(out, format="pdf")
        plt.close()
        print(f" Saved {out}")

print("\n DONE — scoring + exporting + plotting complete.")