#!/usr/bin/env python3
import pandas as pd
import numpy as np
import joblib
import os

os.makedirs("m1A_prediction/results", exist_ok=True)

tgirt = pd.read_csv("data/m1A_prediction/independent/tgirt_independent_processed.tsv", sep="\t")
ss3  = pd.read_csv("data/m1A_prediction/independent/ss3_independent_processed.tsv",  sep="\t")

df = pd.concat([tgirt, ss3], ignore_index=True)
df = df.dropna(subset=["uniq_coord", "norm_rt_end", "A", "T", "C", "G", "coverage"])
df = df[df["coverage"] >= 20]

df = df.sort_values(["Enzyme", "Buffer", "chr", "pos"])

df["norm_rt_roll"] = (
    df.groupby(["Enzyme", "Buffer", "chr"])["norm_rt_end"]
      .transform(lambda x: x.rolling(window=5, center=True, min_periods=1).median())
)
df["norm_rt_shifted"] = (
    df.groupby(["Enzyme", "Buffer", "chr"])["norm_rt_roll"]
      .shift(-13)
)
df = df.drop(columns=["norm_rt_roll"])

df["A_freq"] = df["A"] / df["coverage"]
df["T_freq"] = df["T"] / df["coverage"]
df["C_freq"] = df["C"] / df["coverage"]
df["G_freq"] = df["G"] / df["coverage"]

df = df.dropna(subset=["T_freq", "C_freq", "G_freq", "mis_freq", "norm_rt_shifted"])

# SNP-like sites (very high mis_freq in both enzymes)
snp_candidates = (
    df[df["mis_freq"] > 0.9]
    .groupby(["chr", "pos", "ref_nuc"])
    .filter(lambda g: g["Enzyme"].nunique() >= 2)
)[["chr", "pos", "ref_nuc"]].drop_duplicates()

print(f"🧬 Found {len(snp_candidates)} potential SNP positions to exclude.")

df = df.merge(snp_candidates, on=["chr", "pos", "ref_nuc"], how="left", indicator=True)
df = df[df["_merge"] == "left_only"].drop(columns="_merge")

features = ["T_freq", "C_freq", "G_freq", "mis_freq"]
all_preds = []
all_hits  = []

for enzyme in df["Enzyme"].unique():
    model_path = f"m1A_prediction/models/xgb_m1A_classifier_{enzyme}.joblib"
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model not found for {enzyme}: {model_path}")

    model = joblib.load(model_path)
    df_sub = df[df["Enzyme"] == enzyme].copy()
    X = df_sub[features]

    df_sub["xgb_score"] = model.predict_proba(X)[:, 1]
    df_sub["xgb_pred"]  = model.predict(X)

    df_sub["xgb_boosted_score"] = df_sub["xgb_score"] + (df_sub["norm_rt_shifted"] * 15)

    # dynamic threshold based on background
    background = df_sub[
        (df_sub["ref_nuc"] != "A") | (df_sub["xgb_score"] < 0.5)
    ]["xgb_boosted_score"]
    X_factor = 75.0
    threshold = background.median() * X_factor

    print(f" {enzyme}: boosted threshold = {threshold:.3f}")

    hits = df_sub[
        (df_sub["ref_nuc"] == "A") &
        (df_sub["xgb_pred"] == 1) &
        (df_sub["xgb_boosted_score"] > threshold)
    ].copy()

    all_preds.append(df_sub)
    all_hits.append(hits)

df_all   = pd.concat(all_preds)
hits_all = pd.concat(all_hits)

df_all.to_csv("m1A_prediction/results/predictions_all_independent.tsv", sep="\t", index=False)
hits_all.to_csv("m1A_prediction/results/predicted_m1A_xgb.tsv", sep="\t", index=False)

print(" Prediction complete.")
print(f" Total positions analyzed: {len(df_all)}")
print(f" Total m1A-like hits: {len(hits_all)}")

drop_cols = [
    "A", "T", "C", "G", "Ends", "rtstop", "norm_rt_end", "ins", "del",
    "mean_qual", "median_qual", "Biotype", "Buffer", "Ref", "Pos", "mismatch",
    "norm_cov", "uniq_coord", "xgb_pred"
]
hits_clean = hits_all.drop(columns=[c for c in drop_cols if c in hits_all.columns])
hits_clean.to_csv("m1A_prediction/results/predicted_m1A_xgb_cleaned.tsv", sep="\t", index=False)

print(" Cleaned filtered hits saved.")