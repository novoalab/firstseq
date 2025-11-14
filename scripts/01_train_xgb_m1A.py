#!/usr/bin/env python3
import pandas as pd
import numpy as np
import joblib
import os
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
from xgboost import XGBClassifier

# === Setup ===
os.makedirs("m1A_prediction/models", exist_ok=True)
os.makedirs("m1A_prediction/plots", exist_ok=True)

# === Load training data ===
df = pd.read_csv("data/m1A_prediction/training/combined_training_m1A_unmod.tsv", sep="\t")

# === Preprocessing ===
df = df[df["coverage"] >= 20]
drop_cols = ["A", "T", "C", "G", "Ends", "rtstop", "norm_rt_end"]
df = df.drop(columns=[col for col in drop_cols if col in df.columns])

features = ["T_freq", "C_freq", "G_freq", "mis_freq"]

# === Loop over enzymes ===
for enzyme in df["Enzyme"].unique():
    print(f"\n Training model for {enzyme}...")

    df_enzyme = df[df["Enzyme"] == enzyme].copy()
    X = df_enzyme[features]
    y = (df_enzyme["Mod"] == "m1A").astype(int)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=0.2, random_state=42
    )

    pos_weight = (len(y_train) - sum(y_train)) / sum(y_train)

    clf = XGBClassifier(
        n_estimators=300,
        learning_rate=0.05,
        max_depth=3,
        random_state=42,
        use_label_encoder=False,
        eval_metric="logloss",
        scale_pos_weight=pos_weight
    )
    clf.fit(X_train, y_train)

    model_path = f"m1A_prediction/models/xgb_m1A_classifier_{enzyme}.joblib"
    joblib.dump(clf, model_path)
    print(f" Saved model: {model_path}")

    # Evaluation
    y_probs = clf.predict_proba(X_test)[:, 1]
    y_pred = clf.predict(X_test)

    fpr, tpr, _ = roc_curve(y_test, y_probs)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}", color="crimson")
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve — {enzyme}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"m1A_prediction/plots/roc_m1A_classifier_{enzyme}.pdf")
    plt.close()

    cm = confusion_matrix(y_test, y_pred)
    disp = ConfusionMatrixDisplay(cm, display_labels=["Unm", "m1A"])
    disp.plot()
    plt.title(f"Confusion Matrix — {enzyme}")
    plt.tight_layout()
    plt.savefig(f"m1A_prediction/plots/confusion_m1A_classifier_{enzyme}.pdf")
    plt.close()

    print(f" Train size: {len(y_train)} | Test size: {len(y_test)} | m1A in test: {sum(y_test)}")

print("\n All enzyme-specific XGBoost models trained and saved.")