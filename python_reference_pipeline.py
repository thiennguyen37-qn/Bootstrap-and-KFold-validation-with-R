import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict, train_test_split
from sklearn.preprocessing import StandardScaler


def metric_row(y_true, y_pred):
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "recall1": recall_score(y_true, y_pred, pos_label=1, zero_division=0),
        "precision1": precision_score(y_true, y_pred, pos_label=1, zero_division=0),
        "recall0": recall_score(y_true, y_pred, pos_label=0, zero_division=0),
        "precision0": precision_score(y_true, y_pred, pos_label=0, zero_division=0),
        "macro_f1": f1_score(y_true, y_pred, average="macro"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="Indian Liver Patient Dataset (ILPD).csv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--target-recall", type=float, default=0.70)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    columns = [
        "Age", "Gender", "Total_Bilirubin", "Direct_Bilirubin",
        "Alkaline_Phosphotase", "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
        "Total_Protiens", "Albumin", "Albumin_and_Globulin_Ratio", "Liver cirrhosis",
    ]
    df = pd.read_csv(args.data, header=None, names=columns)
    df["Liver cirrhosis"] = df["Liver cirrhosis"].map({2: 0, 1: 1})

    for col in [
        "Total_Bilirubin", "Direct_Bilirubin", "Alkaline_Phosphotase",
        "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
    ]:
        df[col] = np.log1p(df[col])

    df_model = df.drop(columns=["Direct_Bilirubin", "Gender", "Total_Protiens"])
    X = df_model.drop(columns=["Liver cirrhosis"]).to_numpy()
    y = df_model["Liver cirrhosis"].to_numpy()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )

    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("smote", SMOTE(random_state=42)),
        ("clf", LogisticRegression(max_iter=1000, solver="liblinear", random_state=42)),
    ])

    param_grid = {
        "clf__C": [0.001, 0.01, 0.1, 1, 10, 100],
        "clf__penalty": ["l1", "l2"],
        "smote__k_neighbors": [3, 5, 7],
    }
    cv_inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    gs = GridSearchCV(
        pipe, param_grid, cv=cv_inner, scoring="f1_macro", n_jobs=-1, refit=True
    )
    gs.fit(X_train, y_train)
    best_model = gs.best_estimator_

    cv_df = pd.DataFrame(gs.cv_results_)
    top = cv_df[["rank_test_score", "params", "mean_test_score", "std_test_score"]].copy()
    top = top.sort_values("rank_test_score").head(5)
    top["params"] = top["params"].map(json.dumps)
    top.to_csv(out_dir / "grid_top.csv", index=False)

    oof_proba = cross_val_predict(
        best_model, X_train, y_train, cv=cv_inner, method="predict_proba"
    )[:, 1]

    grid_thr = np.linspace(0.05, 0.95, 181)
    curve = []
    for thr in grid_thr:
        pred_thr = (oof_proba >= thr).astype(int)
        curve.append({
            "threshold": thr,
            "recall1": recall_score(y_train, pred_thr, pos_label=1, zero_division=0),
            "recall0": recall_score(y_train, pred_thr, pos_label=0, zero_division=0),
            "fn": int(((y_train == 1) & (pred_thr == 0)).sum()),
            "fp": int(((y_train == 0) & (pred_thr == 1)).sum()),
        })
    pd.DataFrame(curve).to_csv(out_dir / "threshold_curve.csv", index=False)

    prec, rec, thr = precision_recall_curve(y_train, oof_proba)
    lookup = []
    threshold = 0.5
    for target in [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]:
        ok = rec[:-1] >= target
        if not ok.any():
            continue
        candidates = np.where(ok)[0]
        i = candidates[np.argmax(thr[candidates])]
        pred_thr = (oof_proba >= thr[i]).astype(int)
        row = {
            "target_recall": target,
            "threshold": float(thr[i]),
            "recall1": recall_score(y_train, pred_thr, pos_label=1, zero_division=0),
            "recall0": recall_score(y_train, pred_thr, pos_label=0, zero_division=0),
        }
        lookup.append(row)
        if abs(target - args.target_recall) < 1e-12:
            threshold = float(thr[i])
    pd.DataFrame(lookup).to_csv(out_dir / "threshold_lookup.csv", index=False)

    proba_test = best_model.predict_proba(X_test)[:, 1]
    pred_final = (proba_test >= threshold).astype(int)
    pred_05 = (proba_test >= 0.5).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_test, pred_final).ravel()
    pd.DataFrame({
        "y_test": y_test,
        "proba_test": proba_test,
        "pred_05": pred_05,
        "pred_final": pred_final,
    }).to_csv(out_dir / "test_predictions.csv", index=False)

    final_metrics = pd.DataFrame([metric_row(y_test, pred_final)])
    final_metrics.to_csv(out_dir / "final_metrics.csv", index=False)
    pd.DataFrame([
        {"threshold_label": "0.500 (mac dinh)", **metric_row(y_test, pred_05)},
        {"threshold_label": f"{threshold:.3f} (da chinh)", **metric_row(y_test, pred_final)},
    ]).to_csv(out_dir / "tradeoff.csv", index=False)

    rng = np.random.default_rng(42)
    boot_rows = []
    n_test = len(y_test)
    for _ in range(args.bootstrap):
        idx = rng.integers(0, n_test, n_test)
        if y_test[idx].sum() == 0:
            continue
        boot_rows.append(metric_row(y_test[idx], pred_final[idx]))
    boot_df = pd.DataFrame(boot_rows)
    boot_df.to_csv(out_dir / "bootstrap_samples.csv", index=False)
    ci = pd.DataFrame({
        "metric": final_metrics.columns,
        "diem_uoc_luong": final_metrics.iloc[0].to_numpy(),
        "ci_duoi": boot_df.quantile(0.025).to_numpy(),
        "ci_tren": boot_df.quantile(0.975).to_numpy(),
        "sai_so_chuan": boot_df.std(ddof=1).to_numpy(),
    })
    ci.to_csv(out_dir / "bootstrap_ci.csv", index=False)

    with open(out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump({
            "train_n": int(len(y_train)),
            "test_n": int(len(y_test)),
            "train_pos": int(y_train.sum()),
            "test_pos": int(y_test.sum()),
            "best_params": gs.best_params_,
            "best_score": float(gs.best_score_),
            "n_param_sets": int(len(gs.cv_results_["params"])),
            "threshold": threshold,
            "target_recall": args.target_recall,
            "confusion": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
        }, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
