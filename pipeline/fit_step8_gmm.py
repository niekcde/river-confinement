"""Fit one confinement GMM using the settings of the previous 7-cluster analysis."""

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import xarray as xr
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler


FEATURES = (
    "slope_out_normalized",
    "slope_right_normalized",
    "slope_left_normalized",
    "slope_inn_normalized",
    "slope_out_normalized_smooth",
    "slope_right_normalized_smooth",
    "slope_left_normalized_smooth",
    "slope_inn_normalized_smooth",
)
IDENTIFIERS = ("file", "bendID", "reach_id", "combined_reach_id")


def fit_gmm(input_file: Path, output_dir: Path) -> dict:
    with xr.open_dataset(input_file) as ds:
        missing = set(FEATURES + IDENTIFIERS) - set(ds.data_vars)
        if missing:
            raise ValueError(f"Missing input variables: {sorted(missing)}")
        frame = ds[list(FEATURES + IDENTIFIERS)].to_dataframe().reset_index()
        total_rows = ds.sizes["index"]

    frame = frame.dropna(subset=FEATURES).copy()
    values = frame[list(FEATURES)].to_numpy()
    if len(frame) < 7 or not np.isfinite(values).all():
        raise ValueError("At least seven rows with finite slope features are required")

    scaler = StandardScaler()
    scaled = scaler.fit_transform(values)
    pca = PCA(n_components=len(FEATURES))
    components = pca.fit_transform(scaled)[:, :3]
    model = GaussianMixture(
        n_components=7,
        init_params="kmeans",
        covariance_type="tied",
        max_iter=100,
        n_init=5,
        random_state=42,
    )
    model.fit(components)
    labels = model.predict(components)
    probabilities = model.predict_proba(components)

    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = frame[["index", *IDENTIFIERS]].copy()
    assignments["GMM_7_hard"] = labels
    for cluster in range(7):
        assignments[f"GMM_7_{cluster}_soft"] = probabilities[:, cluster]
    assignments_file = output_dir / "gmm7_hf2_assignments.parquet"
    assignments.to_parquet(assignments_file, index=False)

    model_file = output_dir / "gmm7_hf2_model.joblib"
    joblib.dump({"scaler": scaler, "pca": pca, "gmm": model, "features": FEATURES}, model_file)
    summary = {
        "input_file": str(input_file.resolve()),
        "input_rows": total_rows,
        "eligible_rows": len(values),
        "features": FEATURES,
        "pca_components_used": 3,
        "gmm": {
            "n_components": 7,
            "init_params": "kmeans",
            "covariance_type": "tied",
            "max_iter": 100,
            "n_init": 5,
            "random_state": 42,
            "converged": bool(model.converged_),
            "n_iter": int(model.n_iter_),
            "lower_bound": float(model.lower_bound_),
        },
        "cluster_counts": {str(i): int((labels == i).sum()) for i in range(7)},
        "assignments_file": str(assignments_file),
        "model_file": str(model_file),
    }
    summary_file = output_dir / "gmm7_hf2_summary.json"
    summary_file.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(fit_gmm(args.input_file, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
