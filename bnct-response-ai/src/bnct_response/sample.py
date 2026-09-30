"""Sample a direct-distance response from the chained baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from .baseline import inverse_transport


def _tree_samples(pipeline, frame: pd.DataFrame, count: int, rng: np.random.Generator) -> np.ndarray:
    transformed = pipeline.named_steps["features"].transform(frame)
    forest = pipeline.named_steps["model"]
    indices = rng.integers(0, len(forest.estimators_), size=count)
    return np.asarray([forest.estimators_[index].predict(transformed)[0] for index in indices])


def sample_response(
    artifact: dict,
    particle: str,
    initial_energy_mev: float,
    distance_um: float,
    direction: tuple[float, float, float],
    count: int,
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    norm = float(np.linalg.norm(direction))
    if norm == 0:
        raise ValueError("Direction must be non-zero")
    direction = tuple(value / norm for value in direction)
    base = {
        "Particle": particle,
        "InitialEnergy_MeV": initial_energy_mev,
        "Distance_um": distance_um,
        "InitialDirX": direction[0], "InitialDirY": direction[1], "InitialDirZ": direction[2],
    }
    condition = pd.DataFrame([base])
    range_input = condition[["Particle", "InitialEnergy_MeV", "InitialDirX", "InitialDirY", "InitialDirZ"]]
    support = {
        name: float(model.predict(range_input)[0])
        for name, model in artifact["range_models"].items()
    }
    support["requested_distance_um"] = distance_um
    support["inside_median_range"] = distance_um <= support["q50"]
    support["inside_q95_range"] = distance_um <= support["q95"]

    transport_log = _tree_samples(artifact["transport"], condition, count, rng)
    transport_values = inverse_transport(transport_log)
    samples = pd.DataFrame(transport_values, columns=artifact["transport_targets"])
    expanded = pd.concat([condition] * count, ignore_index=True)
    for target in artifact["transport_targets"]:
        expanded[target] = samples[target]

    direct_condition = pd.concat([condition] * count, ignore_index=True)
    any_probability = artifact["direct_damage_any"].predict_proba(direct_condition)[:, 1]
    dsb_probability = artifact["direct_damage_dsb"].predict_proba(direct_condition)[:, 1]
    count_log = _tree_samples(artifact["direct_damage_counts"], direct_condition, count, rng)
    count_mean = np.expm1(count_log).clip(min=0)
    counts = rng.poisson(count_mean)
    any_damage = rng.random(count) < any_probability
    any_dsb = (rng.random(count) < dsb_probability) & any_damage
    counts[~any_damage, :] = 0
    dsb_index = artifact["damage_count_targets"].index("TotalDSB")
    cdsb_index = artifact["damage_count_targets"].index("TotalCDSB")
    counts[~any_dsb, dsb_index] = 0
    counts[~any_dsb, cdsb_index] = 0
    samples["P_AnyDamage"] = any_probability
    samples["P_AnyDSB"] = dsb_probability
    samples["AnyDamage"] = any_damage.astype(int)
    samples["AnyDSB"] = any_dsb.astype(int)
    for index, target in enumerate(artifact["damage_count_targets"]):
        samples[target] = counts[:, index]
    return samples, support


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--particle", choices=["alpha", "lithium"], required=True)
    parser.add_argument("--initial-energy-mev", type=float, required=True)
    parser.add_argument("--distance-um", type=float, required=True)
    parser.add_argument("--direction", type=float, nargs=3, default=(0.0, 0.0, 1.0))
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = joblib.load(args.model_dir / "models.joblib")
    samples, support = sample_response(
        artifact, args.particle, args.initial_energy_mev, args.distance_um,
        tuple(args.direction), args.samples, args.seed,
    )
    print(json.dumps({"support": support, "sample_mean": samples.mean().to_dict()}, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        samples.to_csv(args.output, index=False)
        print(f"Wrote {len(samples)} samples to {args.output}")


if __name__ == "__main__":
    main()
