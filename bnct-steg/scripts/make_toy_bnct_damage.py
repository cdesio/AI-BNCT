import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from config import BNCT_DAMAGE_FEATURES


def parse_args():
    parser = argparse.ArgumentParser(description="Create a toy BNCT damage ROOT table for StEG smoke tests.")
    parser.add_argument("--output", default="data/molecular_bnct_damage.root")
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--rows", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    n = args.rows

    energy = rng.choice([1.47e6, 1.78e6, 0.84e6, 1.01e6], size=n)
    radius = rng.uniform(3.0, 13.0, size=n)
    theta = rng.uniform(0.0, np.pi, size=n)
    phi = rng.uniform(-np.pi, np.pi, size=n)

    pos_x = radius * np.sin(theta) * np.cos(phi)
    pos_y = radius * np.sin(theta) * np.sin(phi)
    pos_z = radius * np.cos(theta)

    deposited = rng.lognormal(mean=2.4, sigma=0.8, size=n)
    direct_breaks = rng.poisson(np.clip(deposited / 40.0, 0.0, 3.0))
    indirect_breaks = rng.poisson(np.clip(deposited / 70.0, 0.0, 2.0))
    strand_damage = np.clip(direct_breaks + indirect_breaks, 0, 4)
    base_damage = rng.poisson(np.clip(deposited / 55.0, 0.0, 3.0))
    induced_breaks = rng.binomial(1, np.clip(base_damage / 5.0, 0.0, 0.8))

    type_class = np.select(
        [
            strand_damage == 0,
            strand_damage == 1,
            strand_damage == 2,
            strand_damage > 2,
        ],
        [0, 1, 3, 5],
        default=0,
    )
    source_class = np.select(
        [
            direct_breaks > indirect_breaks,
            indirect_breaks > direct_breaks,
            (direct_breaks > 0) & (indirect_breaks > 0),
        ],
        [5, 6, 7],
        default=0,
    )

    arrays = {
        "Energy": energy,
        "Position_x_um": pos_x,
        "Position_y_um": pos_y,
        "Position_z_um": pos_z,
        "Size_nm": rng.gamma(shape=2.0, scale=3.0, size=n),
        "FragmentLength": rng.integers(1, 80, size=n),
        "BaseDamage": base_damage,
        "StrandDamage": strand_damage,
        "DirectBreaks": direct_breaks,
        "IndirectBreaks": indirect_breaks,
        "EaqBaseHits": rng.poisson(np.clip(deposited / 120.0, 0.0, 2.0)),
        "EaqStrandHits": rng.poisson(np.clip(deposited / 160.0, 0.0, 2.0)),
        "OHBaseHits": rng.poisson(np.clip(deposited / 35.0, 0.0, 5.0)),
        "OHStrandHits": rng.poisson(np.clip(deposited / 60.0, 0.0, 3.0)),
        "HBaseHits": rng.poisson(np.clip(deposited / 140.0, 0.0, 2.0)),
        "HStrandHits": rng.poisson(np.clip(deposited / 180.0, 0.0, 2.0)),
        "EnergyDeposited_eV": deposited,
        "InducedBreaks": induced_breaks,
        "TypeClassificationInt": type_class,
        "SourceClassificationInt": source_class,
    }

    arrays = {feature: np.asarray(arrays[feature], dtype=np.float64) for feature in BNCT_DAMAGE_FEATURES}

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with uproot.recreate(args.output) as root_file:
        root_file[args.tree_name] = arrays

    print(f"Wrote {n} rows to {args.output}:{args.tree_name}")


if __name__ == "__main__":
    main()
