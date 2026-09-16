#!/usr/bin/env python3
"""Convert bnctVoxelPS ROOT phase-space rows to DNA-simulation binary PS.

The output format is AlphaGlue/TAT-style 18-double records:

  x, y, z, dir_x, dir_y, dir_z, energy, event_id, particle_id,
  voxel_id, time, source_id, world_x, world_y, world_z, seed_id,
  track_id, parent_id

Geant4 internal units are used in the binary file: mm, MeV, ns.
The input ROOT file stores positions in nm and kinetic energy in eV.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path
from typing import Any


PARTICLE_IDS = {
    "e-": 1,
    "gamma": 2,
    "alpha": 3,
    "alpha+": 3,
    "helium": 3,
    "e+": 11,
    "li7": 53,
    "li-7": 53,
    "lithium": 53,
    "lithium+": 53,
    "lithium++": 53,
    "lithium+++": 53,
}


def _decode(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _particle_id(name: str, pdg: int) -> int:
    key = name.strip().lower()
    if key in PARTICLE_IDS:
        return PARTICLE_IDS[key]
    if pdg == 1000020040:
        return 3
    if pdg == 1000030070:
        return 53
    if pdg == 11:
        return 1
    if pdg == -11:
        return 11
    if pdg == 22:
        return 2
    raise ValueError(f"Unsupported particle for DNA replay: name={name!r}, PDG={pdg}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_bin", type=Path)
    parser.add_argument(
        "--tree",
        default="Ntuples/phase_space",
        help="ROOT tree path. Default: Ntuples/phase_space",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=0,
        help="Convert at most this many rows. Default: all rows.",
    )
    parser.add_argument(
        "--nudge-nm",
        type=float,
        default=0.0,
        help="Move local start position by this many nm along the direction, useful for boundary starts.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "This converter needs uproot. Install it in the Python environment "
            "you use for conversion, then rerun this script."
        ) from exc

    tree = uproot.open(args.input_root)[args.tree]
    branches = [
        "EventID",
        "SeedID",
        "TrackID",
        "ParentID",
        "ParticleType",
        "PDGCode",
        "PositionX_nm",
        "PositionY_nm",
        "PositionZ_nm",
        "LocalPositionX_nm",
        "LocalPositionY_nm",
        "LocalPositionZ_nm",
        "DirectionX",
        "DirectionY",
        "DirectionZ",
        "KineticEnergy_eV",
        "Time_ns",
        "VoxelID",
    ]
    arrays = tree.arrays(branches, library="np")
    n_rows = len(arrays["EventID"])
    if args.max_rows > 0:
        n_rows = min(n_rows, args.max_rows)

    args.output_bin.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.output_bin.open("wb") as output:
        for i in range(n_rows):
            name = _decode(arrays["ParticleType"][i])
            pdg = int(arrays["PDGCode"][i])
            particle_id = _particle_id(name, pdg)

            dir_x = float(arrays["DirectionX"][i])
            dir_y = float(arrays["DirectionY"][i])
            dir_z = float(arrays["DirectionZ"][i])
            nudge = args.nudge_nm

            local_x_nm = float(arrays["LocalPositionX_nm"][i]) + nudge * dir_x
            local_y_nm = float(arrays["LocalPositionY_nm"][i]) + nudge * dir_y
            local_z_nm = float(arrays["LocalPositionZ_nm"][i]) + nudge * dir_z

            record = [
                local_x_nm * 1e-6,
                local_y_nm * 1e-6,
                local_z_nm * 1e-6,
                dir_x,
                dir_y,
                dir_z,
                float(arrays["KineticEnergy_eV"][i]) * 1e-6,
                float(arrays["EventID"][i]),
                float(particle_id),
                float(arrays["VoxelID"][i]),
                float(arrays["Time_ns"][i]),
                float(particle_id),
                float(arrays["PositionX_nm"][i]) * 1e-6,
                float(arrays["PositionY_nm"][i]) * 1e-6,
                float(arrays["PositionZ_nm"][i]) * 1e-6,
                float(arrays["SeedID"][i]),
                float(arrays["TrackID"][i]),
                float(arrays["ParentID"][i]),
            ]
            output.write(struct.pack("<18d", *record))
            written += 1

    print(f"Wrote {written} DNA phase-space records to {args.output_bin}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
