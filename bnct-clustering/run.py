"""Cluster strand breaks for BNCT DNA voxel replay output."""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree
import uproot


ROOT = Path(__file__).resolve().parent
DEFAULT_SUGAR = ROOT.parent / "bnct-dna-simulation/geometryFiles/sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin"
RESULT_NAMES = (
    "direct_sb", "direct_ssb", "direct_cssb", "direct_dsb",
    "indirect_sb", "indirect_ssb", "indirect_cssb", "indirect_dsb",
    "total_sb", "total_ssb", "total_cssb", "total_dsb",
    "direct_cdsb", "indirect_cdsb", "hybrid_cdsb", "mixed_cdsb", "total_cdsb",
)
PROVENANCE = (
    "upstream_eventID", "upstream_voxelID", "upstream_particleID",
    "upstream_primaryID", "upstream_seedID", "upstream_trackID",
    "upstream_parentID",
)


def load_sugars(path):
    data = np.fromfile(path, dtype=np.float32)
    if data.size == 0 or data.size % 12:
        raise ValueError(f"Invalid sugar geometry: {path}")
    data = data.reshape(-1, 12)
    return cKDTree(data[:, :3]), cKDTree(data[:, 3:6]), len(data)


def read_tree(root, name, columns):
    tree = root[f"ntuple/{name}"]
    legacy = {"upstream_eventID": "part1_EventNum",
              "upstream_voxelID": "part1_CopyNum",
              "upstream_particleID": "part1_particleSource"}
    source = {column: column if column in tree.keys() else legacy.get(column, column)
              for column in columns}
    if ("upstream_primaryID" in columns and "upstream_primaryID" not in tree.keys()
            and "part1_particleSource" in tree.keys()):
        source["upstream_primaryID"] = "part1_particleSource"
    missing = set(source.values()) - set(tree.keys())
    if missing:
        raise ValueError(f"{name} is missing: {', '.join(sorted(missing))}")
    arrays = tree.arrays(list(dict.fromkeys(source.values())), library="np")
    return {column: arrays[branch] for column, branch in source.items()}


def nearest_sugars(xyz, trees):
    d0, i0 = trees[0].query(xyz)
    d1, i1 = trees[1].query(xyz)
    first = d0 < d1
    return np.minimum(d0, d1), np.where(first, i0, i1), np.where(first, 0, 1)


def direct_breaks(direct, trees, minimum, maximum, rng):
    if len(direct["EventNum"]) == 0:
        return []
    xyz = np.column_stack([direct[axis] for axis in "xyz"])
    distance, sugar, strand = nearest_sugars(xyz, trees)
    energy = {}
    for event, bp, side, edep, nearby in zip(
        direct["EventNum"], sugar, strand, direct["eDep_eV"], distance < 0.35
    ):
        if nearby:
            key = (int(event), int(bp), int(side))
            energy[key] = energy.get(key, 0.0) + float(edep)
    keys = sorted(energy)
    if minimum == maximum:
        return [key for key in keys if energy[key] >= minimum]
    return [key for key in keys if rng.random() < np.clip((energy[key] - minimum) / (maximum - minimum), 0, 1)]


def indirect_breaks(indirect, trees, probability, rng):
    if len(indirect["EventNum"]) == 0:
        return []
    selected = (indirect["DNAmolecule"] == "Deoxyribose^0") & (indirect["radical"] == "OH^0")
    selected &= rng.random(len(selected)) < probability
    if not np.any(selected):
        return []
    xyz = np.column_stack([indirect[axis][selected] for axis in "xyz"])
    distance, sugar, strand = nearest_sugars(xyz, trees)
    if np.any(distance > 1e-5):
        raise ValueError("An indirect deoxyribose reaction does not match the sugar geometry")
    return sorted({(int(event), int(bp), int(side)) for event, bp, side in zip(
        indirect["EventNum"][selected], sugar, strand
    )})


def run(args):
    sys.path.insert(0, str(ROOT / "build"))
    from clustering import clustering

    trees = load_sugars(args.sugar)
    with uproot.open(args.input) as root:
        event = read_tree(root, "EventEdep", ("EventNum", "Edep_J", *PROVENANCE))
        direct = read_tree(root, "Direct", ("EventNum", "x", "y", "z", "eDep_eV"))
        indirect = read_tree(root, "Indirect", ("EventNum", "x", "y", "z", "DNAmolecule", "radical"))
        info = read_tree(root, "Info", ("ChromatinVolume_m3", "NumBasepairs"))
        replay = read_tree(root, "PS_data", ("EventNum", *PROVENANCE))

    if int(info["NumBasepairs"][0]) != trees[2]:
        raise ValueError("Sugar geometry base-pair count differs from DNA simulation Info")
    if len(np.unique(replay["EventNum"])) != len(replay["EventNum"]):
        raise ValueError("PS_data contains duplicate DNA replay event numbers")
    if len(np.unique(event["EventNum"])) != len(event["EventNum"]):
        raise ValueError("EventEdep contains duplicate DNA replay event numbers")
    replay_provenance = {int(number): tuple(int(replay[name][index]) for name in PROVENANCE)
                         for index, number in enumerate(replay["EventNum"])}
    for index, number in enumerate(event["EventNum"]):
        current = tuple(int(event[name][index]) for name in PROVENANCE)
        if current != replay_provenance.get(int(number)):
            raise ValueError(f"EventEdep provenance differs from PS_data for event {number}")

    known = set(map(int, replay["EventNum"]))
    if not set(map(int, event["EventNum"])) <= known:
        raise ValueError("EventEdep contains events absent from PS_data")
    if not set(map(int, direct["EventNum"])) <= known or not set(map(int, indirect["EventNum"])) <= known:
        raise ValueError("Damage contains events absent from PS_data")
    groups = {}
    replay_group = {}
    for index, number in enumerate(replay["EventNum"]):
        key = tuple(int(replay[name][index]) for name in
                    ("upstream_seedID", "upstream_eventID", "upstream_voxelID"))
        group = groups.setdefault(key, {"dna_events": [], "particle_ids": set(),
                                        "primary_ids": set(), "track_ids": set(),
                                        "parent_ids": set(), "dose_Gy": 0.0})
        group["dna_events"].append(int(number))
        for column, field in (("upstream_particleID", "particle_ids"),
                              ("upstream_primaryID", "primary_ids"),
                              ("upstream_trackID", "track_ids"),
                              ("upstream_parentID", "parent_ids")):
            group[field].add(int(replay[column][index]))
        replay_group[int(number)] = key
    group_ids = {key: index for index, key in enumerate(groups)}
    for index, number in enumerate(event["EventNum"]):
        groups[replay_group[int(number)]]["dose_Gy"] += float(event["Edep_J"][index]) / (
            1000 * info["ChromatinVolume_m3"][0]
        )

    rng = np.random.default_rng(args.seed)
    direct["EventNum"] = np.array([group_ids[replay_group[int(number)]]
                                    for number in direct["EventNum"]], dtype=np.int64)
    indirect["EventNum"] = np.array([group_ids[replay_group[int(number)]]
                                      for number in indirect["EventNum"]], dtype=np.int64)
    d = direct_breaks(direct, trees[:2], args.min_energy, args.max_energy, rng)
    i = indirect_breaks(indirect, trees[:2], args.indirect_probability, rng)
    result = clustering(
        sorted({key[0] for key in d + i}),
        [key[0] for key in d], [key[1] for key in d], [key[2] for key in d],
        [key[0] for key in i], [key[1] for key in i], [key[2] for key in i],
        True,
    )
    clustered = {int(row[0]): row[1:] for row in result[0]}
    fields = ("upstream_seedID", "upstream_eventID", "upstream_voxelID",
              "replay_events", "dna_eventIDs", "upstream_particleIDs", "upstream_primaryIDs",
              "upstream_trackIDs", "upstream_parentIDs", "ix", "iy", "iz",
              "z_center_um", "dose_Gy", *RESULT_NAMES)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    z_totals = {}
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key, group in groups.items():
            seed, event_id, voxel = key
            if not 0 <= voxel < args.nx * args.ny * args.nz:
                raise ValueError(f"Voxel ID {voxel} is outside the configured grid")
            row = {"upstream_seedID": seed, "upstream_eventID": event_id,
                   "upstream_voxelID": voxel, "replay_events": len(group["dna_events"]),
                   "dna_eventIDs": ";".join(map(str, group["dna_events"]))}
            for field, source in (("upstream_particleIDs", "particle_ids"),
                                  ("upstream_primaryIDs", "primary_ids"),
                                  ("upstream_trackIDs", "track_ids"),
                                  ("upstream_parentIDs", "parent_ids")):
                row[field] = ";".join(map(str, sorted(group[source])))
            row["ix"] = voxel % args.nx
            row["iy"] = (voxel // args.nx) % args.ny
            row["iz"] = voxel // (args.nx * args.ny)
            row["z_center_um"] = (row["iz"] + 0.5) * args.voxel_size_nm / 1000
            row["dose_Gy"] = group["dose_Gy"]
            row.update(zip(RESULT_NAMES, clustered.get(group_ids[key], [0] * len(RESULT_NAMES))))
            writer.writerow(row)
            total = z_totals.setdefault(row["iz"], {name: 0 for name in RESULT_NAMES})
            total["voxel_events"] = total.get("voxel_events", 0) + 1
            total["replay_events"] = total.get("replay_events", 0) + len(group["dna_events"])
            total["dose_Gy"] = total.get("dose_Gy", 0.0) + row["dose_Gy"]
            for name in RESULT_NAMES:
                total[name] += row[name]

    z_output = output.with_name(output.stem + "_by_z.csv")
    z_fields = ("iz", "z_center_um", "voxel_events", "replay_events", "dose_Gy_sum", *RESULT_NAMES)
    with z_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=z_fields)
        writer.writeheader()
        for iz, total in sorted(z_totals.items()):
            writer.writerow({"iz": iz, "z_center_um": (iz + 0.5) * args.voxel_size_nm / 1000,
                             "voxel_events": total["voxel_events"], "replay_events": total["replay_events"],
                             "dose_Gy_sum": total["dose_Gy"],
                             **{name: total[name] for name in RESULT_NAMES}})
    print(f"Wrote {len(groups)} upstream voxel events from {len(replay['EventNum'])} DNA replays to {output} and {len(z_totals)} Z bins to {z_output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="DNA replay ROOT output")
    parser.add_argument("--output", type=Path, required=True, help="Per-voxel-event CSV")
    parser.add_argument("--sugar", type=Path, default=DEFAULT_SUGAR)
    parser.add_argument("--nx", type=int, default=4)
    parser.add_argument("--ny", type=int, default=4)
    parser.add_argument("--nz", type=int, default=40)
    parser.add_argument("--voxel-size-nm", type=float, default=300)
    parser.add_argument("--damage-preset", choices=("alphaglue", "molecular-bnct"), default="alphaglue")
    parser.add_argument("--min-energy", type=float)
    parser.add_argument("--max-energy", type=float)
    parser.add_argument("--indirect-probability", type=float)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args()
    default_parameters = {"alphaglue": (5.0, 37.5, 0.405), "molecular-bnct": (17.5, 17.5, 1.0)}
    defaults = default_parameters[args.damage_preset]
    args.min_energy = defaults[0] if args.min_energy is None else args.min_energy
    args.max_energy = defaults[1] if args.max_energy is None else args.max_energy
    args.indirect_probability = defaults[2] if args.indirect_probability is None else args.indirect_probability
    if min(args.nx, args.ny, args.nz, args.voxel_size_nm) <= 0:
        parser.error("Grid dimensions and voxel size must be positive")
    if not 0 <= args.indirect_probability <= 1 or args.max_energy < args.min_energy:
        parser.error("Invalid damage model parameters")
    run(args)


if __name__ == "__main__":
    main()
