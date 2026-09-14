import argparse
import glob
import os
import re

import numpy as np
import uproot


DEFAULT_FEATURES = [
    "Energy",
    "Position_x_um",
    "Position_y_um",
    "Position_z_um",
    "Size_nm",
    "FragmentLength",
    "BaseDamage",
    "StrandDamage",
    "DirectBreaks",
    "IndirectBreaks",
    "EaqBaseHits",
    "EaqStrandHits",
    "OHBaseHits",
    "OHStrandHits",
    "HBaseHits",
    "HStrandHits",
    "EnergyDeposited_eV",
    "InducedBreaks",
    "TypeClassificationInt",
    "SourceClassificationInt",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Flatten MolecularBNCT damage ntuples into one StEG-ready ROOT tree."
    )
    parser.add_argument(
        "--input-glob",
        default="/Users/yw18581/work/MolecularBNCT-main/bnct_campaign/*/*.root",
    )
    parser.add_argument("--output", default="data/molecular_bnct_damage.root")
    parser.add_argument("--tree-name", default="Events")
    parser.add_argument("--max-files", type=int)
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--features", nargs="+", default=DEFAULT_FEATURES)
    return parser.parse_args()


def find_tree(root_file, target_name):
    for key in root_file.keys(recursive=True):
        clean_key = re.sub(r";\d+$", "", key)
        if clean_key.split("/")[-1] == target_name:
            return root_file[key]
    return None


def arrays_from_file(path, features):
    try:
        root_file = uproot.open(path)
    except Exception as exc:
        print(f"skip unreadable: {path} ({exc})")
        return None

    tree = find_tree(root_file, "damage")
    if tree is None:
        return None

    available = set(tree.keys())
    missing = [feature for feature in features if feature not in available]
    if missing:
        raise KeyError(f"{path} is missing required damage branches: {missing}")

    arrays = tree.arrays(features, library="np")
    if not arrays or len(next(iter(arrays.values()))) == 0:
        return None

    return {feature: np.asarray(arrays[feature]) for feature in features}


def append_arrays(target, source, max_rows):
    if source is None:
        return

    current_rows = 0 if not target else sum(len(chunk) for chunk in next(iter(target.values())))
    source_rows = len(next(iter(source.values())))
    if max_rows is not None:
        source_rows = min(source_rows, max_rows - current_rows)
    if source_rows <= 0:
        return

    for feature, values in source.items():
        target.setdefault(feature, []).append(values[:source_rows])


def main():
    args = parse_args()
    paths = sorted(glob.glob(args.input_glob))
    if args.max_files is not None:
        paths = paths[: args.max_files]

    collected = {}
    files_used = 0
    for path in paths:
        arrays = arrays_from_file(path, args.features)
        if arrays is None:
            continue

        append_arrays(collected, arrays, args.max_rows)
        files_used += 1
        total_rows = sum(len(chunk) for chunk in next(iter(collected.values())))
        print(f"used {files_used} files, collected {total_rows} rows", end="\r")

        if args.max_rows is not None and total_rows >= args.max_rows:
            break

    if not collected:
        raise RuntimeError(
            "No damage rows were found. Check that the ROOT files are populated and contain a damage ntuple."
        )

    output_arrays = {
        feature: np.concatenate(chunks).astype(np.float64)
        for feature, chunks in collected.items()
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with uproot.recreate(args.output) as root_file:
        root_file[args.tree_name] = output_arrays

    print()
    print(f"Files used: {files_used}")
    print(f"Rows written: {len(next(iter(output_arrays.values())))}")
    print(f"Output: {args.output}:{args.tree_name}")


if __name__ == "__main__":
    main()
