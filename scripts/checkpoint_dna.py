"""Replay 18-double phase space in restartable batches, then merge DNA ntuples."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import uproot


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate(path, records, save_exits):
    with uproot.open(path) as root:
        for name in ("EventEdep", "PS_data") + (("VoxelExit",) if save_exits else ()):
            rows = root[f"ntuple/{name}"].arrays(library="np")
            ids = rows["EventNum"]
            if not np.array_equal(np.sort(ids), np.arange(len(records))):
                raise ValueError(f"Incomplete or duplicate events in {path}: {name}")
            order = np.argsort(ids)
            for column, index in (("upstream_eventID", 7), ("upstream_voxelID", 9),
                                  ("upstream_seedID", 15), ("upstream_trackID", 16)):
                if not np.array_equal(rows[column][order], records[:, index]):
                    raise ValueError(f"Input mismatch in {path}: {name}/{column}")
        for name in ("Direct", "Indirect", "Info"):
            root[f"ntuple/{name}"]


def merge(parts, output):
    temporary = output.with_suffix(".partial.root")
    with uproot.recreate(temporary) as target:
        trees = {}
        for path, offset in parts:
            with uproot.open(path) as source:
                for name, tree in source["ntuple"].items(cycle=False):
                    key = f"ntuple/{name}"
                    for rows in tree.iterate(step_size="32 MB", library="np"):
                        if "EventNum" in rows:
                            rows["EventNum"] = rows["EventNum"] + offset
                        if key not in trees:
                            schema = {k: ("string" if v.dtype.kind in "OUS" else v.dtype)
                                      for k, v in rows.items()}
                            trees[key] = target.mktree(key, schema)
                        if len(next(iter(rows.values()))):
                            trees[key].extend(rows)
                    # Empty damage trees must remain present for clustering.
                    if key not in trees:
                        rows = tree.arrays(entry_stop=0, library="np")
                        trees[key] = target.mktree(key, {k: ("string" if v.dtype.kind in "OUS" else v.dtype)
                                                        for k, v in rows.items()})
    return temporary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--macro", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--save-exits", action="store_true")
    args = parser.parse_args()
    if args.batch_size <= 0:
        parser.error("Batch size must be positive")
    records = np.fromfile(args.input, dtype="<f8").reshape(-1, 18)
    if not len(records):
        raise ValueError("Empty phase space")
    directory = args.output.parent / (args.output.stem + "_checkpoints")
    directory.mkdir(exist_ok=True)
    identity = {"input": digest(args.input), "executable": digest(args.exe),
                "macro": digest(args.macro), "seed": args.seed,
                "batch_size": args.batch_size, "save_exits": args.save_exits}
    manifest = directory / "manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError("Checkpoint configuration changed; use a new output directory")
    else:
        manifest.write_text(json.dumps(identity, indent=2) + "\n")
    parts = []
    for index, offset in enumerate(range(0, len(records), args.batch_size)):
        batch = records[offset:offset + args.batch_size]
        finished = directory / f"batch_{index:05d}.root"
        if finished.exists():
            validate(finished, batch, args.save_exits)
            print(f"Skipping verified batch {index}", flush=True)
        else:
            binary = directory / f"batch_{index:05d}.bin"
            batch.tofile(binary)
            temporary = directory / f"batch_{index:05d}.partial.root"
            command = [str(args.exe.resolve()), "-mac", str(args.macro.resolve()),
                       "-in", str(binary.resolve()), "-out", str(temporary.resolve()),
                       "-seed", str(args.seed + index)]
            if args.save_exits:
                command.append("--save-exits")
            started = time.monotonic()
            with (directory / f"batch_{index:05d}.log").open("w") as log:
                subprocess.run(command, cwd=args.exe.resolve().parent,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            validate(temporary, batch, args.save_exits)
            temporary.replace(finished)
            seconds = time.monotonic() - started
            (directory / f"batch_{index:05d}.timing.json").write_text(
                json.dumps({"events": len(batch), "seconds": seconds,
                            "seed": args.seed + index}) + "\n")
            print(f"Completed batch {index}: {len(batch)} events, {seconds:.1f}s", flush=True)
        parts.append((finished, offset))
    temporary = merge(parts, args.output)
    validate(temporary, records, args.save_exits)
    temporary.replace(args.output)
    print(f"Verified merged output: {args.output}", flush=True)


if __name__ == "__main__":
    main()
