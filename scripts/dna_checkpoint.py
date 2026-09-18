"""Run resumable DNA record batches and merge their ROOT output."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def batches(args):
    size = args.input.stat().st_size
    width = args.record_doubles * 8
    if size == 0 or size % width:
        raise ValueError("Phase-space file does not contain complete records")
    total = size // width
    return [(start, min(args.checkpoint_events, total - start))
            for start in range(0, total, args.checkpoint_events)]


def paths(args, start):
    stem = args.output.with_suffix("")
    chunk = Path(f"{stem}.part{start:010d}.root")
    return chunk, chunk.with_suffix(".done")


def signature(args, start, count):
    stat = args.input.stat()
    return {"input": str(args.input.resolve()), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "start": start, "count": count,
            "record_doubles": args.record_doubles, "seed": args.seed,
            "macro": str(args.macro.resolve()),
            "macro_mtime_ns": args.macro.stat().st_mtime_ns,
            "executable": str(args.executable.resolve()),
            "executable_mtime_ns": args.executable.stat().st_mtime_ns}


def completed(args, start, count):
    chunk, marker = paths(args, start)
    if not marker.exists() or not chunk.is_file() or chunk.stat().st_size == 0:
        return False
    return json.loads(marker.read_text()) == signature(args, start, count)


def run(args):
    deadline = time.monotonic() + args.budget_seconds
    for start, count in batches(args):
        if completed(args, start, count):
            continue
        remaining = deadline - time.monotonic() - args.reserve_seconds
        if remaining <= 0:
            return 75
        chunk, marker = paths(args, start)
        marker.unlink(missing_ok=True)
        command = [str(args.executable), "-mac", str(args.macro), "-in", str(args.input),
                   "-out", str(chunk), "-seed", str(args.seed),
                   "-start-record", str(start), "-record-count", str(count),
                   "-record-doubles", str(args.record_doubles)]
        print(f"DNA records {start}..{start + count - 1}", flush=True)
        try:
            subprocess.run(command, check=True, timeout=remaining, cwd=args.executable.parent)
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"DNA batch at record {start} exceeded the job budget; reduce --checkpoint-events"
            )
        if not chunk.is_file() or chunk.stat().st_size == 0:
            raise RuntimeError(f"DNA did not produce {chunk}")
        temporary = marker.with_suffix(".done.tmp")
        temporary.write_text(json.dumps(signature(args, start, count)))
        os.replace(temporary, marker)
    return 0


def merge(args):
    chunks = []
    for start, count in batches(args):
        if not completed(args, start, count):
            raise RuntimeError(f"Missing or stale DNA checkpoint at record {start}")
        chunks.append(str(paths(args, start)[0]))
    temporary = args.output.with_name(args.output.stem + ".merging.root")
    subprocess.run(["hadd", "-f", str(temporary), *chunks], check=True)
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError("ROOT merge produced no output")
    os.replace(temporary, args.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "merge"))
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--macro", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--checkpoint-events", type=int, required=True)
    parser.add_argument("--record-doubles", type=int, choices=(16, 18), default=18)
    parser.add_argument("--budget-seconds", type=int, default=23 * 3600)
    parser.add_argument("--reserve-seconds", type=int, default=300)
    args = parser.parse_args()
    if args.checkpoint_events <= 0 or args.budget_seconds <= args.reserve_seconds:
        parser.error("Checkpoint size and time budget must be positive")
    try:
        return run(args) if args.action == "run" else merge(args)
    except (OSError, ValueError, subprocess.CalledProcessError, RuntimeError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    sys.exit(main())
