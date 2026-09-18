"""Generate restartable BC5 jobs for voxel transport, DNA replay, and clustering."""

import argparse
from pathlib import Path
import re
import shlex


SOURCE_ROOT = Path(__file__).resolve().parents[1]
SUGAR_NAME = "sugarPos_n2_300nm_H75nm_R11nm_33hist_50deg.bin"
HISTONE_NAME = "histonePos_n2_300nm_H75nm_R11nm_33hist_50deg.bin"


def shell(value):
    return shlex.quote(str(value))


def bc5_walltime(value):
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):([0-5]\d):([0-5]\d)", value)
    if not match:
        raise argparse.ArgumentTypeError("Use HH:MM:SS or D-HH:MM:SS")
    days, hours, minutes, seconds = (int(part or 0) for part in match.groups())
    duration = days * 86400 + hours * 3600 + minutes * 60 + seconds
    if not 0 < duration <= 86400:
        raise argparse.ArgumentTypeError("BC5 walltime must be at most 24 hours")
    return value


def slurm_header(stage, run_dir, prefix, hours, cpus, memory, email):
    lines = [
        "#!/bin/bash --login",
        f"#SBATCH --job-name={stage}_{prefix}",
        f"#SBATCH --output={run_dir}/logs/{stage}.out.%J",
        f"#SBATCH --error={run_dir}/logs/{stage}.err.%J",
        f"#SBATCH --time={hours}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={memory}",
    ]
    if email:
        lines.extend((f"#SBATCH --mail-user={email}", "#SBATCH --mail-type=FAIL,END,TIME_LIMIT"))
    return "\n".join(lines) + "\n\nset -euo pipefail\n"


def geant4_setup(args):
    setup = f"module use {shell(args.modulefiles_dir)}\nmodule load {shell(args.geant4_module)}\n"
    if args.geant4_sh:
        setup += f"source {shell(args.geant4_sh)}\n"
    return setup


def make_macro(particle, events):
    source = SOURCE_ROOT / "bnct-voxel-ps/macros" / f"{particle}_300nm_4x4x40.mac"
    macro = source.read_text()
    macro, count = re.subn(r"(?m)^/run/beamOn\s+\d+\s*$", f"/run/beamOn {events}", macro)
    if count != 1:
        raise ValueError(f"Expected one /run/beamOn command in {source}")
    return macro.rstrip() + "\n"


def make_dna_macro(threads):
    source = SOURCE_ROOT / "bnct-dna-simulation/rbe_PS.in"
    macro, count = re.subn(r"(?m)^/run/numberOfThreads\s+\d+\s*$",
                           f"/run/numberOfThreads {threads}", source.read_text())
    if count != 1:
        raise ValueError(f"Expected one thread command in {source}")
    return macro.rstrip() + "\n"


def make_scripts(args, run_dir, prefix):
    root = args.project_root
    upstream_exe = args.upstream_exe or root / "bnct-voxel-ps/build/bnctVoxelPS"
    dna_exe = args.dna_exe or root / "bnct-dna-simulation/build/rbe"
    phase_stem = run_dir / f"{prefix}_ps"
    phase_bin = phase_stem.with_suffix(".bin")
    phase_root = phase_stem.with_suffix(".root")
    dna_root = run_dir / f"{prefix}_dna.root"
    damage_csv = run_dir / f"{prefix}_damage.csv"

    upstream = slurm_header("ps", run_dir, prefix, args.upstream_time, 1,
                            args.upstream_mem, args.mail_user)
    upstream += geant4_setup(args)
    upstream += f"\ncd {shell(root)}\n"
    upstream += f"test -x {shell(upstream_exe)}\n"
    upstream += f"time {shell(upstream_exe)} "
    upstream += f"-mac {shell(run_dir / 'upstream.mac')} -out {shell(phase_stem)} -seed {args.seed}\n"
    upstream += f"test -s {shell(phase_bin)}\ntest -s {shell(phase_root)}\n"

    dna = slurm_header("dna", run_dir, prefix, args.dna_time, args.dna_cpus,
                       args.dna_mem, args.mail_user)
    dna += geant4_setup(args)
    dna += f"\ntest -s {shell(phase_bin)}\n"
    dna += f"test -f {shell(root / 'bnct-dna-simulation/geometryFiles' / SUGAR_NAME)}\n"
    dna += f"test -f {shell(root / 'bnct-dna-simulation/geometryFiles' / HISTONE_NAME)}\n"
    dna += f"cd {shell(dna_exe.parent)}\n"
    dna += f"test -x {shell(dna_exe)}\n"
    runner = root / "scripts/dna_checkpoint.py"
    runner_args = (f"--input {shell(phase_bin)} --output {shell(dna_root)} "
                   f"--executable {shell(dna_exe)} --macro {shell(run_dir / 'dna.mac')} "
                   f"--seed {args.seed} --checkpoint-events {args.checkpoint_events}")
    dna += "set +e\n"
    dna += f"python {shell(runner)} run {runner_args} --budget-seconds {args.dna_budget_seconds}\n"
    dna += "status=$?\nset -e\n"
    dna += f"if [[ $status -eq 75 ]]; then sbatch {shell(run_dir / 'dna.sbatch')}; exit 0; fi\n"
    dna += "if [[ $status -ne 0 ]]; then exit \"$status\"; fi\n"
    dna += f"sbatch {shell(run_dir / 'merge.sbatch')}\n"

    merge = slurm_header("merge", run_dir, prefix, args.clustering_time, 1,
                         args.clustering_mem, args.mail_user)
    merge += f"module use {shell(args.modulefiles_dir)}\nmodule load apps/root/6.26.00\n"
    merge += f"python {shell(runner)} merge {runner_args}\n"
    merge += f"test -s {shell(dna_root)}\n"
    merge += f"sbatch {shell(run_dir / 'clustering.sbatch')}\n"

    clustering = slurm_header("cluster", run_dir, prefix, args.clustering_time, 1,
                              args.clustering_mem, args.mail_user)
    clustering += f"module use {shell(args.modulefiles_dir)}\n"
    clustering += "module load clustering/conda\n"
    clustering += f"test -s {shell(dna_root)}\n"
    clustering += f"cd {shell(root)}\n"
    clustering += f"time python {shell(root / 'bnct-clustering/run.py')} {shell(dna_root)} "
    clustering += f"--output {shell(damage_csv)} --seed {args.seed} "
    clustering += f"--damage-preset {args.damage_preset}\n"
    clustering += f"test -s {shell(damage_csv)}\n"
    clustering += f"test -s {shell(run_dir / (prefix + '_damage_by_z.csv'))}\n"

    launcher = """#!/bin/bash
set -euo pipefail
RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
start="${1:-upstream}"
case "$start" in
  upstream)
    ps_job=$(sbatch --parsable "$RUN_DIR/upstream.sbatch")
    ps_job=${ps_job%%;*}
    echo "Submitted upstream: $ps_job"
    ;;
  dna|merge|clustering) ;;
  *) echo "Usage: $0 [upstream|dna|merge|clustering]" >&2; exit 2 ;;
esac
if [[ "$start" == merge ]]; then
  sbatch "$RUN_DIR/merge.sbatch"
elif [[ "$start" != clustering ]]; then
  if [[ "$start" == upstream ]]; then
    dna_job=$(sbatch --parsable --dependency="afterok:$ps_job" "$RUN_DIR/dna.sbatch")
  else
    dna_job=$(sbatch --parsable "$RUN_DIR/dna.sbatch")
  fi
  dna_job=${dna_job%%;*}
  echo "Submitted DNA: $dna_job"
else
  sbatch "$RUN_DIR/clustering.sbatch"
fi
"""
    return {
        "upstream.mac": make_macro(args.particle, args.events),
        "dna.mac": make_dna_macro(args.dna_cpus),
        "upstream.sbatch": upstream,
        "dna.sbatch": dna,
        "merge.sbatch": merge,
        "clustering.sbatch": clustering,
        "submit_all.sh": launcher,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particle", choices=("alpha", "lithium"), required=True)
    parser.add_argument("--events", type=int, required=True, help="Number of upstream primaries")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--name", help="Run label (default: particle name)")
    parser.add_argument("--project-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--upstream-exe", type=Path, help="Upstream executable (default: bnct-voxel-ps/build/bnctVoxelPS)")
    parser.add_argument("--dna-exe", type=Path, help="DNA executable (default: bnct-dna-simulation/build/rbe)")
    parser.add_argument("--output-dir", type=Path, help="Directory for generated jobs and run outputs")
    parser.add_argument("--mail-user", default="yw18581@bristol.ac.uk")
    parser.add_argument("--damage-preset", choices=("alphaglue", "molecular-bnct"), default="alphaglue")
    parser.add_argument("--modulefiles-dir", type=Path,
                        default=Path("/projects/b56v/software/modulefiles"))
    parser.add_argument("--geant4-module", default="geant4/11.3.0-lithium")
    parser.add_argument("--geant4-sh", default="/projects/b56v/software/geant4-v11.3.0-lithium-install/bin/geant4.sh",
                        help="Geant4 setup script; pass an empty string if the module sets the environment")
    parser.add_argument("--upstream-time", type=bc5_walltime, default="24:00:00")
    parser.add_argument("--dna-time", type=bc5_walltime, default="24:00:00")
    parser.add_argument("--clustering-time", type=bc5_walltime, default="08:00:00")
    parser.add_argument("--upstream-mem", default="16GB")
    parser.add_argument("--dna-mem", default="100GB")
    parser.add_argument("--clustering-mem", default="32GB")
    parser.add_argument("--dna-cpus", type=int, default=4)
    parser.add_argument("--checkpoint-events", type=int, default=10000,
                        help="DNA phase-space records per durable ROOT checkpoint")
    parser.add_argument("--dna-budget-seconds", type=int, default=23 * 3600,
                        help="Maximum DNA work time per Slurm job")
    args = parser.parse_args()
    if args.events <= 0 or args.seed <= 0 or args.dna_cpus <= 0 or args.checkpoint_events <= 0:
        parser.error("Events, seed, DNA CPUs, and checkpoint size must be positive")
    if args.dna_budget_seconds <= 300 or args.dna_budget_seconds >= 86400:
        parser.error("DNA budget must be between 300 and 86400 seconds")
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):([0-5]\d):([0-5]\d)", args.dna_time)
    days, hours, minutes, seconds = (int(part or 0) for part in match.groups())
    dna_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
    if args.dna_budget_seconds + 300 >= dna_seconds:
        parser.error("DNA budget needs at least 5 minutes below the DNA walltime")
    name = args.name or args.particle
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        parser.error("Run name must contain only letters, numbers, '_' or '-'")
    prefix = f"{name}_seed{args.seed}"
    args.project_root = args.project_root.resolve()
    if args.upstream_exe:
        args.upstream_exe = args.upstream_exe.resolve()
    if args.dna_exe:
        args.dna_exe = args.dna_exe.resolve()
    run_dir = (args.output_dir or args.project_root / "jobs" / prefix).resolve()
    if any(char.isspace() for char in str(args.project_root) + str(run_dir)):
        parser.error("Project and output paths cannot contain whitespace (SLURM log paths)")
    scripts = make_scripts(args, run_dir, prefix)
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "logs").mkdir()
    for name, content in scripts.items():
        target = run_dir / name
        target.write_text(content)
        if name == "submit_all.sh":
            target.chmod(0o755)
    print(f"Created {run_dir}")
    print(f"Submit: {run_dir / 'submit_all.sh'}")


if __name__ == "__main__":
    main()
