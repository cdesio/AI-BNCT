"""Small end-to-end test of the current DNA replay ROOT schema."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import uproot


ROOT = Path(__file__).resolve().parent


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory)
        sugar = np.zeros((2, 12), dtype=np.float32)
        sugar[0, :3] = (0, 0, 0)
        sugar[0, 3:6] = (0, 0, 1)
        sugar[1, :3] = (0, 0, 2)
        sugar[1, 3:6] = (0, 0, 3)
        sugar.tofile(base / "sugar.bin")

        provenance = {
            "upstream_eventID": np.array([7, 7, 7], dtype=np.int32),
            "upstream_voxelID": np.array([0, 0, 4], dtype=np.int32),
            "upstream_particleID": np.array([3, 1, 3], dtype=np.int32),
            "upstream_primaryID": np.array([3, 3, 3], dtype=np.int32),
            "upstream_seedID": np.array([42, 42, 42], dtype=np.int32),
            "upstream_trackID": np.array([1, 2, 1], dtype=np.int32),
            "upstream_parentID": np.array([0, 1, 0], dtype=np.int32),
        }
        events = np.array([0, 1, 2], dtype=np.int32)
        with uproot.recreate(base / "dna.root") as output:
            output["ntuple/PS_data"] = {"EventNum": events, **provenance}
            output["ntuple/EventEdep"] = {
                "EventNum": events, "Edep_J": np.array([1e-17, 1e-18, 0.0]), **provenance,
            }
            output["ntuple/Direct"] = {
                "EventNum": np.array([0, 1], dtype=np.int32),
                "x": np.array([0.0, 0.0]), "y": np.array([0.0, 0.0]),
                "z": np.array([0.0, 2.0]), "eDep_eV": np.array([40.0, 40.0]),
            }
            output["ntuple/Indirect"] = {
                "EventNum": np.array([0], dtype=np.int32),
                "x": np.array([0.0]), "y": np.array([0.0]), "z": np.array([1.0]),
                "DNAmolecule": np.array(["Deoxyribose^0"]),
                "radical": np.array(["OH^0"]),
            }
            output["ntuple/Info"] = {
                "ChromatinVolume_m3": np.array([27e-21]),
                "NumBasepairs": np.array([2.0]),
            }

        output = base / "damage.csv"
        subprocess.run([
            sys.executable, str(ROOT / "run.py"), str(base / "dna.root"),
            "--sugar", str(base / "sugar.bin"), "--output", str(output),
            "--nx", "2", "--ny", "2", "--nz", "2", "--indirect-probability", "1",
        ], check=True)
        rows = read_rows(output)
        assert len(rows) == 2
        assert rows[0]["total_dsb"] == "1"
        assert rows[0]["total_cdsb"] == "1"
        assert rows[0]["upstream_eventID"] == "7"
        assert rows[0]["replay_events"] == "2"
        assert rows[0]["upstream_trackIDs"] == "1;2"
        assert rows[1]["total_dsb"] == "0"
        assert rows[1]["iz"] == "1"
        z_rows = read_rows(base / "damage_by_z.csv")
        assert len(z_rows) == 2
        assert [row["replay_events"] for row in z_rows] == ["2", "1"]

        bnct_output = base / "bnct_damage.csv"
        subprocess.run([
            sys.executable, str(ROOT / "run.py"), str(base / "dna.root"),
            "--sugar", str(base / "sugar.bin"), "--output", str(bnct_output),
            "--nx", "2", "--ny", "2", "--nz", "2",
            "--damage-preset", "molecular-bnct",
        ], check=True)
        assert read_rows(bnct_output)[0]["total_cdsb"] == "1"


if __name__ == "__main__":
    main()
