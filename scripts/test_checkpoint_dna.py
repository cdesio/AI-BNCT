import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

from checkpoint_dna import merge, validate


class CheckpointTests(unittest.TestCase):
    def test_merge_and_reject_incomplete(self):
        records = np.zeros((4, 18))
        records[:, 7] = [8, 8, 9, 9]
        records[:, 9] = [1, 2, 1, 2]
        records[:, 15] = 1234
        records[:, 16] = 1
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            parts = []
            for index in range(2):
                offset = index * 2
                path = directory / f"batch{index}.root"
                batch = records[offset:offset + 2]
                with uproot.recreate(path) as root:
                    rows = {"EventNum": np.arange(2, dtype="int32")}
                    for column, field in (("upstream_eventID", 7), ("upstream_voxelID", 9),
                                          ("upstream_seedID", 15), ("upstream_trackID", 16)):
                        rows[column] = batch[:, field].astype("int32")
                    for name in ("EventEdep", "PS_data", "VoxelExit"):
                        root.mktree(f"ntuple/{name}", rows)
                    root.mktree("ntuple/Direct", {"EventNum": "int32"})
                    root.mktree("ntuple/Indirect", {"EventNum": "int32", "radical": "string"})
                    info = root.mktree("ntuple/Info", {"G4Version": "string"})
                    info.extend({"G4Version": ["test"]})
                validate(path, batch, True)
                with self.assertRaises(ValueError):
                    validate(path, records, True)
                parts.append((path, offset))
            temporary = merge(parts, directory / "merged.root")
            validate(temporary, records, True)
            with uproot.open(temporary) as root:
                self.assertEqual(root["ntuple/Indirect"].num_entries, 0)
                self.assertEqual(root["ntuple/Info"].num_entries, 2)
                self.assertEqual(root["ntuple/VoxelExit"].num_entries, 4)


if __name__ == "__main__":
    unittest.main()
