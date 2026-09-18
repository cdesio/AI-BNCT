import argparse
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import dna_checkpoint


class CheckpointTests(unittest.TestCase):
    def test_resume_and_merge_require_every_batch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            phase = root / "phase.bin"
            phase.write_bytes(b"\0" * (5 * 18 * 8))
            macro = root / "dna.mac"
            macro.write_text("/run/initialize\n")
            executable = root / "rbe"
            executable.write_bytes(b"exe")
            args = argparse.Namespace(input=phase, output=root / "dna.root",
                                      executable=executable, macro=macro, seed=1234,
                                      checkpoint_events=2, record_doubles=18,
                                      budget_seconds=3600, reserve_seconds=300)
            calls = []

            def fake_run(command, **kwargs):
                calls.append(command)
                if command[0] == "hadd":
                    Path(command[2]).write_bytes(b"merged")
                else:
                    Path(command[command.index("-out") + 1]).write_bytes(b"root")

            with patch.object(dna_checkpoint.subprocess, "run", side_effect=fake_run):
                self.assertEqual(dna_checkpoint.run(args), 0)
                self.assertEqual(len(calls), 3)
                self.assertEqual([call[call.index("-seed") + 1] for call in calls],
                                 ["1234", "1234", "1234"])
                self.assertEqual(dna_checkpoint.run(args), 0)
                self.assertEqual(len(calls), 3)
                dna_checkpoint.merge(args)
                self.assertEqual(args.output.read_bytes(), b"merged")
                self.assertEqual(len(calls[-1]), 6)

                paths = dna_checkpoint.paths(args, 2)
                paths[1].unlink()
                with self.assertRaisesRegex(RuntimeError, "record 2"):
                    dna_checkpoint.merge(args)
                self.assertEqual(dna_checkpoint.run(args), 0)
                self.assertEqual(len(calls), 5)


if __name__ == "__main__":
    unittest.main()
