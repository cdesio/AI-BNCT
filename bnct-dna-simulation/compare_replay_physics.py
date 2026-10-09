"""Compare replay physical records by local event ID, including optional exits."""
import argparse
import json
import numpy as np
import uproot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("candidate")
    args = parser.parse_args()
    reference, candidate = uproot.open(args.reference), uproot.open(args.candidate)
    result = {}
    for path in ["ntuple/EventEdep", "ntuple/PS_data", "ntuple/VoxelExit"]:
        if path not in reference or path not in candidate:
            result[path] = {"available_in_both": False}
            continue
        a, b = reference[path].arrays(library="np"), candidate[path].arrays(library="np")
        ai, bi = np.argsort(a["EventNum"]), np.argsort(b["EventNum"])
        ids, ia, ib = np.intersect1d(a["EventNum"][ai], b["EventNum"][bi], return_indices=True)
        unequal = [key for key in a if key not in b or not np.array_equal(a[key][ai][ia], b[key][bi][ib])]
        result[path] = {"reference_rows": len(ai), "candidate_rows": len(bi),
                        "common_events": len(ids), "unequal_columns": unequal}
    result["versions"] = {
        label: list(np.unique(file["ntuple/Info"]["G4Version"].array(library="np")))
        for label, file in [("reference", reference), ("candidate", candidate)]
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
