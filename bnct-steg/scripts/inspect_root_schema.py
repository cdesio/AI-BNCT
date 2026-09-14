import argparse
import glob

import uproot


def parse_args():
    parser = argparse.ArgumentParser(description="Print ROOT keys and tree branches for candidate files.")
    parser.add_argument("paths", nargs="*", default=["data/molecular_bnct_damage.root"])
    parser.add_argument("--glob", dest="glob_pattern")
    parser.add_argument("--max-files", type=int, default=5)
    return parser.parse_args()


def main():
    args = parse_args()
    paths = sorted(glob.glob(args.glob_pattern)) if args.glob_pattern else args.paths
    for path in paths[: args.max_files]:
        print(path)
        root_file = uproot.open(path)
        keys = root_file.keys(recursive=True)
        if not keys:
            print("  no keys")
            continue
        for key in keys:
            obj = root_file[key]
            if hasattr(obj, "keys"):
                print(f"  {key}: {getattr(obj, 'num_entries', '?')} entries")
                for branch in obj.keys():
                    print(f"    {branch}")
            else:
                print(f"  {key}")


if __name__ == "__main__":
    main()
