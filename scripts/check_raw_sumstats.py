#!/usr/bin/env python3
"""
scripts/check_raw_sumstats.py
-----------------------------
Pre-flight check for Stage A: for every trait x population in a
config_prsxtra_*.yaml, confirm the raw GWAS file exists, is readable, and has
every column named in its merged spec (defaults <- files[pop], as in
Snakefile_prsxtra's file_spec). Run from the repo root, where the config's
relative paths resolve.

Usage:
  python scripts/check_raw_sumstats.py config_prsxtra_respiratory.yaml
"""

import gzip
import re
import sys

import yaml


def file_spec(entry, pop):
    spec = dict(entry.get("defaults", {}))
    over = dict(entry["files"][pop])
    spec["cols"] = {**spec.get("cols", {}), **over.pop("cols", {})}
    spec.update(over)
    return spec


def read_header(path, sep):
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as fh:
        line = fh.readline().rstrip("\n\r")
        first = fh.readline()
    cols = line.split(sep) if sep else re.split(r"\s+", line.strip())
    return cols, bool(first)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    config = yaml.safe_load(open(sys.argv[1]))
    pops, n_bad = config["populations"], 0

    for trait in config["traits"]:
        entry = config["sumstats"][trait]
        if entry.get("source") == "pan_ukbb":
            continue
        for pop in pops:
            spec = file_spec(entry, pop)
            path = spec["path"]
            want = {k: v for k, v in spec["cols"].items() if v}
            try:
                header, has_rows = read_header(path, spec.get("sep"))
            except FileNotFoundError:
                print(f"MISSING  {trait:12s} {pop}  {path}")
                n_bad += 1
                continue
            except (OSError, EOFError, UnicodeDecodeError) as e:
                print(f"UNREADABLE {trait:12s} {pop}  {path}: {e}")
                n_bad += 1
                continue
            absent = {k: v for k, v in want.items() if v not in header}
            if absent or not has_rows:
                n_bad += 1
                why = f"columns not in header {absent}" if absent else "no data rows"
                print(f"BAD      {trait:12s} {pop}  {path}: {why}")
                print(f"         header: {header}")
            else:
                print(f"ok       {trait:12s} {pop}  {path}")

    print(f"\n{n_bad} problem(s)")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
