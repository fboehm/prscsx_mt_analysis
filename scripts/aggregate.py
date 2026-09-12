"""Snakemake script: aggregate per-run JSON evaluations into one CSV."""

import csv
import json
import os

all_jsons = list(snakemake.input.prscsx) + list(snakemake.input.prscsx_mt)

rows = []
for jf in all_jsons:
    with open(jf) as ff:
        data = json.load(ff)

    method = data["method"]
    n_eur = data["n_eur"]
    n_eas = data["n_eas"]
    rg = data["rg"]
    seed = data["seed"]
    trait = data["trait"]

    for pop_label, metrics in data["metrics"].items():
        row = {
            "method": method,
            "n_eur": n_eur,
            "n_eas": n_eas,
            "rg": rg,
            "seed": seed,
            "trait": trait,
            "pop": pop_label,
            "corr": metrics.get("corr"),
            "mse": metrics.get("mse"),
            "r2_effects": metrics.get("r2_effects"),
        }
        rows.append(row)

# Sort for reproducible output
rows.sort(key=lambda r: (r["method"], r["n_eur"], r["n_eas"], r["rg"], r["seed"], r["trait"], r["pop"]))

os.makedirs(os.path.dirname(snakemake.output.csv), exist_ok=True)
fieldnames = ["method", "n_eur", "n_eas", "rg", "seed", "trait", "pop", "corr", "mse", "r2_effects"]
with open(snakemake.output.csv, "w", newline="") as ff:
    writer = csv.DictWriter(ff, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("... aggregated %d rows into %s ..." % (len(rows), snakemake.output.csv))
