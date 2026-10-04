#!/usr/bin/env python3
"""
aou/03_extract_phenotypes.py
----------------------------
Inside the All of Us Researcher Workbench: build the phenotype/covariate table
for the evaluation (04_stack_evaluate.R), one row per participant with
genetic ancestry predictions:

  person_id ancestry age sex PC1..PC16 statin <outcome columns>

  lipids       LDL, HDL, TG: per-person median of EHR lab values (LOINC), in
               mg/dL, within plausible ranges; `statin` = any statin exposure
               on record. Age at the person's most recent lipid measurement.
  respiratory  asthma, copd, lung_cancer (+ asthma_exacerbation,
               copd_exacerbation): 1 = at least --min_code_dates distinct dates
               with a qualifying ICD-9/10-CM code; 0 = has EHR data and no
               code. Age at --ref_date.

Disease codes are He et al.'s (Supplementary Tables 36-40), stored as OMOP
source concept IDs in aou/he2026_codes.tsv. The true ICD codes of those
concepts are looked up in the CDR, and a condition counts if its ICD code
equals one of them or is a child code (prefix match), as selecting a parent
code in the All of Us cohort builder does (e.g. C34 -> C34.90). Not yet
reproduced: He et al.'s SERPINA1 exclusions for COPD.

Ancestry and PCs come from the All of Us genetic ancestry file
(ancestry_preds.tsv: research_id, ancestry_pred, pca_features) of your CDR
version, under .../wgs/short_read/snpindel/aux/ancestry/.

Usage (in a Workbench notebook terminal):
  python aou/03_extract_phenotypes.py --trait_set lipids \\
      --ancestry_preds gs://.../aux/ancestry/ancestry_preds.tsv --out pheno_lipids.tsv
"""

import argparse
import ast
import os

import numpy as np
import pandas as pd

LIPID_LOINC = {
    "LDL": ["13457-7", "18262-6", "2089-1"],   # calculated, direct, LDL mass/vol
    "HDL": ["2085-9"],
    "TG":  ["2571-8"],
}
LIPID_RANGE_MG_DL = {"LDL": (10, 400), "HDL": (5, 200), "TG": (10, 3000)}
STATINS = ["atorvastatin", "fluvastatin", "lovastatin", "pitavastatin",
           "pravastatin", "rosuvastatin", "simvastatin"]

DISEASES = ["asthma", "copd", "lung_cancer", "asthma_exacerbation", "copd_exacerbation"]
CODES_TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "he2026_codes.tsv")


def bq(sql):
    from google.cloud import bigquery
    return bigquery.Client().query(sql).to_dataframe()


def demographics(cdr):
    df = bq(f"""
        SELECT p.person_id, p.birth_datetime, c.concept_name AS sex_at_birth
        FROM `{cdr}.person` p
        LEFT JOIN `{cdr}.concept` c ON p.sex_at_birth_concept_id = c.concept_id
    """)
    df["sex"] = df["sex_at_birth"].map({"Male": 1, "Female": 0})
    df["birth_datetime"] = pd.to_datetime(df["birth_datetime"], utc=True).dt.tz_localize(None)
    return df[["person_id", "birth_datetime", "sex"]]


def ancestry(path, n_pcs):
    a = pd.read_csv(path, sep="\t")
    pcs = np.array([ast.literal_eval(x)[:n_pcs] for x in a["pca_features"]])
    out = pd.DataFrame({"person_id": a["research_id"].astype(np.int64),
                        "ancestry": a["ancestry_pred"].str.upper()})
    for k in range(n_pcs):
        out["PC%d" % (k + 1)] = pcs[:, k]
    return out


def lipids(cdr):
    codes = ",".join("'%s'" % c for v in LIPID_LOINC.values() for c in v)
    m = bq(f"""
        SELECT m.person_id, c.concept_code AS loinc, m.measurement_date, m.value_as_number
        FROM `{cdr}.measurement` m
        JOIN `{cdr}.concept` c ON m.measurement_concept_id = c.concept_id
        WHERE c.vocabulary_id = 'LOINC' AND c.concept_code IN ({codes})
          AND m.value_as_number IS NOT NULL
    """)
    loinc2trait = {c: t for t, cs in LIPID_LOINC.items() for c in cs}
    m["trait"] = m["loinc"].map(loinc2trait)
    lo = m["trait"].map(lambda t: LIPID_RANGE_MG_DL[t][0])
    hi = m["trait"].map(lambda t: LIPID_RANGE_MG_DL[t][1])
    m = m[(m["value_as_number"] >= lo) & (m["value_as_number"] <= hi)]
    m["measurement_date"] = pd.to_datetime(m["measurement_date"])

    vals = m.pivot_table(index="person_id", columns="trait", values="value_as_number",
                         aggfunc="median")
    last = m.groupby("person_id")["measurement_date"].max().rename("age_date")
    statin_list = ",".join("'%s'" % s for s in STATINS)
    st = bq(f"""
        SELECT DISTINCT de.person_id
        FROM `{cdr}.drug_exposure` de
        JOIN `{cdr}.concept_ancestor` ca ON de.drug_concept_id = ca.descendant_concept_id
        JOIN `{cdr}.concept` c ON ca.ancestor_concept_id = c.concept_id
        WHERE c.vocabulary_id = 'RxNorm' AND c.concept_class_id = 'Ingredient'
          AND LOWER(c.concept_name) IN ({statin_list})
    """)
    out = vals.join(last).reset_index()
    out["statin"] = out["person_id"].isin(set(st["person_id"])).astype(int)
    return out


def disease_prefixes(cdr, codes_tsv):
    """{outcome: [(vocabulary_id, icd_code), ...]} for He et al.'s concept IDs."""
    codes = pd.read_csv(codes_tsv, sep="\t", comment="#")
    ids = ",".join(str(i) for i in codes["concept_id"].unique())
    look = bq(f"""
        SELECT concept_id, vocabulary_id, concept_code FROM `{cdr}.concept`
        WHERE concept_id IN ({ids})
    """)
    missing = set(codes["concept_id"]) - set(look["concept_id"])
    if missing:
        print("WARNING: %d He et al. concept IDs not in this CDR: %s" % (len(missing), sorted(missing)))
    codes = codes.merge(look, on="concept_id")
    return {o: list(zip(g["vocabulary_id"], g["concept_code"])) for o, g in codes.groupby("outcome")}


def diseases(cdr, min_dates, ref_date, codes_tsv):
    ehr = bq(f"SELECT person_id FROM `{cdr}.cb_search_person` WHERE has_ehr_data = 1")
    out = ehr.copy()
    prefixes = disease_prefixes(cdr, codes_tsv)
    for name in DISEASES:
        cond = " OR ".join("(c.vocabulary_id = '%s' AND STARTS_WITH(c.concept_code, '%s'))"
                           % (v, code) for v, code in prefixes[name])
        cases = bq(f"""
            SELECT co.person_id, COUNT(DISTINCT co.condition_start_date) AS n_dates
            FROM `{cdr}.condition_occurrence` co
            JOIN `{cdr}.concept` c ON co.condition_source_concept_id = c.concept_id
            WHERE {cond}
            GROUP BY co.person_id
        """)
        case_ids = set(cases.loc[cases["n_dates"] >= min_dates, "person_id"])
        any_code = set(cases["person_id"])
        out[name] = np.where(out["person_id"].isin(case_ids), 1,
                             np.where(out["person_id"].isin(any_code), np.nan, 0))
    out["age_date"] = pd.Timestamp(ref_date)
    out["statin"] = 0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait_set", choices=["lipids", "respiratory"], required=True)
    ap.add_argument("--ancestry_preds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cdr", default=os.environ.get("WORKSPACE_CDR"),
                    help="BigQuery dataset (default $WORKSPACE_CDR)")
    ap.add_argument("--ref_date", default="2023-07-01",
                    help="age reference date for disease outcomes (set to your CDR cutoff)")
    ap.add_argument("--min_code_dates", type=int, default=1,
                    help="distinct code dates needed to call a case (1 or 2 are common)")
    ap.add_argument("--n_pcs", type=int, default=16)
    ap.add_argument("--codes", default=CODES_TSV, help="disease concept IDs (default: He et al.)")
    args = ap.parse_args()
    if not args.cdr:
        raise SystemExit("No CDR dataset: run in the Workbench or pass --cdr")

    demo = demographics(args.cdr)
    anc = ancestry(args.ancestry_preds, args.n_pcs)
    ph = (lipids(args.cdr) if args.trait_set == "lipids"
          else diseases(args.cdr, args.min_code_dates, args.ref_date, args.codes))

    df = anc.merge(demo, on="person_id").merge(ph, on="person_id")
    df["age"] = (df["age_date"] - df["birth_datetime"]).dt.days / 365.25
    df = df.drop(columns=["birth_datetime", "age_date"])
    df = df[df["sex"].notna() & df["age"].between(18, 110)]
    df.to_csv(args.out, sep="\t", index=False)

    outcome_cols = list(LIPID_LOINC) if args.trait_set == "lipids" else DISEASES
    print("%d participants with genetic ancestry and phenotypes -> %s" % (len(df), args.out))
    for o in outcome_cols:
        sub = df[df[o].notna()]
        if args.trait_set == "lipids":
            print("  %-12s n=%d" % (o, len(sub)))
        else:
            print("  %-12s n=%d cases=%d" % (o, len(sub), int(sub[o].sum())))


if __name__ == "__main__":
    main()
