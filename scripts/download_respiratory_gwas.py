#!/usr/bin/env python3
"""
scripts/download_respiratory_gwas.py
------------------------------------
Download the ancestry-specific GWAS summary statistics for the respiratory
trait set (config_prsxtra_respiratory.yaml) into

    data/sumstats/raw/respiratory/{trait}_{POP}.txt.gz

Sources (all public, no login):

  asthma, copd        GBMI (Zhou et al. 2022, Cell Genomics), ancestry-specific
                      inverse-variance meta-analyses, GRCh38
                      s3://gbmi-sumstats  (https://gbmi-sumstats.s3.amazonaws.com)
  fev1, fvc, fev1_fvc EAS: He et al. 2026 (KCPS-II + TWB), GWAS Catalog
                        GCST90705067-69, GRCh37, N = 129,680
                      EUR/AFR/AMR: Shrine et al. 2023 (Nat Genet), GWAS Catalog
                        EUR GCST90292609-11 (all cohorts incl. UK Biobank, z-scores)
                        AFR GCST90292621-23, AMR GCST90292633-35 (UK Biobank)
  smoking, cpd        GSCAN (Saunders et al. 2022, Nature), ancestry-stratified,
                      23andMe removed, GRCh38; https://doi.org/10.13020/przg-dp88
                      Only the SmkInit / CigDay files are pulled out of the
                      per-ancestry zips (HTTP range requests), not whole zips.
  lung_cancer         Ancestry-specific, all public (GWAS Catalog):
                        EUR GCST90475573, AFR GCST90475572, AMR GCST90477192:
                          Million Veteran Program (Verma et al. 2024, Science),
                          PheCode 165.1 "Cancer of bronchus; lung", GRCh38
                        EAS GCST90018655: Biobank Japan (Sakaue et al. 2021,
                          Nat Genet), GRCh37, no rsIDs
                      He et al. instead used Byun et al. 2022 (GCST90134661), a
                      single cross-ancestry meta-analysis, which would count the
                      same participants once per ancestry.

GWAS Catalog files are checked against the catalog's md5sum.txt. Existing
files are skipped; partial downloads (.part) are resumed.

Usage:
  python scripts/download_respiratory_gwas.py                 # everything
  python scripts/download_respiratory_gwas.py --list          # show plan only
  python scripts/download_respiratory_gwas.py --traits fev1,fvc --pops EAS
"""

import argparse
import gzip
import hashlib
import io
import os
import shutil
import sys
import time
import urllib.request
import zipfile

GBMI = "https://gbmi-sumstats.s3.amazonaws.com"
GWASCAT = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
DRUM = "https://conservancy.umn.edu/server/api/core/bitstreams"

# GSCAN per-ancestry zip bitstreams on the UMN Data Repository (DRUM).
GSCAN_ZIP = {
    "EUR": DRUM + "/4e3247b6-f2b8-4d5b-824e-bb88a8e9978a/content",
    "AFR": DRUM + "/37f400b4-be93-40c5-b60e-1fcb0f699471/content",
    "EAS": DRUM + "/8eea02d5-98d8-4232-b6b5-ffe457e42f27/content",
    "AMR": DRUM + "/29ed235c-eafc-4161-b7a6-28f8ca4041ef/content",
}
GSCAN_PHENO = {"smoking": "SmkInit", "cpd": "CigDay"}

LUNG_FUNCTION = {   # trait -> {POP: GWAS Catalog accession}
    "fev1":     {"EAS": "GCST90705067", "EUR": "GCST90292609", "AFR": "GCST90292621", "AMR": "GCST90292633"},
    "fvc":      {"EAS": "GCST90705068", "EUR": "GCST90292610", "AFR": "GCST90292622", "AMR": "GCST90292634"},
    "fev1_fvc": {"EAS": "GCST90705069", "EUR": "GCST90292611", "AFR": "GCST90292623", "AMR": "GCST90292635"},
}
LUNG_CANCER = {"EUR": "GCST90475573", "AFR": "GCST90475572", "AMR": "GCST90477192",
               "EAS": "GCST90018655"}

# File name inside each accession folder when it is not {acc}.tsv.gz
# (He et al.'s EAS files are uncompressed .tsv).
GWASCAT_FILE = {"GCST90705067": "GCST90705067.tsv", "GCST90705068": "GCST90705068.tsv",
                "GCST90705069": "GCST90705069.tsv",
                "GCST90018655": "GCST90018655_buildGRCh37.tsv.gz"}

POPS = ["EUR", "EAS", "AFR", "AMR"]


def gwascat_dir(acc):
    n = int(acc[4:])
    lo = (n - 1) // 1000 * 1000 + 1
    return "%s/GCST%08d-GCST%08d/%s" % (GWASCAT, lo, lo + 999, acc)


def plan():
    """Every (trait, pop) -> a download job dict."""
    jobs = []
    for trait, name in (("asthma", "Asthma"), ("copd", "COPD")):
        for pop in POPS:
            jobs.append(dict(trait=trait, pop=pop, kind="http", source="GBMI",
                             url="%s/%s_Bothsex_%s_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
                                 % (GBMI, name, pop.lower()),
                             gzipped=True))
    for trait, accs in list(LUNG_FUNCTION.items()) + [("lung_cancer", LUNG_CANCER)]:
        for pop, acc in accs.items():
            fname = GWASCAT_FILE.get(acc, acc + ".tsv.gz")
            d = gwascat_dir(acc)
            jobs.append(dict(trait=trait, pop=pop, kind="http", source=acc,
                             url="%s/%s" % (d, fname), md5_url=d + "/md5sum.txt",
                             md5_name=fname, gzipped=fname.endswith(".gz")))
    for trait, pheno in GSCAN_PHENO.items():
        for pop in POPS:
            jobs.append(dict(trait=trait, pop=pop, kind="zip_member", source="GSCAN",
                             url=GSCAN_ZIP[pop],
                             member="GSCAN_%s_2022_GWAS_SUMMARY_STATS_%s.txt.gz" % (pheno, pop)))
    return jobs


# ── HTTP helpers ─────────────────────────────────────────────────────────────
def urlopen(url, headers=None, tries=5):
    for k in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}),
                                          timeout=120)
        except Exception as e:                       # noqa: BLE001 - retry anything transient
            if k == tries - 1:
                raise
            wait = 5 * 2 ** k
            print("    %s; retrying in %ds" % (e, wait), flush=True)
            time.sleep(wait)


def remote_size(url):
    r = urlopen(url, {"Range": "bytes=0-0"})
    cr = r.headers.get("Content-Range", "")
    r.read()
    if "/" in cr:
        return int(cr.rsplit("/", 1)[1])
    raise RuntimeError("server did not report a size for %s" % url)


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only view of a remote file via HTTP range requests."""

    def __init__(self, url):
        self.url, self.pos, self.size = url, 0, remote_size(url)

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        data = urlopen(self.url, {"Range": "bytes=%d-%d" % (self.pos, self.pos + n - 1)}).read()
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def fetch_md5(md5_url, name):
    for line in urlopen(md5_url).read().decode().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == name:
            return parts[0]
    return None


# ── download kinds ───────────────────────────────────────────────────────────
def download_http(job, dest):
    """Stream to dest (gzip-compressing plain-text sources), resuming a .part
    file when the source is already gzipped; verify md5 when published."""
    part = dest + ".part"
    expected = fetch_md5(job["md5_url"], job["md5_name"]) if job.get("md5_url") else None
    if job["gzipped"]:
        total = remote_size(job["url"])
        for attempt in range(20):                 # resume after mid-transfer stalls
            have = os.path.getsize(part) if os.path.exists(part) else 0
            if have >= total:
                break
            headers = {"Range": "bytes=%d-" % have} if have else {}
            try:
                r = urlopen(job["url"], headers)
                if have and r.status != 206:      # server ignored the range: restart
                    have = 0
                with open(part, "ab" if have else "wb") as out:
                    shutil.copyfileobj(r, out, 1 << 22)
            except Exception as e:                # noqa: BLE001 - resume on any stall
                print("    stalled at %.1f MB (%s); resuming"
                      % (os.path.getsize(part) / 1e6 if os.path.exists(part) else 0, e), flush=True)
                time.sleep(min(60, 5 * (attempt + 1)))
        if os.path.getsize(part) != total:
            raise RuntimeError("incomplete download (%d of %d bytes); rerun to resume"
                               % (os.path.getsize(part), total))
        if expected:
            h = hashlib.md5()
            with open(part, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(chunk)
            check(h.hexdigest(), expected, part)
    else:
        h = hashlib.md5()
        r = urlopen(job["url"])
        with gzip.open(part, "wb", compresslevel=4) as out:
            for chunk in iter(lambda: r.read(1 << 22), b""):
                h.update(chunk)
                out.write(chunk)
        if expected:
            check(h.hexdigest(), expected, part)
    os.replace(part, dest)


def check(got, expected, path):
    if got != expected:
        os.remove(path)
        raise RuntimeError("md5 mismatch for %s (got %s, expected %s); removed, rerun to retry"
                           % (path, got, expected))
    print("    md5 ok", flush=True)


def download_zip_member(job, dest):
    part = dest + ".part"
    raw = io.BufferedReader(HTTPRangeFile(job["url"]), buffer_size=8 << 20)
    with zipfile.ZipFile(raw) as z:
        info = z.getinfo(job["member"])
        with z.open(info) as src, open(part, "wb") as out:   # CRC-checked by zipfile
            shutil.copyfileobj(src, out, 8 << 20)
    os.replace(part, dest)


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out_dir", default="data/sumstats/raw/respiratory")
    ap.add_argument("--traits", help="comma-separated subset (default: all)")
    ap.add_argument("--pops", help="comma-separated subset of EUR,EAS,AFR,AMR (default: all)")
    ap.add_argument("--list", action="store_true", help="print the plan and exit")
    args = ap.parse_args()

    jobs = plan()
    if args.traits:
        keep = set(args.traits.split(","))
        jobs = [j for j in jobs if j["trait"] in keep]
    if args.pops:
        keep = set(args.pops.split(","))
        jobs = [j for j in jobs if j["pop"] in keep]
    os.makedirs(args.out_dir, exist_ok=True)

    failed = []
    for j in jobs:
        dest = os.path.join(args.out_dir, "%s_%s.txt.gz" % (j["trait"], j["pop"]))
        where = j.get("member") or j.get("url", "")
        if args.list:
            print("%-12s %-4s %-14s %s" % (j["trait"], j["pop"], j["source"], where))
            continue
        if os.path.exists(dest):
            print("[skip] %s (exists)" % dest)
            continue
        print("[get]  %s  <-  %s" % (dest, where), flush=True)
        try:
            if j["kind"] == "http":
                download_http(j, dest)
            else:
                download_zip_member(j, dest)
            print("    %.1f MB" % (os.path.getsize(dest) / 1e6), flush=True)
        except Exception as e:                       # noqa: BLE001 - report and continue
            print("    FAILED: %s" % e, flush=True)
            failed.append(dest)

    if failed:
        print("\nFailed (rerun to retry/resume):")
        for f in failed:
            print("  " + f)
        sys.exit(1)


if __name__ == "__main__":
    main()
