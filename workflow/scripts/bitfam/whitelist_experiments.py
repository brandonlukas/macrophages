"""Per-experiment provenance for the whitelisted ChIP-Atlas evidence.

For each ChIP-seq experiment that survived the mono/mac/DC whitelist filter and
is present in the BITFAM input network, resolve:

    factor, srx, gsm, gse, cell_type, sample_title, n_targets,
    study_title, pmid, n_pmids, all_pmids,
    doi, paper_title, authors, journal, year

Each row is a (factor, SRX, cell_type) triple. Publication fields (pmid, doi,
paper_title, authors, journal, year) describe the *chosen* paper for the GSE:
GEO's first linked PMID when present, otherwise a Europe PMC search hit on the
GSE accession. NCBI does not document the order of GEO's PubMedIds list, so
"first" is not principled — n_pmids and all_pmids preserve the full set so
multi-publication GSEs can be audited.

Columns:
    n_pmids   - count of PubMed IDs GEO links to this GSE (0 when the EPMC
                fallback supplied the publication, or no paper was found)
    all_pmids - ";"-separated list of all GEO-linked PMIDs

Hops:
    SRX  -> GSM            (NCBI SRA efetch runinfo)
    GSM  -> GSE UID        (NCBI esearch gds)
    GSE  -> title, samples,
            pubmed_ids     (NCBI esummary gds)
    GSE  -> pmid, doi, year (Europe PMC fallback when GEO has no PubMedIds)
    PMID -> doi, year      (NCBI esummary pubmed)

All NCBI/EuropePMC results are cached as JSON under params.cache_dir so reruns
are offline joins.
"""

import ast
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from tqdm import tqdm

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = "macrophages-pipeline/0.1"

# NCBI: 3 req/s anonymous, 10 req/s with API key. EPMC tolerates ~10 req/s.
_HAS_API_KEY = bool(os.environ.get("NCBI_API_KEY", "").strip())
THROTTLE_S = 0.11 if _HAS_API_KEY else 0.34


def _api_key() -> str:
    k = os.environ.get("NCBI_API_KEY", "").strip()
    return f"&api_key={k}" if k else ""


def _http(url: str, retries: int = 4, base_delay: float = 0.5) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    time.sleep(THROTTLE_S)
    for i in range(retries - 1):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode()
        except Exception:
            time.sleep(base_delay * (2**i))
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode()


def _chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i : i + n]


def _cached(path: Path, compute):
    if path.exists():
        return json.loads(path.read_text())
    val = compute()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(val, indent=2))
    return val


def srx_to_gsm(srxs: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    batches = list(_chunks(srxs, 100))
    for batch in tqdm(batches, desc="SRX→GSM"):
        term = urllib.parse.quote(" OR ".join(batch))
        ids_xml = _http(
            f"{EUTILS}/esearch.fcgi?db=sra&term={term}&retmax=10000{_api_key()}"
        )
        ids = re.findall(r"<Id>(\d+)</Id>", ids_xml)
        if not ids:
            continue
        runinfo = _http(
            f"{EUTILS}/efetch.fcgi?db=sra&id={','.join(ids)}"
            f"&rettype=runinfo&retmode=xml{_api_key()}"
        )
        for row in re.findall(r"<Row>.*?</Row>", runinfo, re.DOTALL):
            m_srx = re.search(r"<Experiment>(SRX\d+)</Experiment>", row)
            m_gsm = re.search(r"<SampleName>(GSM\d+)</SampleName>", row)
            if m_srx:
                out[m_srx.group(1)] = m_gsm.group(1) if m_gsm else ""
    for s in srxs:
        out.setdefault(s, "")
    return out


def gsm_to_gse_uid(gsms: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for gsm in tqdm([g for g in gsms if g], desc="GSM→GSE"):
        xml = _http(
            f"{EUTILS}/esearch.fcgi?db=gds"
            f"&term={gsm}[Accession]+AND+gse[ETYP]&retmax=5{_api_key()}"
        )
        ids = re.findall(r"<Id>(\d+)</Id>", xml)
        if ids:
            out[gsm] = ids[0]
    return out


def gse_meta(gse_uids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    uniq = sorted({u for u in gse_uids if u})
    batches = list(_chunks(uniq, 200))
    for batch in tqdm(batches, desc="GSE meta"):
        text = _http(
            f"{EUTILS}/esummary.fcgi?db=gds&id={','.join(batch)}"
            f"&retmode=json{_api_key()}"
        )
        result = json.loads(text).get("result", {})
        for uid in result.get("uids", []):
            rec = result[uid]
            samples = {
                s.get("accession", ""): s.get("title", "")
                for s in rec.get("samples", [])
            }
            out[uid] = {
                "gse": rec.get("accession", ""),
                "title": rec.get("title", ""),
                "pmids": [str(p) for p in rec.get("pubmedids", [])],
                "samples": samples,
            }
    return out


def europepmc_for_gse(gse_acc: str) -> dict | None:
    if not gse_acc:
        return None
    text = _http(
        f"{EPMC}?query={urllib.parse.quote(gse_acc)}"
        "&format=json&resultType=lite&pageSize=5"
    )
    hits = json.loads(text).get("resultList", {}).get("result", [])
    if not hits:
        return None
    hits.sort(key=lambda h: h.get("pubYear", "9999"))
    h = hits[0]
    return {
        "pmid": h.get("pmid", ""),
        "doi": h.get("doi", ""),
        "year": h.get("pubYear", ""),
        "paper_title": h.get("title", ""),
        "authors": h.get("authorString", ""),
        "journal": h.get("journalTitle", ""),
    }


def pubmed_meta(pmids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    uniq = sorted({p for p in pmids if p})
    batches = list(_chunks(uniq, 200))
    for batch in tqdm(batches, desc="PMID meta"):
        text = _http(
            f"{EUTILS}/esummary.fcgi?db=pubmed&id={','.join(batch)}"
            f"&retmode=json{_api_key()}"
        )
        result = json.loads(text).get("result", {})
        for pmid in result.get("uids", []):
            rec = result[pmid]
            doi = ""
            for aid in rec.get("articleids", []):
                if aid.get("idtype", "").lower() == "doi":
                    doi = aid.get("value", "")
                    break
            pub = rec.get("pubdate", "")
            year = pub.split(" ")[0] if pub else ""
            authors = "; ".join(
                a.get("name", "") for a in rec.get("authors", []) if a.get("name")
            )
            out[pmid] = {
                "doi": doi,
                "year": year,
                "paper_title": rec.get("title", ""),
                "authors": authors,
                "journal": rec.get("fulljournalname", "") or rec.get("source", ""),
            }
    return out


PUB_KEYS = ("doi", "paper_title", "authors", "journal", "year")


def whitelist_experiments(input_csv: str, output_tsv: str, cache_dir: str) -> None:
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv, usecols=["source", "target", "evidence"])
    df = df[df["evidence"].str.startswith("[")]
    df["evidence"] = df["evidence"].apply(ast.literal_eval)
    exp = df.explode("evidence")
    exp[["srx", "cell_type"]] = exp["evidence"].str.split("|", n=1, expand=True)
    exp = exp.dropna(subset=["srx"])

    grouped = (
        exp.groupby(["source", "srx", "cell_type"], dropna=False)
        .agg(n_targets=("target", "nunique"))
        .reset_index()
    )

    srxs = sorted(grouped["srx"].unique())
    srx_gsm = _cached(cache / "srx_to_gsm.json", lambda: srx_to_gsm(srxs))

    gsms = sorted({g for g in srx_gsm.values() if g})
    gsm_uid = _cached(cache / "gsm_to_gse_uid.json", lambda: gsm_to_gse_uid(gsms))

    uids = sorted(set(gsm_uid.values()))
    gse = _cached(cache / "gse_meta.json", lambda: gse_meta(uids))

    gse_by_acc: dict[str, dict] = {m["gse"]: m for m in gse.values()}
    gsm_to_acc: dict[str, str] = {
        gsm: acc for acc, m in gse_by_acc.items() for gsm in m["samples"]
    }

    # Resolve each GSE to a single publication record. Prefer GEO-linked PMID,
    # fall back to Europe PMC search on the GSE accession.
    epmc_path = cache / "europepmc_by_gse.json"
    epmc_results: dict[str, dict | None] = (
        json.loads(epmc_path.read_text()) if epmc_path.exists() else {}
    )
    epmc_dirty = False
    gse_pub: dict[str, dict] = {}
    for acc, meta in gse_by_acc.items():
        record = {
            "pmid": "",
            "n_pmids": len(meta["pmids"]),
            "all_pmids": ";".join(meta["pmids"]),
            **{k: "" for k in PUB_KEYS},
        }
        if meta["pmids"]:
            record["pmid"] = meta["pmids"][0]
        else:
            if acc not in epmc_results:
                epmc_results[acc] = europepmc_for_gse(acc)
                epmc_dirty = True
            hit = epmc_results.get(acc) or {}
            record["pmid"] = hit.get("pmid", "")
            for k in PUB_KEYS:
                record[k] = hit.get(k, "")
        gse_pub[acc] = record
    if epmc_dirty:
        epmc_path.write_text(json.dumps(epmc_results, indent=2))

    pmids = sorted({r["pmid"] for r in gse_pub.values() if r["pmid"]})
    pubmed = _cached(cache / "pubmed_meta.json", lambda: pubmed_meta(pmids))

    # Overlay PubMed esummary (authoritative) onto whatever EPMC pre-filled.
    for record in gse_pub.values():
        pm = pubmed.get(record["pmid"], {})
        for k in PUB_KEYS:
            if pm.get(k):
                record[k] = pm[k]

    rows = []
    for r in grouped.itertuples(index=False):
        gsm = srx_gsm.get(r.srx, "")
        gse_acc = gsm_to_acc.get(gsm, "")
        meta = gse_by_acc.get(gse_acc, {})
        pub = gse_pub.get(gse_acc, {})
        rows.append(
            {
                "factor": r.source,
                "srx": r.srx,
                "gsm": gsm,
                "gse": gse_acc,
                "cell_type": r.cell_type,
                "sample_title": meta.get("samples", {}).get(gsm, ""),
                "n_targets": r.n_targets,
                "study_title": meta.get("title", ""),
                "pmid": pub.get("pmid", ""),
                "n_pmids": pub.get("n_pmids", 0),
                "all_pmids": pub.get("all_pmids", ""),
                "doi": pub.get("doi", ""),
                "paper_title": pub.get("paper_title", ""),
                "authors": pub.get("authors", ""),
                "journal": pub.get("journal", ""),
                "year": pub.get("year", ""),
            }
        )

    out = (
        pd.DataFrame(rows)
        .sort_values(["factor", "srx"])
        .reset_index(drop=True)
    )
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_tsv, sep="\t", index=False)


snakemake = snakemake  # type: ignore
whitelist_experiments(
    input_csv=snakemake.input[0],
    output_tsv=snakemake.output[0],
    cache_dir=snakemake.params.cache_dir,
)
