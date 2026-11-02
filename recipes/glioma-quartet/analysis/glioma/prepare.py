"""Extract eligible patient-level measurements from the versioned public sources."""
import argparse
import csv
import gzip
import hashlib
import io
from itertools import takewhile
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
import openpyxl
import pandas as pd
import pdfplumber


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "frozen_panels.tsv"
HASHES = {
    "gse.tsv.gz": "47a6408994e562739286e225ea563aedcd7aaf6992d2d6ec22b40e7ab0eb3a70",
    "gpl.soft.gz": "cf44a7e96244581d7d906858c1a7fd5ced13d1c4cc0497d23477a6f708f6933b",
    "22383048_00085472can092307-sup-stabs_1-6.pdf": "c397b381f0f1573edda4def2b9d691465b09eb3d3703860dcc863f66bbd11d25",
    "huashan.zip": "e906120bb883a865c8c42a60406a9a271b0089d913f8f32cac6823271a45dbd1",
    "huashan_clin.xlsx": "e7ccc5ed9591df13c97e4e2123ee9c2bd30d1ac4b0a9973167323c03a1043771",
}
SUP = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-54352-5/MediaObjects/"
URLS = {
    "gse.tsv.gz":"https://ftp.ncbi.nih.gov/geo/series/GSE16nnn/GSE16011/matrix/GSE16011_series_matrix.txt.gz",
    "gpl.soft.gz":"https://ftp.ncbi.nih.gov/geo/platforms/GPL8nnn/GPL8542/soft/GPL8542_family.soft.gz",
    "22383048_00085472can092307-sup-stabs_1-6.pdf":"https://ndownloader.figshare.com/files/39828447",
    "huashan.zip":SUP+"41467_2024_54352_MOESM11_ESM.zip",
    "huashan_clin.xlsx":SUP+"41467_2024_54352_MOESM2_ESM.xlsx",
}


def erasmus(d, mapping):
    with pdfplumber.open(next(d.glob('*stabs_1-6.pdf'))) as doc:
        repeats = {row[3] for row in doc.pages[18].extract_table()[2:]}
        assert len(repeats) == 7
        ids = []
        for i in range(19, 25):
            table = doc.pages[i].extract_table()
            if i == 19:
                assert table[0][3] == "IDH1 (R132) mutation"
                table = table[1:]
            ids.extend((r[0], r[3]) for r in table)
    calls = dict(ids)
    assert len(ids) == len(calls) == 284 and set(calls.values()) == {"", "no mutation", "mutation"}
    with gzip.open(d / "gse.tsv.gz", "rt") as handle:
        chars = []
        for line in handle:
            if line.startswith("!Sample_"):
                key, *values = next(csv.reader([line], delimiter="\t"))
                if key == "!Sample_title":
                    titles = [v.rsplit(" ", 1)[1] for v in values]
                elif key == "!Sample_accession" or key == "!Sample_geo_accession":
                    gsm = values
                elif key == "!Sample_characteristics_ch1":
                    chars.append(values)
            if line.startswith("!series_matrix_table_begin"):
                break
        mat = pd.read_csv(handle, sep="\t", index_col=0, comment="!").astype(float)
    assert list(mat) == gsm and len(set(titles)) == 284 and set(titles) == set(calls)
    rows = []
    for i, (title, sample) in enumerate(zip(titles, gsm)):
        meta = dict(v.split(": ", 1) for v in (x[i] for x in chars) if ": " in v)
        age = pd.to_numeric(meta.get("age at diagnosis"), errors="coerce")
        keep = (meta.get("histology") == "GBM (grade IV)" and age >= 18 and
                title not in repeats and calls[title] in ("mutation", "no mutation"))
        rows.append((sample, title, calls[title], age, meta.get("histology"), title in repeats, keep))
    audit = pd.DataFrame(rows, columns=["sample","source_id","idh_call","age","histology","repeat","eligible"])
    keep = audit[audit.eligible].set_index("sample")
    # GPL8542's ORF column is the Entrez identifier in this custom CDF.
    with gzip.open(d / "gpl.soft.gz", "rt", errors="replace") as stream:
        for ln in stream:
            if ln.startswith("!platform_table_begin"):
                break
        lines = takewhile(lambda ln: not ln.startswith("!platform_table_end"), stream)
        annot = pd.read_csv(io.StringIO("".join(lines)), sep="\t", usecols=["ID","ORF"], dtype=str)
    annot["ID"] = annot.ID.str.strip()
    genes = annot.set_index("ID", verify_integrity=True).ORF.map(mapping.set_index("entrez_id").gene)
    measured = mat.loc[:,keep.index].groupby(genes.reindex(mat.index)).median().T.reindex(columns=mapping.gene)
    measured.insert(0, "mutant", keep.idh_call.eq("mutation").astype(int))
    return measured, audit


def sheet(bundle, book, name):
    with io.BytesIO(bundle.read("Source data/" + book + ".xlsx")) as buf:
        wb = openpyxl.load_workbook(buf, read_only=True, data_only=True)
        try:
            yield from wb[name].values
        finally:
            wb.close()


def huashan(d, mapping, scale="source"):
    clin = pd.read_excel(d / "huashan_clin.xlsx").set_index("EP_tube_ID", verify_integrity=True)
    eligible = clin[clin.WHO_2021.isin(["Glioblastoma,4", "Astrocytoma,4"]) &
                    (pd.to_numeric(clin.Age,errors="coerce") >= 18) & clin.IDH.isin(["+", "-"])]
    acc = {u:g for g,us in zip(mapping.gene,mapping.uniprot_ids.fillna("")) for u in us.split("|") if u}
    with ZipFile(d / "huashan.zip") as z:
        meta = iter(sheet(z,"Figure2","Figure2A annotation_TN"))
        assert list(next(meta))[2:] == ["EP_tube","Type"]
        tube = {key:(donor,kind) for _,key,donor,kind in meta}
        rows = iter(sheet(z,"FigureS2","FigureS2E"))
        header = next(rows)
        assert header[:3] == ("sampleID","technical_replicates","batch")
        cols = {}
        for i, group in enumerate(header[3:], 3):
            names = {acc.get(u) for u in group.split(";")}
            if len(names) == 1 and None not in names:
                cols[i] = names.pop()
        keep, audit = [], []
        for row in rows:
            sample, rep, batch = row[:3]
            key = f"b{int(batch[1:]):02d}_{int(sample):03d}"
            donor, kind = tube.get(key, ("", "unmatched"))
            include = kind == "Tumor" and donor in eligible.index
            audit.append((key, rep, donor, kind, batch, include))
            if include:
                values = np.array(row[3:], dtype=float)
                if scale == "log2-center":
                    values[values <= 0] = np.nan
                    values = np.log2(values)
                if scale != "source":
                    values -= np.nanmedian(values)
                keep.append([key, donor, *(values[i-3] for i in cols)])
    x = pd.DataFrame(keep, columns=["sample","donor",*cols.values()]).set_index(["donor","sample"])
    # Collapse protein groups, then technical repeats of a specimen, then tumor
    # specimens per patient so that extra runs cannot inflate evidence.
    x = x.T.groupby(level=0).median().T.groupby(level=[0,1]).median().groupby(level=0).median()
    x = x.reindex(columns=mapping.gene)
    x.insert(0, "mutant", eligible.loc[x.index].IDH.eq("+").astype(int))
    return x, pd.DataFrame(audit, columns=["sample","technical_repeat","donor","type","batch","eligible"])


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--inputs", required=True, type=Path)
    cli.add_argument("--out", required=True, type=Path)
    cli.add_argument("--fetch", action="store_true", help="Download missing public source files and verify hashes")
    cli.add_argument("--protein-scale",choices=["source","center","log2-center"],default="source")
    a = cli.parse_args()
    a.inputs.mkdir(parents=True,exist_ok=True)
    for name, digest in HASHES.items():
        path = a.inputs/name
        if a.fetch and not path.exists():
            with tempfile.TemporaryDirectory(dir=a.inputs) as tmp:
                staged = Path(tmp)/name
                with urlopen(URLS[name],timeout=180) as src, staged.open("wb") as dst:
                    shutil.copyfileobj(src,dst)
                with staged.open("rb") as f:
                    if hashlib.file_digest(f,"sha256").hexdigest() != digest:
                        raise ValueError(f"Download differs from reviewed {name}")
                staged.replace(path)
        with path.open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() != digest:
                raise ValueError(f"Source changed: {name}")
    mapping = pd.read_csv(ROOT / "gene_ids.tsv", sep="\t", dtype=str)
    if not mapping.gene.is_unique or set(mapping.gene) != set(pd.read_csv(PANEL,sep="\t").gene):
        raise ValueError("Gene map differs from frozen panels")
    a.out.mkdir(parents=True, exist_ok=False)
    for name in ("RNA","protein"):
        x, audit = erasmus(a.inputs,mapping) if name == "RNA" else huashan(a.inputs,mapping,a.protein_scale)
        counts = x.mutant.value_counts()
        if not x.index.is_unique or set(counts.index) != {0,1} or counts.min() < 5:
            raise ValueError("Invalid patient-level contrast")
        x.to_csv(a.out / (name + ".tsv"), sep="\t", index_label="donor")
        audit.to_csv(a.out / (name + "_eligibility.tsv"),sep="\t",index=False)
        print(name, counts.to_dict(), "mapped", int(x.drop(columns="mutant").notna().any().sum()), flush=True)


if __name__ == "__main__":
    main()
