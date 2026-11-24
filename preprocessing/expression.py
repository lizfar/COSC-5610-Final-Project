"""Read gene-count files and convert them with log2(CPM + 1)."""

from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from config import (
    BEAT_COUNTS_PATH,
    BEAT_MAPPING_PATH,
    CPM_SCALE,
    LITERATURE_GENES,
    PARSE_WORKERS,
)
from preprocessing.harmonize import clean_symbol
from preprocessing.logutil import warn


def attach_and_transform(clinical, counts, library, label, log):
    """Keep clinical patients who have counts, then store log2(CPM + 1)."""
    try:
        counts = counts.copy()
        counts.index = counts.index.astype(str)
        library = library.copy()
        library.index = library.index.astype(str)
        shared = clinical["patient_id"].isin(counts.index) & clinical["patient_id"].isin(library.index)
        missing = clinical.loc[~shared, "patient_id"].tolist()
        log.add(
            f"{label} expression join: {int(shared.sum())} of {len(clinical)} clinical patients have gene counts."
        )
        if missing:
            log.add(
                f"{label} clinical patients without counts: {len(missing)}. Examples: {missing[:12]}."
            )
        model = clinical.loc[shared].copy()
        model["library_size"] = model["patient_id"].map(library).astype(float)
        usable = model["library_size"] > 0
        if (~usable).any():
            log.add_warn(
                f"{label} dropped {int((~usable).sum())} patients with a zero library size."
            )
            model = model.loc[usable].copy()
        expression = counts.loc[model["patient_id"]].fillna(0)
        expressed = expression.sum(axis=0) > 0
        log.add(
            f"{label} genes with a count in at least one joined patient: "
            f"{int(expressed.sum())} of {expression.shape[1]}."
        )
        expression = expression.loc[:, expressed]
        library_aligned = model.set_index("patient_id")["library_size"].astype("float64")
        count_values = expression.to_numpy(dtype="float64")
        library_values = library_aligned.to_numpy(dtype="float64")
        log.add(
            f"{label} raw count range {np.nanmin(count_values):.4g} to {np.nanmax(count_values):.4g}. "
            f"Library size range {np.nanmin(library_values):.4g} to {np.nanmax(library_values):.4g}."
        )
        if library_values.shape[0] != count_values.shape[0]:
            raise RuntimeError(
                f"{label} has {count_values.shape[0]} count rows and {library_values.shape[0]} library sizes."
            )
        per_million = count_values / library_values[:, None] * CPM_SCALE
        log_cpm = pd.DataFrame(
            np.log2(per_million + 1.0),
            index=expression.index,
            columns=expression.columns,
        )
        sparse = int(((expression == 0).mean(axis=0) > 0.75).sum())
        log.add(
            f"{label} log2(CPM + 1) matrix: {log_cpm.shape[0]} patients by {log_cpm.shape[1]} genes. "
            f"Genes that are zero in more than 75% of patients: {sparse}. "
            "Those genes stay in the matrix. Standardization and data-driven gene selection wait until training."
        )
        found = [gene for gene in LITERATURE_GENES if gene in log_cpm.columns]
        absent = [gene for gene in LITERATURE_GENES if gene not in log_cpm.columns]
        log.add(f"{label} literature genes present: {found}. Absent: {absent}.")
        return model, log_cpm
    except Exception as exc:
        raise RuntimeError(f"{label} expression transform failed: {exc}") from exc


def load_beat_counts(clinical, log):
    log.add(f"Beat AML: reading the gene-count matrix {BEAT_COUNTS_PATH.name}.")
    try:
        raw = _read_count_table(BEAT_COUNTS_PATH)
    except Exception as exc:
        raise RuntimeError(f"Could not read Beat AML counts: {exc}") from exc
    log.add(f"Beat AML count table shape before cleanup: {raw.shape[0]} rows by {raw.shape[1]} columns.")
    counts = _collapse_count_matrix(raw, log, "Beat AML")
    wanted = clinical["expression_sample_id"].dropna().astype(str).tolist()
    column_lookup = {str(column): column for column in counts.columns}
    direct = [sample_id for sample_id in wanted if sample_id in column_lookup]
    log.add(f"Beat AML RNA sample ids found directly in the count matrix: {len(direct)} of {len(wanted)}.")
    sample_to_column = {sample_id: sample_id for sample_id in direct}
    if len(direct) < 0.8 * max(len(wanted), 1):
        sample_to_column = _match_with_beat_map(counts.columns, wanted, log)
    if len(sample_to_column) < 0.5 * max(len(wanted), 1):
        raise RuntimeError(
            f"Beat AML count columns matched only {len(sample_to_column)} of {len(wanted)} RNA sample ids."
        )
    _drop_beat_controls(sample_to_column, log)
    patient_by_sample = clinical.set_index("expression_sample_id")["patient_id"].astype(str).to_dict()
    pieces = {}
    libraries = {}
    for sample_id, column in sample_to_column.items():
        patient_id = patient_by_sample.get(sample_id)
        if patient_id is None:
            continue
        values = pd.to_numeric(counts[column], errors="coerce").fillna(0)
        pieces[patient_id] = values
        libraries[patient_id] = float(values.sum())
    patient_counts = pd.DataFrame(pieces).T.fillna(0)
    library = pd.Series(libraries, name="library_size")
    log.add(
        f"Beat AML count matrix aligned to patients: {patient_counts.shape[0]} patients, "
        f"{patient_counts.shape[1]} genes. Library size is the sum of the protein-coding counts."
    )
    return patient_counts, library


def load_star_matrix(manifest, log, label):
    log.add(f"{label}: reading {len(manifest)} STAR gene-count files.")
    tasks = [(row.local_path, row.patient_id) for row in manifest.itertuples(index=False)]
    results = []
    failures = []
    try:
        with ProcessPoolExecutor(max_workers=PARSE_WORKERS) as pool:
            futures = {pool.submit(read_star_task, task): task[1] for task in tasks}
            finished = 0
            for future in as_completed(futures):
                finished += 1
                patient_id = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    failures.append((patient_id, f"{type(exc).__name__}: {exc}"))
                    continue
                if result["error"]:
                    failures.append((patient_id, result["error"]))
                else:
                    results.append(result)
                if finished % 50 == 0 or finished == len(tasks):
                    log.add(f"{label} STAR files read: {finished} of {len(tasks)}.")
    except Exception as exc:
        raise RuntimeError(f"{label} STAR parsing failed: {exc}") from exc
    if failures:
        log.add_warn(f"{label} STAR files that could not be read: {len(failures)}. Examples: {failures[:5]}.")
    if not results:
        raise RuntimeError(f"No {label} STAR files produced gene counts.")
    first = results[0]
    log.add(
        f"{label} first parsed patient {first['patient_id']}: {len(first['counts'])} gene symbols, "
        f"library size {first['library']:.0f}. Example symbols: {list(first['counts'].index[:8])}."
    )
    patient_counts = _patient_frame_from_results(results)
    library = pd.Series({result["patient_id"]: result["library"] for result in results}, name="library_size")
    sample_map = manifest.set_index("patient_id")["sample_submitter_id"].astype(str)
    log.add(
        f"{label} count matrix: {patient_counts.shape[0]} patients by {patient_counts.shape[1]} genes. "
        "Library size is the sum of protein-coding gene counts. "
        "Summary rows beginning with N_ are excluded, and duplicate gene symbols are summed."
    )
    return patient_counts, library, sample_map


def read_star_task(task):
    path, patient_id = task
    try:
        counts, library = read_star_counts(path)
        return {"patient_id": patient_id, "counts": counts, "library": library, "error": None}
    except Exception as exc:
        return {
            "patient_id": patient_id,
            "counts": None,
            "library": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def read_star_counts(path):
    frame = pd.read_csv(path, sep="\t", comment="#", encoding="utf-8")
    columns = {str(column).lower(): column for column in frame.columns}
    for required in ("gene_id", "gene_name", "unstranded"):
        if required not in columns:
            raise KeyError(f"{path} is missing {required}. Columns: {list(frame.columns)}")
    gene_ids = frame[columns["gene_id"]].astype(str)
    ensg = frame.loc[gene_ids.str.startswith("ENSG")].copy()
    unstranded = pd.to_numeric(ensg[columns["unstranded"]], errors="coerce").fillna(0)
    symbols = ensg[columns["gene_name"]].map(clean_symbol)
    keep = symbols.notna()
    if "gene_type" in columns:
        coding = ensg[columns["gene_type"]].astype(str).str.strip().str.lower().eq("protein_coding")
        literature = symbols.isin(LITERATURE_GENES)
        keep = keep & (coding | literature.fillna(False))
    named = unstranded.loc[keep].copy()
    named.index = symbols.loc[keep].to_numpy()
    named = named.groupby(level=0).sum()
    library = float(named.sum())
    return named, library


def _read_count_table(path):
    compression = "gzip" if str(path).endswith(".gz") else None
    try:
        frame = pd.read_csv(
            path,
            sep="\t",
            encoding="utf-8",
            compression=compression,
            low_memory=False,
        )
    except UnicodeDecodeError as exc:
        warn(f"{path.name} was not valid UTF-8 ({exc}). Reading with latin-1.")
        frame = pd.read_csv(
            path,
            sep="\t",
            encoding="latin-1",
            compression=compression,
            low_memory=False,
        )
    if frame.shape[1] == 1:
        frame = pd.read_csv(
            path,
            sep=r"\s+",
            encoding="utf-8",
            compression=compression,
            engine="python",
        )
    frame = frame.set_index(frame.columns[0])
    if _sample_fraction(frame.index) > 0.5 and _sample_fraction(frame.columns) < 0.5:
        frame = frame.T
    return frame


def _collapse_count_matrix(frame, log, label):
    symbol_column = _named_column(frame.columns, ["display_label", "gene_name", "gene_symbol", "hgnc_symbol", "symbol"])
    biotype_column = _named_column(frame.columns, ["biotype", "gene_type"])
    if symbol_column is not None:
        symbols = frame[symbol_column].map(clean_symbol)
        log.add(f"{label} gene symbols read from column {symbol_column}.")
    else:
        symbols = pd.Series(frame.index.map(clean_symbol), index=frame.index)
        log.add(f"{label} gene symbols read from the first column.")
    if biotype_column is not None:
        biotype = frame[biotype_column].astype(str).str.strip().str.lower()
        log.add(f"{label} biotype counts: {biotype.value_counts(dropna=False).head(6).to_dict()}.")
    else:
        biotype = None
    numeric_columns = []
    for column in frame.columns:
        if column in {symbol_column, biotype_column}:
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().mean() > 0.9:
            numeric_columns.append(column)
    log.add(f"{label} numeric sample columns: {len(numeric_columns)}.")
    counts = frame[numeric_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    symbol_values = symbols.to_numpy()
    keep = pd.notna(symbol_values)
    if biotype is not None:
        coding = biotype.eq("protein_coding").to_numpy()
        literature = pd.Series(symbol_values).isin(LITERATURE_GENES).to_numpy()
        set_aside = int((keep & ~coding & ~literature).sum())
        log.add(
            f"{label} protein-coding filter sets aside {set_aside} named non-coding genes. "
            "Literature genes are kept when they are present."
        )
        keep = keep & (coding | literature)
    counts = counts.loc[keep].copy()
    counts.index = pd.Index(symbol_values[keep])
    before = len(counts)
    counts = counts.groupby(level=0).sum()
    log.add(
        f"{label} gene symbols: {counts.shape[0]} after summing duplicate symbols "
        f"({before - counts.shape[0]} duplicate rows combined). Example symbols: {list(counts.index[:8])}."
    )
    return counts


def _named_column(columns, candidates):
    lookup = {str(column).strip().lower(): column for column in columns}
    for candidate in candidates:
        if candidate in lookup:
            return lookup[candidate]
    return None


def _match_with_beat_map(count_columns, wanted, log):
    if not BEAT_MAPPING_PATH.exists():
        raise FileNotFoundError("Beat AML sample map is missing, and the count columns did not match RNA sample ids.")
    mapping = pd.read_excel(BEAT_MAPPING_PATH)
    log.add(f"Beat AML sample-map columns: {list(mapping.columns)}.")
    columns = {str(column).strip() for column in count_columns}
    wanted_set = set(wanted)
    id_columns = []
    for column in mapping.columns:
        values = set(mapping[column].map(lambda value: str(value).strip()))
        overlap_counts = len(values & columns)
        overlap_wanted = len(values & wanted_set)
        if overlap_counts or overlap_wanted:
            id_columns.append((column, overlap_counts, overlap_wanted))
    log.add(f"Beat AML sample-map columns that overlap ids: {id_columns}.")
    count_column = max(id_columns, key=lambda item: item[1])[0] if id_columns else None
    wanted_column = max(id_columns, key=lambda item: item[2])[0] if id_columns else None
    if count_column is None or wanted_column is None or id_columns[0][1] == 0:
        return {}
    pairs = {}
    for _, row in mapping.iterrows():
        sample_id = str(row[wanted_column]).strip()
        column = str(row[count_column]).strip()
        if sample_id in wanted_set and column in columns:
            pairs[sample_id] = column
    log.add(f"Beat AML sample map matched {len(pairs)} RNA sample ids to count columns.")
    return pairs


def _drop_beat_controls(sample_to_column, log):
    if not BEAT_MAPPING_PATH.exists():
        return
    mapping = pd.read_excel(BEAT_MAPPING_PATH)
    control_column = None
    for column in mapping.columns:
        if "control" in str(column).lower():
            control_column = column
            break
    if control_column is None:
        log.add("Beat AML sample map has no control column. Initial-diagnosis filter is the control exclusion.")
        return
    control_values = mapping[control_column].astype(str).str.strip().str.lower()
    control_rows = mapping.loc[control_values.isin({"1", "true", "y", "yes", "control"})]
    log.add(f"Beat AML sample map control rows: {len(control_rows)} using column {control_column}.")
    banned = set()
    for column in mapping.columns:
        banned.update(str(value).strip() for value in control_rows[column].dropna())
    before = len(sample_to_column)
    for sample_id in list(sample_to_column):
        if sample_id in banned or sample_to_column[sample_id] in banned:
            sample_to_column.pop(sample_id, None)
    log.add(f"Beat AML control samples removed from the count alignment: {before - len(sample_to_column)}.")


def _patient_frame_from_results(results):
    chunks = []
    buffer = []
    for result in results:
        buffer.append((result["patient_id"], result["counts"]))
        if len(buffer) == 40:
            chunks.append(pd.DataFrame({patient_id: series for patient_id, series in buffer}))
            buffer = []
    if buffer:
        chunks.append(pd.DataFrame({patient_id: series for patient_id, series in buffer}))
    matrix = pd.concat(chunks, axis=1).fillna(0).astype("float32").T
    matrix.index = matrix.index.astype(str)
    return matrix


def _sample_fraction(values):
    series = pd.Index(values).astype(str)
    return series.str.contains(r"^(?:TCGA-|BA\d+)", case=False).mean()


def _patient_from_tcga_barcode(barcode):
    text = str(barcode).strip().upper()
    if text.startswith("TCGA-") and len(text) >= 12:
        return text[:12]
    return None


def _tcga_sample_rank(barcode):
    parts = str(barcode).split("-")
    sample_type = parts[3][:2] if len(parts) > 3 else ""
    if sample_type == "03":
        return 0
    if sample_type in {"01", "09"}:
        return 1
    return 2
