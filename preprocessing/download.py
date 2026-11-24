"""Download the open clinical tables and gene-count files."""

import json
import shutil
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from config import (
    BEAT_CLINICAL_NAME,
    BEAT_COUNTS_NAME,
    BEAT_DIR,
    BEAT_MAPPING_NAME,
    BEAT_ZIP_BYTES,
    BEAT_ZIP_PATH,
    BEAT_ZENODO_URL,
    DOWNLOAD_WORKERS,
    GDC_DATA_URL,
    GDC_FILES_URL,
    PRIMARY_SAMPLE_PRIORITY,
    RAW_DIR,
    TARGET_DISCOVERY_PATH,
    TARGET_DISCOVERY_URL,
    TARGET_MANIFEST_PATH,
    TARGET_STAR_DIR,
    TARGET_VALIDATION_PATH,
    TARGET_VALIDATION_URL,
    TCGA_CLINICAL_FILE_ID,
    TCGA_CLINICAL_PATH,
    TCGA_MANIFEST_PATH,
    TCGA_STAR_DIR,
    TCGA_SUPPLEMENT_FILE_ID,
    TCGA_SUPPLEMENT_PATH,
    USER_AGENT,
)
from preprocessing.harmonize import clean_id
from preprocessing.logutil import fail, ok, step, warn


def download_reference_files():
    """Download the small clinical tables, the Beat AML archive, and the TCGA count matrix."""
    try:
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        step("Creating the raw data folder.")
        _download_beat_zip()
        _extract_beat_members()
        _download_simple(TARGET_DISCOVERY_URL, TARGET_DISCOVERY_PATH, "TARGET discovery clinical workbook")
        _download_simple(TARGET_VALIDATION_URL, TARGET_VALIDATION_PATH, "TARGET validation clinical workbook")
        _download_gdc_file(TCGA_CLINICAL_FILE_ID, TCGA_CLINICAL_PATH, None)
        _download_gdc_file(TCGA_SUPPLEMENT_FILE_ID, TCGA_SUPPLEMENT_PATH, None)
        ok("Clinical tables and the Beat AML count matrix are on disk.")
    except Exception as exc:
        fail(f"Reference file download failed: {exc}")
        raise


def download_target_expression():
    """Choose one diagnostic TARGET STAR file per clinical patient and download it."""
    try:
        manifest = _build_target_manifest()
        _download_star_files(manifest, "TARGET")
        ok("TARGET RNA-seq downloads are in place.")
        return manifest
    except Exception as exc:
        fail(f"TARGET RNA-seq download failed: {exc}")
        raise


def download_tcga_star(patient_ids, log):
    """Download one open STAR count file for each TCGA-LAML patient."""
    try:
        wanted = {str(patient_id).upper()[:12] for patient_id in patient_ids if pd.notna(patient_id)}
        log.add(
            "TCGA expression uses open GDC STAR unstranded counts. "
            "The 2012 GAF matrix identifies genes by Entrez id, such as ?|100132510_calculated, "
            "and its values are fractional, so it is not the modeling matrix."
        )
        records = _list_star_records("TCGA-LAML")
        chosen = {}
        outside = 0
        for record in records:
            patient_id = str(record["patient_id"]).upper()[:12]
            if patient_id not in wanted:
                outside += 1
                continue
            rank = -int(record["file_size"])
            current = chosen.get(patient_id)
            if current is None or rank < current["rank"]:
                record["patient_id"] = patient_id
                record["rank"] = rank
                chosen[patient_id] = record
        if not chosen:
            raise RuntimeError("No TCGA-LAML STAR count files matched the clinical patient ids.")
        manifest = pd.DataFrame(chosen.values()).drop(columns=["rank"]).sort_values("patient_id")
        TCGA_STAR_DIR.mkdir(parents=True, exist_ok=True)
        manifest["local_path"] = manifest["file_id"].map(lambda file_id: str(TCGA_STAR_DIR / f"{file_id}.tsv"))
        TCGA_MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(TCGA_MANIFEST_PATH, index=False, encoding="utf-8")
        log.add(
            f"TCGA STAR files selected: {len(manifest)} of {len(wanted)} clinical patients. "
            f"{outside} STAR files belonged to ids outside the cleaned clinical table. "
            "Every selected file is primary peripheral blood."
        )
        _download_star_files(manifest, "TCGA")
        return manifest
    except Exception as exc:
        fail(f"TCGA STAR download failed: {exc}")
        raise


def _download_simple(url, dest, label):
    step(f"Downloading {label}.")
    try:
        _stream_download(url, dest, expected_size=None)
        ok(f"{label} saved to {dest.name} ({dest.stat().st_size} bytes).")
    except Exception as exc:
        fail(f"Could not download {label}: {exc}")
        raise


def _download_beat_zip():
    step("Downloading the Beat AML 2.0 Zenodo archive (clinical table, sample map, and gene counts).")
    try:
        _stream_download(BEAT_ZENODO_URL, BEAT_ZIP_PATH, expected_size=BEAT_ZIP_BYTES)
        ok(f"Beat AML archive is {BEAT_ZIP_PATH.stat().st_size} bytes.")
    except Exception as exc:
        fail(f"Could not download the Beat AML archive: {exc}")
        raise


def _extract_beat_members():
    step("Extracting the Beat AML clinical workbook, sample map, and gene-count matrix.")
    wanted = {BEAT_CLINICAL_NAME, BEAT_COUNTS_NAME, BEAT_MAPPING_NAME}
    BEAT_DIR.mkdir(parents=True, exist_ok=True)
    found = set()
    try:
        with zipfile.ZipFile(BEAT_ZIP_PATH) as archive:
            for info in archive.infolist():
                base = Path(info.filename).name
                if base not in wanted:
                    continue
                target = BEAT_DIR / base
                with archive.open(info) as source, open(target, "wb") as output:
                    shutil.copyfileobj(source, output)
                found.add(base)
                ok(f"Extracted {base} ({target.stat().st_size} bytes).")
    except Exception as exc:
        fail(f"Could not extract the Beat AML archive: {exc}")
        raise
    missing = wanted - found
    if missing:
        raise FileNotFoundError(f"Beat AML archive did not contain {sorted(missing)}.")


def _download_gdc_file(file_id, dest, expected_size):
    step(f"Downloading GDC file {file_id} to {dest.name}.")
    try:
        _stream_download(GDC_DATA_URL.format(file_id=file_id), dest, expected_size=expected_size)
        ok(f"Saved {dest.name} ({dest.stat().st_size} bytes).")
    except Exception as exc:
        fail(f"Could not download GDC file {file_id}: {exc}")
        raise


def _stream_download(url, dest, expected_size, quiet=False):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected_size and dest.stat().st_size == expected_size:
        if not quiet:
            ok(f"Using existing {dest.name}; size matches {expected_size} bytes.")
        return
    if dest.exists() and expected_size is None and dest.stat().st_size > 1000:
        if not quiet:
            ok(f"Using existing {dest.name} ({dest.stat().st_size} bytes).")
        return
    partial = dest.with_name(dest.name + ".partial")
    last_error = None
    for attempt in range(1, 6):
        try:
            if dest.exists() and expected_size and dest.stat().st_size == expected_size:
                return
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=300) as response, open(partial, "wb") as output:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
                output.flush()
            if expected_size and partial.stat().st_size != expected_size:
                raise IOError(
                    f"Downloaded {partial.stat().st_size} bytes for {dest.name}; expected {expected_size}."
                )
            _replace_partial(partial, dest)
            return
        except Exception as exc:
            last_error = exc
            if not quiet:
                warn(f"Download attempt {attempt} failed for {dest.name}: {exc}")
            _remove_partial(partial)
            if attempt == 5:
                break
            time.sleep(attempt)
    raise RuntimeError(f"Download failed for {url}: {last_error}")


def _replace_partial(partial, dest):
    last_error = None
    for attempt in range(1, 6):
        try:
            partial.replace(dest)
            return
        except PermissionError as exc:
            last_error = exc
            time.sleep(attempt)
    raise PermissionError(f"Could not move {partial.name} into place: {last_error}")


def _remove_partial(partial):
    try:
        if partial.exists():
            partial.unlink()
    except OSError as exc:
        warn(f"Could not remove the temporary file {partial.name}: {exc}")


def _build_target_manifest():
    step("Reading TARGET clinical IDs so only those patients' RNA-seq files are downloaded.")
    try:
        discovery = pd.read_excel(TARGET_DISCOVERY_PATH, sheet_name="Clinical Data")
        validation = pd.read_excel(TARGET_VALIDATION_PATH, sheet_name="Clinical Data")
    except Exception as exc:
        fail(f"Could not read TARGET clinical workbooks: {exc}")
        raise
    clinical_ids = set()
    for frame in (discovery, validation):
        for value in frame["TARGET USI"]:
            identifier = clean_id(value)
            if identifier:
                clinical_ids.add(identifier.upper())
    step(f"TARGET clinical union contains {len(clinical_ids)} patient IDs.")
    records = _list_target_star_records()
    chosen = {}
    skipped_non_clinical = 0
    skipped_type = 0
    for record in records:
        patient_id = record["patient_id"].upper()
        if patient_id not in clinical_ids:
            skipped_non_clinical += 1
            continue
        sample_type = record["sample_type"].lower()
        if sample_type not in PRIMARY_SAMPLE_PRIORITY:
            skipped_type += 1
            continue
        rank = (PRIMARY_SAMPLE_PRIORITY[sample_type], -int(record["file_size"]))
        current = chosen.get(patient_id)
        if current is None or rank < current["rank"]:
            record["rank"] = rank
            record["patient_id"] = patient_id
            chosen[patient_id] = record
    manifest = pd.DataFrame(chosen.values())
    if manifest.empty:
        raise RuntimeError("No TARGET STAR count files matched the clinical patient IDs.")
    manifest = manifest.drop(columns=["rank"]).sort_values("patient_id")
    manifest["local_path"] = manifest["file_id"].map(lambda file_id: str(TARGET_STAR_DIR / f"{file_id}.tsv"))
    TARGET_MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(TARGET_MANIFEST_PATH, index=False, encoding="utf-8")
    bone_marrow = int(
        manifest["sample_type"].str.lower().eq("primary blood derived cancer - bone marrow").sum()
    )
    blood = int(manifest["sample_type"].str.lower().eq("primary blood derived cancer - peripheral blood").sum())
    missing = len(clinical_ids - set(manifest["patient_id"]))
    ok(
        "TARGET file choice: "
        f"{len(manifest)} diagnostic RNA-seq files "
        f"({bone_marrow} bone marrow, {blood} peripheral blood). "
        f"{missing} clinical patients have no primary STAR file. "
        f"{skipped_non_clinical} primary files belonged to patients outside the clinical workbooks."
    )
    return manifest


def _list_target_star_records():
    return _list_star_records("TARGET-AML")


def _list_star_records(project_id):
    step(f"Asking the GDC API for open {project_id} STAR gene-count files.")
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "access", "value": "open"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.sample_type",
                    "value": list(PRIMARY_SAMPLE_PRIORITY),
                },
            },
        ],
    }
    fields = ",".join(
        [
            "file_id",
            "file_size",
            "cases.submitter_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
        ]
    )
    records = []
    offset = 0
    page_size = 1000
    total = None
    while total is None or offset < total:
        body = json.dumps(
            {"filters": filters, "fields": fields, "size": page_size, "from": offset}
        ).encode("utf-8")
        request = urllib.request.Request(
            GDC_FILES_URL,
            data=body,
            headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = json.load(response)
        except Exception as exc:
            fail(f"GDC file listing failed at offset {offset}: {exc}")
            raise
        hits = payload["data"]["hits"]
        total = payload["data"]["pagination"]["total"]
        for hit in hits:
            parsed = _record_from_hit(hit)
            if parsed is not None:
                records.append(parsed)
        offset += len(hits)
        step(f"{project_id} listing progress: {offset} of {total} primary STAR files.")
        if not hits:
            break
    ok(f"Parsed {len(records)} {project_id} primary STAR file records.")
    return records


def _record_from_hit(hit):
    cases = hit.get("cases") or []
    if not cases:
        return None
    case = cases[0]
    patient_id = clean_id(case.get("submitter_id"))
    if not patient_id:
        return None
    chosen_sample = None
    chosen_rank = None
    for sample in case.get("samples") or []:
        sample_type = str(sample.get("sample_type") or "")
        rank = PRIMARY_SAMPLE_PRIORITY.get(sample_type.lower())
        if rank is None:
            continue
        if chosen_rank is None or rank < chosen_rank:
            chosen_rank = rank
            chosen_sample = sample
    if chosen_sample is None:
        return None
    return {
        "patient_id": patient_id,
        "sample_submitter_id": clean_id(chosen_sample.get("submitter_id")),
        "sample_type": str(chosen_sample.get("sample_type")),
        "file_id": hit.get("file_id") or hit.get("id"),
        "file_size": int(hit.get("file_size") or 0),
    }


def _download_star_files(manifest, label):
    jobs = []
    for row in manifest.itertuples(index=False):
        dest = Path(row.local_path)
        jobs.append((row.file_id, dest, int(row.file_size)))
    pending = [job for job in jobs if not (job[1].exists() and job[1].stat().st_size == job[2])]
    step(
        f"{label} RNA-seq files: {len(jobs) - len(pending)} already downloaded, "
        f"{len(pending)} still to download."
    )
    if not pending:
        ok(f"Every selected {label} STAR file is already on disk.")
        return
    failed = []
    completed = 0
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as pool:
        futures = {
            pool.submit(_stream_download, GDC_DATA_URL.format(file_id=file_id), dest, size, True): file_id
            for file_id, dest, size in pending
        }
        for future in as_completed(futures):
            file_id = futures[future]
            completed += 1
            try:
                future.result()
            except Exception as exc:
                failed.append((file_id, str(exc)))
                fail(f"{label} file {file_id} failed: {exc}")
            if completed % 25 == 0 or completed == len(pending):
                step(f"{label} RNA-seq download progress: {completed} of {len(pending)}.")
    if failed:
        raise RuntimeError(f"{len(failed)} {label} RNA-seq files failed. Re-run to retry the missing files.")
    ok(f"Downloaded {len(pending)} {label} STAR gene-count files.")
