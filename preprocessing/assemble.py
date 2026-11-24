"""Run cohort cleaning, join expression, and write the preprocessing report."""

from datetime import datetime

import pandas as pd

from config import (
    LITERATURE_GENES,
    MIN_BEAT_PATIENTS,
    MIN_TARGET_PATIENTS,
    MIN_TCGA_PATIENTS,
    OS_LANDMARK_DAYS,
    PROCESSED_DIR,
    REPORT_PATH,
    TARGET_MANIFEST_PATH,
)
from preprocessing.clinical import prepare_beat, prepare_target, prepare_tcga
from preprocessing.download import download_tcga_star
from preprocessing.expression import attach_and_transform, load_beat_counts, load_star_matrix
from preprocessing.logutil import StepLog, fail, ok, step

SHARED_RULES = [
    "Shared preprocessing rules",
    "Each modeling table has one row per patient.",
    "A patient is kept only when the sample is from diagnosis, gene counts exist, vital status is Alive or Dead, and overall survival is a positive number of days.",
    "Adult cohorts (Beat AML and TCGA-LAML) and the pediatric cohort (TARGET) are cleaned to the same core columns and saved separately.",
    "Overall survival stays in days. The death indicator is 1 for Dead and 0 for Alive.",
    f"The 2-year label is 1 when death occurs on or before day {OS_LANDMARK_DAYS}.",
    f"The 2-year label is 0 when follow-up reaches day {OS_LANDMARK_DAYS}, including deaths after that day.",
    "A patient who is still alive and last seen before day 730 keeps a survival row and has no 2-year class.",
    "FAB is mapped to M0-M7. Blank, NOS, Not classified, and ambiguous values become UNKNOWN. M6 and M7 stay as their own levels.",
    "Cytogenetic risk is mapped to FAVORABLE, INTERMEDIATE, ADVERSE, or UNKNOWN. Unknown is kept as its own level.",
    "Age, FAB, and risk are not imputed.",
    "Gene ids are converted to uppercase gene symbols. Duplicate symbols are summed. Protein-coding genes are kept.",
    "Expression is log2(CPM + 1). CPM divides by that sample's own library size and multiplies by 1,000,000.",
    "Z-scoring and data-driven gene selection are not done in preprocessing.",
    f"Literature genes checked in every matrix: {', '.join(LITERATURE_GENES)}.",
]


def build_reference_cohorts():
    try:
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        sections = []
        sections.extend(_build_beat())
        sections.extend(_build_tcga())
        return sections
    except Exception as exc:
        fail(f"Could not build the Beat AML or TCGA tables: {exc}")
        raise


def build_target_cohort(manifest=None):
    try:
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        if manifest is None:
            manifest = pd.read_csv(TARGET_MANIFEST_PATH, encoding="utf-8")
        return _build_target(manifest)
    except Exception as exc:
        fail(f"Could not build the TARGET table: {exc}")
        raise


def write_full_report(cohort_sections):
    try:
        _write_report(list(SHARED_RULES) + list(cohort_sections))
        ok(f"Preprocessing report written to {REPORT_PATH}.")
    except Exception as exc:
        fail(f"Could not write the preprocessing report: {exc}")
        raise


def _build_beat():
    log = StepLog("Beat AML")
    clinical = prepare_beat(log)
    _save_table(clinical, "beat_clinical_before_expression.csv")
    counts, library = load_beat_counts(clinical, log)
    model, log_cpm = attach_and_transform(clinical, counts, library, "Beat AML", log)
    _require_patients(model, MIN_BEAT_PATIENTS, "Beat AML")
    _save_outputs("beat", model, log_cpm, log)
    return log.lines


def _build_target(manifest):
    log = StepLog("TARGET-AML")
    clinical = prepare_target(log)
    counts, library, sample_map = load_star_matrix(manifest, log, "TARGET")
    clinical = clinical.copy()
    clinical["expression_sample_id"] = clinical["patient_id"].map(sample_map)
    _save_table(clinical, "target_clinical_before_expression.csv")
    model, log_cpm = attach_and_transform(clinical, counts, library, "TARGET", log)
    bone_marrow = int(
        model["expression_sample_id"].astype(str).str.contains("-09", regex=False).sum()
    )
    log.add(
        f"TARGET joined patients whose selected sample id contains the bone-marrow code 09: {bone_marrow} "
        f"of {len(model)}."
    )
    _require_patients(model, MIN_TARGET_PATIENTS, "TARGET")
    _save_outputs("target", model, log_cpm, log)
    return log.lines


def _build_tcga():
    log = StepLog("TCGA-LAML")
    clinical = prepare_tcga(log)
    manifest = download_tcga_star(clinical["patient_id"].tolist(), log)
    counts, library, sample_map = load_star_matrix(manifest, log, "TCGA")
    clinical = clinical.copy()
    clinical["expression_sample_id"] = clinical["patient_id"].map(sample_map)
    _save_table(clinical, "tcga_clinical_before_expression.csv")
    model, log_cpm = attach_and_transform(clinical, counts, library, "TCGA", log)
    _require_patients(model, MIN_TCGA_PATIENTS, "TCGA")
    _save_outputs("tcga", model, log_cpm, log)
    return log.lines


def _save_outputs(prefix, model, log_cpm, log):
    patient_path = PROCESSED_DIR / f"{prefix}_patients.csv"
    expression_path = PROCESSED_DIR / f"{prefix}_logcpm.csv.gz"
    literature_path = PROCESSED_DIR / f"{prefix}_literature_genes.csv"
    model.to_csv(patient_path, index=False, encoding="utf-8")
    log_cpm.to_csv(expression_path, compression="gzip", encoding="utf-8", float_format="%.5f")
    present = [gene for gene in LITERATURE_GENES if gene in log_cpm.columns]
    log_cpm.loc[:, present].to_csv(literature_path, encoding="utf-8", float_format="%.5f")
    log.add_ok(
        f"Saved {patient_path.name} ({len(model)} patients), {expression_path.name}, and {literature_path.name}."
    )


def _save_table(frame, name):
    path = PROCESSED_DIR / name
    frame.to_csv(path, index=False, encoding="utf-8")
    step(f"Saved {name} ({len(frame)} rows).")


def _require_patients(model, minimum, label):
    if len(model) < minimum:
        raise RuntimeError(
            f"{label} has {len(model)} patients after the expression join. "
            f"The run stops below {minimum} so a bad ID join is not used for modeling."
        )


def _write_report(sections):
    stamp = datetime.now().isoformat(timespec="seconds")
    lines = [
        "AML prognosis preprocessing report",
        f"Generated: {stamp}",
        "",
        *sections,
        "",
    ]
    with open(REPORT_PATH, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
