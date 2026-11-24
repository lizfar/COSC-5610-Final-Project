"""Turn each cohort's clinical export into one row per patient."""

import pandas as pd

from config import (
    BEAT_CLINICAL_PATH,
    TARGET_DISCOVERY_PATH,
    TARGET_LESION_COLUMNS,
    TARGET_VALIDATION_PATH,
    TCGA_CLINICAL_PATH,
    TCGA_SUPPLEMENT_PATH,
)
from preprocessing.harmonize import (
    add_two_year_label,
    age_from_days,
    clean_id,
    complexity_to_number,
    find_column,
    lesion_column_name,
    map_event,
    map_fab,
    map_risk,
    map_yes_no,
    to_number,
    unexpected_risk_token,
)
from preprocessing.logutil import StepLog


def prepare_beat(log: StepLog):
    log.add("Beat AML: reading the clinical summary sheet from the Zenodo workbook.")
    frame = pd.read_excel(BEAT_CLINICAL_PATH, sheet_name="summary")
    required = [
        "dbgap_subject_id",
        "dbgap_rnaseq_sample",
        "rnaSeq",
        "diseaseStageAtSpecimenCollection",
        "vitalStatus",
        "overallSurvival",
        "ageAtDiagnosis",
        "fabBlastMorphology",
        "ELN2017",
        "FLT3-ITD",
        "NPM1",
        "isDenovo",
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise KeyError(f"Beat AML clinical sheet is missing {missing}. Columns: {list(frame.columns)}")
    log.add(f"Beat AML source table: {len(frame)} specimen rows and {frame['dbgap_subject_id'].nunique()} patients.")
    frame = frame.copy()
    frame["patient_id"] = frame["dbgap_subject_id"].map(clean_id)
    frame["expression_sample_id"] = frame["dbgap_rnaseq_sample"].map(clean_id)

    rna = frame["rnaSeq"].astype(str).str.strip().str.lower().eq("y")
    frame = _keep(frame, rna, "Beat AML RNA-seq flag is y", log)
    stage = frame["diseaseStageAtSpecimenCollection"].astype(str).str.strip().eq("Initial Diagnosis")
    frame = _keep(frame, stage, "Beat AML specimen stage is Initial Diagnosis", log)
    vital = frame["vitalStatus"].map(map_event).notna()
    frame = _keep(frame, vital, "Beat AML vital status is Alive or Dead", log)
    survival = to_number(frame["overallSurvival"])
    frame = frame.assign(_os_days=survival)
    frame = _keep(frame, frame["_os_days"] > 0, "Beat AML overall survival is a positive number of days", log)
    has_sample = frame["expression_sample_id"].notna()
    frame = _keep(frame, has_sample, "Beat AML row has an RNA sample id", log)

    before = len(frame)
    frame = frame.sort_values(["patient_id", "expression_sample_id"])
    duplicated = frame["patient_id"].duplicated(keep=False)
    if duplicated.any():
        dropped = frame.loc[frame["patient_id"].duplicated(keep="first"), "expression_sample_id"].tolist()
        log.add(
            "Beat AML patients with two qualifying initial RNA samples: "
            f"{frame.loc[duplicated, 'patient_id'].nunique()}. "
            f"Kept the first sample id in sorted order and dropped {dropped}."
        )
    frame = frame.drop_duplicates("patient_id", keep="first")
    log.add(f"Beat AML one-row-per-patient rule: {len(frame)} patients kept from {before} rows.")

    patients = pd.DataFrame(
        {
            "patient_id": frame["patient_id"].values,
            "cohort": "BEAT",
            "age_years": to_number(frame["ageAtDiagnosis"]).values,
            "fab_raw": frame["fabBlastMorphology"].values,
            "cytogenetic_risk_raw": frame["ELN2017"].values,
            "os_days": frame["_os_days"].values,
            "event": frame["vitalStatus"].map(map_event).astype(int).values,
            "expression_sample_id": frame["expression_sample_id"].values,
            "flt3_itd": frame["FLT3-ITD"].map(map_yes_no).values,
            "npm1": frame["NPM1"].map(map_yes_no).values,
            "is_denovo": frame["isDenovo"].map(map_yes_no).values,
        }
    )
    return _finish_clinical(patients, "Beat AML", log)


def prepare_target(log: StepLog):
    log.add("TARGET: reading the discovery and validation clinical workbooks.")
    discovery = pd.read_excel(TARGET_DISCOVERY_PATH, sheet_name="Clinical Data")
    validation = pd.read_excel(TARGET_VALIDATION_PATH, sheet_name="Clinical Data")
    combined = _deduplicate_target(discovery, validation, log)
    combined["patient_id"] = combined["TARGET USI"].map(clean_id)
    combined = combined[combined["patient_id"].notna()].copy()
    log.add(f"TARGET rows with a patient id: {len(combined)}.")

    event = combined["Vital Status"].map(map_event)
    combined = combined.assign(_event=event)
    combined = _keep(combined, combined["_event"].notna(), "TARGET vital status is Alive or Dead", log)
    survival = to_number(combined["Overall Survival Time in Days"])
    combined = combined.assign(_os_days=survival)
    combined = _keep(combined, combined["_os_days"] > 0, "TARGET overall survival is a positive number of days", log)

    lesion_frame = {}
    for label in TARGET_LESION_COLUMNS:
        if label not in combined.columns:
            log.add_warn(f"TARGET lesion column {label} was not in the workbook.")
            continue
        lesion_frame[lesion_column_name(label)] = combined[label].map(map_yes_no)
    raw_example = combined["t(8;21)"].value_counts(dropna=False).to_dict() if "t(8;21)" in combined.columns else {}
    log.add(f"TARGET raw t(8;21) values before yes/no coding: {raw_example}.")

    patients = pd.DataFrame(
        {
            "patient_id": combined["patient_id"].values,
            "cohort": "TARGET",
            "target_subset": combined["target_subset"].values,
            "age_years": age_from_days(combined["Age at Diagnosis in Days"]).values,
            "fab_raw": combined["FAB Category"].values,
            "cytogenetic_risk_raw": combined["Risk group"].values,
            "os_days": combined["_os_days"].values,
            "event": combined["_event"].astype(int).values,
            "wbc": to_number(combined["WBC at Diagnosis"]).values,
            "cytogenetic_complexity_raw": combined["Cytogenetic Complexity"].astype(str).values,
            "cytogenetic_complexity": combined["Cytogenetic Complexity"].map(complexity_to_number).values,
            "primary_cytogenetic_code": combined["Primary Cytogenetic Code"].astype(str).values,
        }
    )
    for name, values in lesion_frame.items():
        patients[name] = values.values
    patients["expression_sample_id"] = pd.NA
    return _finish_clinical(patients, "TARGET", log)


def prepare_tcga(log: StepLog):
    log.add("TCGA-LAML: reading the open clinical patient table.")
    clinical = _read_tcga_clinical_table(TCGA_CLINICAL_PATH, log)
    log.add(f"TCGA clinical columns ({len(clinical.columns)}): {list(clinical.columns)}.")
    barcode_column = find_column(clinical, ["bcr_patient_barcode", "patient_barcode", "patient_id"])
    if barcode_column is None:
        raise KeyError(f"TCGA clinical table has no patient barcode column. Columns: {list(clinical.columns)}")
    vital_column = find_column(clinical, ["vital_status"])
    death_column = find_column(clinical, ["days_to_death"])
    follow_column = find_column(clinical, ["days_to_last_followup", "days_to_last_follow_up"])
    age_column = find_column(clinical, ["age_at_initial_pathologic_diagnosis", "age_at_diagnosis", "age_at_index"])
    fab_column = find_column(clinical, ["leukemia_french_american_british_morphology_code", "fab"])
    risk_column = find_column(
        clinical,
        ["acute_myeloid_leukemia_calgb_cytogenetics_risk_category", "cytogenetics_risk_category"],
    )
    cyto_column = find_column(clinical, ["cytogenetic_abnormality"])
    gender_column = find_column(clinical, ["gender", "sex"])
    log.add(
        "TCGA fields: "
        f"barcode={barcode_column}, vital={vital_column}, death={death_column}, "
        f"follow-up={follow_column}, age={age_column}, FAB={fab_column}, risk={risk_column}."
    )
    if vital_column is None or (death_column is None and follow_column is None):
        raise KeyError("TCGA clinical table is missing vital status or survival time.")

    patients = pd.DataFrame({"patient_id": clinical[barcode_column].map(clean_id)})
    patients = patients[patients["patient_id"].notna()].copy()
    patients["patient_id"] = patients["patient_id"].str.slice(0, 12)
    patients["cohort"] = "TCGA"
    patients["event"] = clinical.loc[patients.index, vital_column].map(map_event)
    death_days = to_number(clinical.loc[patients.index, death_column]) if death_column else pd.Series(index=patients.index, dtype="float64")
    follow_days = to_number(clinical.loc[patients.index, follow_column]) if follow_column else pd.Series(index=patients.index, dtype="float64")
    survival = follow_days.copy()
    use_death = (patients["event"] == 1) & death_days.notna()
    survival.loc[use_death] = death_days.loc[use_death]
    still_missing = survival.isna() & death_days.notna()
    survival.loc[still_missing] = death_days.loc[still_missing]
    patients["os_days"] = survival
    patients["age_years"] = to_number(clinical.loc[patients.index, age_column]) if age_column else pd.NA
    patients["fab_raw"] = clinical.loc[patients.index, fab_column] if fab_column else pd.NA
    patients["cytogenetic_risk_raw"] = clinical.loc[patients.index, risk_column] if risk_column else pd.NA
    patients["cytogenetic_abnormality"] = clinical.loc[patients.index, cyto_column] if cyto_column else pd.NA
    patients["gender"] = clinical.loc[patients.index, gender_column] if gender_column else pd.NA
    patients = _fill_tcga_from_supplement(patients, log)

    before = len(patients)
    patients = patients[patients["event"].notna() & (patients["os_days"] > 0)].copy()
    log.add(
        f"TCGA Alive/Dead status and positive survival time: {len(patients)} patients kept, "
        f"{before - len(patients)} removed."
    )
    patients["event"] = patients["event"].astype(int)
    patients["expression_sample_id"] = pd.NA
    patients = patients.drop_duplicates("patient_id", keep="first")
    return _finish_clinical(patients, "TCGA", log)


def _finish_clinical(patients, label, log):
    patients["fab"] = patients["fab_raw"].map(map_fab)
    patients["cytogenetic_risk"] = patients["cytogenetic_risk_raw"].map(map_risk)
    unexpected = sorted(
        {
            str(value).strip()
            for value in patients["cytogenetic_risk_raw"]
            if unexpected_risk_token(value)
        }
    )
    if unexpected:
        log.add_warn(f"{label} risk values mapped to UNKNOWN because they were not in the risk map: {unexpected}.")
    patients = add_two_year_label(patients)
    patients = patients.sort_values("patient_id").reset_index(drop=True)
    known_2yr = patients["os_2yr_death"].notna().sum()
    deaths_2yr = int((patients["os_2yr_death"] == 1).sum())
    alive_2yr = int((patients["os_2yr_death"] == 0).sum())
    censored_early = int(patients["os_2yr_death"].isna().sum())
    log.add(
        f"{label} clinical modeling rows before the expression join: {len(patients)}. "
        f"Deaths: {int(patients['event'].sum())}. "
        f"Missing age: {int(patients['age_years'].isna().sum())}. "
        f"FAB unknown: {int((patients['fab'] == 'UNKNOWN').sum())}. "
        f"Risk unknown: {int((patients['cytogenetic_risk'] == 'UNKNOWN').sum())}."
    )
    log.add(
        f"{label} 2-year death label: {deaths_2yr} dead by day 730, {alive_2yr} followed at least 730 days, "
        f"{censored_early} censored before day 730 (kept for survival, left out of classification). "
        f"Known 2-year labels: {known_2yr}."
    )
    return patients


def _keep(frame, mask, reason, log):
    before = len(frame)
    kept = frame.loc[mask].copy()
    log.add(f"{reason}: {len(kept)} rows kept, {before - len(kept)} removed (from {before}).")
    return kept


def _deduplicate_target(discovery, validation, log):
    discovery = discovery.copy()
    validation = validation.copy()
    discovery["TARGET USI"] = discovery["TARGET USI"].map(clean_id)
    validation["TARGET USI"] = validation["TARGET USI"].map(clean_id)
    discovery = discovery[discovery["TARGET USI"].notna()]
    validation = validation[validation["TARGET USI"].notna()]
    overlap = set(discovery["TARGET USI"]) & set(validation["TARGET USI"])
    log.add(
        f"TARGET discovery rows: {len(discovery)}. Validation rows: {len(validation)}. "
        f"Patients in both workbooks: {len(overlap)}."
    )
    key_columns = [
        "Vital Status",
        "Overall Survival Time in Days",
        "Age at Diagnosis in Days",
        "FAB Category",
        "Risk group",
    ]
    merged = discovery.merge(validation, on="TARGET USI", suffixes=("_discovery", "_validation"))
    for column in key_columns:
        same = (
            merged[f"{column}_discovery"].astype(str).fillna("")
            == merged[f"{column}_validation"].astype(str).fillna("")
        ).mean()
        log.add(f"TARGET overlap agreement for {column}: {same:.3f}.")
    discovery_index = discovery.set_index("TARGET USI")
    validation = validation.copy()
    filled = 0
    for column in validation.columns:
        if column == "TARGET USI" or column not in discovery_index.columns:
            continue
        missing = validation[column].isna() & validation["TARGET USI"].isin(overlap)
        if missing.any():
            filled_values = validation.loc[missing, "TARGET USI"].map(discovery_index[column])
            if pd.api.types.is_numeric_dtype(validation[column]):
                filled_values = pd.to_numeric(filled_values, errors="coerce")
            validation.loc[missing, column] = filled_values
            filled += int(missing.sum())
    log.add(
        f"TARGET overlap rows keep the validation value when both workbooks have one. "
        f"Blank validation cells filled from discovery: {filled}."
    )
    validation["target_subset"] = "VALIDATION"
    validation.loc[validation["TARGET USI"].isin(overlap), "target_subset"] = "BOTH"
    discovery_only = discovery[~discovery["TARGET USI"].isin(set(validation["TARGET USI"]))].copy()
    discovery_only["target_subset"] = "DISCOVERY"
    combined = pd.concat([validation, discovery_only], ignore_index=True)
    log.add(
        "TARGET patient table after removing duplicate IDs: "
        f"{len(combined)} patients "
        f"(discovery only {(combined['target_subset'] == 'DISCOVERY').sum()}, "
        f"validation only {(combined['target_subset'] == 'VALIDATION').sum()}, "
        f"in both {(combined['target_subset'] == 'BOTH').sum()})."
    )
    return combined


def _read_tcga_clinical_table(path, log):
    try:
        with open(path, encoding="utf-8") as handle:
            lines = handle.readlines()
    except UnicodeDecodeError as exc:
        log.add_warn(f"TCGA clinical file was not valid UTF-8 ({exc}). Reading with latin-1.")
        with open(path, encoding="latin-1") as handle:
            lines = handle.readlines()
    header_idx = 0
    for index, line in enumerate(lines[:40]):
        lowered = line.lower()
        if "bcr_patient_barcode" in lowered or "vital_status" in lowered:
            header_idx = index
            break
    log.add(f"TCGA clinical header row index: {header_idx}.")
    frame = pd.read_csv(path, sep="\t", encoding="utf-8", skiprows=header_idx, dtype=str)
    barcode_column = find_column(frame, ["bcr_patient_barcode", "patient_barcode"])
    if barcode_column is not None:
        frame = frame[frame[barcode_column].astype(str).str.startswith("TCGA-")].copy()
    log.add(f"TCGA clinical patient rows: {len(frame)}.")
    return frame


def _fill_tcga_from_supplement(patients, log):
    """Fill missing or zero BCR survival times from the NEJM supplemental OS months."""
    if not TCGA_SUPPLEMENT_PATH.exists():
        log.add_warn("TCGA supplemental table was not downloaded. BCR survival times used alone.")
        return patients
    try:
        header = pd.read_excel(TCGA_SUPPLEMENT_PATH, nrows=0)
        columns = list(header.columns)
        id_column = next(column for column in columns if str(column).strip().lower() == "tcga patient id")
        os_column = next(column for column in columns if str(column).strip().lower().startswith("os months"))
        supplement = pd.read_excel(TCGA_SUPPLEMENT_PATH, usecols=[id_column, os_column])
    except Exception as exc:
        log.add_warn(f"TCGA supplemental survival columns could not be read ({exc}). BCR times used alone.")
        return patients
    supplement = supplement.copy()
    supplement["id4"] = supplement[id_column].map(_tcga_short_id)
    supplement = supplement[supplement["id4"].notna()].drop_duplicates("id4")
    months = pd.to_numeric(supplement.set_index("id4")[os_column], errors="coerce")
    days = months * (365.25 / 12.0)
    patients = patients.copy()
    short_id = patients["patient_id"].map(_tcga_short_id)
    supplemental_days = short_id.map(days)
    needs_time = ~(patients["os_days"] > 0)
    can_fill = needs_time & supplemental_days.notna() & (supplemental_days > 0)
    patients.loc[can_fill, "os_days"] = supplemental_days.loc[can_fill]
    log.add(
        "TCGA supplemental table stores overall survival in months as of 2012-03-31, "
        "keyed by the last four digits of the TCGA barcode. "
        f"BCR time was missing or zero for {int(needs_time.sum())} patients. "
        f"{int(can_fill.sum())} of those received a positive time from OS months times 365.25/12. "
        "Patients who already had a positive BCR time kept that time."
    )
    return patients


def _tcga_short_id(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text.upper().startswith("TCGA-") and len(text) >= 4:
        return text[-4:]
    try:
        return f"{int(float(text)):04d}"
    except ValueError:
        return None
