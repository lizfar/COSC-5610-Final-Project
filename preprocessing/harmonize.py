"""Shared clinical cleaning rules for FAB, cytogenetic risk, and the 2-year label."""

import re

import numpy as np
import pandas as pd

from config import DAYS_PER_YEAR, FAB_LEVELS, OS_LANDMARK_DAYS

UNKNOWN_TEXT = {
    "",
    "UNKNOWN",
    "NAN",
    "NA",
    "N/A",
    "NONE",
    "NULL",
    "[NOT AVAILABLE]",
    "[NOT APPLICABLE]",
    "[UNKNOWN]",
    "[DISCREPANCY]",
    "[COMPLETED]",
    "NOT CLASSIFIED",
    "NOTCLASSIFIED",
    "NOS",
    "MISSINGKARYO",
    "MISSINGMUTATIONS",
    "NONAML",
    "NONINITIAL",
    "FAVORABLEORINTERMEDIATE",
    "INTERMEDIATEORADVERSE",
}

FAB_ALIASES = {
    "M0 UNDIFFERENTIATED": "M0",
    "M0 UNDIFFERENTIATED AML": "M0",
    "M5A": "M5",
    "M5B": "M5",
    "M4EO": "M4",
    "M4 EO": "M4",
    "M4 EOS": "M4",
    "M6A": "M6",
    "M6B": "M6",
    "M0/M1": "UNKNOWN",
}

RISK_ALIASES = {
    "FAVORABLE": "FAVORABLE",
    "LOW": "FAVORABLE",
    "LOWRISK": "FAVORABLE",
    "INTERMEDIATE": "INTERMEDIATE",
    "INTERMEDIATE/NORMAL": "INTERMEDIATE",
    "NORMAL": "INTERMEDIATE",
    "STANDARD": "INTERMEDIATE",
    "STANDARDRISK": "INTERMEDIATE",
    "ADVERSE": "ADVERSE",
    "POOR": "ADVERSE",
    "HIGH": "ADVERSE",
    "HIGHRISK": "ADVERSE",
}

YES_TEXT = {"y", "yes", "true", "1", "positive"}
NO_TEXT = {"n", "no", "false", "0", "negative"}


def clean_id(value):
    if pd.isna(value):
        return None
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null"}:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def clean_symbol(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if "|" in text:
        parts = [part.strip() for part in text.split("|") if part.strip() and part.strip() != "?"]
        lettered = [
            part
            for part in parts
            if re.search(r"[A-Za-z]", part) and not part.upper().startswith("ENSG")
        ]
        if not lettered:
            return None
        text = lettered[0]
    text = text.upper()
    if text in UNKNOWN_TEXT or text in {".", "?", "-"} or text.startswith("N_") or text.startswith("ENSG"):
        return None
    return text


def map_fab(value):
    if pd.isna(value):
        return "UNKNOWN"
    text = re.sub(r"\s+", " ", str(value).strip().upper())
    text = FAB_ALIASES.get(text, text)
    if text in UNKNOWN_TEXT:
        return "UNKNOWN"
    if text in FAB_LEVELS:
        return text
    return "UNKNOWN"


def map_risk(value):
    if pd.isna(value):
        return "UNKNOWN"
    text = re.sub(r"\s+", "", str(value).strip().upper())
    unknown = {re.sub(r"\s+", "", item) for item in UNKNOWN_TEXT}
    if text in unknown or _is_number(text):
        return "UNKNOWN"
    mapped = RISK_ALIASES.get(text)
    if mapped is None:
        return "UNKNOWN"
    return mapped


def unexpected_risk_token(value):
    """True when a risk value became UNKNOWN and was not an explicit missing code."""
    if pd.isna(value):
        return False
    text = re.sub(r"\s+", "", str(value).strip().upper())
    unknown = {re.sub(r"\s+", "", item) for item in UNKNOWN_TEXT}
    if text in unknown or _is_number(text):
        return False
    return RISK_ALIASES.get(text) is None


def map_yes_no(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in YES_TEXT:
        return 1.0
    if text in NO_TEXT:
        return 0.0
    return np.nan


def map_event(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper()
    if text in {"DEAD", "DECEASED", "1", "1.0"}:
        return 1
    if text in {"ALIVE", "LIVING", "0", "0.0"}:
        return 0
    return np.nan


def to_number(series):
    if series is None:
        return None
    cleaned = series.astype(str).str.strip()
    cleaned = cleaned.replace(
        {
            "": np.nan,
            "nan": np.nan,
            "NaN": np.nan,
            "None": np.nan,
            "[Not Available]": np.nan,
            "[Not Applicable]": np.nan,
            "[Unknown]": np.nan,
            "[Discrepancy]": np.nan,
            "[Completed]": np.nan,
            "null": np.nan,
        }
    )
    return pd.to_numeric(cleaned, errors="coerce")


def add_two_year_label(frame):
    """Death by 730 days. Earlier censoring stays missing and is kept for survival models."""
    death_by_landmark = (frame["event"] == 1) & (frame["os_days"] <= OS_LANDMARK_DAYS)
    followed_through_landmark = frame["os_days"] >= OS_LANDMARK_DAYS
    label = pd.Series(np.nan, index=frame.index, dtype="float64")
    label.loc[followed_through_landmark] = 0
    label.loc[death_by_landmark] = 1
    frame["os_2yr_death"] = label
    return frame


def age_from_days(series):
    return to_number(series) / DAYS_PER_YEAR


def complexity_to_number(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text == "more than 6":
        return 7
    if text == "4 or more":
        return 4
    return pd.to_numeric(text, errors="coerce")


def lesion_column_name(label):
    cleaned = label.lower()
    cleaned = cleaned.replace("(", "").replace(")", "")
    cleaned = cleaned.replace(";", "_").replace(":", "_").replace(".", "_").replace(" ", "_")
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return f"lesion_{cleaned}"


def find_column(frame, candidates):
    lookup = {str(column).strip().lower(): column for column in frame.columns}
    for name in candidates:
        if name.lower() in lookup:
            return lookup[name.lower()]
    for name in candidates:
        needle = name.lower()
        for lowered, original in lookup.items():
            if needle in lowered:
                return original
    return None


def _is_number(text):
    try:
        float(text)
    except ValueError:
        return False
    return True
