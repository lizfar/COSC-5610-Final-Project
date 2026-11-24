"""Constants for the AML prognosis project. No user input is read here."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
REPORT_PATH = PROCESSED_DIR / "preprocessing_report.txt"

USER_AGENT = "aml-prognosis-preprocessing"
OS_LANDMARK_DAYS = 730
DAYS_PER_YEAR = 365.25
CPM_SCALE = 1_000_000

MIN_BEAT_PATIENTS = 300
MIN_TCGA_PATIENTS = 100
MIN_TARGET_PATIENTS = 100

DOWNLOAD_WORKERS = 4
PARSE_WORKERS = 4

BEAT_ZENODO_URL = (
    "https://zenodo.org/api/records/10654808/files/"
    "biodev/beataml2.0_data-v3.0.zip/content"
)
BEAT_ZIP_BYTES = 160932018
BEAT_ZIP_PATH = RAW_DIR / "beataml2.0_data-v3.0.zip"
BEAT_DIR = RAW_DIR / "beat"
BEAT_CLINICAL_NAME = "beataml_wv1to4_clinical.xlsx"
BEAT_COUNTS_NAME = "beataml_waves1to4_counts_dbgap.txt"
BEAT_MAPPING_NAME = "beataml_waves1to4_sample_mapping.xlsx"
BEAT_CLINICAL_PATH = BEAT_DIR / BEAT_CLINICAL_NAME
BEAT_COUNTS_PATH = BEAT_DIR / BEAT_COUNTS_NAME
BEAT_MAPPING_PATH = BEAT_DIR / BEAT_MAPPING_NAME

TARGET_DISCOVERY_URL = (
    "https://d2l5jy2ao2mx5b.cloudfront.net/target/phs000465/"
    "TARGET_AML_ClinicalData_Discovery_20211201.xlsx"
)
TARGET_VALIDATION_URL = (
    "https://d2l5jy2ao2mx5b.cloudfront.net/target/phs000465/"
    "TARGET_AML_ClinicalData_Validation_20211201.xlsx"
)
TARGET_DISCOVERY_PATH = RAW_DIR / "TARGET_AML_ClinicalData_Discovery_20211201.xlsx"
TARGET_VALIDATION_PATH = RAW_DIR / "TARGET_AML_ClinicalData_Validation_20211201.xlsx"
TARGET_STAR_DIR = RAW_DIR / "target_star"
TARGET_MANIFEST_PATH = RAW_DIR / "target_star_manifest.csv"

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data/{file_id}"
TCGA_CLINICAL_FILE_ID = "1fd61a91-90d5-405a-bcef-b9df8b396a21"
TCGA_COUNTS_FILE_ID = "f60dc63f-4b4f-453d-b3a2-6e327ce7033c"
TCGA_SUPPLEMENT_FILE_ID = "b9196563-a05d-40b8-80dc-640ec712eb06"
TCGA_CLINICAL_PATH = RAW_DIR / "clinical_patient_laml.tsv"
TCGA_SUPPLEMENT_PATH = RAW_DIR / "laml_supplemental_table_01.xlsx"
TCGA_STAR_DIR = RAW_DIR / "tcga_star"
TCGA_MANIFEST_PATH = RAW_DIR / "tcga_star_manifest.csv"

PRIMARY_SAMPLE_PRIORITY = {
    "primary blood derived cancer - bone marrow": 0,
    "primary blood derived cancer - peripheral blood": 1,
}

LITERATURE_GENES = [
    "CD109",
    "CPNE3",
    "CPNE8",
    "INPP4B",
    "SPINK2",
    "PLXNC1",
    "LSP1",
    "DDIT4",
    "SLC40A1",
    "CD79A",
    "CD19",
]

TARGET_LESION_COLUMNS = [
    "t(6;9)",
    "t(8;21)",
    "t(3;5)(q25;q34)",
    "t(6;11)(q27;q23)",
    "t(9;11)(p22;q23)",
    "t(10;11)(p11.2;q23)",
    "t(11:19)(q23:p13.1)",
    "inv(16)",
    "del5q",
    "del7q",
    "del9q",
    "monosomy 5",
    "monosomy 7",
    "trisomy 8",
    "trisomy 21",
    "MLL",
    "Minus Y",
    "Minus X",
]

FAB_LEVELS = ["M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7"]
RISK_LEVELS = ["FAVORABLE", "INTERMEDIATE", "ADVERSE", "UNKNOWN"]
