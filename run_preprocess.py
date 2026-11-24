"""Download the open AML cohorts and build the patient-level modeling tables."""

import traceback

from preprocessing.assemble import build_reference_cohorts, build_target_cohort, write_full_report
from preprocessing.download import download_reference_files, download_target_expression
from preprocessing.logutil import fail, ok, step


def main():
    try:
        step("Starting AML prognosis preprocessing.")
        download_reference_files()
        sections = build_reference_cohorts()
        manifest = download_target_expression()
        sections.extend(build_target_cohort(manifest))
        write_full_report(sections)
        ok("Preprocessing finished.")
    except Exception as exc:
        fail(f"Preprocessing stopped: {exc}")
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
