#!/usr/bin/env python3
"""
Cross-check filled BioSample and SRA templates before upload.

Usage
-----
validate_submission.py BIOSAMPLE.xlsx SRA.xlsx [--sra-sheet SRA_data]
                       [--counts counts.tsv]

``counts.tsv`` has two tab-separated columns: file stem and number of reads.

Notes
-----
Prints a count line, then one "WARN" line per warning and one "ERROR" line
per error. The exit code is 1 if any error was found, otherwise 0.
"""

import argparse
import collections
import re
import sys

import openpyxl

MANDATORY_PREFIX = "*"
SAMPLE_NAME_FIELD = "sample_name"
FILENAME_COLUMN_PREFIX = "filename"
READ_STEM_REGEX = re.compile(r"_R?1(_001)?\.f(ast)?q(\.gz)?$")
AWKWARD_FILENAME_REGEX = re.compile(r"[()\s]")


def clean_field_name(cell_value: object) -> str | None:
    """
    Turn a header cell into a field name.

    Parameters
    ----------
    cell_value : object
        Value of one header cell.

    Returns
    -------
    field_name : str | None
        The cell text without the mandatory "*" prefix and surrounding
        whitespace, or None for an empty cell.
    """
    if not cell_value:
        field_name = None
        return field_name

    field_name = str(cell_value).lstrip(MANDATORY_PREFIX).strip()

    return field_name


def read_template_rows(
                       template_path: str,
                       sheet_name: str | None = None,
                       header_search_rows: int = 60
                       ) -> list[dict]:
    """
    Read the data rows of a filled BioSample or SRA template.

    Parameters
    ----------
    template_path : str
        Path to the .xlsx template.
    sheet_name : str | None
        Sheet to read; None means the first sheet.
    header_search_rows : int
        Number of rows, from the top, to search for the header row.

    Returns
    -------
    template_rows : list[dict]
        One dict per data row, keyed by field name (no "*" prefix), with
        every value as a stripped string. Rows with an empty first cell
        are skipped.

    Raises
    ------
    ValueError
        If no header row starting with "sample_name" is found.
    """
    workbook = openpyxl.load_workbook(
                                      template_path,
                                      data_only=True
                                      )
    worksheet = workbook[sheet_name] if sheet_name else workbook.worksheets[0]
    field_names = None
    first_data_row = None
    for row in worksheet.iter_rows(max_row=header_search_rows):
        first_value = row[0].value
        if not isinstance(first_value, str):
            continue
        if first_value.lstrip(MANDATORY_PREFIX).strip() == SAMPLE_NAME_FIELD:
            field_names = [
                           clean_field_name(cell.value)
                           for cell in row
                           ]
            first_data_row = row[0].row + 1
            break
    if field_names is None:
        message = f"header not found in {template_path}"
        raise ValueError(message)

    template_rows = []
    for values in worksheet.iter_rows(
                                      min_row=first_data_row,
                                      values_only=True
                                      ):
        if values[0] is None or values[0] == "":
            continue
        template_row = {}
        for position, value in enumerate(values):
            if not field_names[position]:
                continue
            template_row[field_names[position]] = "" if value is None else str(value).strip()
        template_rows.append(template_row)

    return template_rows


def find_duplicates(values: list[str]) -> list[str]:
    """
    List the values that occur more than once.

    Parameters
    ----------
    values : list[str]
        Values to check.

    Returns
    -------
    duplicates : list[str]
        Each repeated value once, in order of first appearance.
    """
    value_counts = collections.Counter(values)
    duplicates = []
    for value, count in value_counts.items():
        if count > 1:
            duplicates.append(value)

    return duplicates


def check_biosamples(
                     biosample_rows: list[dict],
                     ignored_columns: tuple[str, ...] = (
                                                         "sample_name",
                                                         "sample_title",
                                                         "description",
                                                         "bioproject_accession"
                                                         ),
                     missing_terms: tuple[str, ...] = (
                                                       "not applicable",
                                                       "not collected",
                                                       "missing",
                                                       "not provided",
                                                       "restricted access"
                                                       ),
                     date_pattern: str = r"^\d{4}(-\d{2}(-\d{2}(T\d{2}:\d{2}(:\d{2})?Z)?)?)?$",
                     lat_lon_pattern: str = r"^\d+(\.\d+)? [NS] \d+(\.\d+)? [EW]$"
                     ) -> tuple[list[str], list[str]]:
    """
    Check the BioSample rows for NCBI rejection causes.

    Parameters
    ----------
    biosample_rows : list[dict]
        Rows from `read_template_rows` for the BioSample template.
    ignored_columns : tuple[str, ...]
        Columns left out when testing whether two samples have identical
        attributes (NCBI rejects such samples).
    missing_terms : tuple[str, ...]
        Controlled "missing value" terms that are accepted in place of a
        date, coordinates or organism.
    date_pattern : str
        Regular expression an ISO 8601 collection_date must match.
    lat_lon_pattern : str
        Regular expression a lat_lon value must match.

    Returns
    -------
    errors : list[str]
        Problems NCBI would reject.
    warnings : list[str]
        Suspicious values that may be fine.
    """
    errors = []
    warnings = []
    sample_names = [
                    row["sample_name"]
                    for row in biosample_rows
                    ]
    for sample_name in find_duplicates(sample_names):
        errors.append(f"BioSample: duplicate sample_name {sample_name}")

    names_by_attributes = collections.defaultdict(list)
    for row in biosample_rows:
        attributes = []
        for column, value in sorted(row.items()):
            if column in ignored_columns or not value:
                continue
            attribute = (
                         column,
                         value
                         )
            attributes.append(attribute)
        names_by_attributes[tuple(attributes)].append(row["sample_name"])
    for sample_group in names_by_attributes.values():
        if len(sample_group) > 1:
            errors.append(f"BioSample: identical attributes (NCBI will reject): {sample_group}")

    date_regex = re.compile(date_pattern)
    lat_lon_regex = re.compile(lat_lon_pattern)
    for row in biosample_rows:
        sample_name = row["sample_name"]
        collection_date = row.get(
                                  "collection_date",
                                  ""
                                  )
        if collection_date and collection_date not in missing_terms:
            if not date_regex.match(collection_date):
                errors.append(
                              f"BioSample {sample_name}: bad collection_date '{collection_date}' "
                              "(ISO 8601, times must end in Z/UTC)"
                              )
        lat_lon = row.get(
                          "lat_lon",
                          ""
                          )
        if lat_lon and lat_lon not in missing_terms:
            if not lat_lon_regex.match(lat_lon):
                errors.append(f"BioSample {sample_name}: bad lat_lon '{lat_lon}'")
        organism = row.get(
                           "organism",
                           ""
                           )
        if organism and "metagenome" not in organism and organism not in missing_terms:
            warnings.append(
                            f"BioSample {sample_name}: organism '{organism}' has no "
                            "'metagenome' (required for MIMARKS/MIMS)"
                            )

    return errors, warnings


def find_filename_columns(sra_rows: list[dict]) -> list[str]:
    """
    List the file-name columns of the SRA template.

    Parameters
    ----------
    sra_rows : list[dict]
        Rows from `read_template_rows` for the SRA template.

    Returns
    -------
    filename_columns : list[str]
        Columns of the first row whose name starts with "filename"; empty
        if there are no rows.
    """
    filename_columns = []
    if not sra_rows:
        return filename_columns

    for column in sra_rows[0].keys():
        if column.startswith(FILENAME_COLUMN_PREFIX):
            filename_columns.append(column)

    return filename_columns


def list_run_files(sra_rows: list[dict]) -> list[str]:
    """
    List every file name in the SRA template.

    Parameters
    ----------
    sra_rows : list[dict]
        Rows from `read_template_rows` for the SRA template.

    Returns
    -------
    run_files : list[str]
        Non-empty values of all "filename*" columns, row by row.
    """
    filename_columns = find_filename_columns(sra_rows)
    run_files = []
    for row in sra_rows:
        for column in filename_columns:
            if row.get(column):
                run_files.append(row[column])

    return run_files


def check_sra_runs(
                   sra_rows: list[dict],
                   biosample_names: list[str],
                   required_columns: tuple[str, ...] = (
                                                        "library_strategy",
                                                        "library_source",
                                                        "library_selection",
                                                        "library_layout",
                                                        "platform",
                                                        "instrument_model",
                                                        "design_description",
                                                        "filetype"
                                                        )
                   ) -> tuple[list[str], list[str]]:
    """
    Check the SRA rows and their links to the BioSample rows.

    Parameters
    ----------
    sra_rows : list[dict]
        Rows from `read_template_rows` for the SRA template.
    biosample_names : list[str]
        The ``sample_name`` of every BioSample row.
    required_columns : tuple[str, ...]
        Columns that must not be empty when present in the template.

    Returns
    -------
    errors : list[str]
        Problems NCBI would reject.
    warnings : list[str]
        Suspicious values that may be fine.
    """
    errors = []
    warnings = []
    library_ids = [row.get(
                           "library_ID",
                           ""
                           ) for row in sra_rows]
    for library_id in find_duplicates(library_ids):
        errors.append(f"SRA: duplicate library_ID {library_id}")

    run_files = list_run_files(sra_rows)
    for run_file in find_duplicates(run_files):
        errors.append(f"SRA: file used more than once {run_file}")
    for run_file in run_files:
        if AWKWARD_FILENAME_REGEX.search(run_file):
            warnings.append(f"SRA: awkward characters in filename {run_file}")

    sra_samples = {
                   row["sample_name"]
                   for row in sra_rows
                   }
    biosample_set = set(biosample_names)
    for sample_name in sorted(sra_samples - biosample_set):
        errors.append(f"SRA sample_name not in BioSample sheet: {sample_name}")
    for sample_name in sorted(biosample_set - sra_samples):
        errors.append(f"BioSample with no SRA run: {sample_name}")

    for row in sra_rows:
        for column in required_columns:
            if column in row and not row[column]:
                errors.append(f"SRA {row.get('library_ID')}: empty {column}")

    return errors, warnings


def read_read_counts(counts_path: str) -> dict[str, int]:
    """
    Read a two-column table of file stems and read counts.

    Parameters
    ----------
    counts_path : str
        Path to the tab-separated counts file.

    Returns
    -------
    read_counts : dict[str, int]
        File stem -> number of reads. Lines whose last column is not a
        whole number (such as a header) are skipped.
    """
    read_counts = {}
    with open(counts_path) as counts_file:
        for line in counts_file:
            columns = line.strip().split("\t")
            if len(columns) >= 2 and columns[-1].isdigit():
                read_counts[columns[0]] = int(columns[-1])

    return read_counts


def check_read_counts(
                      sra_rows: list[dict],
                      read_counts: dict[str, int],
                      low_read_limit: int = 1000
                      ) -> tuple[list[str], list[str]]:
    """
    Flag SRA runs with no reads or very few reads.

    Parameters
    ----------
    sra_rows : list[dict]
        Rows from `read_template_rows` for the SRA template.
    read_counts : dict[str, int]
        File stem -> number of reads, from `read_read_counts`.
    low_read_limit : int
        Runs with fewer reads than this get a warning.

    Returns
    -------
    errors : list[str]
        Runs with 0 reads.
    warnings : list[str]
        Runs with no known read count, or fewer than `low_read_limit`.

    Notes
    -----
    The first "filename*" column of each row gives the file stem; the
    library_ID is the fallback key into `read_counts`.
    """
    errors = []
    warnings = []
    filename_columns = find_filename_columns(sra_rows)
    for row in sra_rows:
        run_file = row.get(
                           filename_columns[0],
                           ""
                           ) if filename_columns else ""
        stem = READ_STEM_REGEX.sub(
                                   "",
                                   run_file
                                   )
        count_by_library = read_counts.get(row.get("library_ID"))
        read_count = read_counts.get(
                                     stem,
                                     count_by_library
                                     )
        if read_count is None:
            warnings.append(f"no read count for {run_file}")
        elif read_count == 0:
            errors.append(
                          f"SRA {run_file}: 0 reads — remove run "
                          "(and BioSample if it has no other runs)"
                          )
        elif read_count < low_read_limit:
            warnings.append(f"SRA {run_file}: only {read_count} reads")

    return errors, warnings


def parse_arguments(argv: list[str] | None = None) -> argparse.Namespace:
    """
    Parse the command line.

    Parameters
    ----------
    argv : list[str] | None
        Arguments to parse; None means ``sys.argv[1:]``.

    Returns
    -------
    args : argparse.Namespace
        Parsed arguments: ``biosample``, ``sra``, ``sra_sheet`` and
        ``counts``.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("biosample")
    parser.add_argument("sra")
    parser.add_argument(
                        "--sra-sheet",
                        default="SRA_data"
                        )
    parser.add_argument("--counts")
    args = parser.parse_args(argv)

    return args


def main(argv: list[str] | None = None) -> int:
    """
    Validate a filled BioSample template and SRA template together.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments; None means ``sys.argv[1:]``.

    Returns
    -------
    exit_code : int
        1 if any error was found or a template has no header row,
        otherwise 0.

    Notes
    -----
    Writes the count line, warnings and errors to standard output.
    """
    args = parse_arguments(argv)
    try:
        biosample_rows = read_template_rows(args.biosample)
        sra_rows = read_template_rows(
                                      args.sra,
                                      args.sra_sheet
                                      )
    except ValueError as error:
        print(
              error,
              file=sys.stderr
              )
        exit_code = 1
        return exit_code

    biosample_names = [
                       row["sample_name"]
                       for row in biosample_rows
                       ]
    biosample_errors, biosample_warnings = check_biosamples(biosample_rows)
    sra_errors, sra_warnings = check_sra_runs(
                                              sra_rows,
                                              biosample_names
                                              )
    errors = biosample_errors + sra_errors
    warnings = biosample_warnings + sra_warnings
    if args.counts:
        read_counts = read_read_counts(args.counts)
        count_errors, count_warnings = check_read_counts(
                                                         sra_rows,
                                                         read_counts
                                                         )
        errors += count_errors
        warnings += count_warnings

    run_files = list_run_files(sra_rows)
    print(
          f"BioSamples: {len(biosample_rows)}  SRA runs: {len(sra_rows)}  "
          f"files: {len(run_files)}"
          )
    for warning in warnings:
        print(f"WARN  {warning}")
    for error in errors:
        print(f"ERROR {error}")
    print("OK" if not errors else f"{len(errors)} error(s)")
    exit_code = 1 if errors else 0

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
