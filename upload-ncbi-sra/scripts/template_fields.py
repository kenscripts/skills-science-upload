#!/usr/bin/env python3
"""
List the fields of an NCBI BioSample or SRA template (.xlsx or .tsv/.txt).

Finds the header row (the first row whose first cell is '*sample_name' or
'sample_name'), then prints each field, whether it is mandatory (a '*'
prefix), and, for .xlsx templates, the first line of its cell comment.

Usage
-----
template_fields.py TEMPLATE [--sheet NAME]

Notes
-----
``--sheet`` only applies to .xlsx templates; the first sheet is used when it
is omitted. A .tsv/.txt template prints one "MANDATORY"/"optional" line per
field. An .xlsx template prints a "# sheet=... header_row=..." line, then one
tab-separated line per field: status, column letter, field name, comment.
"""

import argparse
import csv
import sys

import openpyxl

TSV_SUFFIXES = (
                ".tsv",
                ".txt"
                )
SAMPLE_NAME_FIELD = "sample_name"
MANDATORY_PREFIX = "*"


def find_header_row(
                    worksheet: "openpyxl.worksheet.worksheet.Worksheet",
                    header_search_rows: int = 60
                    ) -> tuple | None:
    """
    Find the header row of a template sheet.

    Parameters
    ----------
    worksheet : openpyxl.worksheet.worksheet.Worksheet
        Sheet to search.
    header_search_rows : int
        Number of rows, from the top, to search.

    Returns
    -------
    header_row : tuple | None
        The cells of the first row whose first cell is "sample_name" (with
        or without a leading "*"), or None if there is no such row.
    """
    header_row = None
    for row in worksheet.iter_rows(max_row=header_search_rows):
        first_value = row[0].value
        if not isinstance(first_value, str):
            continue
        if first_value.lstrip(MANDATORY_PREFIX).strip() == SAMPLE_NAME_FIELD:
            header_row = row
            break

    return header_row


def describe_tsv_fields(template_path: str) -> list[str]:
    """
    Describe the fields of a tab-separated template.

    Parameters
    ----------
    template_path : str
        Path to the .tsv/.txt template.

    Returns
    -------
    field_lines : list[str]
        One "MANDATORY <field>" or "optional  <field>" line per field of
        the header row.

    Raises
    ------
    StopIteration
        If no row starts with "sample_name".
    """
    with open(template_path) as template_file:
        reader = csv.reader(
                            template_file,
                            delimiter="\t"
                            )
        rows = list(reader)

    header_rows = (
                   row
                   for row in rows
                   if row and row[0].lstrip(MANDATORY_PREFIX) == SAMPLE_NAME_FIELD
                   )
    header = next(header_rows)
    field_lines = []
    for field in header:
        status = "MANDATORY " if field.startswith(MANDATORY_PREFIX) else "optional  "
        field_lines.append(status + field.lstrip(MANDATORY_PREFIX))

    return field_lines


def describe_xlsx_fields(
                         template_path: str,
                         sheet_name: str | None = None,
                         max_comment_chars: int = 120
                         ) -> list[str]:
    """
    Describe the fields of an Excel template.

    Parameters
    ----------
    template_path : str
        Path to the .xlsx template.
    sheet_name : str | None
        Sheet to read; None means the first sheet.
    max_comment_chars : int
        Longest comment excerpt to show per field.

    Returns
    -------
    field_lines : list[str]
        A "# sheet=... header_row=... first_data_row=..." line, then one
        tab-separated line per field with a header cell: status, column
        letter, field name, first line of its cell comment.

    Raises
    ------
    ValueError
        If the sheet has no header row.
    """
    workbook = openpyxl.load_workbook(template_path)
    worksheet = workbook[sheet_name] if sheet_name else workbook.worksheets[0]
    header_row = find_header_row(worksheet)
    if header_row is None:
        message = (
                   f"no header row with sample_name found in sheet {worksheet.title}; "
                   f"sheets: {workbook.sheetnames}"
                   )
        raise ValueError(message)

    header_number = header_row[0].row
    first_data_row = header_number + 1
    sheet_line = (
                  f"# sheet={worksheet.title} header_row={header_number} "
                  f"first_data_row={first_data_row}"
                  )
    field_lines = [sheet_line]
    for cell in header_row:
        if cell.value is None:
            continue
        field = str(cell.value)
        comment = ""
        if cell.comment:
            comment = cell.comment.text.strip().splitlines()[0][:max_comment_chars]
        status = "MANDATORY" if field.startswith(MANDATORY_PREFIX) else "optional "
        field_name = field.lstrip(MANDATORY_PREFIX)
        field_lines.append(f"{status}\t{cell.column_letter}\t{field_name}\t{comment}")

    return field_lines


def main(argv: list[str] | None = None) -> int:
    """
    Print the fields of a template given on the command line.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments; None means ``sys.argv[1:]``.

    Returns
    -------
    exit_code : int
        0 on success, 1 when an .xlsx template has no header row.

    Notes
    -----
    Writes the field list to standard output and any error to standard
    error.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("template")
    parser.add_argument("--sheet")
    args = parser.parse_args(argv)
    if args.template.endswith(TSV_SUFFIXES):
        field_lines = describe_tsv_fields(args.template)
    else:
        try:
            field_lines = describe_xlsx_fields(
                                               args.template,
                                               args.sheet
                                               )
        except ValueError as error:
            print(
                  error,
                  file=sys.stderr
                  )
            exit_code = 1
            return exit_code

    for field_line in field_lines:
        print(field_line)
    exit_code = 0

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
