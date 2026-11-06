#!/usr/bin/env python3
"""
Decide whether paired FASTQs must be demultiplexed before SRA submission.

SRA needs one BioSample per run, so a provider file that pools several
libraries has to be split first. This script inspects a folder of paired
FASTQs and labels each pair.

Usage
-----
check_libraries.py FASTQ_DIR MAP.tsv [--n 20000] [--out report.tsv]

MAP.tsv (tab-separated, with header)::

    pool library_id out_name fwd_name fwd_seq rev_name rev_seq

``pool`` is the provider tube/pool/file stem the library was sequenced in,
``out_name`` is the SRA-safe per-library file stem (no spaces/parentheses),
and ``*_seq`` is the inline index at the 5' end of R1 (fwd) and R2 (rev).

Notes
-----
Evidence per file pair:

* name: the stem matches a pool (with >1 library) or a single library.
* content: the first N read pairs are tallied by inline index combination
  (<=1 mismatch); >=2 expected library combinations, each >=1% of the
  reads, means the file is pooled.
* orient: share of reads whose R1 starts with a reverse index (mixed
  orientation means demultiplex with BOTH_ORIENT=1).

The exit code is 0 unless the map or the FASTQ folder is unusable. The
``decision`` column is one of DEMUX, OK_PER_LIBRARY or CHECK.
"""

import argparse
import collections
import collections.abc
import csv
import glob
import gzip
import os
import re
import sys
import typing

COLUMN_ALIASES = {
                  "tube": "pool",
                  "fwd": "fwd_name",
                  "rev": "rev_name"
                  }
REQUIRED_MAP_COLUMNS = {
                        "pool",
                        "library_id",
                        "out_name",
                        "fwd_name",
                        "fwd_seq",
                        "rev_name",
                        "rev_seq"
                        }
REPORT_COLUMNS = [
                  "file",
                  "decision",
                  "reason",
                  "name_evidence",
                  "reads_checked",
                  "frac_index_matched",
                  "frac_reverse_orient",
                  "libraries_detected"
                  ]
R1_NAME_REGEX = re.compile(r"(.+?)[._](R?)1(_001)?\.f(ast)?q(\.gz)?$")
R1_TO_R2_REGEX = re.compile(r"([._]R?)1((_001)?\.f(ast)?q(\.gz)?)$")
R1_TO_R2_REPLACEMENT = r"\g<1>2\2"
LINES_PER_READ = 4
SEQUENCE_LINE_INDEX = 1
REVERSE_ORIENT_WARN_FRACTION = 0.05


def print_stderr(message: str) -> None:
    """
    Print one line to standard error.

    Parameters
    ----------
    message : str
        Text to print.
    """
    print(
          message,
          file=sys.stderr
          )


def normalize_map_row(raw_row: dict) -> dict:
    """
    Strip whitespace from one MAP.tsv row and rename aliased columns.

    Parameters
    ----------
    raw_row : dict
        One row as read by ``csv.DictReader``.

    Returns
    -------
    normalized_row : dict
        Same values, stripped, keyed by canonical column name ("tube"
        becomes "pool", "fwd" becomes "fwd_name", "rev" becomes
        "rev_name"). Empty cells become empty strings.
    """
    normalized_row = {}
    for column, value in raw_row.items():
        column_name = column.strip()
        canonical_name = COLUMN_ALIASES.get(
                                            column_name,
                                            column_name
                                            )
        normalized_row[canonical_name] = (value or "").strip()

    return normalized_row


def read_map_rows(map_path: str) -> list[dict]:
    """
    Read MAP.tsv into one normalized dict per library.

    Parameters
    ----------
    map_path : str
        Path to the tab-separated map file.

    Returns
    -------
    map_rows : list[dict]
        One dict per library, keyed by canonical column name.

    Raises
    ------
    ValueError
        If a required column is missing from the header.
    IndexError
        If the map has a header but no library rows.
    """
    with open(map_path) as map_file:
        reader = csv.DictReader(
                                map_file,
                                delimiter="\t"
                                )
        raw_rows = list(reader)

    map_rows = [
                normalize_map_row(raw_row)
                for raw_row in raw_rows
                ]
    if not REQUIRED_MAP_COLUMNS <= set(map_rows[0]):
        message = f"map needs columns {sorted(REQUIRED_MAP_COLUMNS)}"
        raise ValueError(message)

    return map_rows


def build_index_lookups(map_rows: list[dict]) -> tuple[dict, dict, dict]:
    """
    Build the index-sequence lookups used to read inline indexes.

    Parameters
    ----------
    map_rows : list[dict]
        Normalized rows from `read_map_rows`.

    Returns
    -------
    forward_indexes : dict
        Upper-case forward index sequence -> forward index name.
    reverse_indexes : dict
        Upper-case reverse index sequence -> reverse index name.
    index_pair_to_library : dict
        (forward index name, reverse index name) -> library ``out_name``.
    """
    forward_indexes = {
                       row["fwd_seq"].upper(): row["fwd_name"]
                       for row in map_rows
                       }
    reverse_indexes = {
                       row["rev_seq"].upper(): row["rev_name"]
                       for row in map_rows
                       }
    index_pair_to_library = {}
    for row in map_rows:
        index_pair = (
                      row["fwd_name"],
                      row["rev_name"]
                      )
        index_pair_to_library[index_pair] = row["out_name"]

    return forward_indexes, reverse_indexes, index_pair_to_library


def group_libraries_by_pool(map_rows: list[dict]) -> dict[str, list[dict]]:
    """
    Group the map rows by the pool (provider file stem) they were run in.

    Parameters
    ----------
    map_rows : list[dict]
        Normalized rows from `read_map_rows`.

    Returns
    -------
    pool_libraries : dict[str, list[dict]]
        Pool name -> the map rows of the libraries in that pool, in map
        order.
    """
    pool_libraries = collections.defaultdict(list)
    for row in map_rows:
        pool_libraries[row["pool"]].append(row)

    return pool_libraries


def index_matches_read(
                       read_start: str,
                       index_seq: str,
                       max_mismatch: int = 1
                       ) -> bool:
    """
    Test whether a read begins with an index, allowing a few mismatches.

    Parameters
    ----------
    read_start : str
        Read sequence (only its first ``len(index_seq)`` bases are used).
    index_seq : str
        Index sequence expected at the 5' end of the read.
    max_mismatch : int
        Largest number of mismatching bases still counted as a match.

    Returns
    -------
    is_match : bool
        True if the read is long enough and has at most `max_mismatch`
        mismatches against the index.
    """
    index_len = len(index_seq)
    read_prefix = read_start[:index_len]
    base_pairs = zip(
                     read_prefix,
                     index_seq
                     )
    mismatches = sum(
                     read_base != index_base
                     for read_base, index_base in base_pairs
                     )
    is_match = mismatches <= max_mismatch and len(read_start) >= index_len

    return is_match


def match_index_name(read_start: str, index_names: dict) -> str | None:
    """
    Name the one index in a lookup that a read starts with.

    Parameters
    ----------
    read_start : str
        Read sequence.
    index_names : dict
        Index sequence -> index name.

    Returns
    -------
    index_name : str | None
        The index name if exactly one index matches, otherwise None
        (no match, or an ambiguous match).
    """
    hits = []
    for index_seq, index_name in index_names.items():
        if index_matches_read(read_start, index_seq):
            hits.append(index_name)
    index_name = hits[0] if len(hits) == 1 else None

    return index_name


def open_fastq(fastq_path: str) -> typing.TextIO:
    """
    Open a FASTQ file for text reading, decompressing ``.gz`` files.

    Parameters
    ----------
    fastq_path : str
        Path to a plain or gzip-compressed FASTQ file.

    Returns
    -------
    fastq_file : typing.TextIO
        Open text-mode file handle. The caller closes it.
    """
    if fastq_path.endswith(".gz"):
        fastq_file = gzip.open(
                               fastq_path,
                               "rt"
                               )
        return fastq_file

    fastq_file = open(fastq_path)

    return fastq_file


def iter_read_sequences(
                        fastq_path: str,
                        max_reads: int
                        ) -> collections.abc.Iterator[str]:
    """
    Yield the upper-case sequences of the first reads of a FASTQ file.

    Parameters
    ----------
    fastq_path : str
        Path to a plain or gzip-compressed FASTQ file.
    max_reads : int
        Number of reads to read before stopping.

    Yields
    ------
    sequence : str
        Upper-case sequence line of one read.
    """
    with open_fastq(fastq_path) as fastq_file:
        for line_number, line in enumerate(fastq_file):
            if line_number // LINES_PER_READ >= max_reads:
                break
            if line_number % LINES_PER_READ == SEQUENCE_LINE_INDEX:
                sequence = line.strip().upper()
                yield sequence


def find_fastq_pairs(fastq_dir: str) -> list[tuple[str, str, str | None]]:
    """
    Find the R1 FASTQ files in a folder and their R2 mates.

    Parameters
    ----------
    fastq_dir : str
        Folder holding the provider FASTQ files.

    Returns
    -------
    pairs : list[tuple[str, str, str | None]]
        One (file stem, R1 path, R2 path) tuple per R1 file, sorted by
        path. The R2 path is None when no mate file exists.
    """
    fastq_pattern = os.path.join(
                                 fastq_dir,
                                 "*"
                                 )
    fastq_paths = glob.glob(fastq_pattern)
    pairs = []
    for r1_path in sorted(fastq_paths):
        r1_name = os.path.basename(r1_path)
        name_match = R1_NAME_REGEX.match(r1_name)
        if not name_match:
            continue
        r2_name = R1_TO_R2_REGEX.sub(
                                     R1_TO_R2_REPLACEMENT,
                                     r1_name
                                     )
        r2_path = os.path.join(
                               os.path.dirname(r1_path),
                               r2_name
                               )
        has_mate = os.path.exists(r2_path) and r2_path != r1_path
        mate_path = r2_path if has_mate else None
        pair = (
                name_match.group(1),
                r1_path,
                mate_path
                )
        pairs.append(pair)

    return pairs


def describe_name_evidence(
                           stem: str,
                           pool_libraries: dict[str, list[dict]],
                           library_names: set[str]
                           ) -> str:
    """
    Say what a file stem is, according to the map.

    Parameters
    ----------
    stem : str
        File name without the R1/R2 suffix.
    pool_libraries : dict[str, list[dict]]
        Pool name -> map rows, from `group_libraries_by_pool`.
    library_names : set[str]
        All ``out_name`` and ``library_id`` values in the map.

    Returns
    -------
    name_evidence : str
        "pool with N libraries", "pool with 1 library",
        "matches a library" or "unknown name".
    """
    if stem in pool_libraries:
        pool_size = len(pool_libraries[stem])
        name_evidence = f"pool with {pool_size} libraries"
        if pool_size == 1:
            name_evidence = "pool with 1 library"
        return name_evidence

    if stem in library_names:
        name_evidence = "matches a library"
        return name_evidence

    name_evidence = "unknown name"

    return name_evidence


def tally_inline_indexes(
                         r1_path: str,
                         r2_path: str,
                         forward_indexes: dict,
                         reverse_indexes: dict,
                         max_reads: int
                         ) -> tuple[collections.Counter, int, int, int]:
    """
    Count the inline index combinations in the first read pairs of a file.

    Parameters
    ----------
    r1_path : str
        Path to the R1 FASTQ file.
    r2_path : str
        Path to the R2 FASTQ file.
    forward_indexes : dict
        Forward index sequence -> name.
    reverse_indexes : dict
        Reverse index sequence -> name.
    max_reads : int
        Number of read pairs to inspect.

    Returns
    -------
    index_pair_counts : collections.Counter
        (forward name, reverse name) -> number of read pairs.
    reads_checked : int
        Number of read pairs inspected.
    reads_matched : int
        Read pairs with a forward index on R1 and a reverse index on R2.
    reads_reverse : int
        Read pairs with the indexes the other way round (reverse
        orientation), among those not already matched.
    """
    index_pair_counts = collections.Counter()
    reads_checked = 0
    reads_matched = 0
    reads_reverse = 0
    r1_reads = iter_read_sequences(
                                   r1_path,
                                   max_reads
                                   )
    r2_reads = iter_read_sequences(
                                   r2_path,
                                   max_reads
                                   )
    for r1_seq, r2_seq in zip(
                              r1_reads,
                              r2_reads
                              ):
        reads_checked += 1
        forward_name = match_index_name(
                                        r1_seq,
                                        forward_indexes
                                        )
        reverse_name = match_index_name(
                                        r2_seq,
                                        reverse_indexes
                                        )
        if forward_name and reverse_name:
            index_pair = (
                          forward_name,
                          reverse_name
                          )
            index_pair_counts[index_pair] += 1
            reads_matched += 1
            continue
        swapped_forward = match_index_name(
                                           r1_seq,
                                           reverse_indexes
                                           )
        swapped_reverse = match_index_name(
                                           r2_seq,
                                           forward_indexes
                                           )
        if swapped_forward and swapped_reverse:
            reads_reverse += 1

    return index_pair_counts, reads_checked, reads_matched, reads_reverse


def decide_pair_status(
                       name_evidence: str,
                       pool_size: int,
                       libraries_seen: list[str],
                       matched_fraction: float,
                       reverse_fraction: float,
                       pooled_min_fraction: float = 0.1,
                       single_min_fraction: float = 0.5,
                       unindexed_max_fraction: float = 0.1
                       ) -> tuple[str, str]:
    """
    Decide DEMUX, OK_PER_LIBRARY or CHECK for one FASTQ pair.

    Parameters
    ----------
    name_evidence : str
        Output of `describe_name_evidence`.
    pool_size : int
        Number of libraries in the pool the file stem names; 0 if the stem
        is not a pool in the map.
    libraries_seen : list[str]
        Libraries whose index combination is common in the reads.
    matched_fraction : float
        Share of inspected read pairs carrying a known index pair.
    reverse_fraction : float
        Share of inspected read pairs in reverse orientation.
    pooled_min_fraction : float
        A pool with several libraries is called DEMUX on file-name evidence
        alone when more than this share of reads carry known indexes.
    single_min_fraction : float
        A single detected library is trusted when more than this share of
        reads carry known indexes.
    unindexed_max_fraction : float
        A file named after a library is called OK_PER_LIBRARY without
        inline indexes when fewer than this share of reads carry them.

    Returns
    -------
    decision : str
        "DEMUX", "OK_PER_LIBRARY" or "CHECK".
    reason : str
        One-line explanation shown in the report.
    """
    matches_library = name_evidence == "matches a library"
    pooled_by_name = pool_size > 1 and matched_fraction > pooled_min_fraction
    if len(libraries_seen) >= 2 or pooled_by_name:
        decision = "DEMUX"
        reason = f"{len(libraries_seen)} libraries detected by inline index"
    elif len(libraries_seen) == 1 and matched_fraction > single_min_fraction:
        decision = "OK_PER_LIBRARY" if matches_library else "CHECK"
        reason = f"single library {libraries_seen[0]}"
        if not matches_library:
            reason += " — rename file to its out_name"
    elif matched_fraction < unindexed_max_fraction and matches_library:
        decision = "OK_PER_LIBRARY"
        reason = (
                  "name matches a library; no inline indexes found "
                  "(already trimmed or index-read design?)"
                  )
    else:
        decision = "CHECK"
        reason = "inconclusive — inspect reads / ask user"

    if decision == "DEMUX" and pool_size == 0:
        decision = "CHECK"
        reason += (
                   "; file name is not a pool in MAP — set its pool name "
                   "in MAP so demux_for_sra.sh can find it"
                   )
    if reverse_fraction > REVERSE_ORIENT_WARN_FRACTION and decision == "DEMUX":
        reason += f"; {reverse_fraction:.0%} reverse-orientation reads -> use BOTH_ORIENT=1"

    return decision, reason


def classify_fastq_pair(
                        pair: tuple[str, str, str | None],
                        pool_libraries: dict[str, list[dict]],
                        library_names: set[str],
                        forward_indexes: dict,
                        reverse_indexes: dict,
                        index_pair_to_library: dict,
                        max_reads: int,
                        min_combo_fraction: float = 0.01
                        ) -> dict:
    """
    Build the report row for one FASTQ pair.

    Parameters
    ----------
    pair : tuple[str, str, str | None]
        (file stem, R1 path, R2 path) from `find_fastq_pairs`.
    pool_libraries : dict[str, list[dict]]
        Pool name -> map rows.
    library_names : set[str]
        All ``out_name`` and ``library_id`` values in the map.
    forward_indexes : dict
        Forward index sequence -> name.
    reverse_indexes : dict
        Reverse index sequence -> name.
    index_pair_to_library : dict
        (forward name, reverse name) -> library ``out_name``.
    max_reads : int
        Number of read pairs to inspect.
    min_combo_fraction : float
        An index combination counts as a library in the file when it
        covers at least this share of the inspected reads (and at least
        one read).

    Returns
    -------
    report_row : dict
        Keys from REPORT_COLUMNS. A pair with no R2 mate only has "file",
        "decision" and "reason".
    """
    stem, r1_path, r2_path = pair
    if r2_path is None:
        report_row = {
                      "file": stem,
                      "decision": "CHECK",
                      "reason": "no R2 mate found"
                      }
        return report_row

    name_evidence = describe_name_evidence(
                                           stem,
                                           pool_libraries,
                                           library_names
                                           )
    stem_libraries = pool_libraries.get(
                                        stem,
                                        []
                                        )
    pool_size = len(stem_libraries)
    tallies = tally_inline_indexes(
                                   r1_path,
                                   r2_path,
                                   forward_indexes,
                                   reverse_indexes,
                                   max_reads
                                   )
    index_pair_counts, reads_checked, reads_matched, reads_reverse = tallies
    min_reads = max(
                    1,
                    min_combo_fraction * reads_checked
                    )
    libraries_seen = sorted(
                            index_pair_to_library[index_pair]
                            for index_pair, count in index_pair_counts.items()
                            if index_pair in index_pair_to_library and count >= min_reads
                            )
    matched_fraction = reads_matched / reads_checked if reads_checked else 0
    reverse_fraction = reads_reverse / reads_checked if reads_checked else 0
    decision, reason = decide_pair_status(
                                          name_evidence,
                                          pool_size,
                                          libraries_seen,
                                          matched_fraction,
                                          reverse_fraction
                                          )
    report_row = {
                  "file": stem,
                  "decision": decision,
                  "reason": reason,
                  "name_evidence": name_evidence,
                  "reads_checked": reads_checked,
                  "frac_index_matched": round(
                                              matched_fraction,
                                              3
                                              ),
                  "frac_reverse_orient": round(
                                               reverse_fraction,
                                               3
                                               ),
                  "libraries_detected": ";".join(libraries_seen)
                  }

    return report_row


def write_report_rows(
                      report_rows: list[dict],
                      report_file: typing.TextIO
                      ) -> None:
    """
    Write the report as a tab-separated table.

    Parameters
    ----------
    report_rows : list[dict]
        Rows from `classify_fastq_pair`.
    report_file : typing.TextIO
        Open text file (or stdout) to write to.
    """
    writer = csv.DictWriter(
                            report_file,
                            REPORT_COLUMNS,
                            delimiter="\t",
                            extrasaction="ignore",
                            lineterminator="\n"
                            )
    writer.writeheader()
    for report_row in report_rows:
        writer.writerow(report_row)


def find_absent_pools(
                      pool_libraries: dict[str, list[dict]],
                      stems: set[str]
                      ) -> list[str]:
    """
    List the map pools that have no FASTQ file in the folder.

    Parameters
    ----------
    pool_libraries : dict[str, list[dict]]
        Pool name -> map rows.
    stems : set[str]
        File stems found in the FASTQ folder.

    Returns
    -------
    absent_pools : list[str]
        Sorted names of pools with neither a file named after the pool nor
        a file for every one of their libraries.
    """
    absent_pools = []
    for pool_name, libraries in pool_libraries.items():
        if pool_name in stems:
            continue
        library_files_found = [
                               library["out_name"] in stems or library["library_id"] in stems
                               for library in libraries
                               ]
        if all(library_files_found):
            continue
        absent_pools.append(pool_name)

    return sorted(absent_pools)


def print_summary(
                  report_rows: list[dict],
                  pairs: list[tuple[str, str, str | None]],
                  pool_libraries: dict[str, list[dict]],
                  fastq_dir: str,
                  map_path: str
                  ) -> None:
    """
    Print the decision counts and next-step hints to standard error.

    Parameters
    ----------
    report_rows : list[dict]
        Rows from `classify_fastq_pair`.
    pairs : list[tuple[str, str, str | None]]
        Pairs from `find_fastq_pairs`.
    pool_libraries : dict[str, list[dict]]
        Pool name -> map rows.
    fastq_dir : str
        FASTQ folder, echoed in the suggested demultiplex command.
    map_path : str
        MAP.tsv path, echoed in the suggested demultiplex command.
    """
    decision_counts = collections.Counter(
                                          report_row["decision"]
                                          for report_row in report_rows
                                          )
    print_stderr(f"\n# {dict(decision_counts)}")
    if decision_counts["DEMUX"]:
        demux_rows = [
                      report_row
                      for report_row in report_rows
                      if report_row["decision"] == "DEMUX"
                      ]
        any_reverse = False
        for demux_row in demux_rows:
            if demux_row.get("frac_reverse_orient", 0) > REVERSE_ORIENT_WARN_FRACTION:
                any_reverse = True
        both_orient = "BOTH_ORIENT=1 " if any_reverse else ""
        command = f"{both_orient}bash demux_for_sra.sh {fastq_dir} <OUT_DIR> {map_path}"
        print_stderr(f"# demultiplex: {command}")

    stems = {
             stem
             for stem, _, _ in pairs
             }
    absent_pools = find_absent_pools(
                                     pool_libraries,
                                     stems
                                     )
    if absent_pools:
        print_stderr(f"# pools in map with no FASTQs here: {absent_pools}")


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
        Parsed arguments: ``fastq_dir``, ``map``, ``max_reads`` and
        ``out``.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("fastq_dir")
    parser.add_argument("map")
    parser.add_argument(
                        "--n",
                        dest="max_reads",
                        type=int,
                        default=20000
                        )
    parser.add_argument("--out")
    args = parser.parse_args(argv)

    return args


def main(argv: list[str] | None = None) -> int:
    """
    Check every FASTQ pair in a folder and report which need demultiplexing.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments; None means ``sys.argv[1:]``.

    Returns
    -------
    exit_code : int
        0 when a report was written, 1 when the map or the folder is
        unusable.

    Notes
    -----
    Writes the report table to ``--out`` (or stdout) and the summary to
    standard error.
    """
    args = parse_arguments(argv)
    try:
        map_rows = read_map_rows(args.map)
    except ValueError as error:
        print_stderr(str(error))
        exit_code = 1
        return exit_code

    pairs = find_fastq_pairs(args.fastq_dir)
    if not pairs:
        print_stderr("no *_R1 FASTQ files found")
        exit_code = 1
        return exit_code

    pool_libraries = group_libraries_by_pool(map_rows)
    library_names = {
                     row["out_name"]
                     for row in map_rows
                     }
    library_names |= {
                      row["library_id"]
                      for row in map_rows
                      }
    forward_indexes, reverse_indexes, index_pair_to_library = build_index_lookups(map_rows)
    report_rows = []
    for pair in pairs:
        report_row = classify_fastq_pair(
                                         pair,
                                         pool_libraries,
                                         library_names,
                                         forward_indexes,
                                         reverse_indexes,
                                         index_pair_to_library,
                                         args.max_reads
                                         )
        report_rows.append(report_row)

    if args.out:
        with open(
                  args.out,
                  "w"
                  ) as report_file:
            write_report_rows(
                              report_rows,
                              report_file
                              )
    else:
        write_report_rows(
                          report_rows,
                          sys.stdout
                          )
    print_summary(
                  report_rows,
                  pairs,
                  pool_libraries,
                  args.fastq_dir,
                  args.map
                  )
    exit_code = 0

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
