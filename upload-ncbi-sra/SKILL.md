---
name: upload-ncbi-sra
description: Fill NCBI BioSample and SRA submission templates from a user's curated sample/library metadata file (Excel/TSV), including choosing the correct BioSample package (MIMARKS, MIMS, MIGS, Plant, etc.), mapping fields, validating cross-links, checking a folder of FASTQ libraries for pooled files, and demultiplexing them into per-sample files. Use this skill whenever the user is depositing sequencing data at NCBI — mentions SRA, BioSample, BioProject, the NCBI Submission Portal, MIxS/MIMARKS/MIMS packages, "SRA metadata", accession numbers for a paper, or wants to turn their sample sheet / supplementary metadata table into submission spreadsheets, or asks whether sequencing files need demultiplexing before upload — even if they only ask about one piece (e.g. "which BioSample package should I pick?").
metadata:
  author: Kenneth Acosta
  version: "1.00"
---

# NCBI BioSample + SRA submission from curated metadata

Goal: turn the user's curated metadata workbook into (1) a filled BioSample
template for the correct package and (2) a filled SRA metadata template, both
internally consistent and ready to export as TSV for the Submission Portal.

The user's metadata file is the source of truth. Never invent values; when a
field needs a judgement call (ontology term, organism name, instrument), fill
your best-supported value, highlight the cell, and list it in the summary.

## 0. Check the live NCBI pages first

NCBI revises packages and template versions. At the start of every run, with
web fetch/search if available:

1. Fetch the BioSample template list: https://submit.ncbi.nlm.nih.gov/biosample/template/
   — this is the primary source for current package names, versions, and the
   downloadable template files. Follow its links to get the current template
   for the chosen package (and its field definitions) when the user hasn't
   uploaded one.
2. Fetch the SRA submission overview: https://submit.ncbi.nlm.nih.gov/about/sra/
   — check for changed requirements (file formats, data standards, the SRA
   metadata template, submission steps).
3. If needed, the package documentation (definitions, mandatory attributes):
   https://www.ncbi.nlm.nih.gov/biosample/docs/packages/
4. Compare with `references/biosample_packages.md` (a dated snapshot). If the
   live pages differ (new version number, renamed/retired package, new
   mandatory attributes or SRA rules), trust the live pages and tell the user
   what changed.
5. If web access is unavailable, say so in one line, use the snapshot, and
   rely on the template files the user downloaded from the portal — their
   header rows are the authoritative field lists for that version.

For a walkthrough of the portal submission itself, point the user to NCBI's
video: https://youtu.be/PTg9Ru68fc0 (the skill does not depend on its content;
Claude cannot watch video).

The template the user uploads (downloaded from the portal at submission time)
always wins over anything remembered: read its headers with
`scripts/template_fields.py` rather than assuming column names.

## 1. Inventory the inputs

Read every sheet of the metadata workbook (`pandas.read_excel(sheet_name=None)`).
Identify:
- **sample level**: collection ID, date/time + time zone, site, coordinates
  (note any flagged "approximate"), host, tissue/material, treatment, storage.
- **library level**: library ID, sample ID it belongs to, replicate, primers,
  index sequences, PCR conditions, clean-up, sequencing tube/pool, provider
  project, run.
- **controls** (NTC, extraction blanks, mock communities).
- Column-name quirks (typos, trailing spaces) — strip and map explicitly.

Ask only for what is genuinely missing and needed (e.g. BioProject accession,
time zone of recorded times, file names). Read `references/field_rules.md` for
formatting rules before filling anything.

## 2. Choose the BioSample package

Decide from what was sequenced and how — see the decision guide in
`references/biosample_packages.md`. Key logic:
- Marker-gene/amplicon (16S, 18S, ITS, COI) from environmental or host-associated
  material, uncultured → **MIMARKS: survey**, then the environmental sub-package
  matching the source (plant-associated, host-associated, soil, water, …).
- Shotgun metagenome → **MIMS**; isolate genome → **MIGS**; cultured/vouchered
  specimen marker genes → **MIMARKS: specimen**; the organism itself (plant
  tissue RNA-seq) → Plant / Model organism packages.
State the choice and one-line reason; offer the runner-up if it's close.
If the portal doesn't show the expected package name, map to what it does show
(e.g. "MIMARKS Survey related → plant-associated").

## 3. Fill the BioSample template

One row = one biological sample (not one library). Rules that repeatedly matter:
- Technical/library replicates (A/B/C) are SRA runs under one BioSample.
  Physically different material (e.g. surface-sterilized vs not) = separate BioSamples.
- `organism` must be an NCBI Taxonomy metagenome name for MIMARKS/MIMS
  (e.g. "plant metagenome", "leaf metagenome", "soil metagenome"); it may
  differ per sample. Flag for the user to confirm the exact string exists.
- Real host goes in `host` (Latin binomial); symbionts/parasites/gall inducers
  in `host_symbiont`.
- NCBI rejects samples whose attributes are identical apart from name/title/
  description. Add a distinguishing attribute (`source_material_id`, replicate,
  or a custom field) and verify with `scripts/validate_submission.py`.
- Useful columns from the metadata with no template slot → custom attributes
  appended after the last column (lower_snake_case).
- Controls: include as BioSamples, `neg_cont_type`/`pos_cont_type` filled,
  mandatory environment fields "not applicable".
- Dates ISO 8601; times must be UTC (`2024-07-19T15:00Z`) — convert from the
  recorded local zone (ask if unclear; mind DST, e.g. EDT = UTC−4).
Write into the template workbook starting at the first empty row under the
header, preserving the template's instructions rows; set text cell format
(`@`) so Excel doesn't reformat dates.

## 4. Fill the SRA template

One row = one run (library) linked by `sample_name` to a BioSample.
- `library_ID` unique; `title` like
  "16S rRNA gene V4 amplicon sequencing of plant metagenome: <sample info>".
- `design_description` per library from actual metadata (target region,
  primers + indexes, PCR cycles/protocol, clean-up, provider/service, read
  length) — make it vary where the protocol varied.
- Controlled vocab (strategy/source/selection/layout/platform/instrument) must
  come from the template's "Library and Platform Terms" sheet.
- **Instrument**: take from provider reports; confirm with FASTQ header
  instrument IDs (`references/instrument_ids.md`). If sources conflict, tell
  the user and let them decide; per-project read length can differ.
- **Files**: each run's FASTQs must contain only that BioSample's reads.
  Provider files that pool several samples must be demultiplexed first —
  pooled-tube uploads are not acceptable. Use the library-folder workflow
  below whenever the user has FASTQ files (or a folder of them).

### Library folder check → demultiplex (when FASTQs are provided)

1. Build `MAP.tsv` from the metadata (tab-separated, header row):
   `pool  library_id  out_name  fwd_name  fwd_seq  rev_name  rev_seq`
   - `pool` = the provider file stem each library was sequenced in
     (tube/pool ID, e.g. `KA01` for `KA01_R1.fastq.gz`)
   - `out_name` = SRA-safe library stem (no spaces/parentheses; "(SS)" → "-SS");
     use the same names in the SRA sheet's filename columns
   - `fwd_seq`/`rev_seq` = inline index sequences read at the 5' end of R1/R2
   Aliases `tube`, `fwd`, `rev` are accepted for pool/fwd_name/rev_name.
2. Run `python scripts/check_libraries.py FASTQ_DIR MAP.tsv --out check.tsv`.
   It pairs R1/R2 files, matches file names to pools/libraries, and tallies
   inline index combinations in the first 20,000 read pairs. Decisions:
   - `DEMUX`: a pool containing ≥2 libraries → demultiplex.
   - `OK_PER_LIBRARY`: already one library per file → submit as is.
   - `CHECK`: inconclusive (unknown file name, no inline indexes found, a
     single library under a pool name that needs renaming, missing mate) —
     inspect a few reads (`zcat X_R1.fastq.gz | head -8`) and ask the user.
   It also reports the share of reverse-orientation reads (reverse index on R1).
3. If any file is `DEMUX`, run the command the checker prints:
   `[BOTH_ORIENT=1] bash scripts/demux_for_sra.sh FASTQ_DIR OUT_DIR MAP.tsv`
   (`BOTH_ORIENT=1` only when the checker reports >5% reverse-orientation
   reads; the user's data may be forward-only, e.g. ligation-free amplicons.)
   Output: `OUT_DIR/upload/<out_name>_R1/_R2.fastq.gz` (reads unmodified,
   `--action none`), `OUT_DIR/read_counts.tsv`, logs and unassigned reads
   in `OUT_DIR/demux/<pool>/`. Libraries sharing an index make cutadapt warn
   "specified multiple times" — expected; assignment is still correct.
   If cutadapt isn't installed: `pip install cutadapt` (or conda).
4. Re-run `check_libraries.py OUT_DIR/upload MAP.tsv` — every file should be
   `OK_PER_LIBRARY` — then pass `OUT_DIR/read_counts.tsv` to the validator.
   Mixed folders (some pooled, some already per-library) are fine: submit the
   OK files as they are plus the demultiplexed ones.
5. Use the read counts: drop runs with 0 reads (and any BioSample left with
   no runs); keep low-count libraries unless the user says otherwise, but
   point out suspicious patterns (e.g. failures sharing an index).

If the data are large and the sandbox can't hold them, give the user the exact
commands to run locally instead, and ask them to send back `check.tsv` and
`read_counts.tsv`.

## 5. Validate and deliver

Run `scripts/validate_submission.py <biosample.xlsx> <sra.xlsx> [--counts read_counts.tsv]`
and fix everything it reports. It checks unique names, identical-attribute
collisions, SRA↔BioSample links, date/lat_lon formats, filename uniqueness,
and zero-read files.

Deliver the two filled .xlsx files (uncertain cells highlighted yellow with a
cell comment) and a short summary: package chosen and why, row counts, every
judgement call to verify, and next steps (export each data tab as TSV, upload
FASTQs, set a release/hold date), with a pointer to the walkthrough video
(section 0) for the portal steps. Mention that metadata can be corrected after
submission (SRA: sra@ncbi.nlm.nih.gov; BioSample: portal update or
biosamplehelp@ncbi.nlm.nih.gov) but reads are hard to change, so file
structure matters most.
