# Field formatting rules (BioSample + SRA)

## Missing values
Mandatory fields cannot be blank. Use INSDC terms: "not collected",
"not applicable", "missing", "restricted access" (or "not provided").

## Dates and times
- collection_date: ISO 8601, `YYYY-MM-DD`, `YYYY-MM`, or `YYYY`.
- Times only in UTC: `YYYY-MM-DDThh:mmZ`. Convert from the recorded zone;
  confirm the zone with the user (labels like "EST" in summer are often
  really EDT, UTC−4). Omit time if unknown.
- Store as text in Excel.

## Location
- geo_loc_name: `Country: Region, Locality` e.g. `USA: New Jersey, New Brunswick`.
- lat_lon: `DD.dddddd N DD.dddddd W` (decimal degrees, hemisphere letters,
  single spaces). Convert signed decimals. Note approximate coordinates in
  the description.

## Environment (MIxS)
env_broad_scale / env_local_scale / env_medium: ontology term + ID, e.g.
`temperate broadleaf forest biome [ENVO:01000202]`,
`plant-associated environment [ENVO:01001001]`, `leaf [PO:0025034]`.
Plant anatomy (PO) terms are fine for env_medium. Do not guess IDs silently:
verify in EBI OLS (https://www.ebi.ac.uk/ols4/) when web access exists;
otherwise fill the term, highlight, and ask the user to verify. Consider
urban vs forest biome based on the actual site.

## Organism (MIMARKS / MIMS)
Must be an NCBI Taxonomy name ending in "metagenome". Examples: plant
metagenome, leaf metagenome, phyllosphere metagenome, root metagenome,
rhizosphere metagenome, soil metagenome, gut metagenome, insect metagenome,
metagenome (generic, e.g. for PCR blanks). Confirm exact strings in Taxonomy.

## Uniqueness
- sample_name unique; ≥1 non-name attribute must differ between BioSamples.
- library_ID unique; each FASTQ file name used once.

## SRA controlled vocab (typical amplicon study)
library_strategy AMPLICON · library_source METAGENOMIC · library_selection PCR ·
library_layout paired · platform ILLUMINA · filetype fastq.
Always check against the template's terms sheet.

## Title patterns
BioSample title: "<host> <tissue/condition> (<sample_name>)".
SRA title: "<method> of <organism>: <sample info> (<sample>, replicate <x>)".
