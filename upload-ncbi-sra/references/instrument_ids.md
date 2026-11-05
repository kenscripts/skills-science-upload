# Identifying the sequencer from FASTQ headers

Illumina header: `@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> ...`
Check one read per provider project:
`zcat X_R1.fastq.gz | head -1`

| Instrument ID prefix | Instrument | Confidence |
|---|---|---|
| M | MiSeq | high |
| A | NovaSeq 6000 | high |
| LH | NovaSeq X / X Plus | moderate |
| VH, VL | NextSeq 1000/2000 | moderate |
| NB, NS | NextSeq 500/550 | moderate |
| FS | iSeq 100 | moderate |
| MN | MiniSeq | moderate |
| K, J, E, D | HiSeq 4000/3000/X/2500 families | moderate |
| SH, SL | possibly MiSeq i100 series | low — confirm with provider |
| provider-specific (e.g. GWNJ-1012) | renamed by provider — undecodable | — |

Run/flow-cell fields may contain provider labels (e.g. "...Miseq") that
describe the service, not the instrument. When reports and headers conflict,
ask the provider which instrument each project ran on. The user may choose to
submit with the provider-reported model and correct later.

Read length (2x250, 2x300) can change between provider projects even on the
same service; record it per project in design_description.
