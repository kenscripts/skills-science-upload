# BioSample packages — snapshot + decision guide

Snapshot date: 2026-09 (from memory and the portal list seen during a 16S
submission). ALWAYS compare with the live pages:
- Template list (primary): https://submit.ncbi.nlm.nih.gov/biosample/template/
- Package docs: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/
- SRA submission overview: https://submit.ncbi.nlm.nih.gov/about/sra/
Version numbers (e.g. 6.0) and exact names change; the live page and the
template the user downloads win.

## Decision guide

| What was sequenced | Package family |
|---|---|
| Marker gene (16S/18S/ITS/COI) amplified directly from environment/host material, not cultured | MIMARKS: survey → env sub-package |
| Marker gene from a cultured isolate or voucher specimen | MIMARKS: specimen |
| Shotgun metagenome / metatranscriptome of a community | MIMS: metagenome/environmental → env sub-package |
| Genome of a cultured bacterium/archaeon | MIGS: cultured bacteria/archaea |
| Eukaryote / virus / plasmid / organelle genome | MIGS: eukaryote / virus / plasmid / organelle |
| MAG or SAG | MIMAG / MISAG |
| Clinical or environmental pathogen isolate | Pathogen: clinical or host-associated / environmental |
| The organism itself (tissue RNA-seq, resequencing) | Plant / Model organism or animal / Invertebrate / Human / Microbe |
| Environmental sample, no MIxS reporting wanted | Metagenome or environmental |

## MIxS environmental sub-packages (MIMARKS survey, MIMS, MIGS share these)
air · built environment · host-associated · human-associated · human-gut ·
human-oral · human-skin · human-vaginal · microbial mat/biofilm ·
miscellaneous natural or artificial environment · plant-associated ·
sediment · soil · wastewater/sludge · water

Pick by the material the DNA came from:
- leaves, roots, rhizosphere tissue, galls, seeds, endophytes → plant-associated
- animal gut/skin/tissue (non-human) → host-associated
- bulk soil (incl. rhizosphere soil sampled as soil) → soil
- lab synthetic communities / bioreactors → miscellaneous (or host-associated if in a host)

Portal naming seen: the package list may show "MIMARKS Survey related" and then
ask for the environmental package on a second step.
