#!/usr/bin/env bash
##################################################
# Objective
##################################################
# demultiplex pooled libraries (one R1/R2 pair per sequencing tube/pool)
# into per-library paired fastq files for NCBI SRA

##################################################
# Background
##################################################
# reads are sorted, not trimmed (--action none); SRA receives unmodified reads
# inline forward index at 5' end of R1, reverse index at 5' end of R2
# R1/R2 barcodes are paired by row order (--pair-adapters); works with
#   combinatorial designs where libraries in a pool share one index
# BOTH_ORIENT=1 also captures amplicons in reverse orientation (rev index on R1)
#   and merges them into the same library files; check_libraries.py reports
#   whether this is needed (default: forward orientation only)
# MAP.tsv: pool, library_id, out_name, fwd_name, fwd_seq, rev_name, rev_seq
#   out_name = SRA-safe filename stem, e.g. "(SS)" -> "-SS"
#
# usage: [ERR=0.15] [BOTH_ORIENT=0|1] demux_for_sra.sh FASTQ_DIR OUT_DIR MAP.tsv
# output: OUT_DIR/upload/<out_name>_R1/_R2.fastq.gz, OUT_DIR/read_counts.tsv,
#         per-pool logs + unassigned reads in OUT_DIR/demux/<pool>/
# requires: cutadapt (>=3), gzip, awk
# note: cutadapt warns "Adapter ... specified multiple times" when libraries
#   share an index; expected with --pair-adapters, assignment is still correct

set -euo pipefail
FASTQ_DIR=$(realpath "$1"); OUT=$(realpath -m "$2"); MAP=$(realpath "$3")
ERR=${ERR:-0.15}; BOTH=${BOTH_ORIENT:-0}
DEMUX="$OUT/demux"; UP="$OUT/upload"
mkdir -p "$DEMUX" "$UP"

find_fq(){ for e in .fastq.gz .fastq .fq.gz .fq; do
  [[ -f "$FASTQ_DIR/${1}_$2$e" ]] && { echo "$FASTQ_DIR/${1}_$2$e"; return; }; done; }

##################################################
# De-Multiplexing
##################################################
for POOL in $(tail -n +2 "$MAP" | cut -f1 | sort -u); do
  R1=$(find_fq "$POOL" R1); R2=$(find_fq "$POOL" R2)
  [[ -n $R1 && -n $R2 ]] || { echo "WARN: no FASTQs for $POOL, skipped" >&2; continue; }
  echo "$POOL"
  PD="$DEMUX/$POOL"; mkdir -p "$PD"

  # barcode files, paired by row order
  awk -F'\t' -v p="$POOL" '$1==p{print ">"$3"\n^"$5}' "$MAP" > "$PD/fwd_barcodes.fasta"
  awk -F'\t' -v p="$POOL" '$1==p{print ">"$3"\n^"$7}' "$MAP" > "$PD/rev_barcodes.fasta"
  if [[ $BOTH == 1 ]]; then   # reverse-orientation entries: rev index on R1, fwd on R2
    awk -F'\t' -v p="$POOL" '$1==p{print ">"$3"__rc\n^"$7}' "$MAP" >> "$PD/fwd_barcodes.fasta"
    awk -F'\t' -v p="$POOL" '$1==p{print ">"$3"__rc\n^"$5}' "$MAP" >> "$PD/rev_barcodes.fasta"
  fi

  cutadapt \
    --error-rate "$ERR" --no-indels --action none --pair-adapters \
    -g file:"$PD/fwd_barcodes.fasta" -G file:"$PD/rev_barcodes.fasta" \
    -o "$PD/{name}_R1.fastq.gz" -p "$PD/{name}_R2.fastq.gz" \
    --untrimmed-output "$PD/unknown_R1.fastq.gz" \
    --untrimmed-paired-output "$PD/unknown_R2.fastq.gz" \
    "$R1" "$R2" > "$PD/demux.log"

  # assemble upload files (forward + optional reverse orientation)
  for NAME in $(awk -F'\t' -v p="$POOL" '$1==p{print $3}' "$MAP"); do
    for RD in R1 R2; do
      : > "$UP/${NAME}_${RD}.fastq.gz"
      for S in "$NAME" "${NAME}__rc"; do
        [[ -f "$PD/${S}_${RD}.fastq.gz" ]] && cat "$PD/${S}_${RD}.fastq.gz" >> "$UP/${NAME}_${RD}.fastq.gz"
      done
    done
  done
done

##################################################
# Read Counts
##################################################
for FQ in "$UP"/*_R1.fastq.gz; do
  echo -e "$(basename "$FQ" _R1.fastq.gz)\t$(( $(gzip -cd "$FQ" | wc -l) / 4 ))"
done > "$OUT/read_counts.tsv"
awk -F'\t' '$2==0{z++} END{if(z) print "WARN: " z " libraries with 0 reads - remove before upload (see read_counts.tsv)"}' "$OUT/read_counts.tsv"
echo "done: $UP"
