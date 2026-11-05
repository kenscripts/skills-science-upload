#!/usr/bin/env python3
"""Inspect a folder of paired FASTQs and decide whether each pair must be
demultiplexed before SRA submission (SRA needs one BioSample per run).

usage: check_libraries.py FASTQ_DIR MAP.tsv [--n 20000] [--out report.tsv]
MAP.tsv (tab, header): pool library_id out_name fwd_name fwd_seq rev_name rev_seq
  pool      = provider tube/pool/file stem the library was sequenced in
  out_name  = SRA-safe per-library file stem (no spaces/parentheses)
  *_seq     = inline index at the 5' end of R1 (fwd) and R2 (rev)

Evidence per file pair:
  name    : stem matches a pool (with >1 library) or a single library/out_name
  content : first N read pairs tallied by inline index combination
            (<=1 mismatch); >=2 expected library combos each >=1% => pooled
  orient  : share of reads whose R1 starts with a reverse index (mixed
            orientation => demultiplex with BOTH_ORIENT=1)
Exit 0 always; decision column: DEMUX | OK_PER_LIBRARY | CHECK
"""
import argparse, collections, csv, glob, gzip, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument('fastq_dir'); ap.add_argument('map')
ap.add_argument('--n', type=int, default=20000); ap.add_argument('--out')
a = ap.parse_args()

ALIAS = {'tube': 'pool', 'fwd': 'fwd_name', 'rev': 'rev_name'}
rows = [{ALIAS.get(k.strip(), k.strip()): (v or '').strip() for k, v in r.items()}
        for r in csv.DictReader(open(a.map), delimiter='\t')]
need = {'pool', 'library_id', 'out_name', 'fwd_name', 'fwd_seq', 'rev_name', 'rev_seq'}
if not need <= set(rows[0]): sys.exit(f'map needs columns {sorted(need)}')
pools = collections.defaultdict(list)
for r in rows: pools[r['pool']].append(r)
libnames = {r['out_name']: r for r in rows} | {r['library_id']: r for r in rows}
fwd = {r['fwd_seq'].upper(): r['fwd_name'] for r in rows}
rev = {r['rev_seq'].upper(): r['rev_name'] for r in rows}
combo2lib = {(r['fwd_name'], r['rev_name']): r['out_name'] for r in rows}

def ham1(read, idx):
    return sum(x != y for x, y in zip(read[:len(idx)], idx)) <= 1 and len(read) >= len(idx)

def call(seq, table):
    hits = [n for s, n in table.items() if ham1(seq, s)]
    return hits[0] if len(hits) == 1 else None

def opener(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def seqs(p, n):
    with opener(p) as f:
        for i, line in enumerate(f):
            if i // 4 >= n: break
            if i % 4 == 1: yield line.strip().upper()

pairs = []
for r1 in sorted(glob.glob(os.path.join(a.fastq_dir, '*'))):
    m = re.match(r'(.+?)[._](R?)1(_001)?\.f(ast)?q(\.gz)?$', os.path.basename(r1))
    if not m: continue
    b2 = re.sub(r'([._]R?)1((_001)?\.f(ast)?q(\.gz)?)$', r'\g<1>2\2', os.path.basename(r1))
    r2 = os.path.join(os.path.dirname(r1), b2)
    r2 = r2 if os.path.exists(r2) and r2 != r1 else None
    pairs.append((m.group(1), r1, r2))
if not pairs: sys.exit('no *_R1 FASTQ files found')

out = []
for stem, r1, r2 in pairs:
    if r2 is None:
        out.append(dict(file=stem, decision='CHECK', reason='no R2 mate found')); continue
    if stem in pools:
        name_ev = f'pool with {len(pools[stem])} libraries' if len(pools[stem]) > 1 else 'pool with 1 library'
    elif stem in libnames: name_ev = 'matches a library'
    else: name_ev = 'unknown name'
    tally = collections.Counter(); n = matched = rc = 0
    for s1, s2 in zip(seqs(r1, a.n), seqs(r2, a.n)):
        n += 1
        f, r = call(s1, fwd), call(s2, rev)
        if f and r: tally[(f, r)] += 1; matched += 1; continue
        if call(s1, rev) and call(s2, fwd): rc += 1
    expected = [c for c, k in tally.items() if c in combo2lib and k >= max(1, 0.01 * n)]
    libs_seen = sorted(combo2lib[c] for c in expected)
    frac = matched / n if n else 0
    rc_frac = rc / n if n else 0
    if len(libs_seen) >= 2 or (name_ev.startswith('pool with') and len(pools[stem]) > 1 and frac > 0.1):
        dec, why = 'DEMUX', f'{len(libs_seen)} libraries detected by inline index'
    elif len(libs_seen) == 1 and frac > 0.5:
        dec, why = ('OK_PER_LIBRARY' if name_ev == 'matches a library' else 'CHECK',
                    f'single library {libs_seen[0]}' + ('' if name_ev == 'matches a library' else ' — rename file to its out_name'))
    elif frac < 0.1 and name_ev == 'matches a library':
        dec, why = 'OK_PER_LIBRARY', 'name matches a library; no inline indexes found (already trimmed or index-read design?)'
    else:
        dec, why = 'CHECK', 'inconclusive — inspect reads / ask user'
    if dec == 'DEMUX' and stem not in pools:
        dec, why = 'CHECK', why + '; file name is not a pool in MAP — set its pool name in MAP so demux_for_sra.sh can find it'
    if rc_frac > 0.05 and dec == 'DEMUX': why += f'; {rc_frac:.0%} reverse-orientation reads -> use BOTH_ORIENT=1'
    out.append(dict(file=stem, decision=dec, reason=why, name_evidence=name_ev, reads_checked=n,
                    frac_index_matched=round(frac, 3), frac_reverse_orient=round(rc_frac, 3),
                    libraries_detected=';'.join(libs_seen)))

cols = ['file', 'decision', 'reason', 'name_evidence', 'reads_checked', 'frac_index_matched',
        'frac_reverse_orient', 'libraries_detected']
w = csv.DictWriter(open(a.out, 'w') if a.out else sys.stdout, cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
w.writeheader(); [w.writerow(r) for r in out]
c = collections.Counter(r['decision'] for r in out)
print(f'\n# {dict(c)}', file=sys.stderr)
if c['DEMUX']:
    anyrc = any(r['decision'] == 'DEMUX' and r.get('frac_reverse_orient', 0) > 0.05 for r in out)
    print(f"# demultiplex: {'BOTH_ORIENT=1 ' if anyrc else ''}bash demux_for_sra.sh {a.fastq_dir} <OUT_DIR> {a.map}", file=sys.stderr)
stems = {s for s, _, _ in pairs}
missing = [p for p, libs in pools.items() if p not in stems
           and not all(l['out_name'] in stems or l['library_id'] in stems for l in libs)]
if missing: print(f'# pools in map with no FASTQs here: {sorted(missing)}', file=sys.stderr)
