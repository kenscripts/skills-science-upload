#!/usr/bin/env python3
"""Cross-check filled BioSample and SRA templates before upload.

usage: validate_submission.py BIOSAMPLE.xlsx SRA.xlsx [--sra-sheet SRA_data]
                              [--counts counts.tsv]   # 2 cols: file_stem  n_reads
Exit code 1 if errors were found.
"""
import re, sys, argparse, openpyxl, collections

ap = argparse.ArgumentParser()
ap.add_argument('biosample'); ap.add_argument('sra')
ap.add_argument('--sra-sheet', default='SRA_data'); ap.add_argument('--counts')
a = ap.parse_args()
err, warn = [], []

def table(path, sheet=None):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    for r in ws.iter_rows(max_row=60):
        if isinstance(r[0].value, str) and r[0].value.lstrip('*').strip() == 'sample_name':
            h = [str(c.value).lstrip('*').strip() if c.value else None for c in r]
            start = r[0].row + 1; break
    else:
        sys.exit(f'header not found in {path}')
    rows = []
    for r in ws.iter_rows(min_row=start, values_only=True):
        if r[0] in (None, ''): continue
        rows.append({h[i]: ('' if v is None else str(v).strip()) for i, v in enumerate(r) if h[i]})
    return rows

bs = table(a.biosample); sra = table(a.sra, a.sra_sheet)

# BioSample checks
names = [r['sample_name'] for r in bs]
for n, k in collections.Counter(names).items():
    if k > 1: err.append(f'BioSample: duplicate sample_name {n}')
skip = {'sample_name', 'sample_title', 'description', 'bioproject_accession'}
sig = collections.defaultdict(list)
for r in bs:
    sig[tuple(sorted((k, v) for k, v in r.items() if k not in skip and v))].append(r['sample_name'])
for s in sig.values():
    if len(s) > 1: err.append(f'BioSample: identical attributes (NCBI will reject): {s}')
date_re = re.compile(r'^\d{4}(-\d{2}(-\d{2}(T\d{2}:\d{2}(:\d{2})?Z)?)?)?$')
ll_re = re.compile(r'^\d+(\.\d+)? [NS] \d+(\.\d+)? [EW]$')
missing = {'not applicable', 'not collected', 'missing', 'not provided', 'restricted access'}
for r in bs:
    d = r.get('collection_date', '')
    if d and d not in missing and not date_re.match(d):
        err.append(f"BioSample {r['sample_name']}: bad collection_date '{d}' (ISO 8601, times must end in Z/UTC)")
    ll = r.get('lat_lon', '')
    if ll and ll not in missing and not ll_re.match(ll):
        err.append(f"BioSample {r['sample_name']}: bad lat_lon '{ll}'")
    org = r.get('organism', '')
    if org and 'metagenome' not in org and org not in missing:
        warn.append(f"BioSample {r['sample_name']}: organism '{org}' has no 'metagenome' (required for MIMARKS/MIMS)")

# SRA checks
libs = [r.get('library_ID', '') for r in sra]
for n, k in collections.Counter(libs).items():
    if k > 1: err.append(f'SRA: duplicate library_ID {n}')
fcols = [k for k in (sra[0].keys() if sra else []) if k.startswith('filename')]
files = [r[c] for r in sra for c in fcols if r.get(c)]
for n, k in collections.Counter(files).items():
    if k > 1: err.append(f'SRA: file used more than once {n}')
for f in files:
    if re.search(r'[()\s]', f): warn.append(f'SRA: awkward characters in filename {f}')
sset, bset = {r['sample_name'] for r in sra}, set(names)
for s in sorted(sset - bset): err.append(f'SRA sample_name not in BioSample sheet: {s}')
for s in sorted(bset - sset): err.append(f'BioSample with no SRA run: {s}')
for r in sra:
    for c in ('library_strategy', 'library_source', 'library_selection', 'library_layout',
              'platform', 'instrument_model', 'design_description', 'filetype'):
        if c in r and not r[c]: err.append(f"SRA {r.get('library_ID')}: empty {c}")

# read counts
if a.counts:
    cnt = {}
    for line in open(a.counts):
        p = line.strip().split('\t')
        if len(p) >= 2 and p[-1].isdigit(): cnt[p[0]] = int(p[-1])
    for r in sra:
        f = r.get(fcols[0], '') if fcols else ''
        stem = re.sub(r'_R?1(_001)?\.f(ast)?q(\.gz)?$', '', f)
        n = cnt.get(stem, cnt.get(r.get('library_ID')))
        if n is None: warn.append(f'no read count for {f}')
        elif n == 0: err.append(f'SRA {f}: 0 reads — remove run (and BioSample if it has no other runs)')
        elif n < 1000: warn.append(f'SRA {f}: only {n} reads')

print(f'BioSamples: {len(bs)}  SRA runs: {len(sra)}  files: {len(files)}')
for w in warn: print('WARN ', w)
for e in err: print('ERROR', e)
print('OK' if not err else f'{len(err)} error(s)')
sys.exit(1 if err else 0)
