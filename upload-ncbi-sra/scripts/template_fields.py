#!/usr/bin/env python3
"""List the fields of an NCBI BioSample or SRA template (.xlsx/.tsv).

Finds the header row (first row whose first cell is '*sample_name' or
'sample_name'), prints each field, whether it is mandatory ('*' prefix),
and the first line of its cell comment if any.
usage: template_fields.py TEMPLATE [--sheet NAME]
"""
import sys, argparse, openpyxl, csv

ap = argparse.ArgumentParser(); ap.add_argument('template'); ap.add_argument('--sheet')
a = ap.parse_args()
if a.template.endswith(('.tsv', '.txt')):
    rows = list(csv.reader(open(a.template), delimiter='\t'))
    hdr_i = next(i for i, r in enumerate(rows) if r and r[0].lstrip('*') == 'sample_name')
    for f in rows[hdr_i]:
        print(('MANDATORY ' if f.startswith('*') else 'optional  ') + f.lstrip('*'))
    sys.exit()
wb = openpyxl.load_workbook(a.template)
ws = wb[a.sheet] if a.sheet else wb.worksheets[0]
hdr = None
for r in ws.iter_rows(max_row=60):
    v = r[0].value
    if isinstance(v, str) and v.lstrip('*').strip() == 'sample_name':
        hdr = r; break
if hdr is None:
    sys.exit(f'no header row with sample_name found in sheet {ws.title}; sheets: {wb.sheetnames}')
print(f'# sheet={ws.title} header_row={hdr[0].row} first_data_row={hdr[0].row+1}')
for c in hdr:
    if c.value is None: continue
    f = str(c.value)
    note = (c.comment.text.strip().splitlines()[0][:120] if c.comment else '')
    print(f"{'MANDATORY' if f.startswith('*') else 'optional '}\t{c.column_letter}\t{f.lstrip('*')}\t{note}")
