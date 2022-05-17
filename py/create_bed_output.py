#!/usr/bin/env python

"""
USAGE: ./parse_igdet_csv.py <path_to_predicted_genes> <bed_file_prefix>
"""

import sys 
import csv 

in_dir = sys.argv[1]
out_bed = sys.argv[2]

d_genes = f'{in_dir}/genes_D.tsv'
v_genes = f'{in_dir}/genes_V.tsv'
j_genes = f'{in_dir}/genes_J.tsv'

final_rows = []

def GenerateBedOutput(in_tsv, gene_type, out_list):
    def ref_dict(x):
        d = dict()
        for i in x.split('|'):
            j = i.split(':')
            d[j[0]] = j[1]
        return d
    gene_type_colors = {
        'd': ['D_gene', '0,128,0'],
        'v': ['V_gene', '255,153,153'],
        'j': ['J_gene', '255,255,0'],
    }
    with open(in_tsv, 'rt') as f:
        tbl = csv.reader(f, delimiter = '\t')
        header = dict()
        for i, j in enumerate(next(tbl)):
            header[j] = i
        start_col = header['start of gene']
        end_col = header['end of gene']
        strand_col = header['strand']
        name_col = header.get('best aligned human gene', None)
        for line in tbl:
            d = ref_dict(line[0])
            if name_col:
                name = f'{gene_type_colors[gene_type][0]}|{line[name_col]}'
            else:
                name = gene_type_colors[gene_type][0]
            color = gene_type_colors[gene_type][1]
            out_list.append(
                [
                    d['CONTIG'],
                    int(d['START']) + int(line[start_col]),
                    int(d['START']) + int(line[end_col]) + 1,
                    name,
                    1000,
                    line[strand_col],
                    ".",
                    ".",
                    color
                ]
            )
    return out_list

final_rows = GenerateBedOutput(d_genes, 'd', final_rows)
final_rows = GenerateBedOutput(v_genes, 'v', final_rows)
final_rows = GenerateBedOutput(j_genes, 'j', final_rows)

with open(out_bed, 'wt') as f:
    tbl = csv.writer(f, delimiter = '\t', quotechar = "\xb6",lineterminator = '\n')
    tbl.writerow(['track name="IgDetective Output" itemRgb="On"'])
    for i in final_rows:
        tbl.writerow(i)