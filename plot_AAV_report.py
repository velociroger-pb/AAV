#!/usr/bin/env python3
# Roger Volden

'''
Reproduce Liz's AAV report (originally in R)

Plots:
    - Distribution of mapped reference start position ]
    - Distribution of mapped reference end position   ] *.summary.csv
    - Distribution of mapped reference length         ]

    - Distribution of non-matches by reference position: substitutions ]
    - Distribution of non-matches by reference position: deletions     ] *.nonmatch_stat.csv
    - Distribution of non-matches by reference position: insertions    ]
    
    - Distribution of mapped identity to reference                                         } *.summary.csv
    - Distribution of non-matches by ref position and size of non-match                    ] *.nonmatch_stat.csv
    - Distribution of non-matches by ref position and size of non-match - limit to <100bp  ]

    - RL distribution by assigned AAV type } *.per_read.csv

Tables:
    - Length distribution of different non-matches ] *.nonmatch_stat.csv
        error type, length, count, frequency       ]
    - AAV type breakdown                           ] *.per_read.csv
        assigned type, subtype, count, frequency   ]

Usage:
    python3 plot_AAV_report.py \
        -i bc1001 \
        -o test
'''

import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
from palette import Palette as pbp
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_prefix', '-i', type=str, required=True,
        help='The input file prefix (output from summarize_AAV_alignment.py)'
    )
    parser.add_argument(
        '--output_prefix', '-o', type=str, required=True,
        help='Output file prefix'
    )
    return parser.parse_args()

def check_files(in_prefix):
    for suffix in ['.summary.csv', '.per_read.csv', '.nonmatch_stat.csv']:
        current_file = in_prefix + suffix
        if not os.path.exists(current_file):
            print(f"Couldn't find {current_file}")
            exit(1)

def parse_summary(in_prefix):
    '''
    Gets all of the relevant info out of *.summary.csv
        - map start distribution
        - map end distribution
        - map len distribution
        - map identity distribution

    mapping_dists = {
        starts = [0, 0, ...],
        ends = [2000, 2000, ...],
        lens = [1000, 1000, ...],
        ids = [0.99, 0.99, ...]
    }
    '''
    mapping_dists = {'Starts': [], 'Ends': [], 'Lengths': [], 'Identities': []}
    col_to_idx, first = {}, True
    with open(in_prefix + '.summary.csv') as f:
        for line in f:
            line = line.rstrip().split('\t')
            if first:
                for idx, col in enumerate(line):
                    col_to_idx[col] = idx
                first = False
                continue
            mapping_dists['Starts'].append(int(line[col_to_idx['map_start0']]))
            mapping_dists['Ends'].append(int(line[col_to_idx['map_end1']]))
            mapping_dists['Lengths'].append(int(line[col_to_idx['map_len']]))
            mapping_dists['Identities'].append(float(line[col_to_idx['map_iden']]))
    return mapping_dists

def parse_per_read(in_prefix, out_prefix):
    '''
    RL dist per assigned type
    Table:
        assigned type | assigned subtype | count | frequency (%)

    type_RL = { ssAAV = [2000, 4000, ...], unknown = [5000, ...] } <- to be made into violin plots
    type_table = { ssAAV: {subtype: count}, ... }
    '''
    type_RL, type_table = {}, {}
    col_to_idx, first, nreads = {}, True, 0
    with open(in_prefix + '.per_read.csv') as f:
        for line in f:
            line = line.rstrip().split('\t')
            if first:
                for idx, col in enumerate(line):
                    col_to_idx[col] = idx
                first = False
                continue
            rlen = int(line[col_to_idx['read_len']])
            a_type = line[col_to_idx['assigned_type']]
            a_subtype = line[col_to_idx['assigned_subtype']]

            # numbers for the read lengths
            if a_type not in type_RL:
                type_RL[a_type] = []
            type_RL[a_type].append(rlen)

            # numbers for the table
            if a_type not in type_table:
                type_table[a_type] = {a_subtype: 1}
            elif a_subtype not in type_table[a_type]:
                type_table[a_type][a_subtype] = 1
            else:
                type_table[a_type][a_subtype] += 1
            nreads += 1

    # output the type table
    out_fh = open(out_prefix + '.type_breakdown.tsv', 'w+')
    print('Assigned type\tAssigned subtype\tCount\tFrequency (%)', file=out_fh)
    for a_type, subdict in type_table.items():
        for a_subtype, count in subdict.items():
            print(f'{a_type}\t{a_subtype}\t{count}\t{count/nreads*100:.2f}', file=out_fh)
    out_fh.close()

    return type_RL

def parse_nonmatch(in_prefix):
    pass

def plot_mapping_dists(mapping_dists, ref_range, out):
    for plot_type, dist in mapping_dists.items():
        plt.figure(figsize=(6, 3))
        plt.style.use('clean')
        h = plt.axes([0.15, 0.125, 0.8, 0.8])

        if plot_type != 'Identities':
            xmin, xmax = ref_range
            heights, bins = np.histogram(dist, bins=100, range=(xmin, xmax))
            h_total = sum(heights)
            heights = [x/h_total for x in list(heights)]
        else:
            xmin, xmax = 0.95, 1
            heights, bins = np.histogram(dist, bins=100, range=(xmin, xmax))

        binsize = bins[1] - bins[0]
        for left, height in zip(bins[:-1], heights):
            bar = Rect((left, 0), binsize, height, lw=0, fc='black')
            h.add_patch(bar)

        h.set_xlim(xmin, xmax)
        h.set_ylim(0, max(heights)*1.05)
        if plot_type != 'Identities':
            h.set_xlabel(f'Mapped reference {plot_type.lower()}')
            h.set_ylabel('Fraction of reads')
        else:
            h.set_xlabel('Mapping identity')
            h.set_ylabel('Read count')

        output = out + f'.{plot_type.lower()}_dist.png'
        plt.savefig(output, dpi=600)
        plt.close()

def main(args):
    mapping_dists = parse_summary(args.input_prefix)
    type_RL = parse_per_read(args.input_prefix, args.output_prefix)
    ref_range = (min(mapping_dists['Starts']), max(mapping_dists['Ends']))
    plot_mapping_dists(mapping_dists, ref_range, args.output_prefix)

if __name__ == '__main__':
    args = parse_args()
    check_files(args.input_prefix)
    main(args)
