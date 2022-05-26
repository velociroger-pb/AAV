#!/usr/bin/env python3
# Roger Volden

'''
Reproduce Liz's AAV report (originally in R)

Plots:
    [x] Distribution of mapped reference start position ]
    [x] Distribution of mapped reference end position   ] *.summary.csv
    [x] Distribution of mapped reference length         ]

    [x] Distribution of non-matches by reference position: substitutions ]
    [x] Distribution of non-matches by reference position: deletions     ] *.nonmatch_stat.csv
    [x] Distribution of non-matches by reference position: insertions    ]
    
    [x] Distribution of mapped identity to reference                                         } *.summary.csv
    [x] Distribution of non-matches by ref position and size of non-match                    ] *.nonmatch_stat.csv
    [x] Distribution of non-matches by ref position and size of non-match - limit to <100bp  ]

    [x] RL distribution by assigned AAV type } *.per_read.csv

Tables:
    [x] Length distribution of different non-matches ] *.nonmatch_stat.csv
        error type, length, count, frequency         ]
    [x] AAV type breakdown                           ] *.per_read.csv
        assigned type, subtype, count, frequency     ]

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
    Table: assigned type | assigned subtype | count | frequency (%)

    type_rl = { ssAAV = [2000, 4000, ...], unknown = [5000, ...] } (for violin plots)
    type_table = { ssAAV: {subtype: count}, ... }
    '''
    type_rl, type_table = {}, {}
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
            if a_type not in type_rl:
                type_rl[a_type] = []
            type_rl[a_type].append(rlen)

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

    return type_rl

def parse_nonmatch(in_prefix, out_prefix):
    '''
    Gets all of the relevant info out of *.nonmatch_stat.csv
        - dist of each error mode (per base)
        - position vs error mode and size
            - position vs error mode and size zoom in to <100bp
        - table breakdown of error mode and length
    '''
    error_dict = {} # error type: [(error_pos, error_len), ...]
    col_to_idx, first = {}, True
    total_errs = 0
    with open(in_prefix + '.nonmatch_stat.csv') as f:
        for line in f:
            line = line.rstrip().split('\t')
            if first:
                # columns: read_id, pos0, type, type_len
                for idx, col in enumerate(line):
                    col_to_idx[col] = idx
                first = False
                continue
            position = int(line[col_to_idx['pos0']])
            error_type = line[col_to_idx['type']]
            error_len = int(line[col_to_idx['type_len']])

            if error_type not in error_dict:
                error_dict[error_type] = []
            error_dict[error_type].append((position, error_len))
            total_errs += 1

    out_fh = open(out_prefix + '.nonmatch_lengths.tsv', 'w+')
    print('Error type\tError length\tCount\tFrequency (%)', file=out_fh)
    idx_to_lbin = {0: '0-10', 1: '11-100', 2: '101-500', 3: '>500'}
    for e_type, err_list in error_dict.items():
        tmp_len = [0] * 4 # [0-10, 11-100, 101-500, >500]
        for e_tuple in err_list:
            if e_tuple[1] < 11:
                tmp_len[0] += 1
            elif e_tuple[1] < 101:
                tmp_len[1] += 1
            elif e_tuple[1] < 501:
                tmp_len[2] += 1
            else:
                tmp_len[3] += 1
        for idx, count in enumerate(tmp_len):
            if not count:
                continue
            print(f'{e_type}\t{idx_to_lbin[idx]}\t{count}\t{count/total_errs*100:.2f}', file=out_fh)
    out_fh.close()

    return error_dict

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
            h.set_title(f'Distribution of mapped reference {plot_type.lower()}')
        else:
            h.set_xlabel('Mapping identity')
            h.set_ylabel('Read count')
            h.set_title('Distribution of mapped identity to reference')

        output = out + f'.{plot_type.lower()}_dist.png'
        plt.savefig(output, dpi=600)
        plt.close()

def plot_rl_violins(type_rl, out):
    plt.figure(figsize=(6, 3))
    plt.style.use('clean')
    v = plt.axes([0.15, 0.125, 0.8, 0.8])

    labels, curr_pos = [], 1
    for a_type, lens in type_rl.items():
        # setting points=len(lens) can be dangerous for big datasets
        farts = v.violinplot(lens, positions=[curr_pos], points=len(lens))
        labels.append(a_type)
        curr_pos += 1
        for pc in farts['bodies']:
            pc.set_color('black')
            pc.set_alpha(1)
        farts['cmins'].set_color('black')
        farts['cmaxes'].set_color('black')
        farts['cbars'].set_color('black')

    v.set_xlim(0, curr_pos)
    v.set_xticks(range(1, curr_pos))
    v.set_xticklabels(labels)
    v.set_xlabel('Assigned AAV type')
    v.set_ylabel('Read length')
    v.set_title('Distribution of read lengths by assigned AAV type')
    plt.savefig(out + '_rl_dist.png', dpi=600)

def plot_error_dists(error_dict, ref_range, out):
    # error type: [(error_pos, error_len), ...]
    type_to_name = {'D': 'deletion', 'I': 'insertion', 'X': 'mismatch'}
    type_to_plural = {'D': 'deletions', 'I': 'insertions', 'X': 'mismatches'}
    for err_type, err_list in error_dict.items():
        plt.figure(figsize=(6, 3))
        plt.style.use('clean')
        h = plt.axes([0.15, 0.125, 0.8, 0.8])

        name, plural = type_to_name[err_type], type_to_plural[err_type]

        # err_dist = list(range(ref_range[0], ref_range[1]))
        err_dist = [0] * ref_range[1]
        for e in err_list:
            # index at the err pos - 1
            err_dist[e[0]-1] += 1

        for left, height in enumerate(err_dist):
            # lines are hard to see if the width is only 1 :/
            # another way to get around visibility is to bump dpi to like 1200
            box = Rect((left-0.5, 0), 2, height, lw=1, fc='black')
            h.add_patch(box)

        max_height = max(err_dist)
        h.set_xlim(0-(ref_range[1]*0.05), ref_range[1]*1.05)
        h.set_ylim(0-(max_height*0.05), max_height*1.05)
        h.set_xlabel('Reference position')
        h.set_ylabel('Error count')
        h.set_title(f'Distribution of {plural} by reference position')

        output = out + f'.{name}_dist.png'
        plt.savefig(output, dpi=600)
        plt.close()

def plot_error_lens(error_dict, ref_range, out):
    # error type: [(error_pos, error_len), ...]
    plt.figure(figsize=(6, 3))
    plt.style.use('clean')
    s = plt.axes([0.15, 0.125, 0.8, 0.8])
    zoom = 'full'
    err_to_color = {'D': pbp.pink, 'I': pbp.green, 'X': pbp.blue}

    legend = []
    for e_type, err_list in error_dict.items():
        x, y = list(zip(*err_list))
        color = err_to_color[e_type]
        line, = s.plot(x, y, lw=0, c=color, marker='o', ms=5, mew=0, alpha=0.5, label=e_type)
        legend.append(line)

    s.legend(handles=legend)
    s.set_ylabel('Error length')
    s.set_xlabel('Reference position')
    s.set_title('Distribution of error lengths')
    output = out + f'.{zoom}_err_length_dist.png'
    plt.savefig(output, dpi=600)

    s.set_ylim(0, 100)
    s.set_title('Distribution of error lengths (<100 bp)')
    zoom = 'zoomed'
    output = out + f'.{zoom}_err_length_dist.png'
    plt.savefig(output, dpi=600)

    plt.close()

def main(args):
    mapping_dists = parse_summary(args.input_prefix)
    type_rl = parse_per_read(args.input_prefix, args.output_prefix)
    ref_range = (min(mapping_dists['Starts']), max(mapping_dists['Ends']))
    error_dict = parse_nonmatch(args.input_prefix, args.output_prefix)

    plot_mapping_dists(mapping_dists, ref_range, args.output_prefix)
    plot_rl_violins(type_rl, args.output_prefix)
    plot_error_dists(error_dict, ref_range, args.output_prefix)
    plot_error_lens(error_dict, ref_range, args.output_prefix)

if __name__ == '__main__':
    args = parse_args()
    check_files(args.input_prefix)
    main(args)
