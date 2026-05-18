import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser


def run(args):
    d = {}
    combined_array = np.empty((0, 4))
    all_file = args.in_path
    output = args.out_path
    a = args.threshold
    or_file = args.original_file
    with open(all_file) as f:
        for i in f:
            single_file = i.strip()
            single_array = pd.read_csv(single_file + '.txt', sep=' ', header=None)
            combined_array = np.concatenate([combined_array, single_array], axis=0)
    chr_col = combined_array[:, 0].astype(str)
    pos_col = combined_array[:, 1].astype(int)
    pred_labels_one = combined_array[:, 2:].astype(float)
    print(f"threshold: {a}")
    pred_labels = np.where(pred_labels_one[:, 1] > a, 1, 0)
    fin_file = output + '/vc_boost.vcf'
    for index, x in enumerate(pred_labels):
        if x == 0: d[str(chr_col[index]) + ' ' + str(pos_col[index])] = 1
    add_pos = args.add_file
    with open(add_pos, 'r') as add:
        for line in add:
            atom = line.strip().split(' ')
            d[str(atom[0]) + ' ' + str(atom[1])] = 1
    exclude_xy = args.conting == 'chr1-22'
    try:
        import vcboostrs
        kept = vcboostrs.merge_predictions_to_vcf(or_file, fin_file, d, exclude_xy)
        print(f"[vcboostrs] VCF merge done, kept {kept} records")
    except Exception as e:
        print(f"[vcboostrs] merge failed ({e}), falling back to Python")
        with open(or_file, 'r') as vcf, open(fin_file, 'w') as fl:
            for line in vcf:
                if line[0] == '#': print(line.strip(), file=fl); continue
                atom = line.strip().split('\t')
                if exclude_xy:
                    if atom[0] == 'chrX' or atom[0] == 'chrY': continue
                if atom[0] + ' ' + atom[1] not in d: print(line.strip(), file=fl)


def main():
    parser = ArgumentParser(description="Merge prediction results into VCF file")
    parser.add_argument('--in_path', '-i', type=str, default=None, help="Path to the input file list")
    parser.add_argument('--threshold', '-t', type=float, default=0.02, help="Threshold for filtering")
    parser.add_argument('--out_path', '-o', type=str, default=None, help="Path of the output directory")
    parser.add_argument('--original_file', '-f', type=str, default=None, help="Original VCF file path")
    parser.add_argument('--add_file', '-a', type=str, default=None, help="Additional positions file")
    parser.add_argument('--conting', '-c', type=str, default=None, help="Contig range")
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0: parser.print_help(); sys.exit(1)
    run(args)


if __name__ == '__main__':
    main()
