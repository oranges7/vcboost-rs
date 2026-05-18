import sys
import pysam
import time
import os
import argparse
import numpy as np
from argparse import ArgumentParser

import vcboostrs


def run(cp0, fastafile, samfile, pileup_file, miss_output, hap_reads_0, hap_reads_1):
    mincov = 1
    maxcov = 100

    flag = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
    chrom = cp0.strip().split(' ')[0]
    v_pos = int(cp0.strip().split(' ')[1])
    info1 = cp0.strip().split(' ')[2].replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace(
        'K', '-').replace('M', '-').replace('S', '-').replace('B', '-')
    info2 = cp0.strip().split(' ')[3].replace('N', '-').replace('R', '-').replace('Y', '-').replace('W', '-').replace(
        'K', '-').replace('M', '-').replace('S', '-').replace('B', '-')
    ref = fastafile.fetch(chrom, v_pos - 1, v_pos + 200)
    if ',' in info1:
        info1 = info1.split(',')[0]
    if ',' in info2:
        info2 = info2.split(',')[0]

    ref = info2 + ref[len(info1):]
    ref = ref[:160]

    d = {0: {}, 1: {}}
    in_run = 0
    for pcol in samfile.pileup(chrom, v_pos - 1, v_pos, min_base_quality=12, min_mapping_quality=28, flag_filter=flag,
                               truncate=True):
        for pread in pcol.pileups:
            dt = pread.alignment.query_sequence[max(0, pread.query_position_or_next):pread.query_position_or_next + 160]
            mq = pread.alignment.mapping_quality

            if pread.alignment.qname in hap_reads_0:
                d[0][pread.alignment.qname] = dt
            elif pread.alignment.qname in hap_reads_1:
                d[1][pread.alignment.qname] = dt

        if len(d[0]) == 0:
            ref_data_0, data_0 = np.zeros((72, 5)), np.zeros((72, 5))
        else:
            ref_data_0, data_0 = vcboostrs.fast_msa(d[0], ref, 2, 100)
            ref_data_0 = np.array(ref_data_0, dtype=np.float32)
            data_0 = np.array(data_0, dtype=np.float32)

        if len(d[1]) == 0:
            ref_data_1, data_1 = np.zeros((72, 5)), np.zeros((72, 5))
        else:
            ref_data_1, data_1 = vcboostrs.fast_msa(d[1], ref, 2, 100)
            ref_data_1 = np.array(ref_data_1, dtype=np.float32)
            data_1 = np.array(data_1, dtype=np.float32)

        data = np.dstack([ref_data_0, data_0, ref_data_1, data_1])
        in_run = 1
        s = chrom + ' ' + str(v_pos) + ' ' + ' '.join(str(x) for x in data.reshape(-1).astype(np.int16)) + '\n'
        pileup_file.write(s)

    if in_run == 0:
        s = chrom + ' ' + str(v_pos) + '\n'
        miss_output.write(s)


def gen(args):
    fb = args.pos_path
    sam_path = args.bam_path
    ref_path = args.ref_path

    pileup_file = open(args.out_path + '/' + fb.split('/')[-1].split('\\')[-1] + '.pileup', "w")
    miss_output = open(args.out_path + '/' + 'miss' + '.pos', "a")
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile = pysam.FastaFile(ref_path)
    kf = open(fb, 'r')
    chrom = args.chromosome
    length = fastafile.get_reference_length(chrom)
    hap_dict = {1: [], 2: []}
    for pread in samfile.fetch(chrom, 0, length + 1):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0 = set(hap_dict[1])
    hap_reads_1 = set(hap_dict[2])
    for chr_pos in kf:
        run(chr_pos, fastafile, samfile, pileup_file, miss_output, hap_reads_0, hap_reads_1)


def main():
    parser = ArgumentParser(
        description="VCboost-RS: Rust-accelerated feature generation for variant filtering")
    parser.add_argument('--bam_path', '-b', type=str, default=None,
                        help="Path to the BAM file")
    parser.add_argument('--chromosome', '-c', type=str, default=None,
                        help="Chromosome")
    parser.add_argument('--ref_path', '-r', type=str, default=None,
                        help="Path of the reference file")
    parser.add_argument('--pos_path', '-p', type=str, default=None,
                        help="Path of the input sites file")
    parser.add_argument('--out_path', '-o', type=str, default=None,
                        help="Path of the output file")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    gen(args)


if __name__ == '__main__':
    main()
