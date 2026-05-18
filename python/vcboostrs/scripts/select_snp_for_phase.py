import os
import sys
from argparse import ArgumentParser
from collections import defaultdict


def FiterHeteSnpPhasing(ctgName, args):
    vcf_fn = args.vcf_fn
    split_folder = args.split_folder

    try:
        import vcboostrs
        count = vcboostrs.filter_heterozygous_snps(vcf_fn, split_folder, ctgName)
        print('[INFO] Total heterozygous SNP positions selected: {}: {}'.format(ctgName, count))
        return
    except Exception as e:
        print(f"[vcboostrs] filter failed ({e}), falling back to Python")

    variant_dict = defaultdict(str)
    qual_set = defaultdict(int)
    header = []

    import gzip
    opener = gzip.open if vcf_fn.endswith('.gz') else open
    with opener(vcf_fn, 'rt') as f:
        for row in f:
            row = row.rstrip()
            if row[0] == '#':
                header.append(row + '\n')
                continue
            columns = row.strip().split()
            ctg_name = columns[0]
            if ctgName and ctgName != ctg_name:
                continue
            ref_base = columns[3]
            alt_base = columns[4]
            genotype = columns[9].split(':')[0].replace('|', '/')

            if len(ref_base) == 1 and len(alt_base) == 1:
                if genotype == '0/1' or genotype == '1/0':
                    variant_dict[int(columns[1])] = row

    print('[INFO] Total heterozygous SNP positions selected: {}: {}'.format(ctgName, len(variant_dict)))

    out_path = os.path.join(split_folder, '{}.vcf'.format(ctgName))
    with open(out_path, 'w') as f:
        f.write(''.join(header))
        for key, row in sorted(variant_dict.items(), key=lambda x: x[0]):
            f.write(row + '\n')


def main():
    parser = ArgumentParser(description="Select heterozygous SNP candidates for WhatsHap phasing")
    parser.add_argument('--split_folder', type=str, default=None,
                        help="Path to directory for split VCF files")
    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Path of the input VCF file")
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    chr_list = list(range(1, 23)) + ['X', 'Y']
    for num in chr_list:
        contig_name = 'chr' + str(num)
        FiterHeteSnpPhasing(contig_name, args)


if __name__ == "__main__":
    main()
