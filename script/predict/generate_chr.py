import sys


def main():
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    chromosomes = ["chr" + str(i) for i in range(1, 23)]
    chromosomes_xy = chromosomes + ["chrX", "chrY"]
    with open(f2, 'w') as out:
        if f1 == 'chr1-22':
            for chr in chromosomes: print(chr, file=out)
        if f1 == 'chr1-22XY':
            for chr in chromosomes_xy: print(chr, file=out)
        if f1 in chromosomes_xy: print(f1, file=out)


if __name__ == '__main__':
    main()
