import argparse


def get_options():
    parser = argparse.ArgumentParser(description='Reverse-complement a gff.',
                                     prog='reverse_complement_gff')
    parser.add_argument('--input_gff', help='input gff', required=True)
    parser.add_argument('--output_gff', help='output gff', required=True)
    return parser.parse_args()


def reverse_complement_gff(gff_list):
    """reverse complements a gff i.e. as if annotations refer to reverse-complement
    of sequence region"""
    offset = max([int(entry[4]) for entry in gff_list])
    for entry in gff_list:
        new_start = offset - int(entry[4])
        new_end = offset - int(entry[3])
        new_strand = "-" if entry[6]=="+" else "+"
        entry[3] = new_start
        entry[4] = new_end
        entry[6] = new_strand
    return(gff_list)



def read_gff_to_list(input_gff):
    """reads in a gff to a list"""
    gff_list = []
    with open(input_gff, "r") as f:
        for line in f.readlines():
            if line.startswith("#"):
                    pass
            if line.startswith("##FASTA"): # Don't read in fasta components
                    break
            else:
                line = line.strip("\n").split("\t")
                if len(line)==9:
                    gff_list += [line]
    return(gff_list)

def main():
    args = get_options()
    gff_list = read_gff_to_list(args.input_gff)
    print(gff_list)
    rev_gff_list = reverse_complement_gff(gff_list)
    with open(args.output_gff, "w") as f:
        f.write("##gff-version 3\n")
        for e in rev_gff_list:
            f.write("\t".join([str(x) for x in e])+"\n")


if __name__=="__main__":
    main()

