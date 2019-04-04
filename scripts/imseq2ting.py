import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    argument_parser(parser)
    args = parser.parse_args()
    sample_folder = args.imseq_folder
    samples = [file[:-4] for file in os.listdir(sample_folder)]

    output_folder = args.output
    #parse filename from each sample in sample_sheet
    for sample in samples:
        print(sample)
        with open(f'{output_folder}/{sample}.tsv', 'w') as output_file:
            print(f'CDR3b\tTRBV\tTRBJ\tPatientID\tFrequency', file=output_file)
            #iterate over repertoire and write clonotypes to output_file
            with open(f'{sample_folder}/{sample}.act', 'r') as repertoire:
                for line in repertoire:
                    clonotype, frequency = line.split('\t')
                    TRBV, CDR3b, TRBJ = clonotype.split(':')
                    print(f'{CDR3b}\t{TRBV}\t{TRBJ}\t{sample}\t{frequency}', file=output_file, end='')


def argument_parser(parser):
    parser.add_argument('-i', '--imseq_folder',
                        help='Path to directory containing imseq repertoires as act-file',
                        type=str,
                        required=True)
    parser.add_argument('-o', '--output',
                        help='Path to output directory',
                        type=str,
                        required=True)


if __name__ == '__main__':
    main()
