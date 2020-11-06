import re
import sys
import os.path
import argparse
import warnings
import numpy as np
from collections import Counter
from itertools import combinations
from scipy.stats import fisher_exact


def main():
    parser = argparse.ArgumentParser()
    argument_parser(parser)
    args = parser.parse_args()
    cluster_local = args.no_local
    cluster_global = args.no_global
    use_structural_boundaries = args.use_structural_boundaries
    print(f"Loading sequences from: {args.tcr_sequences}")
    tcr_sequences = (
        load_filtered_tcr_sequences(args.tcr_sequences)
        if args.stringent_filtering
        else load_tcr_sequences(args.tcr_sequences)
    )
    print(f"Unique sequences loaded: {len(tcr_sequences)}")
    if cluster_local:
        if not os.path.exists(args.kmer_file):
            # Process significant kmers
            if args.kmers_gliph:
                kmer_preprocessing_gliph(tcr_sequences, args)
            else:
                kmer_preprocessing(tcr_sequences, args)
        else:
            print("Motif file found")
        # Clustering by local similarity
        print("Local clustering")
        local_struct = local_clustering(
            tcr_sequences, args.kmer_file, use_structural_boundaries
        )
    # Clustering by global similarity
    if cluster_global:
        print("Global clustering")
        global_edges = global_clustering(tcr_sequences, use_structural_boundaries)
        final_struct = local_struct if cluster_local else UnionFind(tcr_sequences)

        for seq1, seq2 in global_edges:
            final_struct.union(seq1, seq2)
        final_clusters = summarize_clusters(tcr_sequences, final_struct)
    else:
        final_clusters = summarize_clusters(tcr_sequences, local_struct)
    output_clusters(args.output, final_clusters)


def kmer_preprocessing(tcr_sequences, args):
    if os.path.isfile(args.reference):
        reference = args.reference
    else:
        script_path = os.path.dirname(os.path.abspath(__file__))
        reference = f"{script_path}/{args.reference}"
    reference_sequences = load_reference_sequences(reference)
    number_seqs_condition = len(tcr_sequences)
    number_seqs_control = len(reference_sequences)
    print("Counting kmers in sample set")
    kmers_condition = count_kmers(tcr_sequences, args.use_structural_boundaries)
    print(f"Unique kmers in sample set: {len(kmers_condition)}")
    print("Counting kmers in control set")
    kmers_control = count_kmers(reference_sequences, args.use_structural_boundaries)
    print(f"Unique kmers in control set: {len(kmers_control)}")
    print("Identifying significant kmers")
    identify_significant_kmers(
        number_seqs_condition,
        kmers_condition,
        number_seqs_control,
        kmers_control,
        args.max_p_value,
        args.kmer_file,
    )


# --------------Optional preprocessing as implemented in gliph-------------------
def kmer_preprocessing_gliph(tcr_sequences, args):
    """
    Non-deterministic determination of significant motifs as performed by bugfixed gliph implemention
    """
    if os.path.isfile(args.reference):
        reference = args.reference
    else:
        script_path = os.path.dirname(os.path.abspath(__file__))
        reference = f"{script_path}/{args.reference}"

    print("Counting kmers in sample set")
    kmers_sample = count_kmers(tcr_sequences, args.use_structural_boundaries)
    most_frequent_kmers = {
        kmer: kmers_sample[kmer] for kmer in kmers_sample if kmers_sample[kmer] >= 3
    }

    print(f"\tkmers obtained: {len(kmers_sample)}")
    print(f"\tminimum frequency >= 3 kmers obtained: {len(most_frequent_kmers)}")
    print(f"Loading reference from {reference}")
    reference_sequences = load_reference_sequences(reference)

    print("Resampling from reference")
    number_sequences = len(tcr_sequences)
    if len(reference_sequences) < number_sequences:
        warnings.warn(
            "Sample size is larger then reference set size.\n"
            "Possible solutions:\n\t- Downsampling of input sequences\n"
            "\t- Use larger reference file",
            Warning,
        )
        sys.exit(1)
    simulated_kmers = [
        simulate_sample_set(
            reference_sequences, number_sequences, args.use_structural_boundaries
        )
        for i in range(1000)
    ]

    print("Analyzing kmers")
    analyze_kmers(
        most_frequent_kmers,
        simulated_kmers,
        number_sequences,
        1000,
        args.gliph_minp,
        args.kmer_file,
    )


def simulate_sample_set(sequences, set_size, use_structural_boundaries):
    random_sequences = np.random.choice(sequences, size=set_size, replace=False)
    kmers = count_kmers(random_sequences, use_structural_boundaries)
    return kmers


def analyze_kmers(
    most_frequent_kmers,
    simulated_kmers,
    set_size,
    number_simulations,
    max_p_value,
    kmer_file,
):
    with open(kmer_file, "w") as output_file:
        print("Motif\tCounts\tavgRef\ttopRef\tOvE\tp-value", file=output_file)
        for kmer in most_frequent_kmers:
            kmer_count_sample = most_frequent_kmers[kmer]
            odds_as_enriched_as_discovery = 0
            highest = 0
            average = 0
            for simulated_set in simulated_kmers:
                kmer_occurrence_control = simulated_set[kmer]
                if kmer_occurrence_control >= kmer_count_sample:
                    odds_as_enriched_as_discovery += 1 / number_simulations
                highest = (
                    kmer_occurrence_control
                    if kmer_occurrence_control > highest
                    else highest
                )
                average += kmer_occurrence_control / number_simulations
            if average > 0:
                ove = kmer_count_sample / average
            else:
                average = 1 / (number_simulations * set_size)
                ove = kmer_count_sample / average  # line contains bug fix
            ove = int(ove * 1000) / 1000  # truncate decimals
            average = int(average * 100) / 100  # truncate decimals
            if odds_as_enriched_as_discovery < max_p_value:
                if odds_as_enriched_as_discovery == 0:
                    odds_as_enriched_as_discovery = 1 / number_simulations
                minfoldchange = get_minfoldchange(kmer_count_sample)
                if ove >= minfoldchange:
                    print(
                        f"{kmer}\t{kmer_count_sample}\t{average}\t{highest}\t{ove}\t{odds_as_enriched_as_discovery}",
                        file=output_file,
                    )


# ----------- Gliph's kmer preprocessing ends here-----------------------


def get_minfoldchange(kmer_count):
    if kmer_count == 1:
        return 10000
    if kmer_count == 2:
        return 1000
    elif kmer_count == 3:
        return 100
    return 10


def identify_significant_kmers(
    seqs_condition,
    kmers_condition,
    seqs_control,
    kmers_control,
    p_value_threshold,
    kmer_file,
):
    seqs_condition += 2  # add pseudo count
    seqs_control += 2  # add pseudo count
    bonferroni_threshold = p_value_threshold / len(kmers_condition)
    with open(kmer_file, "w") as output_file:
        print("Motif\todds\tp-value\tcondition_count\tcontrol_count", file=output_file)
        for condition_kmer, condition_count in kmers_condition.items():
            minfoldchange = get_minfoldchange(condition_count)
            control_count = (
                kmers_control[condition_kmer] + 1
            )  # add one pseudo occurence
            condition_count += 1  # add one pseudo occurence
            table = [[control_count, condition_count], [seqs_control, seqs_condition]]
            oddsratio, p_value = fisher_exact(table, alternative="less")
            if p_value <= bonferroni_threshold and oddsratio < 1 / minfoldchange:
                print(
                    f"{condition_kmer}\t{oddsratio}\t{p_value}\t{condition_count}\t{control_count}",
                    file=output_file,
                )


def load_tcr_sequences(sequence_file):
    tcr_sequences = set()
    with open(sequence_file, "r") as sequences:
        line = sequences.readline()  # check if first line is header
        sequence = line.split("\t")[0]
        if sequence.upper() != "CDR3B":
            tcr_sequences.add(sequence)
        for line in sequences.readlines():
            sequence = line.split("\t")[0]
            sequence = sequence.upper()
            tcr_sequences.add(sequence)
    return list(tcr_sequences)


def load_filtered_tcr_sequences(sequence_file):
    tcr_sequences = set()
    with open(sequence_file, "r") as sequences:
        for line in sequences.readlines():
            sequence = line.split("\t")[0]
            sequence = sequence.upper()
            if re.match("^C[AC-WY][AC-WY][AC-WY][AC-WY]*F", sequence):
                tcr_sequences.add(sequence)
    return list(tcr_sequences)


def load_reference_sequences(reference_file):
    cdr3b_sequences = set()
    with open(reference_file, "r") as reference:
        for i, line in enumerate(reference.readlines()):
            if i % 2 == 1:
                cdr3b_sequences.add(line.strip())
    return list(cdr3b_sequences)


def count_kmers(sequences, use_structural_boundaries):
    # counting all 2-, 3- and 4-mers
    c = Counter()
    offset = 0 if use_structural_boundaries else 3
    for sequence in sequences:
        seqlen = len(sequence)
        for k in (2, 3, 4):
            kmers = [
                sequence[i : i + k] for i in range(offset, seqlen + 1 - k - offset)
            ]
            c.update(kmers)
    return c


def global_clustering(tcr_sequences, use_structural_boundaries):
    max_distance = (
        1 if len(tcr_sequences) >= 125 else 2
    )  # maximum hamming distance defined by gliph
    tcr_sequences_dict, sequences_original = sequences_to_dict(
        tcr_sequences, use_structural_boundaries
    )
    global_edges = set()
    for length in tcr_sequences_dict.keys():
        sequences = np.array(tcr_sequences_dict[length])  # list of sequences
        difference = sequences[:, np.newaxis] - sequences
        nonzeros = np.count_nonzero(difference, axis=2)
        edges = zip(*np.where(nonzeros <= max_distance))
        for x, y in edges:
            global_edges.add(
                (sequences_original[length][x], sequences_original[length][y])
            )
    return global_edges


def sequences_to_dict(tcr_sequences, use_structural_boundaries):
    # Sequences are added to dict by same length
    min_length = 8  # sequences must be at least 8 chars long
    sequences_reduced = dict()
    sequences_original = dict()
    for (i, sequence) in enumerate(tcr_sequences):
        seq_length = len(sequence)
        if seq_length >= min_length:
            sequence_reduced = (
                sequence[3:-3] if not use_structural_boundaries else sequence
            )
            sequence_reduced = np.frombuffer(sequence_reduced.encode(), dtype="uint8")
            if seq_length not in sequences_reduced:
                sequences_reduced[seq_length] = list()
                sequences_original[seq_length] = list()
            sequences_original[seq_length].append(i)
            sequences_reduced[seq_length].append(sequence_reduced)
    return sequences_reduced, sequences_original


def output_clusters(output, clusters_tcr):
    number_sequences = len([seq for cluster in clusters_tcr for seq in cluster])
    number_clusters = len([cluster for cluster in clusters_tcr if len(cluster) > 1])
    print(f"Clusters: {number_clusters}")
    with open(output, "w") as output_file:
        for cluster in clusters_tcr:
            cluster_content = f"{len(cluster)}\t{list(cluster)[0]}\t" + " ".join(
                x for x in cluster
            )
            print(cluster_content, file=output_file)


def local_clustering(tcr_sequences, kmer_file, use_structural_boundaries):
    two_mers, three_mers, four_mers = load_kmers(kmer_file)
    print("\tClustering kmers...")
    distinct_kmers = reduce_kmers(two_mers, three_mers)
    distinct_kmers = reduce_kmers(distinct_kmers, four_mers)
    print("\tClustering CDR3b sequences...")
    clusters_tcr = cluster_sequences(
        distinct_kmers, tcr_sequences, use_structural_boundaries
    )
    return clusters_tcr


def reduce_kmers(distinct_kmers, candidate_kmers):
    redundant_kmers = []
    for (i, candidate_kmer) in enumerate(candidate_kmers):
        for unique_kmer in distinct_kmers:
            if unique_kmer in candidate_kmer:
                redundant_kmers.append(i)
                break
    np.delete(candidate_kmers, redundant_kmers)
    return np.append(distinct_kmers, candidate_kmers)


def cluster_sequences(kmers, sequences, use_structural_boundaries):
    cluster_struct = UnionFind(sequences)
    for kmer in kmers:
        subcluster = []
        for i in range(len(sequences)):
            sequence = (
                sequences[i][3:-3] if not use_structural_boundaries else sequences[i]
            )
            if kmer in sequence:
                subcluster.append(i)
        for i in range(len(subcluster) - 1):
            cluster_struct.union(subcluster[i], subcluster[i + 1])
    return cluster_struct


def summarize_clusters(sequences, cluster_struct):
    clusters = dict()
    for i, sequence in enumerate(sequences):
        cluster_id = cluster_struct.find_representative(i)
        if cluster_id not in clusters:
            clusters[cluster_id] = [sequence]
        else:
            clusters[cluster_id].append(sequence)
    return [set(cluster) for cluster in clusters.values()]


def load_kmers(input_file):
    print(f"Loading kmers from {input_file}")
    two_mers = []
    three_mers = []
    four_mers = []
    with open(input_file, "r") as kmer_file:
        kmer_file.readline()
        for line in kmer_file.readlines():
            kmer = line.split()[0]
            if len(kmer) == 2:
                two_mers.append(kmer)
            elif len(kmer) == 3:
                three_mers.append(kmer)
            elif len(kmer) == 4:
                four_mers.append(kmer)
            else:
                print(kmer)
    return np.array(two_mers), np.array(three_mers), np.array(four_mers)


class UnionFind:
    def __init__(self, sequences):
        self.representatives = np.array(range(len(sequences)))

    def union(self, i, j):
        representative_i = self.find_representative(i)
        representative_j = self.find_representative(j)
        if representative_i < representative_j:
            self.representatives[representative_j] = representative_i
        elif representative_i > representative_j:
            self.representatives[representative_i] = representative_j

    def find_representative(self, i):
        if self.representatives[i] == i:
            return i
        else:
            self.representatives[i] = self.find_representative(self.representatives[i])
        return self.representatives[i]


def argument_parser(parser):
    parser.add_argument(
        "-t", "--tcr_sequences", help="File holding TCRs", type=str, required=True
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference fasta file of naive CDR3 amino acid sequences used for estimation of significant k-mers.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-k",
        "--kmer_file",
        help="tab separated file holding kmers in first row",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="path of output-file", type=str, required=True
    )
    parser.add_argument(
        "-b",
        "--use_structural_boundaries",
        help="First and last three amino acids are included in processing",
        action="store_true",
    )
    parser.add_argument(
        "-ng",
        "--no_global",
        help="If set global clustering is excluded",
        action="store_false",
    )
    parser.add_argument(
        "-nl",
        "--no_local",
        help="If set local clustering is excluded",
        action="store_false",
    )
    parser.add_argument(
        "-p",
        "--max_p_value",
        help="p-value threshold for identifying significant k-mers by fisher exact test",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--gliph_minp",
        help="probability threshold for identifying significant k-mers by gliph test",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "-f",
        "--stringent_filtering",
        help="If used only TCRs starting with a cystein and ending with phenylalanine will be used",
        action="store_true",
    )
    parser.add_argument(
        "-g",
        "--kmers_gliph",
        help="If set kmers are identified by the non-deterministic approach as implemented by gliph",
        action="store_true",
    )


if __name__ == "__main__":
    main()
