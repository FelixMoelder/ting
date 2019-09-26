import re
import sys
import os.path
import argparse
import warnings
import numpy as np
import networkx as nx
from collections import Counter
from multiprocessing import Pool
from scipy.stats import fisher_exact


def main():
    parser = argparse.ArgumentParser()
    argument_parser(parser)
    args = parser.parse_args()
    cluster_local = args.no_local
    cluster_global = args.no_global
    use_structural_boundaries = args.use_structural_boundaries
    print(f"Loading sequences from: {args.tcr_sequences}")
    tcr_sequences = load_filtered_tcr_sequences(args.tcr_sequences) if args.stringent_filtering else load_tcr_sequences(args.tcr_sequences)
    print(f'Unique sequences loaded: {len(tcr_sequences)}')
    final_clusters = nx.Graph()
    final_clusters.add_nodes_from(tcr_sequences)
    if cluster_local:
        if not os.path.exists(args.kmer_file):
            # Process significant kmers
            if args.kmers_gliph:
                kmer_preprocessing_gliph(tcr_sequences, args)
            else:
                kmer_preprocessing(tcr_sequences, args)
        else:
            print('Motif file found')
        # Clustering by local similarity
        print('Local clustering')
        local_clusters = local_clustering(tcr_sequences, args.kmer_file)
        final_clusters.add_edges_from(local_clusters.edges())
    # Clustering by global similarity
    if cluster_global:
        print('Global clustering')
        global_clusters = global_clustering(tcr_sequences, use_structural_boundaries)
        final_clusters.add_edges_from(global_clusters.edges())
    output_clusters(args.output, final_clusters)


def kmer_preprocessing(tcr_sequences, args):
    if os.path.isfile(args.reference):
        reference = args.reference
    else:
        script_path = os.path.dirname(os.path.abspath(__file__))
        reference = f'{script_path}/{args.reference}'
    reference_sequences = load_reference_sequences(reference)
    number_seqs_condition = len(tcr_sequences)
    number_seqs_control = len(reference_sequences)
    print('Counting kmers in sample set')
    kmers_condition = count_kmers(tcr_sequences, args.use_structural_boundaries)
    print('Unique kmers in sample set: {len(kmers_condition)}')
    print('Counting kmers in control set')
    kmers_control = count_kmers(reference_sequences, args.use_structural_boundaries)
    print('Unique kmers in control set: {len(kmers_control)}')
    print('Identifying significant kmers')
    identify_significant_kmers(number_seqs_condition, kmers_condition, number_seqs_control, kmers_control, args.max_p_value, args.kmer_file)


# --------------Optional preprocessing as implemented in gliph-------------------
def kmer_preprocessing_gliph(tcr_sequences, args):
    """
    Non-deterministic determination of significant motifs as performed by bugfixed gliph implemention
    """
    if os.path.isfile(args.reference):
        reference = args.reference
    else:
        script_path = os.path.dirname(os.path.abspath(__file__))
        reference = f'{script_path}/{args.reference}'

    print('Counting kmers in sample set')
    kmers_sample = count_kmers(tcr_sequences, args.use_structural_boundaries)
    most_frequent_kmers = {kmer: kmers_sample[kmer] for kmer in kmers_sample
                           if kmers_sample[kmer] >= 3}

    print(f'\tkmers obtained: {len(kmers_sample)}')
    print(f'\tminimum frequency >= 3 kmers obtained: {len(most_frequent_kmers)}')
    print(f'Loading reference from {reference}')
    reference_sequences = load_reference_sequences(reference)

    print('Resampling from reference')
    number_sequences = len(tcr_sequences)
    if len(reference_sequences) < number_sequences:
        warnings.warn("Sample size is larger then reference set size.\n"
        "Possible solutions:\n\t- Downsampling of input sequences\n"
        "\t- Use larger reference file", Warning)
        sys.exit(1)
    with Pool(1) as p:
        simulation_params = [(reference_sequences, number_sequences, args.use_structural_boundaries)
                            for i in range(1000)]
        simulated_kmers = p.starmap(simulate_sample_set, simulation_params)

    print('Analyzing kmers')
    analyze_kmers(most_frequent_kmers, simulated_kmers, number_sequences, 1000, args.gliph_minp,
                  args.kmer_file)


def simulate_sample_set(sequences, set_size, use_structural_boundaries):
    random_sequences = np.random.choice(sequences, size=set_size, replace=False)
    kmers = count_kmers(random_sequences, use_structural_boundaries)
    # total_counted_k_mers.append(kmers)
    return kmers


def analyze_kmers(most_frequent_kmers, simulated_kmers, set_size, number_simulations, max_p_value, kmer_file):
    with open(kmer_file, 'w') as output_file:
        print('Motif\tCounts\tavgRef\ttopRef\tOvE\tp-value', file=output_file)
        for kmer in most_frequent_kmers:
            kmer_count_sample = most_frequent_kmers[kmer]
            odds_as_enriched_as_discovery = 0
            highest = 0
            average = 0
            for simulated_set in simulated_kmers:
                kmer_occurrence_sim = simulated_set[kmer]
                if kmer_occurrence_control >= kmer_count_sample:
                    odds_as_enriched_as_discovery += 1/number_simulations
                highest = kmer_occurrence_control if kmer_occurrence_control > highest else highest
                average += kmer_occurrence_control/number_simulations
            if average > 0:
                ove = kmer_count_sample/average
            else:
                average = 1/(number_simulations*set_size)
                ove = kmer_count_sample/average # line contains bug fix
            ove = int(ove * 1000)/1000  # truncate decimals
            average = int(average*100)/100  # truncate decimals
            if odds_as_enriched_as_discovery < max_p_value:
                if odds_as_enriched_as_discovery == 0:
                    odds_as_enriched_as_discovery = 1/number_simulations
                minfoldchange = get_minfoldchange(kmer_count_sample)
                if ove >= minfoldchange:
                    print(f'{kmer}\t{kmer_count_sample}\t{average}\t{highest}\t{ove}\t{odds_as_enriched_as_discovery}',
                          file=output_file)

# ----------- Gliph's kmer preprocessing ends here-----------------------

def get_minfoldchange(kmer_count):
    if kmer_count == 2:
        return 1000
    elif kmer_count == 3:
        return 100
    elif kmer_count >= 4:
        return 10



def identify_significant_kmers(seqs_condition, kmers_condition, seqs_control, kmers_control, p_value_threshold, kmer_file):
    kmers_condition_count = sum(kmers_condition.values())
    kmers_control_count = sum(kmers_control.values())
    seqs_condition += 2 # add pseudo count
    seqs_control += 2 # add pseudo count
    bonferroni_threshold = p_value_threshold/len(kmers_condition)
    with open(kmer_file, 'w') as output_file:
         print('Motif\todds\tp-value\tcondition_count\tcontrol_count', file=output_file)
         for condition_kmer, condition_count in kmers_condition.items():
            minfoldchange = get_minfoldchange(condition_count)
            control_count = kmers_control[condition_kmer] + 1 # add one pseudo occurence
            condition_count += 1 # add one pseudo occurence
            table = [[control_count, condition_count], [seqs_control, seqs_condition]]
            oddsratio, p_value = fisher_exact(table, alternative='less')
            if p_value <= bonferroni_threshold and oddsratio < 1/minfoldchange:
                print(f'{condition_kmer}\t{oddsratio}\t{p_value}\t{condition_count}\t{control_count}', file=output_file)


def load_tcr_sequences(sequence_file):
    tcr_sequences = set()
    with open(sequence_file, 'r') as sequences:
        line = sequences.readline()  # check if first line is header
        sequence = line.split('\t')[0]
        if sequence.upper() != 'CDR3B':
            tcr_sequences.add(sequence)
        for line in sequences.readlines():
            sequence = line.split('\t')[0]
            sequence = sequence.upper()
            tcr_sequences.add(sequence)
    return tcr_sequences


def load_filtered_tcr_sequences(sequence_file):
    tcr_sequences = set()
    with open(sequence_file, 'r') as sequences:
        for line in sequences.readlines():
            sequence = line.split('\t')[0]
            sequence = sequence.upper()
            #if re.match('Ĉ[AC-WY]{3,}F$', sequence):
            if re.match('^C[AC-WY][AC-WY][AC-WY][AC-WY]*F', sequence):
                tcr_sequences.add(sequence)
    return tcr_sequences


def load_reference_sequences(reference_file):
    cdr3b_sequences = set()
    with open(reference_file, 'r') as reference:
        for line in reference.readlines():
            line = line.split(';')
            cdr3b_sequences.add(line[-1])
    return list(cdr3b_sequences)


def count_kmers(sequences, use_structural_boundaries):
    # counting all 2-, 3- and 4-mers
    c = Counter()
    offset = 0 if use_structural_boundaries else 3
    for sequence in sequences:
        seqlen = len(sequence)
        for k in (2, 3, 4):
            kmers = [sequence[i:i+k] for i in range(offset, seqlen+1-k-offset)]
            c.update(kmers)
    return c


def global_clustering(tcr_sequences, use_structural_boundaries):
    max_distance = 1 if len(tcr_sequences) >= 125 else 2  # maximum hamming distance defined by gliph
    tcr_sequences_dict, sequences_original = sequences_to_dict(tcr_sequences, use_structural_boundaries)
    clusters = nx.Graph()
    clusters.add_nodes_from(tcr_sequences)
    for length in tcr_sequences_dict.keys():
        sequences = np.array(tcr_sequences_dict[length])  # list of sequences
        difference = sequences[:, np.newaxis] - sequences
        nonzeros = np.count_nonzero(difference, axis=2)
        edges = zip(*np.where(nonzeros <= max_distance))
        edges = [(sequences_original[length][x], sequences_original[length][y]) for x, y in edges]
        clusters.add_edges_from(edges)
    return clusters


def sequences_to_dict(tcr_sequences, use_structural_boundaries):
    # Sequences are added to dict by same length
    min_length = 8  # sequences must be at least 8 chars long
    sequences_reduced = dict()
    sequences_original = dict()
    for sequence in tcr_sequences:
        seq_length = len(sequence)
        if seq_length >= min_length:
            sequence_reduced = sequence[3:-3] if not use_structural_boundaries else sequence
            sequence_reduced = np.fromstring(sequence_reduced, dtype='uint8')
            if seq_length not in sequences_reduced:
                sequences_reduced[seq_length] = list()
                sequences_original[seq_length] = list()
            sequences_original[seq_length].append(sequence)
            sequences_reduced[seq_length].append(sequence_reduced)
    return sequences_reduced, sequences_original


def calculate_distance(sequence1, sequence2, use_structural_boundaries):
    # Hamming distance of two sequences
    if not use_structural_boundaries:
        sequence1 = sequence1[3:-3]
        sequence2 = sequence2[3:-3]
    sequence1 = np.fromstring(sequence1, dtype='uint8')
    sequence2 = np.fromstring(sequence2, dtype='uint8')
    distance = np.count_nonzero(sequence1-sequence2)
    return distance


def output_clusters(output, clusters_tcr):
    nx.write_gml(clusters_tcr, f'{output[:-4]}.gml')
    number_sequences = len(clusters_tcr.nodes())
    clusters_tcr = list(nx.connected_components(clusters_tcr))
    number_clusters = len([cluster for cluster in clusters_tcr if len(cluster) > 1])
    print(f'Clusters: {number_clusters}')
    print(f'Unique Sequences: {number_sequences}')
    with open(output, 'w') as output_file:
        for cluster in clusters_tcr:
            cluster_content = f'{len(cluster)}\t{list(cluster)[0]}\t'+' '.join(x for x in cluster)
            print(cluster_content, file=output_file)


def local_clustering(tcr_sequences, kmer_file):
    two_mers, three_mers, kmer_clusters = load_kmers(kmer_file)
    print('\tClustering kmers...')
    kmer_clusters = cluster_kmers(three_mers, kmer_clusters)
    kmer_clusters = cluster_kmers(two_mers, kmer_clusters)
    kmer_clusters = list(nx.connected_components(kmer_clusters))  # separate clusters
    print('\tClustering CDR3b sequences...')
    clusters_tcr = cluster_cdr3b(kmer_clusters, tcr_sequences)
    return clusters_tcr


def cluster_kmers(kmers, clusters):
    for kmer in kmers:
        for node in list(clusters.nodes):
            if kmer in node:
                clusters.add_edge(kmer, node)
            else:
                clusters.add_node(kmer)
    return clusters


def cluster_cdr3b(kmer_clusters, tcr_sequences):
    clusters_tcr = [nx.Graph() for _ in range(len(kmer_clusters))]
    # add sequences to clusters
    for sequence in tcr_sequences:
        sequence_seen = False
        for i, kmer_cluster in enumerate(kmer_clusters):
            for kmer in kmer_cluster:
                if kmer in sequence[3:-3]:
                    clusters_tcr[i].add_node(sequence)
                    sequence_seen = True
                    break
        # add sequence if not in cluster
        if not sequence_seen:
            new_cluster = nx.Graph()
            new_cluster.add_node(sequence)
            clusters_tcr.append(new_cluster)
    # connect nodes
    for i in range(len(clusters_tcr)):
        clusters_tcr[i].add_path(clusters_tcr[i].nodes())
    # join clusters
    clusters_tcr = nx.compose_all(clusters_tcr)
    return clusters_tcr


def union_clusters(clusters):
    joined_clusters = nx.Graph()
    for cluster in clusters:
        joined_clusters = nx.compose(joined_clusters, cluster)
    return joined_clusters


def load_kmers(input_file):
    print(f"Loading kmers from {input_file}")
    two_mers = []
    three_mers = []
    clusters = nx.Graph()
    with open(input_file, 'r') as kmer_file:
        kmer_file.readline()
        for line in kmer_file.readlines():
            kmer = line.split()[0]
            if len(kmer) == 2:
                two_mers.append(kmer)
            elif len(kmer) == 3:
                three_mers.append(kmer)
            elif len(kmer) == 4:
                clusters.add_node(kmer)
            else:
                print(kmer)
    return two_mers, three_mers, clusters


def argument_parser(parser):
    parser.add_argument('-t', '--tcr_sequences',
                        help='File holding TCRs',
                        type=str,
                        required=True)
    parser.add_argument('-r', '--reference',
                        help='Reference fasta file of naive CDR3 amino acid sequences used for simulation of sample sets. \
                        Example references are available on https://bitbucket.org/genomeinformatics/ting/',
                        type=str,
                        required=True)
    parser.add_argument('-k', '--kmer_file',
                        help='tab separated file holding kmers in first row',
                        type=str,
                        required=True)
    parser.add_argument('-o', '--output',
                        help='path of output-file',
                        type=str,
                        required=True)
    parser.add_argument('-b', '--use_structural_boundaries',
                        help='First and last three amino acids are included in processing',
                        action='store_true')
    parser.add_argument('-ng', '--no_global',
                        help='If set global clustering is excluded',
                        action='store_false')
    parser.add_argument('-nl', '--no_local',
                        help='If set local clustering is excluded',
                        action='store_false')
    parser.add_argument('-p', '--max_p_value',
                        help='p-value threshold for identifying significant motifs by fisher exact test',
                        type=float,
                        default=0.05)
    parser.add_argument('--gliph_minp',
                        help='probability threshold for identifying significant motifs by gliph test',
                        type=float,
                        default=0.001)
    parser.add_argument('-f', '--stringent_filtering',
                        help='If used only TCRs starting with a cystein and ending with phenylalanine will be used',
                        action='store_true')
    parser.add_argument('-c', '--cores',
                        help='Number of cores used for resampling. Default: 1',
                        type=int,
                        default=1)
    parser.add_argument('-g', '--kmers_gliph',
                        help='If set kmers are identified by the non-deterministic approach as implemented by gliph',
                        action='store_true')


if __name__ == '__main__':
    main()
