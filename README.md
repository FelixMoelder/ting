# ting

ting - T cell receptor interaction grouping

Synopsis
--------

ting [options]

Options
-------

Required Input
~~~~~~~~~~~~~~

The user must provide a list of CDR3b sequences.
For compatibility reasons the tab seperated table of TCR sequences required for gliph is supported, too.

    --tcr_sequences tcr_sequences   The format of the table is tab delimited, expecting only the first
                                    column. The header is optional, but if included only use column
                                    names as shown in the example.

    --kmer_file K-MER_FILE          The k-mer file holds all 2-, 3- and 4-mers considered for local
                                    clustering. If file does not exist it will automatically be
                                    generated.

    --reference                     Reference file used for simulation of sample sets.

Example:

CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	Sample-ID
CAADTSSGANVLTF	TRBV30	TRBJ2-6	CALSDEDTGRRALTF	TRAV19		TRAJ5	09/02171
CAATGGDRAYEQYF	TRBV2	TRBJ2-7	CAASSGANSKLTF	TRAV13-1	TRAJ56	03/04922
CAATQQGETQYF	TRBV2	TRBJ2-5	CAASYGGSARQLTF	TRAV13-1	TRAJ22	02/02591
CACVSNTEAFF	TRBV28	TRBJ1-1	CAGDLNGAGSYQLTF	TRAV25		TRAJ28	PBMC8631
CAGGKGNSPLHF	TRBV2	TRBJ1-6	CVVLRGGSQGNLIF	TRAV12-1	TRAJ42	02/02071
CAGQILAGSDTQYF	TRBV6-4	TRBJ2-3	CATASGNTPLVF	TRAV17		TRAJ29	09/00181
CAGRTGVSTDTQYF	TRBV5-1	TRBJ2-3	CAVTPGGGADGLTF	TRAV41		TRAJ45	02/02591
CAGYTGRANYGYTF	TRBV2	TRBJ1-2	CVVNGGFGNVLHC	TRAV12-1	TRAJ35	01/08733

Optional Input



    --use_structural_boundaries     If set, the first and last three amino acids will be included
                                    in kmer counting and global clustering.

    --no_global                     No global clustering will be performed.

    --no_local                      No local clustering will be performed.

    --min_kmer_occurence            Only kmers which occure at least min_kmer_occurences times in the
                                    sequence sample set will be taken in account. Default is 3.
    
    --max_p_value                   p-value threshold for identifying significant motifs by fisher exact test
    
    --gliph_minp                    probability threshold for identifying significant motifs by gliph test

    --stringent_filtering           Only TCRs starting with a cystein and ending with phenylalanine will be
                                    used (IGMT definition of CDR3 region). Default: False
                                    
    --kmers_gliph                   If set kmers are identified by the non-deterministic approach as implemented by gliph

