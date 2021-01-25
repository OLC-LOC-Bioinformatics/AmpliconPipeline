import os
import logging
import multiprocessing
import qiime2
import pandas as pd

from qiime2 import Metadata
from qiime2.plugins import feature_table, \
    dada2, \
    demux, \
    metadata, \
    alignment, \
    phylogeny, \
    diversity, \
    emperor, \
    feature_classifier, \
    taxa


def load_data_artifact(filepath):
    """
    :param filepath: Path to qiime2 artifact created with helper_functions.create_sampledata_artifact
    :return: QIIME2 object
    """
    data_artifact = qiime2.Artifact.load(filepath)
    logging.info('Loaded {}'.format(filepath))
    return data_artifact


def load_sample_metadata(filepath):
    """
    :param filepath: Path to the sample metadata file
    :return: QIIME2 metadata object
    """
    metadata_object = qiime2.Metadata.load(filepath)
    logging.info('Loaded {}'.format(filepath))
    return metadata_object


def visualize_metadata(base_dir, metadata_object):
    """
    :param base_dir: Main working directory filepath
    :param metadata_object: QIIME2 metadata object
    :return: QIIME2 metadata visualization object
    """
    # Path setup
    export_path = os.path.join(base_dir, 'sample-metadata-tabulate')

    # Prepare and save metadata visualization
    metadata_viz = metadata.visualizers.tabulate(metadata_object)
    metadata_viz.visualization.save(export_path)
    logging.info('Saved {}'.format(export_path))

    return metadata_viz


def visualize_demux(base_dir, data_artifact):
    """
    :param base_dir: Main working directory filepath
    :param data_artifact: QIIME2 data artifact object
    :return: QIIME2 demux visualization object
    """
    logging.info('Visualizing demux...')

    # Path setup
    export_path = os.path.join(base_dir, 'demux_summary.qzv')

    # Prepare and save demux summary visualization
    demux_viz = demux.visualizers.summarize(data=data_artifact)
    demux_viz.visualization.save(export_path)
    logging.info('Saved {}'.format(export_path))

    return demux_viz


def dada2_qc(base_dir, demultiplexed_seqs, trim_left_f, trim_left_r, trunc_len_f, trunc_len_r, max_ee=3,
             chimera_method='consensus', cpu_count=None):
    """
    :param base_dir: Main working directory filepath
    :param demultiplexed_seqs: QIIME2 object
    :param trim_left_f: Number of bases to trim from 5' of forward read
    :param trim_left_r: Number of bases to trim from 5' of reverse read
    :param trunc_len_f: Number of bases for forward read truncation
    :param trunc_len_r: Number of bases for reverse read truncation
    :param max_ee: number of errors allowed before rejecting a read
    :param chimera_method: Method for chimera detection
    :param cpu_count: Number of CPUs to use for DADA2
    :return: QIIME2/DADA2 filtered table and representative sequences objects
    """
    logging.info('Running DADA2 (this could take awhile)...')

    # Grab all CPUs if parameter is not specified
    if cpu_count is None:
        cpu_count = multiprocessing.cpu_count()
        logging.info('Set CPU count to {}'.format(cpu_count))

    logging.info('DADA2 trimming/truncation parameters:')
    logging.info('trim_left_f: {}'.format(trim_left_f))
    logging.info('trim_left_r: {}'.format(trim_left_r))
    logging.info('trunc_len_f: {}'.format(trunc_len_f))
    logging.info('trunc_len_r: {}'.format(trunc_len_r))

    # Run dada2
    (dada2_filtered_table, dada2_filtered_rep_seqs, denoising_stats) = dada2.methods.denoise_paired(
        demultiplexed_seqs=demultiplexed_seqs, trim_left_f=trim_left_f, trim_left_r=trim_left_r,
        trunc_len_f=trunc_len_f, trunc_len_r=trunc_len_r, max_ee=max_ee, chimera_method=chimera_method,
        n_threads=cpu_count)

    # Save artifacts
    dada2_filtered_table.save(os.path.join(base_dir, 'table-dada2.qza'))
    dada2_filtered_rep_seqs.save(os.path.join(base_dir, 'rep-seqs-dada2.qza'))
    denoising_stats.save(os.path.join(base_dir, 'denoising-stats.qza'))
    logging.info('Completed running DADA2')

    return dada2_filtered_table, dada2_filtered_rep_seqs, denoising_stats


def visualize_dada2(base_dir, dada2_filtered_table, dada2_filtered_rep_seqs, denoising_stats, metadata_object):
    """
    :param base_dir: Main working directory filepath
    :param dada2_filtered_table: DADA2 filtered table object
    :param dada2_filtered_rep_seqs: DADA2 representative sequences object
    :param metadata_object: QIIME2 metadata object
    :return: QIIME2 feature table summary object
    """
    logging.info('Visualizing DADA2 results...')

    # Prepare feature table
    feature_table_summary = feature_table.visualizers.summarize(table=dada2_filtered_table,
                                                                sample_metadata=metadata_object)

    # Prepare sequence table
    feature_table_seqs = feature_table.visualizers.tabulate_seqs(data=dada2_filtered_rep_seqs)

    # Prepare denoising stats table
    denoising_stats_table = feature_table.visualizers.tabulate_seqs(data=denoising_stats)

    # Path setup
    table_dada2_path = os.path.join(base_dir, 'table-dada2-summary.qzv')
    rep_seqs_path = os.path.join(base_dir, 'rep-seqs-summary.qzv')
    denoising_stats_path = os.path.join(base_dir, 'denoising-stats.qzv')

    # Save visualizations
    feature_table_summary.visualization.save(table_dada2_path)
    feature_table_seqs.visualization.save(rep_seqs_path)
    denoising_stats_table.visualization.save(denoising_stats_path)
    logging.info('Saved {}'.format(table_dada2_path))
    logging.info('Saved {}'.format(rep_seqs_path))
    logging.info('Saved {}'.format(denoising_stats_path))
    return feature_table_summary


def seq_alignment_mask(base_dir, dada2_filtered_rep_seqs, cpu_count=None):
    """
    :param base_dir: Main working directory filepath
    :param dada2_filtered_rep_seqs: DADA2 filtered representative sequence object
    :param cpu_count: Number of CPUs to use for analysis
    :return: QIIME2 sequence mask, sequence alignment objects
    """

    # Threading setup
    if cpu_count is None:
        cpu_count = multiprocessing.cpu_count()

    # Path setup
    aligned_export_path = os.path.join(base_dir, 'aligned-rep-seqs.qza')
    mask_export_path = os.path.join(base_dir, 'masked-aligned-rep-seqs.qza')

    # Perform and save sequence alignment
    logging.info('Running sequence alignment...')
    seq_alignment = alignment.methods.mafft(sequences=dada2_filtered_rep_seqs, n_threads=cpu_count)
    seq_alignment.alignment.save(aligned_export_path)
    logging.info('Saved {}'.format(aligned_export_path))

    # Perform and save alignment mask
    logging.info('Running alignment mask...')
    seq_mask = alignment.methods.mask(alignment=seq_alignment.alignment)
    seq_mask.masked_alignment.save(mask_export_path)
    logging.info('Saved {}'.format(mask_export_path))

    return seq_mask, seq_alignment


def phylo_tree(base_dir, seq_mask):
    """
    :param base_dir: Main working directory filepath
    :param seq_mask: QIIME2 sequence mask object
    :return: QIIME2 unrooted, rooted tree objects
    """
    # Path setup
    unrooted_export_path = os.path.join(base_dir, 'unrooted-tree.qza')
    rooted_export_path = os.path.join(base_dir, 'rooted-tree.qza')

    # Run and save unrooted tree
    logging.info('Generating unrooted tree...')
    phylo_unrooted_tree = phylogeny.methods.fasttree(alignment=seq_mask.masked_alignment)
    phylo_unrooted_tree.tree.save(unrooted_export_path)
    logging.info('Saved {}'.format(unrooted_export_path))

    # Run and save rooted tree
    logging.info('Generating rooted tree...')
    phylo_rooted_tree = phylogeny.methods.midpoint_root(tree=phylo_unrooted_tree.tree)
    phylo_rooted_tree.rooted_tree.save(rooted_export_path)
    logging.info('Saved {}'.format(rooted_export_path))

    return phylo_unrooted_tree, phylo_rooted_tree


def export_newick(base_dir, tree):
    """
    :param base_dir: Main working directory filepath
    :param tree: QIIME2 tree object
    :return: Path to tree file in newick format exported from the tree object
    """
    # Path setup
    export_path = os.path.join(base_dir, 'newick.tree')

    # Export data
    tree.rooted_tree.export_data(export_path)
    logging.info('Exported tree: {}'.format(export_path))
    return export_path


def load_artifact(artifact_path):
    """
    Generic loading of a QIIME2 artifact (.qza)

    :param artifact_path: Path to a QIIME .qza artifact
    :return: QIIME2 object
    """
    # Load existing artifact
    artifact = qiime2.Artifact.load(artifact_path)
    logging.info('Loaded {}'.format(artifact_path))
    return artifact


def calculate_maximum_depth(dada2_table):
    """
    Extracts the maximum observed read depth from post-filtering sequence object

    :param dada2_table: QIIME2 DADA2 table object
    :return: Maximum depth retrieved from the QIIME2/DADA2 table object
    """
    # Read in dada2 table
    df = dada2_table.view(pd.DataFrame)
    df = df.transpose()
    value_dict = {}

    # Calculate the sum of each column (sample)
    for column in df:
        value_dict[column] = df[column].sum()

    max_depth = max(value_dict.values())
    logging.info('Maximum depth in DADA2 table: {}'.format(str(max_depth)))
    return max_depth


def alpha_rarefaction_visualization(base_dir, dada2_filtered_table, max_depth=None):
    """
    Produces rarefaction visualization object

    :param base_dir: Main working directory filepath
    :param dada2_filtered_table: QIIME2 DADA2 filtered table object
    :param max_depth: Maximum depth value (integer)
    :return: QIIME2 alpha rarefaction visualization object
    """
    logging.info('Generating rarefaction curves...')

    # Path setup
    alpha_rarefaction_export_path = os.path.join(base_dir, 'alpha-rarefaction.qzv')

    # Max depth calculation. This (arbitrarily) sets the maximum depth to 80% of the highest value found.
    if max_depth is None:
        max_depth = int(calculate_maximum_depth(dada2_filtered_table) * 0.8)

    # Produce rarefaction curve
    alpha_rarefaction_viz = diversity.visualizers.alpha_rarefaction(table=dada2_filtered_table, max_depth=max_depth)

    # Save
    alpha_rarefaction_viz.visualization.save(alpha_rarefaction_export_path)
    logging.info('Saved {}'.format(alpha_rarefaction_export_path))

    return alpha_rarefaction_viz


def classify_taxonomy(base_dir, dada2_filtered_rep_seqs, classifier, cpu_count=None):
    """
    Uses a provided pre-trained classifier object to classify reads by taxonomy

    :param base_dir: Main working directory filepath
    :param dada2_filtered_rep_seqs: DADA2 filtered representative sequences object
    :param classifier: QIIME2 classifier object
    :param cpu_count: Number of CPUs to use for taxonomy classification
    :return: QIIME2 post-classification taxonomy object
    """
    logging.info('Classifying reads...')

    # Path setup
    export_path = os.path.join(base_dir, 'taxonomy.qza')

    # Threading setup
    if cpu_count is None:
        cpu_count = multiprocessing.cpu_count()
        logging.info('Set CPU count to {}'.format(cpu_count))

    # Classify reads
    taxonomy_analysis = feature_classifier.methods.classify_sklearn(reads=dada2_filtered_rep_seqs,
                                                                    classifier=classifier,
                                                                    n_jobs=cpu_count)
    # Save the resulting artifact
    taxonomy_analysis.classification.save(export_path)
    logging.info('Saved {}'.format(export_path))

    return taxonomy_analysis


def visualize_taxonomy(base_dir, metadata_object, taxonomy_analysis, dada2_filtered_table):
    """
    Generates .qzv visualization files (taxonomy_barplot, taxonomy) from a QIIME2 taxonomy object

    :param base_dir: Main working directory filepath
    :param metadata_object: QIIME2 metadata object
    :param taxonomy_analysis: QIIME2 taxonomy object
    :param dada2_filtered_table: DADA2 filtered table object
    :return: QIIME2 taxonomy metadata object
    """
    logging.info('Visualizing taxonomy...')

    # Path setup
    tax_export_path = os.path.join(base_dir, 'taxonomy.qzv')
    barplot_export_path = os.path.join(base_dir, 'taxonomy_barplot.qzv')

    # Load metadata
    taxonomy_metadata = taxonomy_analysis.classification.view(Metadata)

    # Create taxonomy visualization
    taxonomy_visualization = metadata.visualizers.tabulate(taxonomy_metadata)

    # Save taxonomy visualization
    taxonomy_visualization.visualization.save(tax_export_path)
    logging.info('Saved {}'.format(tax_export_path))

    # Create and save barplot visualization
    taxonomy_barplot = taxa.visualizers.barplot(table=dada2_filtered_table, taxonomy=taxonomy_analysis.classification,
                                                metadata=metadata_object)
    taxonomy_barplot.visualization.save(barplot_export_path)
    logging.info('Saved {}'.format(barplot_export_path))

    return taxonomy_metadata


def run_diversity_metrics(base_dir, dada2_filtered_table, phylo_rooted_tree, metadata_object,
                          sampling_depth=None, beta_column='sample_subsubtype'):
    """
    TODO: Allow beta_column (for beta diversity calculation) to be set through CLI

    :param base_dir: Main working directory filepath
    :param dada2_filtered_table: 
    :param phylo_rooted_tree: 
    :param metadata_object: 
    :param sampling_depth:
    :param beta_column: Column name to use for the beta group significance step of the pipeline
    :return: QIIME2 diversity core metrics object
    """
    logging.info('Running diversity metrics...')

    # Set sampling_depth to 10% of the maximum if no value is provided. Should probably rework this.
    if sampling_depth is None:
        sampling_depth = int(calculate_maximum_depth(dada2_filtered_table) * 0.1)

    # Path setup
    bray_curtis_path = os.path.join(base_dir, 'bray_curtis_emperor.qzv')
    jaccard_emperor_path = os.path.join(base_dir, 'jaccard_emperor.qzv')
    unweighted_unifrac_emperor_path = os.path.join(base_dir, 'unweighted_unifrac_emperor.qzv')
    weighted_unifrac_emperor_path = os.path.join(base_dir, 'weighted_unifrac_emperor.qzv')

    faith_visualization_path = os.path.join(base_dir, 'faith-pd-group-significance.qzv')
    evenness_visualization_path = os.path.join(base_dir, 'evenness-group-significance.qzv')
    beta_visualization_path = os.path.join(base_dir, 'unweighted-unifrac-sample-type-significance.qzv')

    # Retrieve diversity metrics
    diversity_metrics = diversity.pipelines.core_metrics_phylogenetic(table=dada2_filtered_table,
                                                                      phylogeny=phylo_rooted_tree.rooted_tree,
                                                                      sampling_depth=sampling_depth,
                                                                      metadata=metadata_object)

    # Save
    diversity_metrics.bray_curtis_emperor.save(bray_curtis_path)
    logging.info('Saved {}'.format(bray_curtis_path))

    diversity_metrics.jaccard_emperor.save(jaccard_emperor_path)
    logging.info('Saved {}'.format(jaccard_emperor_path))

    diversity_metrics.unweighted_unifrac_emperor.save(unweighted_unifrac_emperor_path)
    logging.info('Saved {}'.format(unweighted_unifrac_emperor_path))

    diversity_metrics.weighted_unifrac_emperor.save(weighted_unifrac_emperor_path)
    logging.info('Saved {}'.format(weighted_unifrac_emperor_path))

    # Alpha group significance
    alpha_group_faith = diversity.visualizers.alpha_group_significance(
        alpha_diversity=diversity_metrics.faith_pd_vector,
        metadata=metadata_object)

    alpha_group_evenness = diversity.visualizers.alpha_group_significance(
        alpha_diversity=diversity_metrics.evenness_vector,
        metadata=metadata_object)

    # Save
    alpha_group_faith.visualization.save(faith_visualization_path)
    logging.info('Saved {}'.format(faith_visualization_path))
    alpha_group_evenness.visualization.save(evenness_visualization_path)
    logging.info('Saved {}'.format(evenness_visualization_path))

    # Beta group significance
    try:
        beta_group = diversity.visualizers.beta_group_significance(
            distance_matrix=diversity_metrics.unweighted_unifrac_distance_matrix,
            metadata=metadata_object.get_column(beta_column),
            pairwise=True)
        beta_group.visualization.save(beta_visualization_path)
    except:
        logging.info('Could not calculate beta group significance with metadata feature {}\n'.format(beta_column))

    return diversity_metrics


# TODO: Implement this function to allow for training using custom primer sets/reference taxonomies
def train_feature_classifier(base_dir, otu_filepath, reference_taxonomy_filepath, f_primer=None, r_primer=None):
    """
    Trains a Naive Bayes classifier based on a reference database/taxonomy

    Primers for V3-V4 region:
    F: S-D-Bact-0341-b-S-17, 5′-CCTACGGGNGGCWGCAG-3′,
    R: S-D-Bact-0785-a-A-21, 5′-GACTACHVGGGTATCTAATCC-3
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592464/

    :param base_dir: Main working directory filepath
    :param otu_filepath: File path to reference OTU .qza file
    :param reference_taxonomy_filepath: File path to reference taxonomy .qza file
    :param f_primer: String containing forward primer sequence. Default V3-V4 regions.
    :param r_primer: String containing reverse primer sequence Default V3-V4 regions.
    :return: Returns the trained feature classifier
    """
    if f_primer is None and r_primer is None:
        f_primer = 'CCTACGGGNGGCWGCAG'  # V3-V4
        r_primer = 'GACTACHVGGGTATCTAATCC'  # V3-V4

    # Path setup
    ref_seqs_filepath = os.path.join(base_dir, 'ref-seqs.qza')

    otus = qiime2.Artifact.load(otu_filepath)
    ref_taxonomy = qiime2.Artifact.load(reference_taxonomy_filepath)
    reference_seqs = feature_classifier.methods.extract_reads(sequences=otus, f_primer=f_primer, r_primer=r_primer)
    reference_seqs.reads.save(ref_seqs_filepath)
    naive_bayes_classifier = feature_classifier.methods.fit_classifier_naive_bayes(reference_reads=reference_seqs.reads,
                                                                                   reference_taxonomy=ref_taxonomy)
    return naive_bayes_classifier


def run_qc_pipeline(base_dir, data_artifact_path, sample_metadata_path):
    """
    1. Loads qiime2 data artifact and sample metadata file
    2. Generates a visualization of pre-filtering sequence quality

    :param base_dir: Main working directory filepath
    :param data_artifact_path: Artifact generated via helper_functions.create_sampledata_artifact()
    :param sample_metadata_path: Path to sample metadata tsv file
    """
    # Load seed objects
    data_artifact = load_data_artifact(data_artifact_path)
    metadata_object = load_sample_metadata(sample_metadata_path)

    # Visualize metadata
    visualize_metadata(base_dir=base_dir, metadata_object=metadata_object)

    # Demux
    visualize_demux(base_dir=base_dir, data_artifact=data_artifact)


def read_metadata_df(sample_metadata_path):
    """
    :param sample_metadata_path: Path to .tsv metadata file
    :return: Pandas DataFrame of metadata
    """
    df = pd.read_csv(sample_metadata_path, delimiter='\t')
    return df


def validate_sample_id(sample_id):
    """
    This is a validation function to make sure _00 was properly appended to the sample IDs.
    Accomodates the QIIME2 import functionality via qiime tools import (CasavaOneEightSingleLanePerSampleDirFmt)
    Automatically corrects the IDs if necessary.

    :param sample_id: OLC Sample ID
    :return: Validated/corrected OLC Sample ID
    """
    if not sample_id.endswith('_00'):
        sample_id += '_00'
        logging.debug('Sample ID in metadata not correctly named --> corrected to: {}'.format(sample_id))
    return sample_id


def write_new_metadata(df, sample_metadata_path):
    """
    Creates a new metadata .tsv file where the Sample IDs have been validated

    :param df: Pandas DataFrame of the validated metadata
    :param sample_metadata_path: Path to the original metadata .tsv file
    :return: Path to the new metadata file
    """
    new_metadata_path = os.path.join(os.path.dirname(sample_metadata_path),
                                     os.path.basename(sample_metadata_path).replace('.tsv', '_Validated.tsv'))
    df.to_csv(new_metadata_path, sep='\t', index=None)
    return new_metadata_path


def validate_metadata(sample_metadata_path):
    """
    Validates the Sample IDs provided in the .tsv metadata file

    :param sample_metadata_path: Path to .tsv sample metadata file
    :return: Path to new validated .tsv sample metadata file
    """
    logging.info('Validating metadata file: {}'.format(sample_metadata_path))
    df = read_metadata_df(sample_metadata_path)
    df['#SampleID'] = df['#SampleID'].apply(validate_sample_id)  # Assumption that first column is the SampleID column
    new_metadata_path = write_new_metadata(df, sample_metadata_path)
    return new_metadata_path


def run_pipeline(base_dir, data_artifact_path, sample_metadata_path, classifier_artifact_path,
                 trim_left_f, trim_left_r, trunc_len_f, trunc_len_r, filtering_flag=False, include_chloroplast=False, include_mitochondria=False):
    """
    1. Load sequence data and sample metadata file into a QIIME 2 Artifact
    2. Filter, denoise reads with dada2
    3. Multiple sequence alignment and masking of highly variable regions
    4. Generate a phylogenetic tree
    5. Generate alpha rarefaction curves
    6. Conduct taxonomic analysis
    7. Generate taxonomy barplots
    8. Run diversity metrics

    :param base_dir: Main working directory filepath
    :param data_artifact_path: Artifact generated via helper_functions.create_sampledata_artifact()
    :param sample_metadata_path: Path to .tsv sample metadata file
    :param classifier_artifact_path: Path to the .qza classifer for assigning reads to taxonomy
    :param filtering_flag: Flag to determine which steps of the pipeline to execute
    :param trim_left_f: Number of bases to trim from 5' of forward read
    :param trim_left_r: Number of bases to trim from 5' of reverse read
    :param trunc_len_f: Number of bases for forward read truncation
    :param trunc_len_r: Number of bases for reverse read truncation
    :param include_chloroplast: If True, include reads that map to chloroplast in the analysis.
    :param include_mitochondria: If True, include reads that map to mitochondria in the analysis.
    """

    # Load seed object
    data_artifact = load_data_artifact(data_artifact_path)

    # Validate and correct metadata (currently only adds _00 if the SampleID doesn't end with it already)
    new_metadata_path = validate_metadata(sample_metadata_path)

    # Load metadata
    metadata_object = load_sample_metadata(new_metadata_path)

    # Visualize metadata
    visualize_metadata(base_dir=base_dir, metadata_object=metadata_object)

    # Demux
    visualize_demux(base_dir=base_dir, data_artifact=data_artifact)

    # Filter & denoise w/dada2
    (dada2_filtered_table, dada2_filtered_rep_seqs, denoising_stats) = dada2_qc(base_dir=base_dir, demultiplexed_seqs=data_artifact,
                                                               trim_left_f=trim_left_f, trim_left_r=trim_left_r,
                                                               trunc_len_f=trunc_len_f, trunc_len_r=trunc_len_r)
    # Visualize dada2
    visualize_dada2(base_dir=base_dir, dada2_filtered_table=dada2_filtered_table,
                    dada2_filtered_rep_seqs=dada2_filtered_rep_seqs, denoising_stats=denoising_stats, metadata_object=metadata_object)

    # Only do these steps if the filtering_flag is false
    if filtering_flag is False:
        # Mask and alignment
        (seq_mask, seq_alignment) = seq_alignment_mask(base_dir=base_dir,
                                                       dada2_filtered_rep_seqs=dada2_filtered_rep_seqs)

        # Phylogenetic tree
        (phylo_unrooted_tree, phylo_rooted_tree) = phylo_tree(base_dir=base_dir, seq_mask=seq_mask)

        # Export tree
        export_newick(base_dir=base_dir, tree=phylo_rooted_tree)

        # Load classifier
        classifier = load_artifact(artifact_path=classifier_artifact_path)

        # Produce rarefaction visualization
        alpha_rarefaction_visualization(base_dir=base_dir, dada2_filtered_table=dada2_filtered_table)

        # Run taxonomic analysis
        taxonomy_analysis = classify_taxonomy(base_dir=base_dir, dada2_filtered_rep_seqs=dada2_filtered_rep_seqs,
                                              classifier=classifier)
        # filter out chloroplasts and/or mitochondria
        # note that this step will also remove samples with zero hits
        if include_chloroplast is False or include_mitochondria is False:
            if include_chloroplast is False:
                if include_mitochondria is False:
                    filterflag = "mitochondria,chloroplast "
                else:
                    filterflag = "chloroplast "
            else:
                filterflag = "mitochondria "
            cmd = 'qiime taxa filter-table --i-table {input_table} --i-taxonomy {input_tax} ' \
                  '--p-exclude {filter_flag} ' \
                  '--o-filtered-table {output_table}'.format(input_table=os.path.join(base_dir, 'table-dada2.qza'),
                                                             input_tax=os.path.join(base_dir, 'taxonomy.qza'),
                                                             output_table=os.path.join(base_dir, 'table-filtered.qza'), filter_flag=filterflag)
            os.system(cmd)
            dada2_filtered_table = load_artifact(os.path.join(base_dir, 'table-filtered.qza'))

        # Visualize taxonomy
        visualize_taxonomy(base_dir=base_dir, metadata_object=metadata_object,
                           taxonomy_analysis=taxonomy_analysis, dada2_filtered_table=dada2_filtered_table)

        # Alpha and beta diversity
        # TODO: requires metadata object with some sort of sample information (sample type)
        run_diversity_metrics(base_dir=base_dir, dada2_filtered_table=dada2_filtered_table,
                              phylo_rooted_tree=phylo_rooted_tree, metadata_object=metadata_object)

        with open(os.path.join(base_dir, 'trimming_parameters.txt'), 'w') as f:
            f.write('The following trim parameters were passed to dada2 for this run.\n')
            f.write('Forward read, 5\' trim: {} bases\n'.format(trim_left_f))
            f.write('Forward read, 3\' trim: {} bases\n'.format(trunc_len_f))
            f.write('Reverse read, 5\' trim: {} bases\n'.format(trim_left_r))
            f.write('Reverse read, 3\' trim: {} bases\n'.format(trunc_len_r))
