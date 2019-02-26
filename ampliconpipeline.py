#!/usr/bin/env python3

import logging
import click
import os

from bin import helper_functions
from bin import qiime2_pipeline


@click.command()
@click.option('-i', '--inputdir',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing your raw MiSeq output (files must be *.fastq.gz)')
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              help='Base directory for all output from AmpliconPipeline. '
                   'Note that this directory must NOT already exist.')
@click.option('-m', '--metadata',
              type=click.Path(exists=True),
              required=True,
              help='Path to QIIME2 tab-separated metadata file. This must be a *.tsv file.')
@click.option('-c', '--classifier',
              type=click.Path(exists=True),
              default='./classifiers/99_V3V4_Silva_naive_bayes_classifier.qza',
              required=False,
              help='Path to a QIIME2 Classifier Artifact. By default this will point to a previously trained '
                   'V3-V4 classifier using SILVA taxonomy.')
@click.option('-f', '--filtering_flag',
              is_flag=True,
              default=False,
              help='Set flag to only proceed to the filtering step of analysis. This is useful for testing/optimizing '
                   'trimming parameters for a full run, or for generating files to be merged for later analysis.')
@click.option('-eq', '--evaluate_quality',
              is_flag=True,
              default=False,
              help='Setting this flag will only run the pipeline up until generating the demux_summary.qzv file. '
                   'This is important to do before running the pipeline to establish acceptable trimming/truncation '
                   'parameters to pass to dada2.')
@click.option('-tf', '--trim_left_f',
              default=10,
              help='Trim n bases from the 5\' end of the forward reads. Defaults to 10.')
@click.option('-tr', '--trim_left_r',
              default=5,
              help='Trim n bases from the 5\' end of the reverse reads. Defaults to 5.')
@click.option('-trf', '--trunc_len_f',
              default=280,
              help='Truncate the forward reads to n bases. Defaults to 280.')
@click.option('-trr', '--trunc_len_r',
              default=280,
              help='Truncate the reverse reads to n bases. Defaults to 280.')
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,
              help='Set this flag to enable more verbose output.')
@click.pass_context
def cli(ctx, inputdir, outdir, metadata, classifier, evaluate_quality, filtering_flag,
        trim_left_f, trim_left_r, trunc_len_f, trunc_len_r, verbose):
    # Logging setup
    if verbose:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')

    if evaluate_quality:
        logging.info('Starting QIIME2-QC Pipeline with output routing to {}'.format(outdir))
        data_artifact_path = helper_functions.project_setup(outdir=outdir, inputdir=inputdir)
        qiime2_pipeline.run_qc_pipeline(base_dir=os.path.join(outdir, 'qiime2'),
                                        data_artifact_path=data_artifact_path,
                                        sample_metadata_path=metadata)
        logging.info('QIIME2-QC Pipeline Completed')
        ctx.exit()

    # Input validation
    if os.path.isdir(outdir):
        click.echo(ctx.get_help(), err=True)
        click.echo('\nERROR: Specified output directory already exists. '
                   'Please provide a new path that does not already exist.', err=True)
        ctx.exit()

    # Classifier check
    if not os.path.isfile(classifier):
        click.echo(ctx.get_help(), err=True)
        click.echo('\nERROR: Classifier path is not valid. Please point to an existing classifier .qza file.', err=True)
        ctx.exit()
    else:
        logging.debug('Classifier path found at {}'.format(os.path.abspath(classifier)))

    # Project setup + get path to data artifact
    data_artifact_path = helper_functions.project_setup(outdir=outdir, inputdir=inputdir)

    # Filtering flag
    if filtering_flag:
        logging.info('FILTERING_FLAG SET. Pipeline will only proceed to DADA2 filtering step.')

    # Run the full pipeline
    qiime2_pipeline.run_pipeline(base_dir=os.path.join(outdir, 'qiime2'),
                                 data_artifact_path=data_artifact_path,
                                 sample_metadata_path=metadata,
                                 classifier_artifact_path=classifier,
                                 filtering_flag=filtering_flag,
                                 trim_left_f=trim_left_f, trim_left_r=trim_left_r,
                                 trunc_len_f=trunc_len_f, trunc_len_r=trunc_len_r)
    logging.info('QIIME2 Pipeline Completed')
    ctx.exit()


if __name__ == '__main__':
    cli()
