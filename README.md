# AmpliconPipeline

### Description

Based on: https://github.com/forestdussault/AmpliconPipeline

This program will take a 16S MiSeq amplicon sequencing run and
process the output with QIIME 2.

Currently, this is only  intended to be used internally within the OLC R&D division.

This program was built with the
[QIIME 2 Artifact API](https://docs.qiime2.org/2018.2/interfaces/artifact-api/).
Documentation on how to interact with various plugins
for QIIME 2 via Python can be found in the [QIIME 2 plugin docs](https://docs.qiime2.org/2018.2/plugins/).

### Installation Instructions

The AmpliconPipeline can be installed and run via conda. Follow the install instructions for 
QIIME2 2019.4 at https://docs.qiime2.org/2019.4/install/native/#install-qiime-2-within-a-conda-environment

1. Install AmpliconPipeline and activate the environment, get scikit-learn to the same version as the classifier.
```
git clone https://github.com/lowandrew/AmpliconPipeline.git
cd AmpliconPipeline
conda activate your_environment_name
pip install scikit-learn==0.19.1
```

2. Retrieve the classifier trained on the V3-V4 region and drop it into `.../AmpliconPipeline/classifiers`
```
mkdir classifiers
curl -L https://ndownloader.figshare.com/files/10970087 -o classifiers/99_V3V4_Silva_naive_bayes_classifier.qza
```

### Command Line Arguments
NOTE: `ampliconpipeline.py` must be used from the root of this directory.

```
Usage: ampliconpipeline.py [OPTIONS]

Options:
  -i, --inputdir PATH          Directory containing your raw MiSeq output
                               (files must be *.fastq.gz)  [required]
  -o, --outdir PATH            Base directory for all output from
                               AmpliconPipeline. Note that this directory must
                               NOT already exist.  [required]
  -m, --metadata PATH          Path to QIIME2 tab-separated metadata file.
                               This must be a *.tsv file.  [required]
  -c, --classifier PATH        Path to a QIIME2 Classifier Artifact. By
                               default this will point to a previously trained
                               V3-V4 classifier using SILVA taxonomy.
  -f, --filtering_flag         Set flag to only proceed to the filtering step
                               of analysis. This is useful for
                               testing/optimizing trimming parameters for a
                               full run, or for generating files to be merged
                               for later analysis.
  -eq, --evaluate_quality      Setting this flag will only run the pipeline up
                               until generating the demux_summary.qzv file.
                               This is important to do before running the
                               pipeline to establish acceptable
                               trimming/truncation parameters to pass to
                               dada2.
  -tf,  --trim_left_f INTEGER  Trim n bases from the 5' end of the forward
                               reads. Defaults to 10.
  -tr,  --trim_left_r INTEGER  Trim n bases from the 5' end of the reverse
                               reads. Defaults to 5.
  -trf, --trunc_len_f INTEGER  Truncate the forward reads to n bases. Defaults
                               to 280.
  -trr, --trunc_len_r INTEGER  Truncate the reverse reads to n bases. Defaults
                               to 280.
  -v, --verbose                Set this flag to enable more verbose output.
  --help                       Show this message and exit.
```

### Usage notes
#### Output
- Symlinks to the provided raw data will be created at `outdir/data`
- All QIIME 2 output will be available in `outdir/qiime2`

#### Metadata
A valid tab delimited metadata file must be provided to run `ampliconpipeline.py`.
An example of one using the standard OLC format is provided in the root
of this repository (**sample_metadata_example.tsv**).

The only fields that are
absolutely **required** are `#SampleID` (containing OLC Seq IDs) and `sample_annotation` (containing secondary IDs), though downstream analysis is
dependent on having detailed metadata available.

#### Flags & Output Details
There are three separate paths the pipeline can take depending on the
flag provided to `ampliconpipeline.py`. The relevant output file is listed at the end of each step.
1. Standard
    1. Load sequence data and sample metadata file into a QIIME 2 Artifact (`paired-sample-data.qza`)
    2. Filter, denoise reads with dada2 (`table-dada2-summary.qzv`)
    3. Multiple sequence alignment and masking of highly variable regions (`masked-aligned-rep-seqs.qza`)
    4. Generate a phylogenetic tree (`rooted-tree.qza`, `unrooted-tree.qza`)
    5. Generate alpha rarefaction curves (`alpha-rarefaction.qzv`)
    6. Conduct taxonomic analysis (`taxonomy.qzv`, `taxonomy.qza`)
    7. Generate taxonomy barplot (`taxonomy_barplot.qzv`)
    8. Run diversity metrics (`bray_curtis_emperor.qzv`,
    `evenness-group-significance.qzv`,
    `faith-pd-group-significance.qzv`,
    `jaccard_emperor.qzv`)
2. Evaluate quality (`-eq`, `--evaluate_quality`)
    1. Load sequence data and sample metadata file into a QIIME 2 Artifact (`paired-sample-data.qza`)
    2. Produce data quality visualization (`demux_summary.qzv`)
3. Filtering flag (`-f`, `--filtering_flag`)
    1. Load sequence data and sample metadata file into a QIIME 2 Artifact (`paired-sample-data.qza`)
    2. Filter, denoise reads with dada2 (`table-dada2-summary.qzv`)

#### Viewing data
All `.qzv` and `.qza` output files are viewable at https://view.qiime2.org/

An alternative pie chart view of the `taxonomic_barplot.qzv` files generated by the pipeline can be
created with `qiimegraph.py`. This has also been implemented on the [OLC Redmine](https://redmine.biodiversity.agr.gc.ca/projects/cfia/wiki#211-QIIMEGraph).

![QIIMEGraph](misc/qiimegraph_example.png?raw=true "QIIMEGraph")

**Note** that this script expects the metadata
file used to generate this .qzv visualization to contain a secondary
ID column called `sample_annotation`. If this column isn't present, the script won't work.
See the following for usage:
```
Usage: qiimegraph.py [OPTIONS]

Options:
  -i, --input_file PATH       CSV file exported from taxonomy_barplot
                              visualization (*.qzv). You can also just point
                              to the *.qzv file, in which case the taxonomy
                              level specified will be exported. Defaults to
                              family-level.  [required]
  -o, --out_dir PATH          Folder to save output file into  [required]
  -s, --samples TEXT          List of samples to provide. Must be delimited by
                              commas, e.g. -s SAMPLE1,SAMPLE2,SAMPLE3
                              [required]
  -t, --taxonomic_level TEXT  Taxonomic level to generate pie charts with.
                              Defaults to "family". Options: ["kingdom",
                              "phylum", "class", "order", "family", "genus",
                              "species"]
  -f, --filtering TEXT        Filter dataset to a single group (e.g.
                              Enterobacteriaceae)
  --help                      Show this message and exit.

```


#### Merging data
It is possible to merge previous runs together for combined analysis.
This also allows for filtering to samples of interest.

Steps:
  1. Loads and merges dada2 results from two previous analyses (requires
  rep-seqs-dada2 and table-dada2 for each run)

  2. Runs the full pipeline as usual with the merged run

  NOTE: A metadata file containing information for BOTH runs is required.

  Optionally, you may also provide a 'filtering list' which is a text file
  containing a sample ID on each new line which will only run the pipeline
  on those provided samples.

```
Usage: merge_runs.py [OPTIONS]

Options:
  -b, --base_dir PATH             Base directory for all output from
                                  QIIME2-MERGE Pipeline.  [required]
  -m, --sample_metadata_path PATH
                                  Path to QIIME2 tab-separated metadata file
                                  [required]
  -c, --classifier_artifact_path PATH
                                  Path to QIIME2 Classifier Artifact
  -t1, --table1_artifact_path PATH
                                  Path to first table artifact generated by
                                  DADA2 for merging  [required]
  -t2, --table2_artifact_path PATH
                                  Path to second table artifact generated by
                                  DADA2 for merging  [required]
  -rs1, --repseqs1_artifact_path PATH
                                  Path to first representative sequences
                                  artifact generated by DADA2 for merging
                                  [required]
  -rs2, --repseqs2_artifact_path PATH
                                  Path to second representative sequences
                                  artifact generated by DADA2 for merging
                                  [required]
  -f, --filtering_list PATH       Path to text file containing sample IDs that
                                  you wish to keep for the analysis
  --help                          Show this message and exit.
```

#### Classifier
By default, this pipeline uses a pre-trained classifier using the V3-V4 region.

The classifier and some additional details on how it was trained can be retrieved here:
https://figshare.com/articles/99_V3V4_Silva_naive_bayes_classifier_qza/6087197

### Tests and Example Data

Basic tests can be found in `tests/`.

To run these tests, use: `pytest -v tests/` from the root of the repository.