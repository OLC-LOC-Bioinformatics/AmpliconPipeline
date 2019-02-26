import os
import glob
import click
import qiime2
import shutil
import pandas as pd


TAXONOMIC_DICT = {
    'kingdom': ('level-1', 'D_0__'),
    'phylum': ('level-2', 'D_1__'),
    'class': ('level-3', 'D_2__'),
    'order': ('level-4', 'D_3__'),
    'family': ('level-5', 'D_4__'),
    'genus': ('level-6', 'D_5__'),
    'species': ('level-7', 'D_6__'),
}


# NOTE: Note that the index_col is by default 'sample_annotation'. This assumes that the METADATA file used to generate
# this run has a column called 'sample_annotation' which acts as a secondary ID alongside the Seq-IDs provided.
def prepare_df(filepath, taxonomic_level, sample, index_col='sample_annotation', filtering=None, cutoff=None):
    df = pd.read_csv(filepath, index_col=index_col)

    # Remove all extraneous metadata columns
    df = df.drop((x for x in df.columns.tolist()
                  if (x.startswith('D_0__') is False) and
                  (x.startswith('Unassigned;') is False)),
                 axis=1)

    # Remove columns that don't have target filtering keyword; e.g. remove everything that isn't Bacteroidales
    if filtering is not None:
        df = df.drop((x for x in df.columns.tolist() if filtering not in x), axis=1)

    # Transpose
    df = df.transpose()
    df = df.reset_index()

    # Drop columns that aren't our sample of interest
    df = df.drop((x for x in df.columns.tolist() if sample not in x and 'index' not in x), axis=1)

    # Columns to target for conversion to percentage
    columns_to_target = [x for x in df.columns.tolist() if x not in ['index', taxonomic_level]]

    # Convert
    df = convert_to_percentages(df, columns_to_target)

    # Remove rows where the value for the sample is 0
    try:
        df = df[df.iloc[:, 1] != 0]
    except IndexError:
        print('The specified sample {} could not be found. Quitting.'.format(sample))
        quit()

    if cutoff is not None:
        # Remove rows where the value for the sample is < cutoff
        df = df[df.iloc[:, 1] >= cutoff]

    # Fix names of index
    df['index'] = df['index'].map(extract_taxonomy)

    return df


def convert_to_percentages(df, cols):
    """
    :param df:
    :param cols:
    :return:
    """
    df[cols] = df[cols].div(df[cols].sum(axis=0), axis=1).multiply(100)
    return df


def extract_csv_files(input_path, out_dir):
    """
    :param input_path:
    :param out_dir:
    :return:
    """
    # Load visualization file
    try:
        qzv = qiime2.Visualization.load(input_path)
    except:
        print('Could not load .qzv file. Quitting.')
        return None

    # Create temporary directory to dump contents into
    temp_dir = os.path.join(out_dir, 'temporary_qiime2')
    try:
        os.makedirs(temp_dir)
    except OSError:
        shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)

    # Export contents of QZV
    qzv.export_data(temp_dir)

    # Make list of files to keep
    all_files = glob.glob(os.path.join(temp_dir, 'level-*.csv'))
    files_to_move = []
    for file in all_files:
        for key, value in TAXONOMIC_DICT.items():
            if value[0]+'.' in os.path.basename(file):
                files_to_move.append(file)

    # Move .csv files of interest
    out_files = []
    for file in sorted(files_to_move):
        out_files.append(os.path.join(out_dir, os.path.basename(file)))
        shutil.move(file, os.path.join(out_dir, os.path.basename(file)))

    # Cleanup
    shutil.rmtree(temp_dir)

    return out_files


def extract_taxonomy(value):
    if 'Unassigned;_' in value:
        return 'Unassigned'

    junk_list = [
        'D_0__', 'D_1__', 'D_2__', 'D_3__',
        'D_4__', 'D_5__', 'D_6__', 'D_7__',
    ]

    # Initial cleanup
    for junk in junk_list:
        value = value.replace(junk, '')

    tax_string = value.split(';')[-1:][0]  # last item in list

    if tax_string == '':
        tax_list = value.split(';')
        tax_string = [x for x in reversed(tax_list) if x != ''][0] + ' (closest known classification)'
    elif tax_string == '__':
        tax_list = value.split(';')
        for value in tax_list:
            if value != '__':
                tax_string = value + ' (closest known classification)'
    elif 'uncultured' in tax_string.lower():
        tax_list = value.split(';')
        for value in tax_list:
            if 'uncultured' not in value:
                tax_string = value + ' (closest known classification)'

    return tax_string


@click.command()
@click.option('-i', '--input_file',
              type=click.Path(exists=True),
              required=True,
              help='Path to taxonomy barplot *.qzv file')
@click.option('-o', '--out_dir',
              type=click.Path(exists=True),
              required=True,
              help='Folder to save output file into')
@click.option('-s', '--sample',
              default=None,
              required=True,
              help='Sample name to prepare data for')
@click.option('-t', '--taxonomic_level',
              required=True,
              help='Taxonomic level to generate report for. Options: '
                   '["kingdom", "phylum", "class", "order", "family", "genus", "species"]')
@click.option('-c', '--cutoff',
              required=False,
              default=0.0,
              help='Filter dataset to a specified cutoff level. For example, setting this to 5.5 will only show '
                   'rows with values >= 5.5%')
def taxonomy_report_generator(input_file, out_dir, sample, taxonomic_level, cutoff):
    csv_files = extract_csv_files(input_file, out_dir)

    taxonomic_level = taxonomic_level.lower()

    target_file = None
    for file in csv_files:
        if TAXONOMIC_DICT[taxonomic_level][0] in file:
            target_file = file

    df = prepare_df(filepath=target_file, taxonomic_level=taxonomic_level, sample=sample, cutoff=cutoff)
    csv_out_path = os.path.join(out_dir, 'taxonomy_report_{}_{}.csv'.format(taxonomic_level, sample))
    df.to_csv(csv_out_path, index=False)


if __name__ == '__main__':
    taxonomy_report_generator()
