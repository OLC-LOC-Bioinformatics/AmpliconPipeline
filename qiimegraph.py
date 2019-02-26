#!/usr/bin/env python3

from collections import OrderedDict

import os
import re
import click
import pickle
import qiime2
import shutil
import random

import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def extract_taxonomy(value):
    """
    :param value:
    :return:
    """
    if 'Unassigned;_' in value:
        return 'Unassigned'

    try:
        tax_string = value.split(TAXONOMIC_DICT[TAXONOMIC_LEVEL][1])[1]
        if tax_string == '':
            tax_string = value
    except:
        tax_string = value

    # Final cleanup. If there are remaining taxonomy characters it means a call to the specified level couldn't be made.
    if ';' in tax_string and '__' in tax_string:
        tax_string = 'Unclassified'

    return tax_string


def convert_to_percentages(df, cols):
    """
    :param df:
    :param cols:
    :return:
    """
    df[cols] = df[cols].div(df[cols].sum(axis=0), axis=1).multiply(100)
    return df


def get_column_pair(df, col1, col2):
    """
    :param df:
    :param col1:
    :param col2:
    :return:
    """
    # Filter columns
    df = df.filter([col1, col2], axis=1)

    # Sort
    df = df.sort_values(col1)
    return df


def prepare_df(filepath, index_col, filtering=None):
    """
    :param filepath:
    :param index_col:
    :param filtering:
    :return:
    """
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

    # Create taxonomic basename column
    df[TAXONOMIC_LEVEL] = df['index'].map(extract_taxonomy)

    # Columns to target for conversion to percentage
    columns_to_target = [x for x in df.columns.tolist() if x not in ['index', TAXONOMIC_LEVEL]]

    # Convert
    df = convert_to_percentages(df, columns_to_target)

    return df


def fixed_df(filename, index='sample_annotation', filtering=None):
    """
    :param filename:
    :param index:
    :param filtering:
    :return:
    """
    df = prepare_df(filepath=filename, index_col=index, filtering=filtering)
    new_filename = filename.replace('.csv', '_temp.csv')

    # Stupid hack
    df.to_csv(new_filename, index=None)
    df = pd.read_csv(new_filename, index_col=TAXONOMIC_LEVEL).fillna('NA')

    # Cleanup
    os.remove(new_filename)
    return df


def load_visualization(filepath):
    """
    :param filepath: path to qiime2 visualization
    :return: qiime2 object containing all information on viz
    """
    data_visualization = qiime2.Visualization.load(filepath)
    return data_visualization


def prepare_plot(df, sampleid):
    """
    :param df:
    :param sampleid:
    :return:
    """

    ordered_dict = df.to_dict(into=OrderedDict)[sampleid]
    ordered_dict = OrderedDict(sorted(ordered_dict.items(), key=lambda x: x[1], reverse=True))

    # Set explodes (i.e. separation from the other wedges). This will take the top 3 wedges and explode them.
    explode = [0 for x in range(len(ordered_dict))]

    explode[0] = 0.1

    try:
        explode[1] = 0.1
    except IndexError:
        pass

    try:
        explode[2] = 0.1
    except IndexError:
        pass

    explode = tuple(explode)

    values = []
    labels = []
    for key, value in ordered_dict.items():
        values.append(value)
        labels.append(key)

    return values, labels, explode


def style_wedges(wedges, colordict):
    """
    :param wedges:
    :param colordict:
    :return:
    """
    # Wedges
    for wedge in wedges:
        # wedge.set_color('black') # This puts outlines around the wedges.
        try:
            wedge.set_facecolor(colordict[wedge.get_label()])
        except:
            wedge.set_facecolor(colordict[random.sample(list(colordict), 1)[0]])


def generate_pct_labels(values, labels):
    """
    This function takes the labels (i.e. taxonomy) and values (i.e. 50%) and creates new labels
    with both pieces of information combined.
    :param values:
    :param labels:
    :return:
    """
    labels_values = zip(labels, values)
    pct_labels = []
    for label in labels_values:
        if label[0] != '':
            # This determined if there is an integer or a float or whatever in the label. Change %i to %.2f for floats.
            pct_labels.append(label[0] + ' (%i' % label[1] + '%)')
    return pct_labels


def paired_multi_pie_charts(samples, out_dir, filtering):
    """
    :param samples:
    :param out_dir:
    :param filtering:
    :return:
    """
    # Style setup
    plt.style.use('fivethirtyeight')

    # Consistent colouring across taxonomy e.g. Listeria will always be red
    colordict = read_color_pickle()

    # Font size
    mpl.rcParams['font.size'] = 8

    # Setup figure canvas
    plt.figure(figsize=(24, 14))

    # Regex setup -- grabs numeric value from string
    reg_pattern = re.compile(r'^\D*(\d+(?:\.\d+)?)\D*$')

    # Dictionary to store all labels for each sample
    pct_labels_dict = {}

    # Number of samples to figure out grid size
    n_samples = len(samples)

    # These counters control the positioning of pie charts on the grid
    x_counter = 1
    y_counter = 1

    # Create a pie chart for every sample
    for sample, attributes in samples.items():
        # print('Generating plot for {}'.format(sample))
        pct_labels_dict[sample] = generate_pct_labels(attributes[0], attributes[1])

        # Grid position for plot depending on n_samples
        if n_samples == 4:
            ax = plt.subplot2grid((4, 4), (y_counter, x_counter))
        elif n_samples == 6:
            ax = plt.subplot2grid((3, 4), (y_counter, x_counter))
        elif n_samples < 4:
            ax = plt.subplot2grid((3, 4), (y_counter, x_counter))

        # Create raw pie chart
        wedges, labels = ax.pie(attributes[0], labels=attributes[1], explode=attributes[2],
                                startangle=0, shadow=False)

        # Fix labels on the pie chart
        for label, pct_label in zip(labels, pct_labels_dict[sample]):

            # This will allow for the regex to properly extract the value
            pct_label_fix = pct_label
            for x in range(6):
                pct_label_fix = pct_label_fix.replace('D_{}__'.format(str(x)), '')

            # Grab percentage value from string
            try:
                pct_value = float(re.findall(reg_pattern, pct_label_fix)[0])
            except:
                pct_value = 0

            # Only show the label if it's value is >= 2%
            if pct_value >= 2:
                label.set_text(pct_label)
            else:
                label.set_text('')

        # Make the wedges look nice with a fixed color
        style_wedges(wedges=wedges, colordict=colordict)

        # Make sure pie chart is a circle
        ax.axis('equal')

        # Add title
        plt.title(sample)

        # Increment grid logic for pie chart placement
        if n_samples == 6:
            if x_counter == 3:
                x_counter = 1
                y_counter += 1
            else:
                x_counter += 1

            if y_counter == 3:
                y_counter = 1

        elif n_samples == 4:
            if x_counter == 2:
                x_counter = 1
                y_counter += 1
            else:
                x_counter += 1
        elif n_samples < 4:
            x_counter += 1

    # File naming
    prepend = ''
    for key in samples:
        prepend += key
        prepend += '_'

    # Set the filename
    if filtering is not None:
        outfile = os.path.join(out_dir, '{}_[{}]_{}_plot.png'.format(prepend, filtering, TAXONOMIC_LEVEL.capitalize()))
    else:
        outfile = os.path.join(out_dir, '{}_{}_plot.png'.format(prepend, TAXONOMIC_LEVEL.capitalize()))

    # Save the file (set transparent=True if you want to eliminate the background)
    plt.savefig(outfile, bbox_inches='tight', transparent=True)

    return outfile


def create_paired_pie_wrapper(filename, out_dir, samples, filtering):
    """
    :param filename:
    :param out_dir:
    :param samples:
    :param filtering:
    :return:
    """
    df = fixed_df(filename=filename, filtering=filtering)

    sample_dict = OrderedDict()
    for sample in samples:
        (values, labels, explode) = prepare_plot(df, sample)
        sample_dict[sample] = (values, labels, explode)

    filename = paired_multi_pie_charts(sample_dict, out_dir, filtering)

    return filename


def supress_autopct(pct):
    return ''


def my_autopct(pct):
    """
    :param pct:
    :return:
    """
    # Only return a label if it is > 2% else return an empty string
    return (('%.2f' % pct) + '%') if pct > 2 else ''


def get_spaced_colors(n):
    """
    :param n:
    :return: list of n colors in RGB
    """
    max_value = 16581375  # 255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return [((int(i[:2], 16)) / 255, (int(i[2:4], 16)) / 255, (int(i[4:], 16) / 255)) for i in colors]


def generate_color_pickle():
    """
    Takes the SILVA QIIME formatted taxonomy text file and assigns a unique colour to every single entry.
    Drops a pickled dictionary of {organism:color} associations onto the NAS for usage whenever necessary.
    :return:
    """
    taxonomy_file = '/mnt/nas/Databases/16S/Silva/qiime2/SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.txt'
    f = open(taxonomy_file, 'r')
    curated_tax = []

    # Extract text from each taxon
    levels_range = [x for x in range(15)]
    for line in f.readlines():
        # Strip ID
        line = line.split('\t', 1)[1]
        # Strip all taxonomic level delineations (i.e. D_0__)
        for level in levels_range:
            line = line.replace('D_{}__'.format(str(level)), '')
        # Cleanup
        all_tax = line.split(';')
        for tax in all_tax:
            tax = tax.strip()
            if tax is not '':
                curated_tax.append(tax)
        # Reduce list down to only unique entries
        curated_tax = list(set(curated_tax))

    # Associate a unique color with each entry
    colors = get_spaced_colors(len(curated_tax))
    colordict = {}
    for l, c in zip(curated_tax, colors):
        colordict[l] = c

    # Manual additions
    colordict['Unassigned'] = 'grey'
    colordict['Unclassified'] = 'grey'

    # Save the dictionary as a pickle for later reference
    pickle.dump(colordict,
                open('/mnt/nas/Redmine/QIIME2_CondaEnv/qiimegraph_taxonomic_color_dictionary_V2.pickle', 'wb'))


def read_color_pickle():
    colordict = pickle.load(open('./bin/qiimegraph_taxonomic_color_dictionary.pickle', 'rb'))
    return colordict


def extract_viz_csv(input_path, out_dir):
    """
    :param input_path:
    :param out_dir:
    :return:
    """
    # Load visualization file
    try:
        qzv = load_visualization(input_path)
    except:
        print('Could not load .qzv file. Quitting.')
        return None

    # Create temporary directory to dump contents into
    temp_dir = os.path.join(os.path.dirname(out_dir), 'temporary_qiime2_extraction')

    # Outfile path
    out_file = os.path.join(out_dir, 'qiime2_data_extract.csv')

    try:
        os.mkdir(temp_dir)
    except:
        shutil.rmtree(temp_dir)
        os.mkdir(temp_dir)

    # Grab CSV
    qzv.export_data(temp_dir)
    taxonomic_csv_path = os.path.join(temp_dir, TAXONOMIC_DICT[TAXONOMIC_LEVEL][0]+'.csv')

    # Move file
    os.rename(taxonomic_csv_path, out_file)

    # Cleanup
    shutil.rmtree(temp_dir)

    return out_file


@click.command()
@click.option('-i', '--input_file',
              type=click.Path(exists=True),
              required=True,
              help='CSV file exported from taxonomy_barplot visualization (*.qzv). '
                   'You can also just point to the *.qzv file, in which case the '
                   'taxonomy level specified will be exported. Defaults to family-level.')
@click.option('-o', '--out_dir',
              type=click.Path(exists=True),
              required=True,
              help='Folder to save output file into')
@click.option('-s', '--samples',
              default=None,
              help='List of samples to provide. Must be delimited by commas, e.g. -s SAMPLE1,SAMPLE2,SAMPLE3',
              required=True)
@click.option('-t', '--taxonomic_level',
              required=False,
              default="family",
              help='Taxonomic level to generate pie charts with. Defaults to "family". Options: '
                   '["kingdom", "phylum", "class", "order", "family", "genus", "species"]')
@click.option('-f', '--filtering',
              required=False,
              help='Filter dataset to a single group (e.g. Enterobacteriaceae)')
def cli(input_file, out_dir, samples, taxonomic_level, filtering):
    # generate_color_pickle()

    if samples is not None:
        samples = tuple(samples.split(','))

        if len(samples) > 9:
            print('ERROR: Cannot process more than 9 samples. Quitting.')
            quit()

    # Quick validation
    if not os.path.isdir(out_dir):
        click.echo('ERROR: Provided parameter to [-o, --out_dir] is not a valid directory. Try again.')
        quit()

    # Global variables. This is a hacky way of accomodating a few functions.
    global TAXONOMIC_LEVEL
    TAXONOMIC_LEVEL = taxonomic_level

    global TAXONOMIC_DICT
    TAXONOMIC_DICT = {
        'kingdom': ('level-1', 'D_0__'),
        'phylum': ('level-2', 'D_1__'),
        'class': ('level-3', 'D_2__'),
        'order': ('level-4', 'D_3__'),
        'family': ('level-5', 'D_4__'),
        'genus': ('level-6', 'D_5__'),
        'species': ('level-7', 'D_6__'),
    }

    filename = None

    # Input file handling
    if input_file.endswith('.csv'):
        filename = create_paired_pie_wrapper(input_file, out_dir, samples, filtering)
    elif input_file.endswith('.qzv'):
        input_file = extract_viz_csv(input_path=input_file, out_dir=out_dir)
        if input_file is None:
            quit()
        else:
            filename = create_paired_pie_wrapper(input_file, out_dir, samples, filtering)
    else:
        click.echo('ERROR: Invalid input_file provided. Please ensure file is .csv or .qzv.')
        quit()

    click.echo('Created chart at {} successfully'.format(filename))


if __name__ == '__main__':
    cli()
