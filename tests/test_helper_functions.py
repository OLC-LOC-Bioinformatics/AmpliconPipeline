import os
import shutil
import pytest

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from bin.helper_functions import *

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

# Clear out sample_outdir before testing
outdir = os.path.join(parentdir, 'tests', 'sample_outdir')
try:
    shutil.rmtree(outdir)
except OSError:
    pass


def test_retrieve_fastqgz():
    assert len(retrieve_fastqgz('tests/sample_miseq')) == 2


def test_execute_command():
    (empty_test_out, empty_test_err) = execute_command('')
    assert (empty_test_out, empty_test_err) == ('', '')

    (echo_test_out, echo_test_err) = execute_command('echo test')
    assert (echo_test_out, echo_test_err) == ('test\n', '')


def test_retrieve_unique_sampleids():
    fastq_list = retrieve_fastqgz('tests/sample_miseq')
    assert len(retrieve_unique_sampleids(fastq_list)) == 1


def test_get_readpair():
    fastq_list = retrieve_fastqgz('tests/sample_miseq')
    sample_id_list = retrieve_unique_sampleids(fastq_list)
    for sample_id in sample_id_list:
        assert len(get_readpair(sample_id, fastq_list)) == 2


def test_project_setup():
    outdir = os.path.join(parentdir, 'tests', 'sample_outdir')
    inputdir = os.path.join(parentdir, 'tests', 'sample_miseq')

    # Test
    project_setup(outdir, inputdir)
    assert os.path.isdir(outdir)


def test_populate_sample_dictionary():
    id_list = retrieve_unique_sampleids(retrieve_fastqgz('tests/sample_miseq'))
    sample_dictionary = populate_sample_dictionary(id_list, retrieve_fastqgz('tests/sample_miseq'))
    assert len(sample_dictionary) == len(id_list)


def test_get_sample_dictionary():
    test_dict = get_sample_dictionary('tests/sample_miseq')
    assert len(test_dict) == 1
    assert '2017-SEQ-1114' in test_dict


def test_valid_olc_id():
    test_dict = get_sample_dictionary('tests/sample_miseq')
    for key, value in test_dict.items():
        assert valid_olc_id(key) is True


def test_create_sampledata_artifact():
    datadir = os.path.join(parentdir, 'tests', 'sample_miseq')
    qiimedir = os.path.join(parentdir, 'tests', 'sample_outdir', 'qiime2')
    path = create_sampledata_artifact(datadir, qiimedir)
    assert os.path.isfile(path)

