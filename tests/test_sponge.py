### Imports ###
import os
import pytest

import pandas as pd

### Fixtures ###
# Core config fixture
from sponge.config_manager import ConfigManager

@pytest.fixture
def core_config():
    return ConfigManager()

### Unit tests ###
# Data retrieval functions
import sponge.modules.utils.data_retrieval as data_f

@pytest.mark.network
@pytest.mark.parametrize('input, compare_to', [
    (('https://raw.githubusercontent.com/kuijjerlab/sponge/main/LICENSE',
        'LICENSE'), 'LICENSE'),
    (('https://raw.githubusercontent.com/kuijjerlab/sponge/main/LICENSE',
        None), 'LICENSE'),
])
def test_download_with_progress(input, compare_to, tmp_path):
    if input[1] == None:
        data = data_f.download_with_progress(*input).read().decode()
    else:
        file_path = os.path.join(tmp_path, input[1])
        data_f.download_with_progress(input[0], file_path)
        data = open(file_path, 'r').read()

    comp_data = open(compare_to, 'r').read()

    assert data == comp_data


@pytest.mark.parametrize('input', [
    ['test_dataset', ['field1', 'field2', 'field3']],
    ['test_dataset', []],
])
def test_create_xml_query(input):
    xml_query = data_f.create_xml_query(*input)

    assert xml_query[:38].lower() == "<?xml version='1.0' encoding='utf-8'?>"
    assert xml_query.count('Attribute') == len(input[1])


@pytest.mark.network
@pytest.mark.parametrize('input', [
    ['hsapiens_gene_ensembl', ['ensembl_transcript_id', 'ensembl_gene_id']],
    ['hsapiens_gene_ensembl', ['ensembl_transcript_id']],
])
def test_retrieve_ensembl_data(input, core_config):
    input.append(core_config['url']['region']['xml'])
    df = pd.read_csv(data_f.retrieve_ensembl_data(*input), sep='\t')

    assert type(df) == pd.DataFrame
    assert len(df.columns) == len(input[1])


@pytest.mark.network
def test_get_ensembl_version(core_config):
    version_string = data_f.get_ensembl_version(
        core_config['url']['region']['rest'])
    split_version = version_string.split('.')

    assert len(split_version) == 2
    assert split_version[0] == 'GRCh38'


@pytest.mark.network
@pytest.mark.parametrize('input, expected_type', [
    ('hg38', pd.Series),
    ('random_assembly', type(None)),
    ('hg1992', type(None)),
])
def test_get_chromosome_mapping(input, expected_type, core_config):
    mapping = data_f.get_chromosome_mapping(input,
        core_config['url']['chrom_mapping'])

    assert type(mapping) == expected_type


# ConfigReader class


# VersionLogger class


# FileRetriever class


# TFBSRetriever class


# RegionRetriever class


# DataRetriever class