### Imports ###
import bioframe
import datetime
import os
import pytest
import yaml

import pandas as pd

from Bio.motifs.jaspar import Motif
from pathlib import Path
from pyjaspar import jaspardb
from typing import Tuple

from sponge.config_manager import ConfigManager

### Fixtures ###
# Core config fixture
@pytest.fixture
def core_config():
    return ConfigManager()

# A motif without any information
@pytest.fixture
def no_info_motif():
    no_info_row = [0.25] * 4
    no_info_counts = [no_info_row] * 6
    no_info_pwm = pd.DataFrame(no_info_counts, columns=['A', 'C', 'G', 'T'])
    no_info_motif = Motif(matrix_id='XXX', name='XXX', counts=no_info_pwm)

    yield no_info_motif

# A motif with perfect information
@pytest.fixture
def all_A_motif():
    all_A_row = [1] + [0] * 3
    all_A_counts = [all_A_row] * 6
    all_A_pwm = pd.DataFrame(all_A_counts, columns=['A', 'C', 'G', 'T'])
    all_A_motif = Motif(matrix_id='XXX', name='XXX', counts=all_A_pwm)

    yield all_A_motif

# A real motif for SOX2
@pytest.fixture
def SOX2_motif():
    jdb_obj = jaspardb(release='JASPAR2024')
    SOX2_motif = jdb_obj.fetch_motif_by_id('MA0143.1')

    yield SOX2_motif

# A subset of promoters on chromosome 19
@pytest.fixture
def chr19_promoters():
    path_to_file = os.path.join('tests', 'sponge', 'chr19_subset.tsv')
    df = bioframe.read_table(path_to_file, header=0)
    df.set_index('Transcript stable ID', inplace=True)

    yield df

# Part of the FOXF2 track for chromosome 19
@pytest.fixture
def foxf2_chr19():
    path_to_file = os.path.join('tests', 'sponge', 'foxf2_chr19_subset.tsv')
    df = pd.read_csv(path_to_file, sep='\t')

    yield df

### Unit tests ###
# Analysis functions
import sponge.modules.analysis as anal_f

@pytest.mark.parametrize('input, n_tfs, n_genes, n_edges', [
    (os.path.join('tests', 'sponge', 'comp_motif_prior_1.tsv'), 3, 3, 5),
    (os.path.join('tests', 'sponge', 'comp_motif_prior_2.tsv'), 4, 4, 5),
])
def test_load_prior(input, n_tfs, n_genes, n_edges):
    prior_df = anal_f.load_prior(input)

    assert prior_df['tf'].nunique() == n_tfs
    assert prior_df['gene'].nunique() == n_genes
    assert len(prior_df) == n_edges


@pytest.mark.parametrize('input, n_tfs, n_genes, n_edges', [
    (os.path.join('tests', 'sponge', 'comp_motif_prior_1.tsv'), 3, 3, 5),
    (os.path.join('tests', 'sponge', 'comp_motif_prior_2.tsv'), 4, 4, 5),
])
def test_describe_prior(input, n_tfs, n_genes, n_edges, capsys):
    prior_df = anal_f.load_prior(input)
    anal_f.describe_prior(prior_df)

    captured = capsys.readouterr()
    lines = captured.out.splitlines()

    assert lines[0] == f'Number of unique TFs: {n_tfs}'
    assert lines[1] == f'Number of unique genes: {n_genes}'
    assert lines[2] == f'Number of edges: {n_edges}'


def test_plot_confusion_matrix():
    df_1 = anal_f.load_prior(os.path.join('tests', 'sponge',
        'comp_motif_prior_1.tsv'))
    df_2 = anal_f.load_prior(os.path.join('tests', 'sponge',
        'comp_motif_prior_2.tsv'))

    common_tfs = set(df_1['tf'].unique()).intersection(
        df_2['tf'].unique())
    common_genes = set(df_1['gene'].unique()).intersection(
        df_2['gene'].unique())

    common_index = pd.MultiIndex.from_product([sorted(common_tfs),
        sorted(common_genes)])
    prior_1_mod = df_1.set_index(['tf', 'gene']).reindex(
        common_index, fill_value=0)
    prior_2_mod = df_2.set_index(['tf', 'gene']).reindex(
        common_index, fill_value=0)
    comp_df = prior_1_mod.join(prior_2_mod, lsuffix='_1', rsuffix='_2')

    cm = anal_f.confusion_matrix(comp_df['edge_1'], comp_df['edge_2'])

    ax = anal_f.plot_confusion_matrix(cm)

    assert type(ax.figure) == anal_f.plt.Figure


def test_compare_priors(capsys):
    df_1 = anal_f.load_prior(os.path.join('tests', 'sponge',
        'comp_motif_prior_1.tsv'))
    df_2 = anal_f.load_prior(os.path.join('tests', 'sponge',
        'comp_motif_prior_2.tsv'))

    _ = anal_f.compare_priors(df_1, df_2)

    captured = capsys.readouterr()
    lines = captured.out.splitlines()

    assert lines[12] == 'Number of common TFs: 3'
    assert lines[13] == 'Number of common genes: 3'

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

# Dictionary update functions
import sponge.modules.utils.dictionary_update as dict_f

@pytest.mark.parametrize('input, expected_output', [
    (({}, {'a': {'b': 1}}), {'a': {'b': 1}}),
    (({'a': 1, 'b': {'c': 2}}, {'b': {'d': 3}}), 
        {'a': 1, 'b': {'c': 2, 'd': 3}}),
])
def test_recursive_update(input, expected_output):
    assert dict_f.recursive_update(*input) == expected_output

# JASPAR versioning functions
import sponge.modules.utils.jaspar_versioning as jaspar_f

@pytest.mark.parametrize('input, expected_output', [
    (None, 'JASPAR2024'),
    ('JASPAR2022', 'JASPAR2022'),
    ('2024', 'JASPAR2024'),
])
def test_process_jaspar_version(input, expected_output):
    assert jaspar_f.process_jaspar_version(input) == expected_output

# Motif information functions
import sponge.modules.utils.motif_information as motif_f

@pytest.mark.parametrize('input, expected_output', [
    (0, 0),
    (0.5, -0.5),
    (1, 0),
])
def test_plogp(input, expected_output):
    assert motif_f.plogp(input) == expected_output


def test_calculate_ic_no_info(no_info_motif):
    assert motif_f.calculate_ic(no_info_motif) == 0


def test_calculate_ic_all_the_same(all_A_motif):
    # Length of the test motif is 6, so expected value is 2 * 6 = 12
    assert motif_f.calculate_ic(all_A_motif) == 12


def test_calculate_ic_SOX2(SOX2_motif):
    assert (motif_f.calculate_ic(SOX2_motif) ==
        pytest.approx(12.95, abs=0.01))

# String manipulation functions
import sponge.modules.utils.string_manipulation as string_f

@pytest.mark.parametrize('input, expected_output', [
    ('CAB', 'Cab'),
    ('SOX2', 'SOx2'),
    ('ARHGAP21', 'ARHGAP21'),
    ('ABC2DE', 'ABC2de'),
    ('ABCDE::FGHIJ', 'ABCde::FGHij'),
])
def test_adjust_gene_name(input, expected_output):
    assert string_f.adjust_gene_name(input) == expected_output


@pytest.mark.parametrize('input, expected_output', [
    ('a_string', 'a_string'),
    (datetime.datetime(1992, 5, 29, 23, 15), '29/05/1992, 23:15:00'),
])
def test_parse_datetime(input, expected_output):
    assert string_f.parse_datetime(input) == expected_output

# TFBS filtering functions
import sponge.modules.utils.tfbs_filtering as filter_f

@pytest.mark.parametrize('input, expected_length', [
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'),
        'chr19', ['MA0036.4', 'MA0030.2', 'MA0147.4'], 0, 2_000_000), 62),
])
def test_filter_edges(input, expected_length, chr19_promoters):
    df = filter_f.filter_edges(input[0], chr19_promoters, *input[1:])

    assert type(df) == pd.DataFrame
    assert len(df) == expected_length


@pytest.mark.parametrize('input, expected_length', [
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'),
        ['chr1', 'chr19'], ['MA0036.4', 'MA0030.2', 'MA0147.4']), 62),
])
def test_iterate_chromosomes(input, expected_length, chr19_promoters):
    df_list = filter_f.iterate_chromosomes(input[0], chr19_promoters,
        *input[1:])

    assert sum(len(df) for df in df_list) == expected_length


def test_process_chromosome(chr19_promoters, foxf2_chr19):
    df = filter_f.process_chromosome(foxf2_chr19, chr19_promoters)

    assert len(df) == 38


def test_process_motif(chr19_promoters, foxf2_chr19):
    df = filter_f.process_motif(foxf2_chr19, chr19_promoters)

    assert len(df) == 38


@pytest.mark.network
@pytest.mark.parametrize('input, expected_length', [
    ((['chr1', 'chr19'], ['MA0030.2'], ['FOXF2'], 'hg38', 'JASPAR2024'), 38),
])
def test_iterate_motifs(input, expected_length, core_config, chr19_promoters):
    df_list = filter_f.iterate_motifs(core_config['url']['motif']['by_tf'],
        chr19_promoters, *input)

    assert sum(len(df) for df in df_list) == expected_length

# TODO: Add tests for individual classes
# ConfigReader class


# VersionLogger class


# FileRetriever class


# TFBSRetriever class


# RegionRetriever class


# DataRetriever class


# ProteinIDMapper class


# JasparRetriever class


# HomologyRetriever class


# MotifSelector class


# MatchFilter class


# PPIRetriever class


# MatchAggregator class


# FileWriter class


### Integration tests ###
from sponge.sponge import Sponge

def run_integration_test_common(
    tmp_path: Path,
    config_file: Path,
) -> Tuple[Path, Path]:
    """
    Run the common part of the integration tests, which includes
    modifying the output file paths to be in the temporary directory and
    checking they exist after running SPONGE.

    Parameters
    ----------
    tmp_path : Path
        Path to the temporary directory generated by pytest
    config_file : Path
        Path to the config file

    Returns
    -------
    Tuple[Path, Path]
        Paths to the generated motif and PPI priors
    """

    motif_output = os.path.join(tmp_path, 'motif_prior.tsv')
    ppi_output = os.path.join(tmp_path, 'ppi_prior.tsv')

    settings = yaml.safe_load(open(config_file, 'r'))
    settings['motif_output']['file_name'] = motif_output
    settings['ppi_output']['file_name'] = ppi_output

    # Using the default user config file
    _ = Sponge(
        config=settings,
        temp_folder=os.path.join(tmp_path, '.sponge_temp'),
    )

    assert os.path.exists(motif_output)
    assert os.path.exists(ppi_output)

    return (motif_output, ppi_output)


# The test is marked as slow because the download of the bigbed file takes
# a lot of time and the filtering is also time consuming unless parallelised
@pytest.mark.integration
@pytest.mark.network
@pytest.mark.slow
def test_full_default_workflow(tmp_path):
    _,_ = run_integration_test_common(
        tmp_path,
        # Default config file
        os.path.join('sponge', 'user_config.yaml'),
    )


@pytest.mark.integration
@pytest.mark.network
def test_small_workflow(tmp_path):
    motif_output,ppi_output = run_integration_test_common(
        tmp_path,
        os.path.join('tests', 'sponge', 'test_user_config.yaml'),
    )

    motif_df = pd.read_csv(motif_output, sep='\t', header=None)
    motif_df_t = pd.read_csv(os.path.join('tests', 'sponge',
        'test_motif_prior.tsv'), sep='\t', header=None)

    pd.testing.assert_frame_equal(motif_df, motif_df_t)

    ppi_df = pd.read_csv(ppi_output, sep='\t', header=None)
    ppi_df_t = pd.read_csv(os.path.join('tests', 'sponge',
        'test_ppi_prior.tsv'), sep='\t', header=None)

    pd.testing.assert_frame_equal(ppi_df, ppi_df_t)