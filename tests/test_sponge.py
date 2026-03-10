### Imports ###
import bioframe
import datetime
import os
import pytest
import yaml

import pandas as pd

from Bio.motifs.jaspar import Motif
from pathlib import Path
from pyjaspar import jaspardb, JASPAR_LATEST_RELEASE
from typing import Any, Iterable, Tuple

from sponge.config_manager import ConfigManager
from sponge.modules.match_aggregator import MatchAggregator
from sponge.modules.ppi_retriever import PPIRetriever

### Fixtures ###
# Core config fixture
@pytest.fixture
def core_config():
    yield ConfigManager()

# Default user config fixture
@pytest.fixture
def default_user_config():
    yield ConfigManager({})

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
    df = pd.read_table(path_to_file)

    yield df

# PPI retriever instance
@pytest.fixture
def ppi_retriever(core_config):
    yield PPIRetriever(
        core_config['url']['ppi'],
        core_config['url']['protein'],
    )

# Match aggregator instance with initial edges
@pytest.fixture
def match_aggregator():
    yield MatchAggregator(
        pd.DataFrame(
            [
                ['aaa', 1000, 'MA0001.1', 'ENST00000530711'],
                ['aaa', 900, 'MA0001.1', 'ENST00000530711'],
                ['aaa', 800, 'MA0001.1', 'ENST00000530711'],
                ['BBB', 400, 'MA0002.1', 'ENST00000530711'],
                ['CCC', 600, 'MA0003.1', 'ENST00000585993'],
                ['DDD', 700, 'MA0004.1' ,'ENST00000585993'],
                ['BBB', 500, 'MA0002.1', 'ENST00000633742'],
                ['CCC', 600, 'MA0003.1', 'ENST00000633742'],
                ['BBB', 550, 'MA0002.1', 'ENST00000633703'],
            ],
            columns=['TFName', 'score', 'name', 'transcript'],
        ),
        os.path.join('tests', 'sponge', 'chr19_subset.tsv'),
        {'aaa': 'AAA', 'xxx': 'XXX'},
    )

# A mock prior frame
@pytest.fixture
def prior_frame():
    yield pd.DataFrame(
        [
            ['SOX2', 'FOXF1', 'BRCA1', 1, 4.0],
            ['SOX2', 'SOX2', 'MHCA2', 1, 0.4],
            ['ABC', 'DEF', 'XYZ', 0, 0.0],
            ['AAA', 'AAA', 'AAA', 4, 3.14],
        ],
        columns=['tf1', 'tf2', 'gene', 'edge', 'score'],
    )

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
    ('test_dataset', ['field1', 'field2', 'field3']),
    ('test_dataset', []),
    ('test_dataset', [], {}),
    ('test_dataset', [], {'field1': ['val1', 'val2', 'val3']}),
])
def test_create_xml_query(input):
    xml_query = data_f.create_xml_query(*input)

    assert xml_query[:38].lower() == "<?xml version='1.0' encoding='utf-8'?>"
    assert xml_query.count('Attribute') == len(input[1])


@pytest.mark.network
@pytest.mark.parametrize('dataset, fields, filters', [
    ('hsapiens_gene_ensembl', ['ensembl_transcript_id', 'ensembl_gene_id'],
        None),
    ('hsapiens_gene_ensembl', ['ensembl_transcript_id'], None),
    ('hsapiens_gene_ensembl', ['ensembl_transcript_id'],
        {'chromosome_name': ['chr19', 'chrX']}),
])
def test_retrieve_ensembl_data(dataset, fields, filters, core_config):
    ensembl_url = core_config['url']['region']['xml']
    df = pd.read_table(data_f.retrieve_ensembl_data(dataset, fields,
        ensembl_url, filters))

    assert type(df) == pd.DataFrame
    assert len(df.columns) == len(fields)


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
    (None, JASPAR_LATEST_RELEASE),
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
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'), 'chr19',
        ['MA0036.4', 'MA0030.2', 'MA0147.4'], 0, 2_000_000), 62),
])
def test_filter_edges(input, expected_length, chr19_promoters):
    df = filter_f.filter_edges(input[0], chr19_promoters, *input[1:])

    assert type(df) == pd.DataFrame
    assert len(df) == expected_length


@pytest.mark.parametrize('input, expected_length', [
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'), ['chr1', 'chr19'],
        ['MA0036.4', 'MA0030.2', 'MA0147.4']), 62),
])
def test_iterate_chromosomes(input, expected_length, chr19_promoters):
    df_list = filter_f.iterate_chromosomes(input[0], chr19_promoters,
        *input[1:])

    assert sum(len(df) for df in df_list) == expected_length


def test_process_chromosome(chr19_promoters, foxf2_chr19):
    df = filter_f.process_chromosome(foxf2_chr19, chr19_promoters)

    assert len(df) == 128


def test_process_motif(chr19_promoters, foxf2_chr19):
    df = filter_f.process_motif(foxf2_chr19, chr19_promoters)

    assert len(df) == 128


@pytest.mark.network
@pytest.mark.parametrize('input, expected_length', [
    ((['chr1', 'chr19'], ['MA0030.2'], ['FOXF2'], 'hg38', 'JASPAR2024'), 38),
])
def test_iterate_motifs(input, expected_length, core_config, chr19_promoters):
    df_list = filter_f.iterate_motifs(core_config['url']['motif']['by_tf'],
        chr19_promoters, *input)

    assert sum(len(df) for df in df_list) == expected_length

### Class unit tests ###
def call_class_methods(
    class_object: Any,
    input: Iterable[Tuple[str, Iterable[Any]]],
) -> None:
    """
    Calls the specified methods with given arguments on a class object.

    Parameters
    ----------
    class_object : Any
        Class object on which methods should be called
    input : Iterable[Tuple[str, Iterable[Any]]]
        Iterable with method names and their arguments
    """

    for func,args in input:
        f_call = getattr(class_object, func)
        f_call(*args)

# ConfigManager class
@pytest.mark.parametrize('c_path, input, expected_length', [
    (None, (
    ), 4),
    (os.path.join('src', 'sponge', 'user_config.yaml'), (
    ), 9),
    (os.path.join('src', 'sponge', 'user_config.yaml'), (
        ('get_value', ('genome_assembly',)),
        ('set_value', ('X', 2)),
        ('deep_update', ({'Y': {'Z': 'def'}},)),
        ('set_value', (['Y', 'Z'], 'abc')),
    ), 11),
])
def test_config_manager(c_path, input, expected_length, tmp_path):
    config_manager = ConfigManager(c_path, tmp_path)
    call_class_methods(config_manager, input)
    del config_manager
    config_file = os.path.join(tmp_path, 'last_user_config.yaml')

    assert os.path.exists(config_file)

    data = yaml.safe_load(open(config_file, 'r', encoding='utf-8'))

    assert len(data) == expected_length

# VersionLogger class
from sponge.modules.version_logger import VersionLogger

@pytest.mark.parametrize('input, expected_length', [
    ((
        ('write_provided', ('X',)),
    ), 1),
    ((
        ('write_default', ('X',)),
        ('write_provided', ('X',)),
        ('update_cached', ('X',)),
        ('write_retrieved', ('X', '1.0')),
    ), 1),
    ((
        ('write_default', ('X',)),
        ('write_provided', ('Y',)),
        ('write_retrieved', ('Z', '1.0')),
    ), 3),
])
def test_version_logger(input, expected_length, tmp_path):
    version_logger = VersionLogger(tmp_path)
    call_class_methods(version_logger, input)
    del version_logger  # Need to delete both references
    fp_file = os.path.join(tmp_path, 'fingerprint.yaml')

    assert os.path.exists(fp_file)

    data = yaml.safe_load(open(fp_file, 'r', encoding='utf-8'))

    assert len(data) == expected_length

# FileRetriever class
from sponge.modules.data_retriever.file_retriever import FileRetriever

@pytest.mark.parametrize('key, tmp_file, file_path', [
    ('A', 'tmp.file', os.path.join('src', 'sponge', 'user_config.yaml')),
    ('B', 'tmp.file', None),
    ('C', 'preexisting.file', None),
])
def test_file_retriever(key, tmp_file, file_path, tmp_path):
    Path(os.path.join(tmp_path, 'preexisting.file')).touch()
    tmp_file = os.path.join(tmp_path, tmp_file)
    file_retriever = FileRetriever(key, tmp_file, file_path)
    file_retriever.retrieve_file(lambda: Path(tmp_file).touch())

    if file_path is not None:
        assert os.path.exists(file_path)
        assert file_retriever.actual_path == file_path
    else:
        assert os.path.exists(tmp_file)
        assert file_retriever.actual_path == tmp_file

# TFBSRetriever class
from sponge.modules.data_retriever.tfbs_retriever import TFBSRetriever

@pytest.mark.parametrize('settings, assembly, on_the_fly', [
    ({}, 'hg38', True),
    ({}, 'hg19', False),
    ({}, 'hg38', False),
    ({'tfbs_file': 'LICENSE'}, 'hg38', False),
])
def test_tfbs_retriever(settings, assembly, on_the_fly, default_user_config,
    tmp_path):
    # The full bigbed file is way too big, just a placeholder
    test_file = ('https://raw.githubusercontent.com/kuijjerlab/sponge/main/'
        'LICENSE')
    default_user_config.deep_update({'motif': settings})
    motif_settings = default_user_config['motif']
    if motif_settings['jaspar_release'] is None:
        motif_settings['jaspar_release'] = 'JASPAR2024'
    tfbs_retriever = TFBSRetriever(
        tmp_path,
        test_file,
        motif_settings,
        assembly,
        on_the_fly
    )
    tfbs_retriever.retrieve_file()

    if not on_the_fly:
        if motif_settings['tfbs_file'] is None:
            assert os.path.exists(os.path.join(tmp_path, 'tfbs.bb'))
        else:
            assert os.path.exists(motif_settings['tfbs_file'])

# RegionRetriever class
from sponge.modules.data_retriever.region_retriever import RegionRetriever

@pytest.mark.parametrize('settings, assembly', [
    ({'chromosomes': ['chr19']}, 'hg38'),
    ({'chromosomes': ['chr19', 'chrX']}, 'hg19'),
    ({'region_file': 'LICENSE'}, 'hg38'),
    ({'chromosomes': ['chr19']}, 'random_assembly'),
])
def test_region_retriever(settings, assembly, core_config, default_user_config,
    tmp_path):
    default_user_config.deep_update({'region': settings})
    region_settings = default_user_config['region']
    region_retriever = RegionRetriever(
        tmp_path,
        core_config['url']['region']['xml'],
        core_config['url']['region']['rest'],
        core_config['url']['chrom_mapping'],
        region_settings,
        assembly,
        core_config['default_mapping'],
        core_config['default_chromosomes'],
    )
    region_retriever.retrieve_file()

    if region_settings['region_file'] is None:
        assert os.path.exists(os.path.join(tmp_path, 'regions.tsv'))
    else:
        assert os.path.exists(region_settings['region_file'])

# DataRetriever class
from sponge.modules.data_retriever import DataRetriever

def test_data_retriever():
    pass

# ProteinIDMapper class


# JasparRetriever class


# HomologyRetriever class


# MotifSelector class


# MatchFilter class


# PPIRetriever class
@pytest.mark.network
@pytest.mark.parametrize('input, expected_length', [
    (([],), 0),
    ((['GATA2', 'FOXF2', 'JUN'],), 4),
    ((['GATA2', 'FOXF2', 'JUN'], 400, False), 4),
    ((['GATA2', 'FOXF2', 'JUN'], 1000, False), 3),
])
def test_ppi_retriever(input, expected_length, ppi_retriever):
    ppi_retriever.retrieve_ppi(*input)

    assert len(ppi_retriever.get_ppi_frame()) == expected_length

# MatchAggregator class
@pytest.mark.parametrize('input, expected_length', [
    ((), 6), # Defaults are True, False
    ((False, False), 7),
    ((False, True), 4),
    ((True, True), 4),
])
def test_match_aggregator(input, expected_length, match_aggregator):
    match_aggregator.aggregate_matches(*input)

    assert len(match_aggregator.get_edges()) == expected_length

# FileWriter class
from sponge.modules.file_writer import FileWriter

@pytest.mark.parametrize('input', [
    (['tf1', 'tf2'], 'edge'),
    (['tf1', 'tf2'], 'score'),
    (['tf1', 'gene'], 'edge'),
    (['tf1', 'tf2', 'gene'], 'edge'),
    (['tf1'], 'edge'),
    ('tf1', 'edge'),
])
def test_file_writer(input, prior_frame, tmp_path):
    file_writer = FileWriter()
    write_path = os.path.join(tmp_path, 'test_frame.tsv')
    file_writer.write_network_file(prior_frame, *input, write_path)

    assert os.path.exists(write_path)

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
        os.path.join('src', 'sponge', 'user_config.yaml'),
    )


@pytest.mark.integration
@pytest.mark.network
def test_small_workflow(tmp_path):
    motif_output,ppi_output = run_integration_test_common(
        tmp_path,
        os.path.join('tests', 'sponge', 'test_user_config.yaml'),
    )

    motif_df = pd.read_table(motif_output, header=None)
    motif_df_t = pd.read_table(os.path.join('tests', 'sponge',
        'test_motif_prior.tsv'), header=None)

    pd.testing.assert_frame_equal(motif_df, motif_df_t)

    ppi_df = pd.read_table(ppi_output, header=None)
    ppi_df_t = pd.read_table(os.path.join('tests', 'sponge',
        'test_ppi_prior.tsv'), header=None)

    pd.testing.assert_frame_equal(ppi_df, ppi_df_t)