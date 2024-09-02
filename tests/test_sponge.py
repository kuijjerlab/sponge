import datetime
import os
import pytest

import pandas as pd

from sponge.sponge import Sponge
from sponge.test_fixtures import *

### Unit tests ###
# Helper functions
import sponge.helper_functions as helper_f

@pytest.mark.parametrize('input, expected_output', [
    (0, 0),
    (0.5, -0.5),
    (1, 0),
])
def test_plogp(input, expected_output):
    assert helper_f.plogp(input) == expected_output


def test_calculate_ic_no_info(no_info_motif):
    assert helper_f.calculate_ic(no_info_motif) == 0


def test_calculate_ic_all_the_same(all_A_motif):
    # Length of the test motif is 6, so expected value is 2 * 6 = 12
    assert helper_f.calculate_ic(all_A_motif) == 12


def test_calculate_ic_SOX2(SOX2_motif):
    assert (helper_f.calculate_ic(SOX2_motif) ==
        pytest.approx(12.95, abs=0.01))


@pytest.mark.parametrize('input, expected_output', [
    ('CAB', 'Cab'),
    ('SOX2', 'SOx2'),
    ('ARHGAP21', 'ARHGAP21'),
    ('ABC2DE', 'ABC2de'),
])
def test_adjust_gene_name(input, expected_output):
    assert helper_f.adjust_gene_name(input) == expected_output


@pytest.mark.parametrize('input, expected_output', [
    ('a_string', 'a_string'),
    (datetime.datetime(1992, 5, 29, 23, 15), '29/05/1992, 23:15:00'),
])
def test_parse_datetime(input, expected_output):
    assert helper_f.parse_datetime(input) == expected_output


# Filtering functions
import sponge.filtering as filter_f

@pytest.mark.parametrize('input, expected_length', [
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'),
        ['MA0036.4', 'MA0030.2', 'MA0147.4'], 'chr19', 0, 2_000_000), 55),
])
def test_filter_edges(input, expected_length, chr19_promoters):
    df = filter_f.filter_edges(input[0], chr19_promoters, *input[1:])

    assert type(df) == pd.DataFrame
    assert len(df) == expected_length


@pytest.mark.parametrize('input, expected_length', [
    ((os.path.join('tests', 'sponge', 'chr19_subset.bb'),
        ['MA0036.4', 'MA0030.2', 'MA0147.4'], ['chr1', 'chr19']), 55),
])
def test_iterate_chromosomes(input, expected_length, chr19_promoters):
    df_list = filter_f.iterate_chromosomes(input[0], chr19_promoters,
        *input[1:])

    assert sum(len(df) for df in df_list) == expected_length


def test_process_chromosome(chr19_promoters, foxf2_chr19):
    df = filter_f.process_chromosome(foxf2_chr19, chr19_promoters)

    assert len(df) == 35


def test_process_motif(chr19_promoters, foxf2_chr19):
    df = filter_f.process_motif(foxf2_chr19, chr19_promoters)

    assert len(df) == 35


@pytest.mark.network
@pytest.mark.parametrize('input, expected_length', [
    ((['FOXF2'], ['MA0030.2'], ['chr1', 'chr19'], 'JASPAR2024', 'hg38'), 51),
])
def test_iterate_motifs(input, expected_length, chr19_promoters):
    df_list = filter_f.iterate_motifs(chr19_promoters, *input)

    assert sum(len(df) for df in df_list) == expected_length


# File retrieval functions
import sponge.data_retrieval as data_f

@pytest.mark.parametrize('cl_input, expected_output', [
    (['yes'], True),
    (['n'], False),
    (['what?', 'hell yeah'], True),
    (['why?', 'no'], False),
    (['1', '2', '3', '4', 'yes', 'irrelevant', 'no'], True),
])
def test_prompt_to_confirm(cl_input, expected_output, monkeypatch):

    def mock_input():
        for i in cl_input:
            yield i
    mi = mock_input()

    monkeypatch.setattr('builtins.input', lambda _: mi.__next__())
    reply = data_f.prompt_to_confirm('Yes or no?')

    assert reply == expected_output


@pytest.mark.parametrize('input, expected_output', [
    (('jaspar_bigbed', '.sponge_temp'),
        os.path.join('.sponge_temp', 'JASPAR.bb')),
    (('nothing', '.sponge_temp'), None),
])
def test_description_to_path(input, expected_output):
    assert data_f.description_to_path(*input) == expected_output


@pytest.mark.parametrize('input, expected_output', [
    (('jaspar_bigbed', os.path.join('tests', 'sponge')), True),
    (('jaspar_bigbed', 'nonexistent_dir'), False),
    (('nothing', '.sponge_temp'), False),
])
def test_check_file_exists(input, expected_output):
    assert data_f.check_file_exists(*input) == expected_output


@pytest.mark.network
def test_load_promoters_from_biomart(tmp_path):
    file_path = os.path.join(tmp_path, 'promoters.bed')
    data = data_f.load_promoters_from_biomart(file_path)

    assert os.path.exists(file_path)
    assert type(data['version']) == str
    assert type(data['ensembl']) == pd.DataFrame


@pytest.mark.network
def test_load_ensembl_from_biomart(tmp_path):
    file_path = os.path.join(tmp_path, 'ensembl.tsv')
    data = data_f.load_ensembl_from_biomart(file_path)

    assert os.path.exists(file_path)
    assert type(data['version']) == str
    assert type(data['ensembl']) == pd.DataFrame


@pytest.mark.network
@pytest.mark.parametrize('input, compare_to', [
    (('https://raw.githubusercontent.com/ladislav-hovan/sponge/main/LICENSE',
        'LICENSE'), 'LICENSE'),
    (('https://raw.githubusercontent.com/ladislav-hovan/sponge/main/LICENSE',
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


@pytest.mark.network
@pytest.mark.parametrize('input, expected_length', [
    (('UniProtKB_AC-ID', 'GeneID', ['O95905', 'P12312', 'P04439']), 2),
    (('Gene_Name', 'UniProtKB', 'BRCA1'), 1),
])
def test_get_uniprot_mapping(input, expected_length):
    mapping = data_f.get_uniprot_mapping(*input)

    assert len(mapping) == expected_length


@pytest.mark.network
def test_get_ensembl_assembly():
    assert data_f.get_ensembl_assembly() == 'hg38'


@pytest.mark.network
@pytest.mark.parametrize('input', [
    'hg38',
    'random_assembly',
    'hg1992',
])
def test_get_chromosome_mapping(input):
    mapping = data_f.get_chromosome_mapping(input)

    for chr in [str(i) for i in range(1, 23)] + ['X', 'Y']:
        assert chr in mapping[0].index
        assert chr in mapping[1].values


# Analysis functions
import sponge.analysis as anal_f

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


### Integration tests ###
# The test is marked as slow because the download of the bigbed file takes
# a lot of time and the filtering is also time consuming unless parallelised
@pytest.mark.integration
@pytest.mark.network
@pytest.mark.slow
def test_full_default_workflow(tmp_path):
    ppi_output = os.path.join(tmp_path, 'ppi_prior.tsv')
    motif_output = os.path.join(tmp_path, 'motif_prior.tsv')

    _ = Sponge(
        run_default=True,
        prompt=False,
        temp_folder=tmp_path,
        ppi_outfile=ppi_output,
        motif_outfile=motif_output,
    )

    assert os.path.exists(ppi_output)
    assert os.path.exists(motif_output)


@pytest.mark.integration
@pytest.mark.network
def test_small_workflow(tmp_path):
    ppi_output = os.path.join(tmp_path, 'ppi_prior.tsv')
    motif_output = os.path.join(tmp_path, 'motif_prior.tsv')

    _ = Sponge(
        jaspar_release='JASPAR2024',
        run_default=True,
        prompt=False,
        on_the_fly_processing=True,
        tf_names=['GATA2', 'FOXF2', 'MYC'],
        chromosomes=['chr21'],
        temp_folder=tmp_path,
        ppi_outfile=ppi_output,
        motif_outfile=motif_output,
    )

    assert os.path.exists(ppi_output)
    assert os.path.exists(motif_output)

    ppi_df = pd.read_csv(ppi_output, sep='\t', header=None)
    ppi_df_t = pd.read_csv(os.path.join('tests', 'sponge',
        'test_ppi_prior.tsv'), sep='\t', header=None)

    pd.testing.assert_frame_equal(ppi_df, ppi_df_t)

    motif_df = pd.read_csv(motif_output, sep='\t', header=None)
    motif_df_t = pd.read_csv(os.path.join('tests', 'sponge',
        'test_motif_prior.tsv'), sep='\t', header=None)

    pd.testing.assert_frame_equal(motif_df, motif_df_t)