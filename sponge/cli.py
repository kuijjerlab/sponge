### Imports ###
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from sponge import Sponge

### Functions ###
def cli(
) -> None:
    """
    Command line interface to the SPONGE class. All options can be
    specified, except run_default which is fixed to True. Execute
    with the --help option for more details.
    """

    DESCRIPTION = """
    SPONGE - Simple Prior Omics Network GEnerator.
    Generates prior motif and PPI networks, usable by other NetZoo tools
    (most notably PANDA).
    Uses the Ensembl, JASPAR, NCBI and STRING databases.
    Developed by Ladislav Hovan (ladislav.hovan@ncmbm.uio.no).
    """
    EPILOG = """
    Code available under GPL-3.0 license:
    https://github.com/kuijjerlab/sponge
    """

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
        description=DESCRIPTION, epilog=EPILOG)

    parser.add_argument('-t', '--temp-folder', dest='temp_folder',
        help='folder to save temporary files to',
        default='.sponge_temp', metavar='DIR')
    parser.add_argument('-jr', '--jaspar-release', dest='jaspar_release',
        help='JASPAR release to be used, newest if not specified')
    parser.add_argument('-ga', '--genome-assembly', dest='genome_assembly',
        help='genome assembly to be used, from Ensembl if not specified')
    parser.add_argument('-np', '--n-processes', dest='n_processes',
        help='number of processes to run in parallel',
        type=int, default=1, metavar='INT')
    parser.add_argument('-pf', '--paths-to-files', dest='paths_to_files',
        help='paths to files in the format of file_decriptor=path/to/file, '
            'any required file not provided will be downloaded',
        nargs='*', default=[], metavar='DESC=PATH')
    parser.add_argument('-tfn', '--tf-names', dest='tf_names',
        help='TF names to select, no filtering if empty',
        nargs='*', default=[], metavar='NAME')
    parser.add_argument('-mid', '--matrix-ids', dest='matrix_ids',
        help='TF JASPAR matrix IDs to select, no filtering if empty',
        nargs='*', default=[], metavar='ID')
    parser.add_argument('-kh', '--keep-heterodimers', dest='keep_hd',
        help='whether to keep TF heterodimers',
        action='store_true')
    parser.add_argument('-c', '--chromosomes', dest='chromosomes',
        help='chromosomes to select',
        nargs='*', default=None, metavar='CHR')
    parser.add_argument('-to', '--tss-offset', dest='tss_offset',
        help='offset from TSS for promoters or regions of interest in which '
            'to search for TF binding sites',
        type=int, nargs=2, default=[-750, 250], metavar='INT')
    parser.add_argument('-st', '--score-threshold', dest='score_threshold',
        help='score threshold for filtering TF binding sites',
        type=float, default=400, metavar='FLOAT')
    parser.add_argument('-otf', '--on-the-fly', dest='on_the_fly',
        help='whether to perform on the fly download of individual TF tracks '
            'instead of using the whole bigbed file',
        action='store_true')
    parser.add_argument('-pco', '--protein-coding-only', dest='pco',
        help='whether to only include protein coding genes in the prior',
        action='store_true')
    parser.add_argument('-ugi', '--use-gene-ids', dest='use_gene_ids',
        help='whether to use gene IDs instead of gene names',
        action='store_true')
    parser.add_argument('-w', '--weighted', dest='weighted',
        help='whether the motif prior should use edge weights based on score '
            'rather than binary',
        action='store_true')
    parser.add_argument('-mo', '--motif-outfile', dest='motif_outfile',
        help='file where the motif prior will be saved',
        default='motif_prior.tsv', metavar='FILE')
    parser.add_argument('-po', '--ppi-outfile', dest='ppi_outfile',
        help='file where the PPI prior will be saved',
        default='ppi_prior.tsv', metavar='FILE')
    parser.add_argument('-y', '--yes', dest='yes',
        help='whether to skip input prompts',
        action='store_true')

    args = parser.parse_args()

    sponge_obj = Sponge()
    sponge_obj.show_fingerprint()