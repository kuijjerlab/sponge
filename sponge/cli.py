### Imports ###
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from sponge import Sponge
from sponge.config import DEFAULT_CHROMOSOMES, FILE_DF
from typing import Dict, List
from warnings import warn

### Functions ###
def process_file_paths(
    fp_list: List[str],
) -> Dict[str, str]:
    """
    Converts a list of descriptor to path assignments to a dictionary,
    invalid assignments result in a warning and are skipped.
    Valid form: valid_descriptor=non-empty/path/to/file
    Only one assignment operator (=) can be present.
    Valid descriptors are: promoter, jaspar_bigbed, ensembl

    Parameters
    ----------
    fp_list : List[str]
        List of descriptor to path assignments to parse

    Returns
    -------
    Dict[str, str]
        Dictionary linking valid file descriptors to paths
    """

    processed = {}

    for fp in fp_list:
        split = fp.split('=')
        if len(split) == 1:
            warn('Path assignment does not contain an equal sign: '
                f'{fp}, skipping')
            continue
        if len(split) > 2:
            warn('Path assignment contains more than one equal sign: '
                f'{fp}, skipping')
            continue
        if split[0] not in FILE_DF.index:
            warn(f'File description not recognised: {split[0]}, skipping')
            continue
        if len(split[1]) == 0:
            warn(f'Path to {split[0]} is empty, skipping')
            continue
        processed[split[0]] = split[1]

    return processed


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
    Developed by Ladislav Hovan (ladislav.hovan@ncmm.uio.no).
    """
    EPILOG = """
    Code available under GPL-3.0 license:
    https://github.com/ladislav-hovan/sponge
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
        nargs='*', default=DEFAULT_CHROMOSOMES, metavar='CHR')
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

    sponge_obj = Sponge(
        temp_folder=args.temp_folder,
        run_default=True,
        jaspar_release=args.jaspar_release,
        genome_assembly=args.genome_assembly,
        n_processes=args.n_processes,
        paths_to_files=process_file_paths(args.paths_to_files),
        tf_names=args.tf_names,
        matrix_ids=args.matrix_ids,
        drop_heterodimers=(not args.keep_hd),
        chromosomes=args.chromosomes,
        tss_offset=tuple(args.tss_offset),
        score_threshold=args.score_threshold,
        on_the_fly_processing=args.on_the_fly,
        protein_coding_only=args.pco,
        use_gene_names=(not args.use_gene_ids),
        weighted=args.weighted,
        motif_outfile=args.motif_outfile,
        ppi_outfile=args.ppi_outfile,
        prompt=(not args.yes),
    )
    sponge_obj.show_fingerprint()