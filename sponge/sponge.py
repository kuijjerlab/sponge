### Imports ###
import pickle

import numpy as np

from datetime import datetime

from pyjaspar import jaspardb

from typing import Mapping

from multiprocessing import Pool

from math import ceil, sqrt

from sponge.file_retrieval import *
from sponge.helper_functions import *
from sponge.filtering import *

from shutil import rmtree

FINGERPRINTS = defaultdict(Dict[str, Union[str, datetime, bool]])

### Class definition ###
class Sponge:
    """
    Sponge class can process the data necessary for creating a prior
    TF-gene regulatory network, along with a prior protein-protein 
    interaction network. It also contains tools to download these data,
    so that minimal input from the user is required. The networks are 
    provided in a format compatible with PANDA/LIONESS and other
    NetZoo tools.

    The run_default option is implemented in the constructor which will
    run the entire pipeline after class instance initialisation using 
    the provided options and defaults where appropriate.

    Usage:
    sponge_obj = Sponge(run_default=True)
    OR
    sponge_obj = Sponge()
    sponge_obj.select_tfs()
    sponge_obj.find_human_homologs()
    sponge_obj.filter_matches()
    sponge_obj.retrieve_ppi()
    sponge_obj.write_ppi_prior()
    sponge_obj.aggregate_matches()
    sponge_obj.write_motif_prior()

    Most functions have options which can usually also be provided 
    to the constructor, for more details refer to the documentation of 
    the individual functions.
    """

    def __init__(
        self,
        temp_folder: PATH = '.sponge_temp',
        run_default: bool = False,
        jaspar_release: Optional[str] = None,
        genome_assembly: Optional[str] = None,
        n_processes: int = 1,
        paths_to_files: Mapping[str, FILE_LIKE] = {},
        drop_heterodimers: bool = True,
        chromosomes: Iterable[str] = [f'chr{i}' for i in [j for j in 
            range(1, 23)] + ['M', 'X', 'Y']],
        tss_offset: Tuple[int, int] = (-750, 250),
        score_threshold: float = 400,
        protein_coding_only: bool = False,
        use_gene_names: bool = True,
        weighted: bool = False,
        motif_outfile: FILE_LIKE = 'motif_prior.tsv',
        ppi_outfile: FILE_LIKE = 'ppi_prior.tsv',
        prompt: bool = True
    ):
        """
        Initialises an instance of the Sponge class.

        Parameters
        ----------
        temp_folder : PATH, optional
            The temporary folder for saving downloaded files, 
            by default '.sponge_temp'
        run_default : bool, optional
            Whether to run the default pipeline automatically after
            class instance initialisation, by default False
        jaspar_release : Optional[str], optional
            Which JASPAR release to use or None to select the newest, 
            by default None
        genome_assembly : Optional[str], optional
            The genome assembly used in the provided files or None to
            deduce automatically during file download, by default None
            If this option is not provided it is assumed that the files
            match the default assembly of Ensembl
        n_processes : int, optional
            The number of processes to run in parallel for the filtering
            of the bigbed file, by default 1
        paths_to_files : Mapping[str, FILE_LIKE], optional
            A dictionary of paths to required files, keyed by their 
            descriptions; if a path is not provided the files will be 
            looked for in the temp folder or downloaded, by default {}
            Available descriptions: 'homologene', 'promoter', 
                'jaspar_bigbed', 'ensembl'
        drop_heterodimers : bool, optional
            Whether to drop the heterodimer motifs from consideration,
            by default True
        chromosomes : Iterable[str], optional
            Which chromosomes to get the promoters from, by default
            [f'chr{i}' for i in [j for j in range(1, 23)] + 
            ['M', 'X', 'Y']] 
        tss_offset : Tuple[int, int], optional
            The offset from the transcription start site to use for the
            assignment of transcription factors to promoters, 
            by default (-750, 250)
        score_threshold : float, optional
            The minimal score of a match for it to be included in the 
            prior, by default 400
        protein_coding_only : bool, optional
            Whether to restrict the gene selection to only protein 
            coding genes, by default False
        use_gene_names : bool, optional
            Whether to use gene names instead of Ensemble IDs in the 
            output, by default True
        weighted : bool, optional
            Whether to use weighted rather than binary prior networks,
            by default False
        motif_outfile : FILE_LIKE, optional
            The path to save the motif prior into, by default 
            'motif_prior.tsv'
        ppi_outfile : FILE_LIKE, optional
            The path to save the PPI prior into, by default
            'ppi_prior.tsv'
        prompt : bool, optional
            Whether to prompt to confirm the file downloads, 
            by default True
        """

        # Initialise important variables
        self.temp_folder = temp_folder
        self.ensembl = None
        self.ppi_frame = None
        self.motif_frame = None
        self.n_processes = n_processes
        self.fingerprint = defaultdict(dict)
        self.provided_paths = paths_to_files
        self.jaspar_release = None

        # Store the provided options
        self.drop_heterodimers = drop_heterodimers
        self.chromosomes = chromosomes
        self.tss_offset = tss_offset
        self.score_threshold = score_threshold
        self.protein_coding_only = protein_coding_only
        self.use_gene_names = use_gene_names
        self.weighted = weighted
        self.motif_outfile = motif_outfile
        self.ppi_outfile = ppi_outfile
        self.assembly = genome_assembly

        # Retrieve genome assembly if necessary
        if self.assembly is None:
            self.assembly = get_ensembl_assembly()

        # Create the temporary folder
        if not os.path.exists(self.temp_folder):
            os.mkdir(self.temp_folder)

        # Initialise the JASPAR database object
        self.initialise_jaspar(jaspar_release)
        if self.jaspar_release is None:
            return
        self.log_fingerprint('JASPAR', self.jaspar_release)

        # Locate or download all the required files
        self.files_ready = self.prepare_files(prompt)

        # Run default workflow if selected
        if run_default and self.files_ready:
            self.run_default_workflow()


    def initialise_jaspar(
        self,
        jaspar_release: Optional[str] = None
    ) -> None:
        """
        Initialises the JASPAR database object to be used for motif
        selection.

        Parameters
        ----------
        jaspar_release : Optional[str], optional
            Which JASPAR release to use or None to select the newest, 
            by default None
        """

        # Initialise a database object to interact with
        self.jdb_obj = jaspardb()
        if jaspar_release is None:
            # Just keep the current object, log the release version
            from pyjaspar import JASPAR_LATEST_RELEASE
            self.jaspar_release = JASPAR_LATEST_RELEASE
        else:
            jaspar_available = self.jdb_obj.get_releases()
            if jaspar_release in jaspar_available:
                # Release found as specified
                self.jaspar_release = jaspar_release
                self.jdb_obj = jaspardb(self.jaspar_release)
            elif 'JASPAR' + jaspar_release in jaspar_available:
                # Try adding JASPAR to the provided release
                # Converts e.g. 2022 to JASPAR2022 (actual release name)
                print (f'Found {"JASPAR" + jaspar_release} in available '
                    'releases, assuming this matches your choice')
                self.jaspar_release = 'JASPAR' + jaspar_release
                self.jdb_obj = jaspardb(self.jaspar_release)
            else:
                # Default self.jaspar_release is None, so this will be
                # registered as failure downstream
                print ('The specified version of the JASPAR release '
                    f'({jaspar_release}) is not available')
                print ('Available versions:')
                print (', '.join(jaspar_available))
                print ('The downstream pipeline will not be run')


    ### File retrieval functions ###
    def prepare_files(
        self,
        prompt: bool = True
    ) -> bool:
        """
        Locates or downloads the files required for the running of the
        prior creation pipeline.

        Parameters
        ----------
        prompt : bool, optional
            Whether to prompt to confirm the file downloads, 
            by default True

        Returns
        -------
        bool
            Whether the files were located or downloaded successfully
        """

        print ()
        print ('--- Running prepare_files() ---')

        # Retrieve and save the chromosome name mappings
        chrom_mappings = get_chromosome_mapping(self.assembly)
        self.ens_to_ucsc = chrom_mappings[0]
        self.ucsc_to_ens = chrom_mappings[1]

        # Attemp to first use the provided paths
        provided = []
        for k,v in self.provided_paths.items():
            if k not in FILE_DF.index:
                print (f'Unrecognised file type: {k}, ignoring provided path')
                continue
            if not os.path.exists(v):
                print (f'The path provided for file type {k} ({v}) is invalid '
                    'and will be ignored')
            else:
                provided.append(k)
                self.log_fingerprint(k.upper(), '', provided=True)
        self.provided_paths = {k: v for k,v in self.provided_paths.items()
            if k in provided}

        # Check for the files in the temp folder
        to_check = [file for file in FILE_DF.index if file not in provided]
        to_retrieve = {}
        for file in to_check:
            if not check_file_exists(file, self.temp_folder):
                to_retrieve[file] = FILE_DF.loc[file, 'name']
            else:
                self.log_fingerprint(file.upper(), '', cached=True)

        if len(to_retrieve) == 0:
            print ('All the required files were provided or located in the '
                f'{self.temp_folder} directory')
            return True
        else:
            # Download anything that's missing
            print ('The following files were not found:')
            print (', '.join(to_retrieve.values()))
            if prompt:
                reply = prompt_to_confirm('Do you want to download them '
                    'automatically?')
            else:
                print ('They will be downloaded automatically')
                print ()
            if not prompt or reply:
                for description in to_retrieve.keys():
                    file_path = self.retrieve_file(description, prompt=False)
                    # This shouldn't really happen as we're only using valid
                    # descriptions, but better safe than sorry
                    if file_path is None:
                        print (f'File {to_retrieve[description]} could not be '
                            'retrieved')
                        return False
            else:
                # Prompt was refused
                return False


    def retrieve_file(
        self,
        description: str,
        prompt: bool = True
    ) -> Optional[str]:
        """
        Attempts to retrieve a file corresponding to a given 
        description. If the file is found in the temporary directory or 
        has been flagged as provided it is not downloaded, but the path
        to it will be returned in any case. If the file cannot be
        retrieved (e.g. because the description is invalid), None is
        returned.

        Parameters
        ----------
        description : str
            The description of a file to be retrieved
            Available descriptions: 'homologene', 'promoter', 
                'jaspar_bigbed', 'ensembl'
        prompt : bool, optional
            Whether to prompt to confirm the file downloads, 
            by default True

        Returns
        -------
        Optional[str]
            The path to the retrieved file or None if the file was not
            retrieved

        Raises
        ------
        ValueError
            If self.jaspar_release is not specified prior to attempting
            to retrieve the jaspar_bigbed file
        """
        
        # Just use the provided path
        if description in self.provided_paths:
            print ('Using provided file', self.provided_paths[description])
            print ()
            return self.provided_paths[description]
        # Figure out where in the temp folder the file should be
        file_path = description_to_path(description, self.temp_folder)
        # Redundant now as the class initialisation creates the temp folder
        if not os.path.exists(self.temp_folder):
            os.mkdir(self.temp_folder)
        # Use the cached file if it exists
        if os.path.exists(file_path):
            print ('Using cached file', file_path)
            print ()
        else:
            print (f'File {FILE_DF.loc[description, "name"]} not found ' 
                f'in directory {self.temp_folder}')
            if prompt:
                reply = prompt_to_confirm('Do you want to download it?')
                if not reply:
                    return None
            # The URL to retrieve the file from
            to_request = FILE_DF.loc[description, 'url']
            if to_request is None:
                # No URL in this case, retrieve via function call
                # These options are passed to the function as kwargs
                options = {'file_path': file_path}
                if description == 'promoter':
                    # Extra options for the promoter file
                    options['tss_offset'] = self.tss_offset
                    if self.chromosomes is None:
                        options['chromosomes'] = None
                    else:
                        options['chromosomes'] = [self.ucsc_to_ens[x] for x in
                            self.chromosomes]
                    options['chromosome_mapping'] = self.ens_to_ucsc
                # Call the responsible function
                result = eval(FILE_DF.loc[description, 'eval'])
                # Log the fingerprint for the retrieved file
                self.log_fingerprint(description.upper(), result['version'])
                # Save the Ensembl DataFrame
                if 'ensembl' in result:
                    self.ensembl = result['ensembl']
            else:
                if description == 'jaspar_bigbed':
                    # The file to download depends on the JASPAR release
                    if self.jaspar_release is None:
                        raise ValueError('The release of jaspar has to be '
                            'specified in order to retrieve the bigbed file')
                    to_request = to_request.format(
                        year=self.jaspar_release[-4:],
                        genome_assembly=self.assembly)
                    version = self.jaspar_release
                elif description == 'homologene':
                    # The version needs to be determined via a separate request
                    version_url = '/'.join(to_request.split('/')[:-1] + 
                        ['RELEASE_NUMBER'])
                    version_request = requests.get(version_url, stream=True)
                    version = version_request.text.strip()
                # Download the specified file
                print (f'Downloading data into {file_path}...')
                download_with_progress(to_request, file_path)
                print ()
                # Log the fingerprint for the retrieved file
                self.log_fingerprint(description.upper(), version)

        return file_path


    def update_label_in_cache(
        self,
        temp_fingerprint: FINGERPRINTS,
        label: str
    ) -> None:
        """
        Updates a given label in the cached fingerprint file.

        Parameters
        ----------
        temp_fingerprint : FINGERPRINTS
            The fingerprint to be updated, can be empty or initially
            loaded from the cached fingerprint
        label : str
            The label to update
        """

        # Update the label
        temp_fingerprint[label] = self.fingerprint[label]
        # Dump into the fingerprint file
        fingerprint_file = os.path.join(self.temp_folder, '.fingerprint')
        pickle.dump(temp_fingerprint, open(fingerprint_file, 'wb'))


    def log_fingerprint(
        self,
        label: str,
        version: str,
        provided: bool = False,
        cached: bool = False
    ) -> None:
        """
        Logs the fingerprint (label, version, retrieval time) for a
        given file. The fingerprints are stored within the class
        instance and also as a hidden file in the temp folder.

        Parameters
        ----------
        label : str
            The label corresponding to the file
        version : str
            The version of the database used to make the file
        provided : bool, optional
            Whether the path to the file was provided by user, 
            by default False
        cached : bool, optional
            Whether the file was retrieved from cache, by default False
        """
        
        # Update the details in the fingerprint based on the provided data
        self.fingerprint[label]['version'] = version
        self.fingerprint[label]['datetime'] = datetime.fromtimestamp(
            time.time())
        self.fingerprint[label]['provided'] = provided
        self.fingerprint[label]['cached'] = cached

        # Update the fingerprint file as well
        fingerprint_file = os.path.join(self.temp_folder, '.fingerprint')
        # Load the current fingerprint or make a new one
        if os.path.exists(fingerprint_file):
            temp_fingerprint = pickle.load(open(fingerprint_file, 'rb'))
        else:
            temp_fingerprint = {}
        if not provided and not cached:
            # No prior data, just rewrite
            self.update_label_in_cache(temp_fingerprint, label)
        if cached:
            if label in temp_fingerprint:
                # Previously cached, so keep the cached data (more precise)
                self.fingerprint[label] = temp_fingerprint[label]
                self.fingerprint[label]['cached'] = True
            else:
                # File is cached but not in the cached fingerprint
                self.fingerprint[label]['version'] = 'unknown version'
                self.fingerprint[label]['datetime'] = 'at unknown time'
                self.update_label_in_cache(temp_fingerprint, label)


    ### Main workflow functions ###
    def run_default_workflow(  
        self
    ) -> None:
        """
        Runs the default workflow by running the following functions
        in sequence using the options provided during class instance
        initialisation:

        self.select_tfs()
        self.find_human_homologs(prompt=False)
        self.filter_matches(prompt=False)
        self.retrieve_ppi()
        self.write_ppi_prior()
        self.aggregate_matches(prompt=False)
        self.write_motif_prior()

        This will result in the generation of priors as specified and
        should be used unless there is a desire to change the provided
        settings after initialisation.
        """

        self.select_tfs()
        self.find_human_homologs(prompt=False)
        self.filter_matches(prompt=False)
        self.retrieve_ppi()
        self.write_ppi_prior()
        self.aggregate_matches(prompt=False)
        self.write_motif_prior()


    def select_tfs(
        self,
        drop_heterodimers: Optional[bool] = None
    ) -> None:
        """
        Selects transcription factors from the newest version of the
        JASPAR database and stores them in the class instance.

        Parameters
        ----------
        drop_heterodimers : Optional[bool], optional
            Whether to drop heterodimer motifs or None to follow the
            option from the initialisation, by default None
        """

        print ()
        print ('--- Running select_tfs() ---')
        print (f'Using: {self.jaspar_release}')

        if drop_heterodimers is None:
            drop_heterodimers = self.drop_heterodimers

        # All vertebrate motifs
        motifs = self.jdb_obj.fetch_motifs(collection='CORE', 
            tax_group='vertebrates', all_versions=True)
        print ()
        print ('All motif versions:', len(motifs))
        print ('Motif base IDs:', len(set([i.base_id for i in motifs])))

        # Select latest, preferring human ones
        latest = {}
        for i in motifs:
            if i.base_id not in latest:
                latest[i.base_id] = [i.matrix_id, i.species]
            else:
                # Replace with newer version if the new one is human or the old 
                # one isn't       
                if (('9606' in i.species) or 
                    ('9606' not in latest[i.base_id][1])):
                    # This could be added to the logical condition above but 
                    # this is more readable
                    if int(i.matrix_id[-1]) > int(latest[i.base_id][0][-1]):
                        latest[i.base_id] = [i.matrix_id, i.species]
        motifs_latest = [i for i in motifs if 
            i.matrix_id == latest[i.base_id][0]]
        
        # Keep only one motif per TF
        # Consider dropping this requirement maybe
        tf_to_motif = {}
        for i in motifs_latest:
            if i.name not in tf_to_motif:
                tf_to_motif[i.name] = {i.matrix_id: calculate_ic(i)}
            else:
                tf_to_motif[i.name][i.matrix_id] = calculate_ic(i)
        self.tf_to_motif = tf_to_motif
        motifs_unique = [i for i in motifs_latest if 
            (tf_to_motif[i.name][i.matrix_id] == 
            max(tf_to_motif[i.name].values()))]
        print ('Unique motifs:', len(motifs_unique))

        # Drop heterodimers
        if drop_heterodimers:
            motifs_nohd = [i for i in motifs_unique if '::' not in i.name]
            print ('Motifs without heterodimers:', len(motifs_nohd))
            self.motifs = motifs_nohd
        else:
            self.motifs = motifs_unique


    def find_human_homologs(
        self, 
        homologene_file: Optional[FILE_LIKE] = None,
        prompt: bool = True
    ) -> None:
        """
        Attempts to map all initially selected non-human transcription
        factors to their human homologs.

        Parameters
        ----------
        homologene_file : Optional[FILE_LIKE], 
            optional
            The path to a homologene file or None to use cache or
            download it, by default None
        prompt : bool, optional
            Whether to prompt before downloading, by default True
        """
        
        print ()
        print ('--- Running find_human_homologs() ---')

        if homologene_file is None:
            homologene_file = self.retrieve_file('homologene', prompt=prompt)
            if homologene_file is None:
                print ('Unable to find or retrieve the homologene file, ' 
                    'exiting')
                return

        # Get the non-human motifs
        non_human_motifs = [i for i in self.motifs if '9606' not in i.species]
        print ('Non-human motifs:', len(non_human_motifs))

        # Read the homologene database
        hg_df = pd.read_csv(homologene_file, sep='\t', header=None, 
                    names=['HG Group ID', 'TaxID', 'Gene ID', 
                           'Gene Symbol', 'Protein GI', 'Protein Accession'])

        # Get the non-human motif names
        non_human_motif_names = [i.name for i in non_human_motifs]
        # Compare against homologene
        found_names = hg_df[hg_df['Gene Symbol'].isin([adjust_gene_name(i) for 
            i in non_human_motif_names])]['Gene Symbol'].unique()
        # Find the missing ones
        missing = (set([adjust_gene_name(i) for i in non_human_motif_names]) - 
            set(found_names))
        print ()
        print ('Names missing from the homologene database:')
        for i in [(i.name, i.acc) for i in non_human_motifs if 
            i.name in missing]:
            print (i[0], *i[1])

        # Get the missing IDs from Uniprot API
        print ()
        print ('Retrieving matches from UniProt...')
        mapping = get_uniprot_mapping('UniProtKB_AC-ID', 'RefSeq_Protein',
            [i.acc[0] for i in non_human_motifs if i.name in missing])
        mapping.columns = ['Uniprot', 'Accession']

        # Create a DataFrame for matching missing entries
        hg_df[hg_df['Protein Accession'].isin(mapping['Accession'])]
        missing_df = pd.DataFrame([(i.name, i.acc[0]) for i in non_human_motifs 
            if i.name in missing], columns = ['Gene Symbol', 'Uniprot'])        
        matching_df = missing_df.join(mapping.set_index('Uniprot'), 
            on='Uniprot').join(hg_df.set_index('Protein Accession'), 
            on='Accession', rsuffix='_HG')

        def corresponding_id(
            name: str
        ) -> np.array:
            """
            Retrieves the corresponding group ID for a given gene name 
            from the homologene DataFrame, or matching DataFrame if it
            is not found.

            Parameters
            ----------
            name : str
                The gene name of interest

            Returns
            -------
            np.array
                The array containing the corresponding group ID
            """

            # Choose the values of the group ID from the homologene DataFrame
            # which correspond to the gene of interest
            values = hg_df[hg_df['Gene Symbol'] == name]['HG Group ID'].values
            if len(values) == 0:
                # If none are found, use the matching DataFrame instead
                matches = matching_df[matching_df['Gene Symbol'] == name]
                values = matches['HG Group ID'].values

            return values

        # Create a DataFrame of corresponding names
        corr_df = pd.DataFrame(non_human_motif_names, 
            columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(
            adjust_gene_name)
        corr_df['Group ID'] = corr_df['Adjusted Name'].apply(corresponding_id)
        corr_df['Group ID'] = corr_df['Group ID'].apply(lambda x: 
            x[0] if len(x) > 0 else np.nan)
        corr_df['Human Name'] = corr_df['Group ID'].apply(lambda x: 
            hg_df[(hg_df['HG Group ID'] == x) & 
            (hg_df['TaxID'] == 9606)]['Gene Symbol'].values)
        corr_df['Human Name'] = corr_df['Human Name'].apply(lambda x: 
            x[0] if len(x) > 0 else '')
        corr_df['Trivial'] = corr_df['Original Name'].apply(lambda x: 
            x.upper()) == corr_df['Human Name']

        # Find duplicates
        duplicated = corr_df[corr_df['Human Name'].duplicated(keep=False) & 
            (corr_df['Human Name'] != '')].copy()
        to_print = duplicated.groupby('Human Name')['Original Name'].unique(
            ).apply(lambda x: ' '.join(x))
        print ()
        print ('Duplicate names:')
        for i in to_print.index:
            print (f'{i}:', to_print.loc[i])

        # Calculate the information content for duplicates
        duplicated['IC'] = duplicated['Original Name'].apply(lambda x: 
            max(self.tf_to_motif[x].values()))
        # Keep the highest IC amongst the duplicates
        to_drop = duplicated['Original Name'][duplicated.sort_values(
            'IC').duplicated('Human Name', keep='last')]

        # Exlude the IDs which are already present among the human ones
        human_motif_names = [i.name for i in self.motifs if 
            '9606' in i.species]
        corr_df['Duplicate'] = corr_df['Human Name'].isin(human_motif_names)

        # Perform the final filtering - discard all duplicates and TFs without 
        # homologs
        corr_df_final = corr_df[(corr_df['Duplicate'] == False) & 
            (corr_df['Human Name'] != '') & 
            (corr_df['Original Name'].isin(to_drop) == False)]

        # The mapping of original to human names and the matrix IDs to be kept
        animal_to_human = {animal_name: human_name for animal_name, human_name 
            in zip(corr_df_final['Original Name'], 
            corr_df_final['Human Name'])}
        print ()
        print ('Final number of IDs which will be replaced by human homologs:', 
               len(animal_to_human))
        matrix_ids = [motif.matrix_id for motif in self.motifs if 
            (motif.name in human_motif_names or motif.name in animal_to_human)]
        print ('Final number of total matrix IDs:', len(matrix_ids))

        self.animal_to_human = animal_to_human
        self.matrix_ids = matrix_ids
        self.tf_names = human_motif_names + list(animal_to_human.keys())


    def filter_matches(
        self, 
        promoter_file: Optional[FILE_LIKE] = None, 
        bigbed_file: Optional[FILE_LIKE] = None,
        score_threshold: Optional[float] = None,
        chromosomes: Optional[Iterable[str]] = None,
        n_processes: Optional[int] = None,
        prompt: bool = True
    ) -> None:
        """
        Filters all the binding sites in the JASPAR bigbed file to
        select only the ones in the promoter regions of genes on given
        chromosomes, subject to a score threshold. Stores the result
        internally.

        Parameters
        ----------
        promoter_file : Optional[FILE_LIKE], optional
            The path to a promoter file or None to use cache or
            download it, by default None
        bigbed_file : Optional[FILE_LIKE], optional
            The path to a JASPAR bigbed file or None to use cache or
            download it, by default None
        score_threshold : Optional[float], optional
            The minimal score of a match for it to be included in the 
            prior or None to follow the option from the initialisation, 
            by default None
        chromosomes : Optional[Iterable[str]], optional
            Which chromosomes to get the promoters from or None to 
            follow the option from the initialisation, by default None
        n_processes : Optional[int], optional
            The number of processes to run in parallel or None to 
            follow the option from the initialisation, by default None
        prompt : bool, optional
            Whether to prompt before downloading, by default True
        """
        
        print ()
        print ('--- Running filter_matches() ---')

        if n_processes is None:
            n_processes = self.n_processes
        if chromosomes is None:
            chromosomes = self.chromosomes
            if chromosomes is None:
                chromosomes = self.ucsc_to_ens.index
        if score_threshold is None:
            score_threshold = self.score_threshold

        if promoter_file is None:
            promoter_file = self.retrieve_file('promoter', prompt=prompt)
            if promoter_file is None:
                print ('Unable to find or retrieve the promoter file, exiting')
                return
        
        if bigbed_file is None:
            bigbed_file = self.retrieve_file('jaspar_bigbed', prompt=prompt)
            if bigbed_file is None:
                print ('Unable to find or retrieve the JASPAR bigbed file, '
                    'exiting')
                return

        print ('Loading the promoter bed file...')
        df_full = bioframe.read_table(promoter_file, schema='bed')
        df_full['name'] = df_full['name'].apply(lambda x: x.split('.')[0])
        # Some parts of the bed file are unnecessary
        df_full.drop(columns=['score', 'strand'], inplace=True)
        df_full.set_index('name', inplace=True)

        results_list = []
        p = Pool(n_processes)

        print ()
        print ('Iterating over the chromosomes...')
        start_time = time.time()
        for chrom in chromosomes:
            st_chr = time.time()
            df_chrom = df_full[df_full['chrom'] == chrom]
            if len(df_chrom) == 0:
                suffix = 'no transcripts'
            elif len(df_chrom) == 1:
                suffix = '1 transcript'
            else:
                suffix = f'{len(df_chrom)} transcripts'
            print (f'Chromosome {chrom[3:]} with ' + suffix)
            if len(df_chrom) == 0:
                continue
            # This is a heuristic approximation of the ideal chunk size
            # Based off of performance benchmarking
            chunk_size = ceil(sqrt(len(df_chrom) / n_processes))
            chunk_divisions = [i for i in range(0, len(df_chrom), chunk_size)]
            input_tuples = [(bigbed_file, df_chrom, self.tf_names, chrom, i, 
                i+chunk_size, score_threshold) for i in chunk_divisions]
            # Run the calculations in parallel
            result = p.starmap_async(filter_edges, input_tuples, 
                chunksize=n_processes)
            edges_chrom_list = result.get()
            results_list += edges_chrom_list
            elapsed_chr = time.time() - st_chr
            print (f'Done in: {elapsed_chr // 60:n} m '
                f'{elapsed_chr % 60:.2f} s')

        elapsed = time.time() - start_time
        print ()
        print (f'Total time: {elapsed // 60:n} m {elapsed % 60:.2f} s')

        # Save the final results, ignoring the index makes this fast
        # The index is irrelevant
        self.all_edges = pd.concat(results_list, ignore_index=True)

    
    def load_matches(
        self,
        file_path: FILE_LIKE
    ):
        """
        Loads the filtered matches from a file, allows the use of
        the downstream SPONGE functions without running the steps up to 
        and including filter_matches 

        Parameters
        ----------
        file_path : FILE_LIKE
            The path to a file that contains the filtered matches 
            in a format compatible with what filter_matches generates
        """

        print ()
        print ('--- Running load_matches() ---')
        
        self.all_edges = pd.read_csv(file_path, sep='\t')


    def retrieve_ppi(
        self
    ) -> None:
        """
        Retrieves the protein-protein interaction data from the STRING
        database for the previously identified transcription factors.
        Stores the resulting network internally.
        """
        
        print ()
        print ('--- Running retrieve_ppi() ---')
        
        # Use the human names for the TFs
        filtered_tfs = self.all_edges['TFName'].unique()
        humanised_tfs = [self.animal_to_human[x] if x in self.animal_to_human
            else x for x in filtered_tfs]
        query_string = '%0d'.join(humanised_tfs)

        print ('Retrieving mapping from STRING...')
        mapping_request = requests.get(f'{STRING_URL}get_string_ids?'
            f'identifiers={query_string}&species=9606')
        mapping_df = pd.read_csv(BytesIO(mapping_request.content), sep='\t')
        mapping_df['queryName'] = mapping_df['queryIndex'].apply(
            lambda i: humanised_tfs[i])
        # Check where the preferred name doesn't match the query
        diff_df = mapping_df[mapping_df['queryName'] != 
            mapping_df['preferredName']]
        ids_to_check = np.concatenate((diff_df['queryName'], 
            diff_df['preferredName']))
        matching_ids = list(mapping_df[mapping_df['queryName'] == 
            mapping_df['preferredName']]['preferredName'])
        # Log the STRING version in the fingerprint
        version_request = requests.get(f'{STRING_URL}version')
        version_df = pd.read_csv(BytesIO(version_request.content), sep='\t')
        self.log_fingerprint('STRING', version_df['string_version'])
        
        if len(ids_to_check) > 0:
            # Retrieve UniProt identifiers for the genes with differing names
            print ('Checking the conflicts in the UniProt database...')
            uniprot_df = get_uniprot_mapping('Gene_Name', 'UniProtKB', 
                ids_to_check).set_index('from')
            p_to_q = {p: q for q,p in zip(diff_df['queryName'], 
                diff_df['preferredName'])}
            # Keep the conflicts where there is a match or where one or both
            # of the names doesn't find an identifier
            for p,q in p_to_q.items():
                if (p not in uniprot_df.index or q not in uniprot_df.index
                    or uniprot_df.loc[p, 'to'] == uniprot_df.loc[q, 'to']):
                    matching_ids.append(p)
        query_string_filt = '%0d'.join(matching_ids)

        print ('Retrieving the network from STRING...')
        request = requests.get(f'{STRING_URL}network?'
            f'identifiers={query_string_filt}&species=9606')
        ppi_df = pd.read_csv(BytesIO(request.content), sep='\t')

        print ('Processing the results...')
        ppi_df.drop(['stringId_A', 'stringId_B', 'ncbiTaxonId', 'nscore', 
            'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore'], 
            axis=1, inplace=True)
        ppi_df.rename(columns={'preferredName_A': 'tf1', 
            'preferredName_B': 'tf2'}, inplace=True)
        if len(ids_to_check) > 0:
            # Replace with names that have been queried (as used by JASPAR)
            ppi_df['tf1'].replace(p_to_q, inplace=True)
            ppi_df['tf2'].replace(p_to_q, inplace=True)
        ppi_df.sort_values(by=['tf1', 'tf2'], inplace=True)

        print ()
        print ('Final number of TFs in the PPI network: '
            f'{len(set(ppi_df["tf1"]).union(set(ppi_df["tf2"])))}')
        print (f'Final number of edges: {len(ppi_df)}')

        self.ppi_frame = ppi_df


    def write_ppi_prior(
        self,
        output_path: Optional[FILE_LIKE] = None,
        weighted: Optional[bool] = None
    ) -> None:
        """
        Writes the protein-protein interaction prior network into a 
        file.

        Parameters
        ----------
        output_path : Optional[FILE_LIKE], optional
            The path to write the prior into or None to follow the 
            option from the initialisation, by default None
        weighted : Optional[bool], optional
            Whether to use weights for the edges as opposed to making 
            them binary or None to follow the option from the 
            initialisation, by default None
        """
        
        print ()
        print ('--- Running write_ppi_prior() ---')
        
        if output_path is None:
            output_path = self.ppi_outfile
        if weighted is None:
            weighted = self.weighted

        if self.ppi_frame is None:
            print ('No motif prior has been generated yet, please run '
                'retrieve_ppi() first')
            return
        
        if weighted:
            self.ppi_frame[['tf1', 'tf2', 'score']].to_csv(output_path, 
                sep='\t', index=False, header=False)
        else:
            self.ppi_frame['edge'] = 1
            self.ppi_frame[['tf1', 'tf2', 'edge']].to_csv(output_path, 
                sep='\t', index=False, header=False)
            

    def aggregate_matches(
        self,
        ensembl_file: Optional[FILE_LIKE] = None,
        prompt: bool = True,
        use_gene_names: Optional[bool] = None,
        protein_coding_only: Optional[bool] = None
    ) -> None:
        """
        Aggregates all the matches corresponding to individual 
        transcripts into genes, creating a transcription factor - gene
        matrix. Stores the result internally.

        Parameters
        ----------
        ensembl_file : Optional[FILE_LIKE], optional
            The path to an Ensembl file or None to use cache or
            download it, by default None
        prompt : bool, optional
            Whether to prompt before downloading, by default True
        use_gene_names : Optional[bool], optional
            Whether to use gene names instead of Ensembl IDs or None
            to follow the option from the initialisation, 
            by default None
        protein_coding_only : Optional[bool], optional
            Whether to restrict the selection to only protein coding 
            genes or None to follow the option from the initialisation, 
            by default None
        """
        
        print ()
        print ('--- Running aggregate_matches() ---')

        if self.ensembl is None or ensembl_file is not None:
            if ensembl_file is None:
                ensembl_file = self.retrieve_file('ensembl', prompt=prompt)
                if ensembl_file is None:
                    print ('Unable to find or retrieve the ensembl file, ' 
                        'exiting')
                    return
            self.ensembl = pd.read_csv(ensembl_file, sep='\t')
        
        if use_gene_names is None:
            use_gene_names = self.use_gene_names
        if protein_coding_only is None:
            protein_coding_only = self.protein_coding_only

        # Add the Ensembl data (gene names) to the edges previously found
        motif_df = self.all_edges.join(other=self.ensembl.set_index(
            'Transcript stable ID'), on='transcript')
        print ('Number of TF - transcript edges:', len(motif_df))
        if protein_coding_only:
            motif_df = motif_df[motif_df['Gene type'] == 
                'protein_coding'].copy()
        # Drop columns that are not required anymore
        motif_df.drop(columns=['Gene type', 'name'], inplace=True)
        # Humanise the TF names
        motif_df['TFName'] = motif_df['TFName'].apply(lambda x: 
            self.animal_to_human[x] if x in self.animal_to_human else x)
        # Ignore genes without identifiers
        motif_df.dropna(subset=['Gene stable ID'], inplace=True)
        motif_df.sort_values('score', ascending=False, inplace=True)
        # Sometimes edges are identified from multiple transcripts
        motif_df.drop_duplicates(subset=['TFName', 'Gene stable ID'],
            inplace=True)
        print ('Number of TF - gene edges:', len(motif_df))
        if use_gene_names:
            # Names are not unique - filtering needed
            # Fill empty gene names with IDs
            motif_df['Gene name'] = motif_df.apply(lambda x: x['Gene name'] if 
                type(x['Gene name']) == str else x['Gene stable ID'], axis=1)
            # Count the number of edges for every name/ID pair
            name_id_matching = motif_df.groupby(
                ['Gene name', 'Gene stable ID'])['Gene name'].count()
            # Use the name for the ID that has the most edges
            id_to_name = {i[1]: i[0] for i in name_id_matching.groupby(
                level=0).idxmax().values}
            # Convert selected gene IDs to names
            motif_df['Gene name'] = motif_df['Gene stable ID'].apply(
                lambda x: id_to_name[x] if x in id_to_name else np.nan)
            # Drop the rest
            motif_df.dropna(subset='Gene name', inplace=True)
            print ('Number of TF - gene edges after name conversion:',
                len(motif_df))
        
        self.motif_frame = motif_df
        if (self.use_gene_names is not None and 
            self.use_gene_names != use_gene_names):
            # Notify that the provided setting for gene name use has changed
            print ('Changing the use_gene_names setting to ', use_gene_names)
        # Update the setting in any case, to prevent unwanted saving issues
        self.use_gene_names = use_gene_names
    

    def write_motif_prior(
        self,
        output_path: Optional[FILE_LIKE] = None,
        use_gene_names: Optional[bool] = None,
        weighted: Optional[bool] = None
    ) -> None:
        """
        Writes the motif (transcription factor - gene) prior network 
        into a file.

        Parameters
        ----------
        output_path : Optional[FILE_LIKE], optional
            The path to write the prior into or None to follow the 
            option from the initialisation, by default None
        use_gene_names : Optional[bool], optional
            Whether to use gene names instead of Ensembl IDs or None
            to follow the option from the initialisation, 
            by default None
        weighted : Optional[bool], optional
            Whether to use weights for the edges as opposed to making 
            them binary or None to follow the option from the 
            initialisation, by default None
        """
        
        print ()
        print ('--- Running write_motif_prior() ---')

        if self.motif_frame is None:
            print ('No motif prior has been generated yet, please run '
                'aggregate_matches() first')
            return
        
        if output_path is None:
            output_path = self.motif_outfile
        if weighted is None:
            weighted = self.weighted
        if use_gene_names is None:
            use_gene_names = self.use_gene_names
            
        if use_gene_names:
            column = 'Gene name'
        else:
            column = 'Gene stable ID'

        self.motif_frame.sort_values(by=['TFName', column], inplace=True)

        if weighted:
            # The scores are converted to floats < 10 (typically)
            self.motif_frame['weight'] = self.motif_frame['score'] / 100
            self.motif_frame[['TFName', column, 'weight']].to_csv(
                output_path, sep='\t', index=False, header=False)
        else:
            self.motif_frame['edge'] = 1
            self.motif_frame[['TFName', column, 'edge']].to_csv(
                output_path, sep='\t', index=False, header=False)


    def show_fingerprint(
        self
    ) -> None:
        """
        Prints the fingerprint for the files and databases used by this
        instance of the Sponge class.
        """
        
        for k,v in self.fingerprint.items():
            if v['provided']:
                # No information is known about this file as it was provided
                print (f'{k}: provided by user')
            elif v['cached']:
                # Some information may be available in the cache
                if v['version'] == '':
                    # Not available
                    print (f'{k}: retrieved from cache')
                else:
                    # Available
                    print (f'{k}: {v["version"]}, retrieved from cache,',
                        'originally retrieved', 
                        parse_datetime(v['datetime']))                      
            else:
                # File was retrieved by this instance, all info available
                print (f'{k}: {v["version"]}, retrieved',
                    parse_datetime(v['datetime']))


    def clear_cache(
        self
    ) -> None:
        """
        Removes the temporary folder and everything in it.
        """
        
        if os.path.exists(self.temp_folder):
            rmtree(self.temp_folder)