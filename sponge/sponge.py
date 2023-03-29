### Imports ###
import pandas as pd
import numpy as np
import bioframe
import time

from pyjaspar import jaspardb

from multiprocessing import Pool

from math import ceil, sqrt

from sponge.functions import *

### Class definition ###
class Sponge:
    def __init__(
        self,
        temp_folder: str = '.sponge_temp/'
    ):
        """
        Usual order of operations:
        sponge_obj = Sponge()
        sponge_obj.select_tfs()
        sponge_obj.find_human_homologs()
        sponge_obj.filter_matches()
        sponge_obj.retrieve_ppi()
        sponge_obj.write_ppi_prior()
        sponge_obj.aggregate_matches()    
        sponge_obj.write_motif_prior()

        Parameters
        ----------
        temp_folder : str, optional
            _description_, by default '.sponge_temp/'
        """

        self.temp_folder = temp_folder
        self.ensembl = None


    def select_tfs(
        self,
        drop_heterodimers: bool = True
    ) -> None:

        # Database object
        jdb_obj = jaspardb()
        self.release = jdb_obj.release
        print ('Using:', self.release)

        # All vertebrate motifs
        motifs = jdb_obj.fetch_motifs(collection='CORE', 
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


    def load_promoters_from_biomart(
        self,
        file_path: str,
        filter_basic: bool = True,
        chromosomes: Iterable[str] = [str(i) for i in range(1,23)] + 
            ['MT', 'X', 'Y'],
        keep_ensembl: bool = True
    ) -> None:

        bm_server = BiomartServer(ENSEMBL_URL)
        ensembl = bm_server.datasets['hsapiens_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'transcript_gencode_basic', 
            'chromosome_name', 'transcription_start_site', 'strand']
        if keep_ensembl:
            attributes += ['ensembl_gene_id', 'external_gene_name', 'gene_biotype']
        print ('Retrieving response to query...')
        response = ensembl.search({'attributes': attributes}, header=1)
        buffer = download_with_progress(response)
        dtype_dict = defaultdict(lambda: str)
        dtype_dict['Transcription start site (TSS)'] = int
        dtype_dict['Strand'] = int
        df = pd.read_csv(buffer, sep='\t', dtype=dtype_dict)
        print ('Filtering and modifying dataframe...')
        if filter_basic:
            df = df[df['GENCODE basic annotation'] == 'GENCODE basic'].copy()
            df.drop(columns='GENCODE basic annotation', inplace=True)
        if chromosomes is not None:
            df = df[df['Chromosome/scaffold name'].isin(chromosomes)]
        df['Chromosome'] = df['Chromosome/scaffold name'].apply(lambda x: 'chrM' 
            if x == 'MT' else f'chr{x}')
        df['Strand'] = df['Strand'].apply(lambda x: '+' if x > 0 else '-')
        df['Start'] = df.apply(lambda row: 
            row['Transcription start site (TSS)'] - 750 if row['Strand'] == '+' 
            else row['Transcription start site (TSS)'] - 250, axis=1)
        df['End'] = df['Start'] + 1000
        df['Score'] = 0
        df.sort_values(['Chromosome', 'Start'], inplace=True)
        columns = ['Chromosome', 'Start', 'End', 'Transcript stable ID', 'Score', 
            'Strand']
        print (f'Saving data to {file_path}...')
        df[columns].to_csv(file_path, sep='\t', header=False, index=False)
        print ()
        if keep_ensembl:
            self.ensembl = df[['Gene stable ID', 'Transcript stable ID', 
                'Gene name', 'Gene type']]


    def load_ensembl_from_biomart(
        self,
        file_path: str
    ) -> None:
        
        bm_server = BiomartServer(ENSEMBL_URL)
        ensembl = bm_server.datasets['hsapiens_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 
            'external_gene_name', 'gene_biotype']
        print ('Retrieving response to query...')
        response = ensembl.search({'attributes': attributes}, header=1)
        buffer = download_with_progress(response)
        df = pd.read_csv(buffer, sep='\t')     
        df.to_csv(file_path, sep='\t', index=False)
        self.ensembl = df


    def retrieve_file(
        self,
        description: str,
        temp_folder: str,
        prompt: bool = True,
        jaspar_release: Optional[str] = None
    ) -> str:

        
        if description not in FILE_DF.index:
            print (f'File description not recognised: {description}')
            return
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)
        file_name = FILE_DF.loc[description, 'name']
        file_path = os.path.join(temp_folder, file_name)
        if os.path.exists(file_path):
            print (f'Using cached file {file_path}')
            print ()
        else:
            print (f'File {file_name} not found in directory {temp_folder}')
            if prompt:
                key = None
                positive = ['y', 'yes', 'hell yeah']
                negative = ['n', 'no', 'nope']
                while key is None or key.lower() not in positive + negative:
                    if key is not None:
                        print (f'Input not recognised: {key}')
                    print ('Do you want to download it? Y/N', flush=True)
                    key = input()
                    print (key)
                print ()
                if key.lower() in negative:
                    return None
            to_request = FILE_DF.loc[description, 'url']
            if to_request is None:
                # These options are not unused: they are passed to the
                # evaluated function call as kwargs
                options = {'file_path': file_path}
                eval(FILE_DF.loc[description, 'eval'])
            else:
                if description == 'jaspar_bigbed':
                    if jaspar_release is None:
                        raise ValueError('The release of jaspar has to be '
                            'specified in order to retrieve the bigbed file')
                    to_request = to_request.format(year=jaspar_release[-4:])
                print (f'Downloading data into {file_path}...')
                download_with_progress(to_request, file_path)    

        return file_path


    def find_human_homologs(
        self, 
        homologene_file: Optional[str] = None,
        prompt: bool = True
    ) -> None:
        
        if homologene_file is None:
            homologene_file = self.retrieve_file('homologene', 
                self.temp_folder, prompt=prompt)
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
        found_names = hg_df[hg_df['Gene Symbol'].isin([adjust_gene_name(i) for i in 
            non_human_motif_names])]['Gene Symbol'].unique()
        # Find the missing ones
        missing = (set([adjust_gene_name(i) for i in non_human_motif_names]) - 
            set(found_names))
        print ()
        print ('Names missing from the homologene database:')
        for i in [(i.name,i.acc) for i in non_human_motifs if 
            i.name in missing]:
            print (i[0], *i[1])

        # Get the missing IDs from Uniprot API
        print ()
        print ('Retrieving matches from UniProt...')
        mapping = get_uniprot_mapping("UniProtKB_AC-ID", "RefSeq_Protein",
            [i.acc[0] for i in non_human_motifs if i.name in missing])
        mapping.columns = ['Uniprot', 'Accession']

        # Create a DataFrame for matching missing entries
        hg_df[hg_df['Protein Accession'].isin(mapping['Accession'])]
        missing_df = pd.DataFrame([(i.name, i.acc[0]) for i in non_human_motifs 
            if i.name in missing], columns = ['Gene Symbol', 'Uniprot'])        
        matching_df = missing_df.join(mapping.set_index('Uniprot'), 
            on='Uniprot').join(hg_df.set_index('Protein Accession'), 
            on='Accession', rsuffix='_HG')

        def corresponding_id(name):
            values = hg_df[hg_df['Gene Symbol'] == name]['HG Group ID'].values
            if len(values) == 0:
                matches = matching_df[matching_df['Gene Symbol'] == name]
                values = matches['HG Group ID'].values
            return values

        # Create a DataFrame of corresponding names
        corr_df = pd.DataFrame(non_human_motif_names, 
            columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(adjust_gene_name)
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

        self.matrix_ids = matrix_ids
        self.tf_names = human_motif_names + list(animal_to_human.keys())
        self.human_names = human_motif_names + list(animal_to_human.values())


    def filter_matches(
        self, 
        promoter_file: Optional[str] = None, 
        bigbed_file: Optional[str] = None,
        n_processes: int = 1,
        prompt: bool = True
    ) -> None:
        
        if promoter_file is None:
            promoter_file = self.retrieve_file('promoter', self.temp_folder, 
                prompt=prompt)
            if promoter_file is None:
                print ('Unable to find or retrieve the promoter file, exiting')
                return
        
        if bigbed_file is None:
            bigbed_file = self.retrieve_file('jaspar_bigbed', self.temp_folder, 
                prompt=prompt, jaspar_release=self.release)
            if bigbed_file is None:
                print ('Unable to find or retrieve the JASPAR bigbed file, '
                    'exiting')
                return

        print ('Loading the promoter bed file...')
        df_full = bioframe.read_table(promoter_file, schema='bed')
        df_full['name'] = df_full['name'].apply(lambda x: x.split('.')[0])
        df_full.drop(columns=['score', 'strand'], inplace=True)
        df_full.set_index('name', inplace=True)

        chr_list = ['chr{}'.format(i) for i in [j for j in range(1,23)] + 
            ['M', 'X', 'Y']]
        results_list = []
        p = Pool(n_processes)

        print ()
        print ('Iterating over the chromosomes...')
        start_time = time.time()
        for chrom in chr_list:
            st_chr = time.time()
            df_chrom = df_full[df_full['chrom'] == chrom]
            if len(df_chrom) == 0:
                suffix = 'no transcripts'
            elif len(df_chrom) == 1:
                suffix = '1 transcript'
            else:
                suffix = f'{len(df_chrom)} transcripts'
            print (f'Chromosome {chrom[3:]} with ' + suffix)
            chunk_size = ceil(sqrt(len(df_chrom) / n_processes))
            chunk_divisions = [i for i in range(0, len(df_chrom), chunk_size)]
            input_tuples = [(bigbed_file, df_chrom, self.tf_names, chrom, i, 
                i+chunk_size) for i in chunk_divisions]
            result = p.map_async(filter_edges_helper, input_tuples, 
                chunksize=n_processes)
            edges_chrom_list = result.get()
            results_list += edges_chrom_list
            elapsed_chr = time.time() - st_chr
            print (f'Done in: {elapsed_chr // 60:n} m '
                f'{elapsed_chr % 60:.2f} s')

        elapsed = time.time() - start_time
        print ()
        print (f'Total time: {elapsed // 60:n} m {elapsed % 60:.2f} s')

        self.all_edges = pd.concat(results_list, ignore_index=True)

    
    def retrieve_ppi(
        self
    ) -> None:
        
        query_string = '%0d'.join(self.human_names)

        print ('Retrieving mapping from STRING...')
        mapping_request = requests.get(f'{STRING_URL}get_string_ids?'
            f'identifiers={query_string}&species=9606')
        mapping_df = pd.read_csv(BytesIO(mapping_request.content), sep='\t')
        mapping_df['queryName'] = mapping_df['queryIndex'].apply(
            lambda i: self.human_names[i])
        diff_df = mapping_df[mapping_df['queryName'] != 
            mapping_df['preferredName']]
        ids_to_check = np.concatenate((diff_df['queryName'], 
            diff_df['preferredName']))
        
        print ('Checking the conflicts in the UniProt database...')
        uniprot_df = get_uniprot_mapping('Gene_Name', 'UniProtKB', 
            ids_to_check).set_index('from')
        matching_ids = list(mapping_df[mapping_df['queryName'] == 
            mapping_df['preferredName']]['preferredName'])
        p_to_q = {p: q for q,p in zip(diff_df['queryName'], 
            diff_df['preferredName'])}
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
        output_path: str = 'ppi_prior.tsv',
        weighted: bool = False
    ) -> None:
        
        if weighted:
            self.ppi_frame[['tf1', 'tf2', 'score']].to_csv(output_path, 
                sep='\t', index=False, header=False)
        else:
            self.ppi_frame['edge'] = 1
            self.ppi_frame[['tf1', 'tf2', 'edge']].to_csv(output_path, 
                sep='\t', index=False, header=False)
            

    def aggregate_matches(
        self,
        ensembl_file: Optional[str] = None,
        prompt: bool = True,
        use_gene_names: bool = True
    ) -> None:
        
        if self.ensembl is None or ensembl_file is not None:
            if ensembl_file is None:
                ensembl_file = self.retrieve_file('ensembl', 
                    self.temp_folder, prompt=prompt)
                if ensembl_file is None:
                    print ('Unable to find or retrieve the ensembl file, ' 
                        'exiting')
                    return
            self.ensembl = pd.read_csv(ensembl_file, sep='\t')
        
        


    def write_motif_prior(
        self
    ) -> None:
        
        pass


    def clear_cache(
        self
    ) -> None:
        
        if os.path.exists(self.temp_folder):
            for file in os.listdir(self.temp_folder):
                os.remove(os.path.join(self.temp_folder, file))
            os.rmdir(self.temp_folder)