### Imports ###
import os

from pathlib import Path

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever.file_retriever import FileRetriever
from sponge.modules.utils import download_with_progress

### Class definition ###
class TFBSRetriever(FileRetriever):
    _default_filename = 'tfbs.bb'

    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigManager,
        user_config: ConfigManager,
    ):

        self.on_the_fly = user_config['on_the_fly_processing']

        path_to_file = None
        if user_config.exists(['motif', 'tfbs_file']):
            path_to_file = user_config.get_value(['motif', 'tfbs_file'])
        temp_filename = os.path.join(temp_folder, self._default_filename)

        self.jaspar_release = user_config.get_value(
            ['motif', 'jaspar_release'])
        self.genome_assembly = user_config.get_value('genome_assembly')
        self.urls = core_config['url']['motif']['full']

        super().__init__(
            key='tfbs_file',
            temp_filename=temp_filename,
            path_to_file=path_to_file,
        )


    def _retrieve_tfbs(
        self,
    ) -> str:

        year = self.jaspar_release[-4:]
        urls_to_try = [url.format(year=year,
            genome_assembly=self.genome_assembly) for url in self.urls]
        download_with_progress(urls_to_try, self.temp_filename)

        return self.jaspar_release


    def retrieve_file(
        self,
    ) -> None:

        if self.on_the_fly:
            print ('Retrieval of tfbs_file is skipped as on the fly '
                'processing was requested.')
            self.actual_path = None
        else:
            super().retrieve_file(self._retrieve_tfbs)