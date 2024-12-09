### Imports ###
import os

from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.utils import download_with_progress
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class MotifRetriever:
    _default_filename = 'jaspar.bb'

    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigReader,
        user_config: ConfigReader,
        version_logger: VersionLogger,
    ):
        
        # TODO: Announce what is being retrieved

        # Check if on the fly download was requested 
        # (in which case we don't do anything)
        if user_config.is_true('on_the_fly_processing'):
            # TODO: Report
            return
        temp_filename = os.path.join(temp_folder, self._default_filename)
        # Check existence of user-provided bigbed file
        if user_config.exists(['motif', 'tfbs_file']):
            # TODO: Report
            version_logger.write_provided('jaspar_bibged')
        # Check for a cached bigbed file
        elif os.path.exists(temp_filename):
            # TODO: Report
            version_logger.update_cached('jaspar_bigbed')
        # Retrieve a bigbed file from JASPAR
        else:
            # TODO: Report on file retrieval
            jaspar_release = user_config.get_value(['motif', 'jaspar_release'])
            year = jaspar_release[-4:]
            assembly = user_config.get_value('genome_assembly')
            urls_to_try = [url.format(year=year, genome_assembly=assembly)
                for url in core_config['url']['motif']['full']]
            download_with_progress(urls_to_try, temp_filename)
            version_logger.write_retrieved('jaspar_bigbed', jaspar_release)