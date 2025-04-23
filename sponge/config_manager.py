### Imports ###
import os
import yaml

from pathlib import Path
from typing import List, Optional, Union

from sponge.modules.utils import recursive_update

### Class definition ###
class ConfigManager:
    # Variables
    _default_core_config = 'config.yaml'
    _default_user_config = 'user_config.yaml'
    _log_filename = 'last_user_config.yaml'

    # Methods
    def __init__(
        self,
        config: Union[Path, dict, None] = None,
        temp_folder: Optional[Path] = None,
    ):
        """
        

        Parameters
        ----------
        config : Union[Path, dict, None], optional
            _description_, by default None
        temp_folder : Optional[Path], optional
            _description_, by default None
        """

        self.temp_folder = temp_folder
        
        self.file_dir = Path(__file__).parents[0]
        # Handle the special case of core config
        if config is None:
            # Find the core config file in the module directory
            config = os.path.join(self.file_dir, self._default_core_config)
            self.config = {}
        else:
            # Load the default user configuration
            self._load_default()
        # Now handle based on the update type and content
        if type(config) == dict:
            # Directly use the provided dictionary
            self._deep_update(config)
        elif not os.path.isfile(config):
            # Use the default user config file - don't change any defaults
            print ('Using the default settings.')
        else:
            # Update with provided config (core or user)
            with open(config, 'r') as f:
                update = yaml.safe_load(f)
            self._deep_update(update)

    # Overwritten built-in methods
    def __del__(
        self,
    ):
        """
        When the object is deleted, if a temporary folder was specified,
        the managed configuration is saved there.
        """

        if self.temp_folder is not None:
            log_file = os.path.join(self.temp_folder, self._log_filename)
            yaml.safe_dump(self.config, open(log_file, 'w', encoding='utf-8'))


    def __getitem__(
        self,
        key: str,
    ) -> Union[dict, str]:
        """
        Accesses the given key in the managed configuration.

        Parameters
        ----------
        key : str
            Key to access

        Returns
        -------
        Union[dict, str]
            Value referred to by the given key
        """

        return self.config[key]


    def __setitem__(
        self,
        key: str,
        val: Union[dict, str],
    ) -> None:
        """
        Sets the value for the given key in the managed configuration.

        Parameters
        ----------
        key : str
            Key to change
        val : Union[dict, str]
            Value to set
        """

        self.config[key] = val


    def __delitem__(
        self,
        key: str,
    ) -> None:
        """
        Deletes the given key from the managed configuration.

        Parameters
        ----------
        key : str
            Key to delete
        """

        del self.config[key]


    def __contains__(
        self,
        key: str,
    ) -> bool:
        """
        Checks if the given key exists in the managed configuration.

        Parameters
        ----------
        key : str
            Key to verify

        Returns
        -------
        bool
            Whether the given key exists in the managed configuration
        """

        return key in self.config

    # Rest of the methods
    def _load_default(
        self,
    ) -> None:
        """
        Loads the default user configuration into the managed
        configuration.
        """
        
        config = os.path.join(self.file_dir, self._default_user_config)
        with open(config, 'r') as f:
            self.config = yaml.safe_load(f)


    def _deep_update(
        self,
        config: dict,
    ) -> None:
        """
        Updates the managed configuration with a dictionary, recursively
        so that multiple level updates are possible.

        Parameters
        ----------
        config : dict
            Configuration update
        """

        # Call the above function to update the internal config
        self.config = recursive_update(self.config, config)


    def _retrieve_level(
        self,
        keys: List[str],
    ) -> Union[dict, str]:

        curr = self.config
        for k in keys:
            curr = curr[k]

        return curr


    def exists(
        self,
        key: Union[str, List[str]],
    ) -> bool:

        if type(key) == str:
            key = [key]

        try:
            last_level = self._retrieve_level(key[:-1])
        except KeyError:
            return False

        return key[-1] in last_level


    def get_value(
        self,
        key: Union[str, List[str]],
    ) -> str:

        if type(key) == str:
            key = [key]

        # Can throw KeyErrors on purpose
        last_level = self._retrieve_level(key[:-1])

        return last_level[key[-1]]