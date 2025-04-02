### Imports ###
import os
import yaml

from pathlib import Path
from typing import Any, List, Optional, Union

### Class definition ###
class ConfigManager:
    # Variables
    _default_core_config = 'config.yaml'
    _default_user_config = 'user_config.yaml'
    _log_filename = 'last_user_config.yaml'

    # Functions
    def __init__(
        self,
        config_path: Optional[Path] = None,
        temp_folder: Optional[Path] = None,
    ):

        self.temp_folder = temp_folder

        file_dir = Path(__file__).parents[0]
        if config_path is None:
            # Find the config file in the module directory
            config_path = os.path.join(file_dir, self._default_core_config)
        if not os.path.isfile(config_path):
            # Use the default user config file
            config_path = os.path.join(file_dir, self._default_user_config)

        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)


    def __del__(
        self,
    ):

        if self.temp_folder is not None:
            log_file = os.path.join(self.temp_folder, self._log_filename)
            yaml.safe_dump(self.config, open(log_file, 'w', encoding='utf-8'))


    def __getitem__(
        self,
        key: str,
    ) -> dict:

        return self.config[key]


    def __setitem__(
        self,
        key: str,
        val: Any,
    ) -> None:

        self.config[key] = val


    def __delitem__(
        self,
        key: str,
    ) -> None:

        del self.config[key]


    def __contains__(
        self,
        key: str,
    ) -> bool:

        return key in self.config


    def _retrieve_level(
        self,
        keys: List[str],
    ) -> dict:

        curr = self.config
        for k in keys:
            curr = curr[k]

        return curr


    def is_true(
        self,
        key: Union[str, List[str]],
    ) -> bool:
        """
        Return if the key is present and evaluates to True. It is mostly
        useful for parameters which are implied False by default.

        Parameters
        ----------
        key : Union[str, List[str]]
            Key to be investigated, nested key can be provided as a list

        Returns
        -------
        bool
            Whether the key exists and is True
        """

        return self.exists(key) and bool(self.get_value(key))


    def is_false(
        self,
        key: Union[str, List[str]],
    ) -> bool:
        """
        Return if the key is present and evaluates to False. It is
        mostly useful for parameters which are implied True by default.

        Parameters
        ----------
        key : Union[str, List[str]]
            Key to be investigated, nested key can be provided as a list

        Returns
        -------
        bool
            Whether the key exists and is False
        """

        return self.exists(key) and not bool(self.get_value(key))


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