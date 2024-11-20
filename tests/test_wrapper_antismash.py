import unittest
from unittest.mock import patch
from pyBioinfo_modules.wrappers.antismash import log_antismash_version
import subprocess


class Test_antiSmash(unittest.TestCase):

    @patch("pyBioinfo_modules.wrappers.antismash.subprocess.run")
    @patch("pyBioinfo_modules.wrappers.antismash.logger")
    def test_log_antismash_version(self, mock_logger, mock_run):
        mock_run.return_value.stdout = b"antiSMASH 6.0.1\n"
        version = log_antismash_version(
            condaexe="conda", antismash_env="env", shell="bash"
        )
        self.assertEqual(version, "antiSMASH 6.0.1\n")
        mock_run.assert_called_once()
        mock_logger.info.assert_any_call("antiSMASH version: 6.0.1")

    @patch("pyBioinfo_modules.wrappers.antismash.subprocess.run")
    def test_log_antismash_version_command_failed(self, mock_run):
        mock_run.side_effect = subprocess.CalledProcessError(1, 'cmd')
        with self.assertRaises(subprocess.CalledProcessError):
            log_antismash_version(
                condaexe="conda", antismash_env="env", shell="bash"
            )
        mock_run.assert_called_once()


if __name__ == "__main__":
    unittest.main()
