import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
from pyBioinfo_modules.wrappers.bigscape import get_bigscape_version


class TestBigscape(unittest.TestCase):

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_with_date(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = (
            "BiG-SCAPE 1.2.3 (2023-01-01)\n"
        )
        mock_run.return_value = mock_process

        version, date = get_bigscape_version()
        self.assertEqual(version, "1.2.3")
        self.assertEqual(date, "(2023-01-01)")

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_without_date(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = "BiG-SCAPE 1.2.3\n"
        mock_run.return_value = mock_process

        version, date = get_bigscape_version()
        self.assertEqual(version, "1.2.3")
        self.assertIsNone(date)

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_unexpected_format(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = "unexpected format\n"
        mock_run.return_value = mock_process

        version, date = get_bigscape_version()
        self.assertEqual(version, "unexpected format")
        self.assertIsNone(date)

        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = (
            "bigscape 1.2.3 (2023-01-01) unexpected\n"
        )
        mock_run.return_value = mock_process

        version, date = get_bigscape_version()
        self.assertEqual(version, "1.2.3")
        self.assertIsNone(date)

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    @patch("pyBioinfo_modules.wrappers.bigscape.logger")
    def test_get_bigscape_version_logging(self, mock_logger, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = (
            "BiG-SCAPE 1.2.3 (2023-01-01)\n"
        )
        mock_run.return_value = mock_process

        version, date = get_bigscape_version()
        mock_logger.info.assert_any_call(
            "BiG-SCAPE version: 1.2.3, (2023-01-01)"
        )


if __name__ == "__main__":
    unittest.main()
