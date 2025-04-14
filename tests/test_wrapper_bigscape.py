import os
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from pyBioinfo_modules.wrappers.bigscape import BigscapeWrapper


class TestBigscapeUnit(unittest.TestCase):
    """Unit tests that don't require an actual BigScape executable."""

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_with_date(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = (
            "BiG-SCAPE 1.2.3 (2023-01-01)\n"
        )
        mock_run.return_value = mock_process

        wrapper = BigscapeWrapper(bigscape_exe=Path("/mock/path/bigscape"))
        version, date = wrapper.get_version()

        self.assertEqual(version, "1.2.3")
        self.assertEqual(date, "(2023-01-01)")

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_without_date(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = "BiG-SCAPE 1.2.3\n"
        mock_run.return_value = mock_process

        wrapper = BigscapeWrapper(bigscape_exe=Path("/mock/path/bigscape"))
        version, date = wrapper.get_version()

        self.assertEqual(version, "1.2.3")
        self.assertIsNone(date)

    @patch("pyBioinfo_modules.wrappers.bigscape.subprocess.run")
    def test_get_bigscape_version_unexpected_format(self, mock_run):
        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = "unexpected format\n"
        mock_run.return_value = mock_process

        wrapper = BigscapeWrapper(bigscape_exe=Path("/mock/path/bigscape"))
        version, date = wrapper.get_version()

        self.assertEqual(version, "unexpected format")
        self.assertIsNone(date)

        mock_process = MagicMock()
        mock_process.stdout.decode.return_value = (
            "bigscape 1.2.3 (2023-01-01) unexpected\n"
        )
        mock_run.return_value = mock_process

        wrapper = BigscapeWrapper(bigscape_exe=Path("/mock/path/bigscape"))
        version, date = wrapper.get_version()

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

        wrapper = BigscapeWrapper(bigscape_exe=Path("/mock/path/bigscape"))
        version, date = wrapper.get_version()

        mock_logger.info.assert_any_call(
            "BiG-SCAPE version: 1.2.3, (2023-01-01)"
        )


@unittest.skipIf(
    os.environ.get("SKIP_INTEGRATION_TESTS"), "Skipping integration tests"
)
class TestBigscapeIntegration(unittest.TestCase):
    """Integration tests that require an actual BigScape executable."""

    def setUp(self):
        self.wrapper = BigscapeWrapper()
        try:
            # Try to access the executable
            self.wrapper.get_bigscape_exe()
            self.bigscape_available = True
        except FileNotFoundError:
            self.bigscape_available = False

    def test_get_bigscape_exe(self):
        if not self.bigscape_available:
            self.skipTest("BigScape executable not available")

        exe_path = self.wrapper.get_bigscape_exe()
        self.assertTrue(exe_path.exists())
        self.assertTrue(
            str(exe_path).endswith("bigscape")
            or str(exe_path).endswith("bigscape.py")
        )

    def test_get_version(self):
        if not self.bigscape_available:
            self.skipTest("BigScape executable not available")

        version, date = self.wrapper.get_version()
        # Just verify we get something back without errors
        self.assertIsInstance(version, str)


if __name__ == "__main__":
    unittest.main()
