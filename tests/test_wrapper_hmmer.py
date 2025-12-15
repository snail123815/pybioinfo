import subprocess
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

from pyBioinfo_modules.wrappers.hmmer import run_jackhmmer
from pyBioinfo_modules.wrappers.hmmer_config import (
    BaseHmmerTblFilters,
    HmmerHomologousProtFillters,
)


class TestRunJackhmmer(unittest.TestCase):
    """Test run_jackhmmer function with different configurations."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parent / "test_data"
        self.query_fasta = self.test_data_dir / "scp1_5prots.faa"
        self.target_fasta = self.test_data_dir / "scp1.faa"
        self.output_path = Path("/tmp/test_output.domtblout")

    @patch("pyBioinfo_modules.wrappers.hmmer.subprocess.run")
    def test_run_jackhmmer_with_default_config(self, mock_run):
        """Test jackhmmer with default BaseHmmerTblFilters configuration."""
        # Setup mock to simulate successful jackhmmer execution
        mock_run.return_value = MagicMock(returncode=0, stderr=b"")

        # Run with default config
        result = run_jackhmmer(
            self.query_fasta,
            self.target_fasta,
            self.output_path,
            jackhmmer_cfg=BaseHmmerTblFilters(),
            cpus=4,
        )

        # Verify subprocess was called
        self.assertTrue(mock_run.called)
        call_args = mock_run.call_args

        # Extract the command from the call
        cmd = call_args[0][0]

        # Verify the command contains the expected parameters from default config
        self.assertIn("-E 10", cmd)
        self.assertIn("--incE 0.001", cmd)
        self.assertIn("--domE 10", cmd)
        self.assertIn("--incdomE 0.001", cmd)
        self.assertIn("--cpu 4", cmd)
        self.assertIn(f"--domtblout {self.output_path}", cmd)
        self.assertIn(str(self.target_fasta), cmd)

        # Verify result
        self.assertEqual(result, 1)

    @patch("pyBioinfo_modules.wrappers.hmmer.subprocess.run")
    def test_run_jackhmmer_with_homologous_prot_filters(self, mock_run):
        """Test jackhmmer with HmmerHomologousProtFillters configuration."""
        # Setup mock to simulate successful jackhmmer execution
        mock_run.return_value = MagicMock(returncode=0, stderr=b"")

        # Run with homologous protein filters
        result = run_jackhmmer(
            self.query_fasta,
            self.target_fasta,
            self.output_path,
            jackhmmer_cfg=HmmerHomologousProtFillters,
            cpus=8,
        )

        # Verify subprocess was called
        self.assertTrue(mock_run.called)
        call_args = mock_run.call_args

        # Extract the command from the call
        cmd = call_args[0][0]

        # Verify the command contains the expected parameters from HmmerHomologousProtFillters
        self.assertIn("-E 1e-05", cmd)
        self.assertIn("--incE 1e-10", cmd)
        self.assertIn("--domE 1e-05", cmd)
        self.assertIn("--incdomE 1e-10", cmd)
        self.assertIn("--cpu 8", cmd)

        # Verify result
        self.assertEqual(result, 1)

    @patch("pyBioinfo_modules.wrappers.hmmer.subprocess.run")
    def test_run_jackhmmer_with_stringio_query(self, mock_run):
        """Test jackhmmer with StringIO query instead of Path."""
        # Setup mock
        mock_run.return_value = MagicMock(returncode=0, stderr=b"")

        # Create a StringIO query
        query_stringio = StringIO(">test_protein\nMKTAYIAKQRQISFVKSHFSRQLE\n")

        # Run with StringIO input
        result = run_jackhmmer(
            query_stringio,
            self.target_fasta,
            self.output_path,
            jackhmmer_cfg=BaseHmmerTblFilters(),
        )

        # Verify subprocess was called
        self.assertTrue(mock_run.called)

        # Verify the input was passed correctly (as encoded bytes)
        call_args = mock_run.call_args
        self.assertIn("input", call_args[1])
        self.assertIsInstance(call_args[1]["input"], bytes)

        # Verify result
        self.assertEqual(result, 1)

    @patch("pyBioinfo_modules.wrappers.hmmer.subprocess.run")
    def test_run_jackhmmer_failure(self, mock_run):
        """Test jackhmmer raises RuntimeError on failure."""
        # Setup mock to simulate failed execution
        mock_run.return_value = MagicMock(
            returncode=1, stderr=b"Error: some error"
        )

        # Verify that RuntimeError is raised
        with self.assertRaises(RuntimeError) as context:
            run_jackhmmer(
                self.query_fasta,
                self.target_fasta,
                self.output_path,
            )

        self.assertIn("Jackhmmer failed", str(context.exception))

    @patch("pyBioinfo_modules.wrappers.hmmer.subprocess.run")
    def test_run_jackhmmer_custom_config(self, mock_run):
        """Test jackhmmer with custom configuration values."""
        # Setup mock
        mock_run.return_value = MagicMock(returncode=0, stderr=b"")

        # Create custom config
        custom_config = BaseHmmerTblFilters(
            T_E=0.01,
            T_INCE=0.0001,
            T_DOME=0.05,
            T_INCDOME=0.0005,
        )

        # Run with custom config
        result = run_jackhmmer(
            self.query_fasta,
            self.target_fasta,
            self.output_path,
            jackhmmer_cfg=custom_config,
            cpus=16,
        )

        # Verify the custom parameters are in the command
        call_args = mock_run.call_args
        cmd = call_args[0][0]

        self.assertIn("-E 0.01", cmd)
        self.assertIn("--incE 0.0001", cmd)
        self.assertIn("--domE 0.05", cmd)
        self.assertIn("--incdomE 0.0005", cmd)
        self.assertIn("--cpu 16", cmd)

        # Verify result
        self.assertEqual(result, 1)


if __name__ == "__main__":
    unittest.main()
