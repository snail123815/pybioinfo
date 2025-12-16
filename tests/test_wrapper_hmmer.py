import subprocess
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

from pyBioinfo_modules.wrappers.hmmer import (
    calculate_domain_coverage,
    full_proteome_jackhmmer,
    parse_dom_table_mt,
    read_domtbl,
    run_jackhmmer,
)
from pyBioinfo_modules.wrappers.hmmer_config import (
    BaseHmmerTblFilters,
    HmmerHomologousProtFillters,
)


class TestRunJackhmmer(unittest.TestCase):
    """Test run_jackhmmer function with different configurations."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parent / "test_data" / "hmmer"
        self.query_fasta = self.test_data_dir / "queries.faa"
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


class TestFullProteomeJackhmmer(unittest.TestCase):
    """Test batching and merge behavior of full_proteome_jackhmmer."""

    def setUp(self):
        """Set up test fixtures."""
        self.hmmer_dir = Path(__file__).parent / "test_data" / "hmmer"
        self.query_fasta = self.hmmer_dir / "queries.faa"
        self.target_fasta = self.hmmer_dir / "scp1.faa"

    @patch("pyBioinfo_modules.wrappers.hmmer.run_jackhmmer")
    def test_split_by_two_merges_three_chunks(self, mock_run):
        """Split queries into groups of 2 and merge three domtbl chunks."""

        fixtures = [
            self.hmmer_dir / "jackhmmer_test_split2-1.domtblout",
            self.hmmer_dir / "jackhmmer_test_split2-2.domtblout",
            self.hmmer_dir / "jackhmmer_test_split2-3.domtblout",
        ]

        def fake_run(group_io, target, domtblout_path, jackhmmer_cfg, cpus):
            fixture_path = fixtures.pop(0)
            domtblout_path.write_text(fixture_path.read_text())
            return 1

        mock_run.side_effect = fake_run

        with tempfile.NamedTemporaryFile(delete=False) as tmp_out:
            out_path = Path(tmp_out.name)
        try:
            full_proteome_jackhmmer(
                self.query_fasta,
                self.target_fasta,
                out_path,
                prots_per_group=2,
                cpus=2,
            )

            self.assertEqual(mock_run.call_count, 3)

            merged = out_path.read_text()
            expected = (self.hmmer_dir / "merged_split2.domtblout").read_text()
            self.assertEqual(merged, expected)
        finally:
            out_path.unlink(missing_ok=True)

    @patch("pyBioinfo_modules.wrappers.hmmer.run_jackhmmer")
    def test_split_by_four_merges_two_chunks(self, mock_run):
        """Split queries into groups of 4 and merge two domtbl chunks."""

        fixtures = [
            self.hmmer_dir / "jackhmmer_test_split4-1.domtblout",
            self.hmmer_dir / "jackhmmer_test_split4-2.domtblout",
        ]

        def fake_run(group_io, target, domtblout_path, jackhmmer_cfg, cpus):
            fixture_path = fixtures.pop(0)
            domtblout_path.write_text(fixture_path.read_text())
            return 1

        mock_run.side_effect = fake_run

        with tempfile.NamedTemporaryFile(delete=False) as tmp_out:
            out_path = Path(tmp_out.name)
        try:
            full_proteome_jackhmmer(
                self.query_fasta,
                self.target_fasta,
                out_path,
                prots_per_group=4,
                cpus=2,
            )

            self.assertEqual(mock_run.call_count, 2)

            merged = out_path.read_text()
            expected = (self.hmmer_dir / "merged_split4.domtblout").read_text()
            self.assertEqual(merged, expected)
        finally:
            out_path.unlink(missing_ok=True)


class TestCalculateDomainCoverage(unittest.TestCase):
    """Test calculate_domain_coverage function."""

    def setUp(self):
        self.hmmer_filters = BaseHmmerTblFilters()

    def test_single_domain_above_threshold(self):
        """Test single domain with E-value above threshold (passes filter)."""
        line = {
            "dom_total": 1,
            "dom_i_E": 1e-20,
            "ali_from": 50,
            "ali_to": 100,
            "qlen": 200,
            "tlen": 150,
            "dom_n": 1,
        }
        expected_is_end = True
        expected_locations = list(range(50, 101))
        expected_covq = 51 / 200  # (100 - 50 + 1) / 200 = 0.255
        expected_covt = 51 / 150  # (100 - 50 + 1) / 150 = 0.34

        is_end, locations, covq, covt = calculate_domain_coverage(
            line, [], self.hmmer_filters
        )

        self.assertEqual(is_end, expected_is_end)
        self.assertEqual(locations, expected_locations)
        self.assertAlmostEqual(covq, expected_covq)
        self.assertAlmostEqual(covt, expected_covt)

    def test_single_domain_below_threshold(self):
        """Test single domain with E-value below threshold (fails filter)."""
        line = {
            "dom_total": 1,
            "dom_i_E": 1.1,  # GATHER_T_DOME is 1 in BaseHmmerTblFilters()
            "ali_from": 50,
            "ali_to": 100,
            "qlen": 200,
            "tlen": 150,
            "dom_n": 1,
        }
        expected_is_end = True
        expected_locations = []
        expected_covq = 0.0
        expected_covt = 0.0

        is_end, locations, covq, covt = calculate_domain_coverage(
            line, [], self.hmmer_filters
        )

        self.assertEqual(is_end, expected_is_end)
        self.assertEqual(locations, expected_locations)
        self.assertEqual(covq, expected_covq)
        self.assertEqual(covt, expected_covt)

    def test_multiple_domains(self):
        """Test multiple domains with varying E-values."""
        # Domain 1: dom_total=3, dom_n=1, dom_i_E=1.1 (fails filter, E > 1)
        line1 = {
            "dom_total": 3,
            "dom_i_E": 1.1,
            "ali_from": 10,
            "ali_to": 50,
            "qlen": 300,
            "tlen": 250,
            "dom_n": 1,
        }
        is_end1, locations1, covq1, covt1 = calculate_domain_coverage(
            line1, [], self.hmmer_filters
        )
        self.assertFalse(is_end1)
        self.assertEqual(locations1, [])
        self.assertEqual(covq1, 0.0)
        self.assertEqual(covt1, 0.0)

        # Domain 2: dom_total=3, dom_n=2, dom_i_E=0.2 (passes filter, E <= 1)
        line2 = {
            "dom_total": 3,
            "dom_i_E": 0.2,
            "ali_from": 60,
            "ali_to": 100,
            "qlen": 300,
            "tlen": 250,
            "dom_n": 2,
        }
        is_end2, locations2, covq2, covt2 = calculate_domain_coverage(
            line2, locations1, self.hmmer_filters
        )
        self.assertFalse(is_end2)
        self.assertEqual(locations2, list(range(60, 101)))
        self.assertEqual(covq2, 0.0)
        self.assertEqual(covt2, 0.0)

        # Domain 3: dom_total=3, dom_n=3, dom_i_E=0.1 (passes filter, is last domain)
        line3 = {
            "dom_total": 3,
            "dom_i_E": 0.1,
            "ali_from": 90,
            "ali_to": 180,
            "qlen": 300,
            "tlen": 250,
            "dom_n": 3,
        }
        is_end3, locations3, covq3, covt3 = calculate_domain_coverage(
            line3, locations2, self.hmmer_filters
        )
        self.assertTrue(is_end3)
        # Combined locations from domain 2 and 3: 60-100 (41 positions) + 90-180 (91 positions)
        expected_locations3 = list(range(60, 101)) + list(range(90, 181))
        self.assertEqual(locations3, expected_locations3)
        # Coverage: 121 unique positions / 300 query length
        self.assertAlmostEqual(covq3, 121 / 300)
        # Coverage: 121 unique positions / 250 target length
        self.assertAlmostEqual(covt3, 121 / 250)


class TestReadDomtbl(unittest.TestCase):
    """Test read_domtbl function."""

    def setUp(self):
        self.hmmer_dir = Path(__file__).parent / "test_data" / "hmmer"
        self.domtbl_file = self.hmmer_dir / "merged_split2.domtblout"

    @patch("pyBioinfo_modules.wrappers.hmmer.logger")
    def test_read_domtbl_logs_line_count(self, mock_logger):
        """Test that read_domtbl logs the correct line count."""
        read_domtbl(self.domtbl_file, BaseHmmerTblFilters())

        # Verify logger.info was called
        self.assertTrue(mock_logger.info.called)

        # Find the call that logs the line count
        line_count_logged = False
        for call in mock_logger.info.call_args_list:
            args, kwargs = call
            if args and "10 lines in total" in str(args[0]):
                line_count_logged = True
                break

        self.assertTrue(
            line_count_logged,
            "Expected log message '13 lines in total' not found",
        )

    def test_read_domtbl_filters(self):
        """Test reading domtbl file with default filters (no filtering)."""
        filters = BaseHmmerTblFilters(
            T_E=1e-50,
            LEN_DIFF=0.3,
        )
        df = read_domtbl(self.domtbl_file, filters)

        # Number of lines pass the intermediate filters
        self.assertEqual(df.shape[0], 3)

        # Check that all expected columns are present
        expected_columns = [
            "tp",
            "qp",
            "tlen",
            "qlen",
            "ali_from",
            "ali_to",
            "dom_i_E",
            "dom_n",
            "dom_total",
            "full_E",
            "anno",
        ]
        for col in expected_columns:
            self.assertIn(col, df.columns)


class TestParseDomTableMt(unittest.TestCase):
    """Test parse_dom_table_mt function."""

    def setUp(self):
        self.hmmer_dir = Path(__file__).parent / "test_data" / "hmmer"
        self.domtbl_file = self.hmmer_dir / "merged_split2.domtblout"

    def test_parse_dom_table_mt_with_filters(self):
        """Test parsing domtbl with specific filters."""
        filters = BaseHmmerTblFilters(
            T_E=1e-50,
            LEN_DIFF=0.3,
        )

        # Read and filter the domtbl
        df = read_domtbl(self.domtbl_file, filters)
        self.assertEqual(len(df), 3)  # Should have 4 rows after filtering

        # Parse the domtbl
        result = parse_dom_table_mt(df, cpus=2, hmmer_filters=filters)

        # Check number of results
        self.assertEqual(len(result), 3)

        # Check expected columns
        expected_columns = [
            "Query",
            "Target protein",
            "Coverage on Query",
            "Coverage on Target",
            "Expect protein",
            "Target description",
        ]
        for col in expected_columns:
            self.assertIn(col, result.columns)

        # Check unique queries
        unique_queries = result["Query"].unique()
        self.assertEqual(len(unique_queries), 2)
        self.assertIn("WP_011026846.1", unique_queries)
        self.assertIn("WP_011029008.1", unique_queries)

        # Check coverage values are between 0 and 1
        self.assertTrue((result["Coverage on Query"] >= 0).all())
        self.assertTrue((result["Coverage on Query"] <= 1).all())
        self.assertTrue((result["Coverage on Target"] >= 0).all())
        self.assertTrue((result["Coverage on Target"] <= 1).all())

        # Check specific row for WP_011029008.1 (self-match with 100% coverage)
        self_match = result[result["Query"] == "WP_011029008.1"]
        self.assertEqual(len(self_match), 1)
        self.assertEqual(self_match.iloc[0]["Target protein"], "WP_011029008.1")
        self.assertAlmostEqual(self_match.iloc[0]["Coverage on Query"], 1.0)
        self.assertAlmostEqual(self_match.iloc[0]["Coverage on Target"], 1.0)
        self.assertEqual(self_match.iloc[0]["Expect protein"], 0.0)


if __name__ == "__main__":
    unittest.main()
