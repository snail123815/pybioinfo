import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    BeforePosition,
    AfterPosition,
    ExactPosition,
)

from pyBioinfo_modules.bio_sequences.bio_features import (
    truncate_feat_translation,
    find_truncated_features,
    slice_sequence,
    reverse_complement_seqrecord_with_features,
    reverse_complement_position,
    reverse_complement_location,
)


class TestBioFeatures(unittest.TestCase):
    # def test_truncate_feat_translation_left(self):
    #     feat = SeqFeature(FeatureLocation(0, 9, strand=1), type="CDS", qualifiers={"translation":["MKTWYIK"]})
    #     truncated = truncate_feat_translation(feat, "left")
    #     self.assertIn("translation", truncated.qualifiers)
    #     self.assertEqual(truncated.qualifiers["translation"][0], "TWYIK")

    # def test_find_truncated_features_basic(self):
    #     seq = SeqRecord(Seq("ATGCATGCATGC"), id="test")
    #     feat = SeqFeature(FeatureLocation(0, 12, strand=1), type="gene")
    #     seq.features.append(feat)
    #     found = find_truncated_features(seq, (0, 6))
    #     self.assertEqual(len(found), 1)
    #     self.assertLess(found[0].location.end, 12)

    # def test_slice_sequence_bounds(self):
    #     seq = SeqRecord(Seq("ATGCATGCATGC"), id="test")
    #     feat = SeqFeature(FeatureLocation(0, 12, strand=1), type="CDS")
    #     seq.features.append(feat)
    #     sliced = slice_sequence(seq, FeatureLocation(3, 9, strand=1), with_features=True)
    #     self.assertEqual(str(sliced.seq), "CATGCA")
    #     self.assertTrue(len(sliced.features) > 0)

    # def test_reverse_complement_seqrecord_with_features(self):
    #     seq = SeqRecord(Seq("ATGCATGC"), id="test")
    #     feat = SeqFeature(FeatureLocation(0, 8, strand=1), type="misc")
    #     seq.features.append(feat)
    #     rc = reverse_complement_seqrecord_with_features(seq)
    #     self.assertEqual(str(rc.seq), "GCATGCAT")
    #     self.assertEqual(rc.features[0].location.strand, -1)

    def test_reverse_complement_position(self):
        pos = reverse_complement_position(12, 5)
        self.assertEqual(pos, 7)
        pos = reverse_complement_position(12, 0)
        self.assertEqual(pos, 12)
        pos = reverse_complement_position(12, BeforePosition(5))
        self.assertIsInstance(pos, AfterPosition)
        self.assertEqual(pos, AfterPosition(7))
        pos = reverse_complement_position(12, AfterPosition(5))
        self.assertIsInstance(pos, BeforePosition)
        self.assertEqual(pos, BeforePosition(7))
        pos = reverse_complement_position(12, ExactPosition(5))
        self.assertIsInstance(pos, ExactPosition)
        self.assertEqual(pos, ExactPosition(7))

    def test_reverse_complement_location(self):
        loc = reverse_complement_location(12, (3, 5))
        self.assertEqual(loc.start, 7)
        self.assertEqual(loc.end, 9)
        self.assertIs(loc.strand, None)
        loc = reverse_complement_location(12, FeatureLocation(3, 5, strand=1))
        self.assertEqual(loc.start, 7)
        self.assertEqual(loc.end, 9)
        self.assertEqual(loc.strand, -1)
        loc = reverse_complement_location(
            12, FeatureLocation(BeforePosition(3), AfterPosition(5), strand=1)
        )
        self.assertIsInstance(loc.start, BeforePosition)
        self.assertEqual(loc.start, BeforePosition(7))
        self.assertIsInstance(loc.end, AfterPosition)
        self.assertEqual(loc.end, AfterPosition(9))
        self.assertEqual(loc.strand, -1)
        with self.assertRaises(ValueError):
            reverse_complement_location(12, (3, 15, 1, 4))


if __name__ == "__main__":
    unittest.main()
