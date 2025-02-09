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
    slice_sequence_keep_truncated_features,
    seqFeature_to_tuple,
    add_seq_to_SeqRecord_as_feature,
    reverse_complement_SeqRecord_with_features,
    reverse_complement_position,
    reverse_complement_location,
)


class TestBioFeatures(unittest.TestCase):

    def test_seqFeature_to_tuple(self):
        feat1 = SeqFeature(
            FeatureLocation(3, 5, strand=-1),
            type="gene",
            qualifiers={"gene": ["test_gene"]},
        )
        feat2 = SeqFeature(
            FeatureLocation(3, 5, strand=1),
            type="gene",
            qualifiers={"gene": ["test_gene"]},
        )
        feat3 = SeqFeature(
            FeatureLocation(3, 5, 1),
            type="gene",
            qualifiers={"gene": ["test_gene"]},
        )
        self.assertNotEqual(
            seqFeature_to_tuple(feat1), seqFeature_to_tuple(feat2)
        )
        self.assertEqual(seqFeature_to_tuple(feat2), seqFeature_to_tuple(feat3))

    def test_add_seq_to_SeqRecord_as_feature(self):
        seq = SeqRecord(
            Seq("ATGCaTGCgaaatc"),
            id="test",
            features=[
                SeqFeature(FeatureLocation(0, 8), type="source"),
            ],
        )
        new_seq = add_seq_to_SeqRecord_as_feature(
            seq,
            "aTTTc",
            "gene",
            {"id": ["gene_id"]},
        )
        new_seq = add_seq_to_SeqRecord_as_feature(
            new_seq,
            "CCCGGcA",
            "primer_bind",
            {"id": ["pb1"], "overhang": ["CcCG"]},
        )
        new_seq = add_seq_to_SeqRecord_as_feature(
            new_seq,
            "CCCGttcgc",
            "primer_bind",
            {"id": ["rev"], "overhang": ["4"]},
        )
        # original seq should not be modified
        self.assertEqual(str(seq.seq), "ATGCaTGCgaaatc")
        self.assertEqual(str(new_seq.seq), "ATGCaTGCgaaatc")
        self.assertEqual(seq.id, "test")
        self.assertEqual(new_seq.id, "test")
        self.assertEqual(len(seq.features), 1)

        # new_seq should have old feature
        self.assertEqual(new_seq.features[0].type, "source")
        # new_seq should have new features
        self.assertEqual(len(new_seq.features), 4)
        gene_feat = new_seq.features[1]
        self.assertEqual(gene_feat.type, "gene")
        self.assertEqual(gene_feat.qualifiers["id"], ["gene_id"])
        self.assertEqual(gene_feat.location.start, 8)
        self.assertEqual(gene_feat.location.end, 13)
        self.assertEqual(gene_feat.location.strand, -1)
        pb_feat1 = new_seq.features[2]
        self.assertEqual(pb_feat1.type, "primer_bind")
        self.assertEqual(pb_feat1.location.start, 2)
        self.assertEqual(pb_feat1.location.end, 5)
        self.assertEqual(pb_feat1.location.strand, 1)
        self.assertEqual(pb_feat1.qualifiers["id"], ["pb1"])
        self.assertEqual(pb_feat1.qualifiers["overhang"], ["CcCG"])
        pb_feat2 = new_seq.features[3]
        self.assertEqual(pb_feat2.location.start, 6)
        self.assertEqual(pb_feat2.location.end, 11)
        self.assertEqual(pb_feat2.location.strand, -1)
        self.assertEqual(pb_feat2.qualifiers["id"], ["rev"])
        self.assertEqual(pb_feat2.qualifiers["overhang"], ["CCCG"])

    def test_reverse_complement_seqrecord_with_features(self):
        seq = SeqRecord(
            Seq("ATGCATGC"),
            id="test",
            features=[
                SeqFeature(FeatureLocation(0, 8), type="source"),
                SeqFeature(
                    FeatureLocation(2, AfterPosition(4), strand=-1),
                    type="gene",
                    qualifiers={"gene": ["test_gene"]},
                ),
            ],
        )
        rc = reverse_complement_SeqRecord_with_features(seq)
        self.assertEqual(str(rc.seq), "GCATGCAT")
        feat1, feat2 = rc.features

        self.assertEqual(feat1.location.start, 0)
        self.assertEqual(feat1.location.end, 8)
        self.assertIs(rc.features[0].location.strand, None)
        self.assertEqual(feat1.type, "source")

        self.assertIsInstance(feat2.location.start, BeforePosition)
        self.assertEqual(feat2.location.start, BeforePosition(4))
        self.assertIsInstance(feat2.location.end, ExactPosition)
        self.assertEqual(feat2.location.end, ExactPosition(6))
        self.assertEqual(feat2.location.strand, 1)
        self.assertEqual(feat2.type, "gene")
        self.assertEqual(feat2.qualifiers["gene"], ["test_gene"])

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
