import unittest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from predict_recombination_products import  BridgeRNAReaction
import find_brna

class TestBridgeRNAReaction(unittest.TestCase):
    
    def test_find_matches_single_strand(self):
        seq = "AAACCCCTGGGTTTAAACCCCTGGGTTT"
        site = "CCCCTGGG"
        matches = BridgeRNAReaction.find_matches(site, seq)
        # print(f"MATCHES\n\n\t{matches}")
        self.assertEqual(len(matches), 2)
        self.assertTrue(all(isinstance(m, tuple) for m in matches))
        self.assertIn(1, [m[0] for m in matches])  # strand +1
        self.assertTrue(all(m[1] < m[2] for m in matches))  # Start < End

    def test_find_matches_reverse_strand(self):
        seq = "AAACCCGGGTTT"
        site = "AAACCC"
        # Reverse complement of site is GGGTTT
        matches = BridgeRNAReaction.find_matches(site, seq)
        self.assertEqual(len(matches), 2)
        strands = [m[0] for m in matches]
        self.assertIn(1, strands)
        self.assertIn(-1, strands)

    def test_annotate_site_creates_features(self):
        seq = Seq("ATGCGTACTCGTAGCTAGCTAGCTAGCTAGCTAGCTGACT")
        record = SeqRecord(seq, id="test")
        result = BridgeRNAReaction.annotate_site(
            record=record,
            site_start=5,
            site_seq="CGTACGTA",
            le_seq="CGTA",
            re_seq="CGTA",
            strand=1,
            mode="Target",
            inplace=False
        )
        self.assertEqual(len(result.features), 3)
        labels = [f.qualifiers["label"] for f in result.features]
        self.assertIn("Target LE", labels)
        self.assertIn("Target RE", labels)
        self.assertIn("Target Core", labels)

    def test_split_feature(self):
        feat = SeqFeature(FeatureLocation(5, 15, strand=1), type="gene", qualifiers={"label": ["TestFeature"]})
        up, down = BridgeRNAReaction.split_feature(feat, 10)
        self.assertEqual(int(up.location.start), 5)
        self.assertEqual(int(up.location.end), 10)
        self.assertEqual(int(down.location.start), 10)
        self.assertEqual(int(down.location.end), 15)
        self.assertEqual(up.qualifiers["partial"], ["true"])
        self.assertEqual(down.qualifiers["partial"], ["true"])

    def test_find_partial_site_on_positive_strand(self):
        part = "CGTA"
        site = "CGTACGTA"
        contig_pos = 10
        strand = 1
        start, end = BridgeRNAReaction.find_partial_site(part, site, contig_pos, strand)
        self.assertEqual(start, 10)
        self.assertEqual(end, 14)

    def test_find_partial_site_on_negative_strand(self):
        part = "CGTA"
        site = "CGTACGTA"
        contig_pos = 10
        strand = -1
        start, end = BridgeRNAReaction.find_partial_site(part, site, contig_pos, strand)
        self.assertTrue(isinstance(start, int))
        self.assertTrue(isinstance(end, int))
        self.assertTrue(start < end)

    def test_verify_unique_site_success(self):
        mock_self = MagicMock()
        matches = [[(1, 100, 110)], [], []]
        try:
            BridgeRNAReaction._verify_unique_site(mock_self, matches, "target")
        except ValueError:
            self.fail("Unexpected ValueError raised for unique site.")

    def test_verify_unique_site_failure_multiple_plasmids(self):
        mock_self = MagicMock()
        matches = [[(1, 100, 110)], [(1, 200, 210)]]
        with self.assertRaises(ValueError):
            BridgeRNAReaction._verify_unique_site(mock_self, matches, "target")

    def test_verify_unique_site_failure_multiple_sites(self):
        mock_self = MagicMock()
        matches = [[(1, 100, 110), (1, 200, 210)]]
        with self.assertRaises(ValueError):
            BridgeRNAReaction._verify_unique_site(mock_self, matches, "target")

    @patch("find_brna.find_target_donor_sites")
    @patch("Bio.SeqIO.read")
    def test_process_brna_sets_sites(self, mock_read, mock_find):
        mock_seq = SeqRecord(Seq("ATGCGTACGTAGCTAGCTAGCTA"))
        mock_read.return_value = mock_seq
        # Donor core should be CT, donor RE should be TCG, LE is TA
        mock_find.return_value = (
            ("CGCTGTA", "CG", "CT", "GTA"), 
            ("TACTTCG", "TA", "CT", "TCG")
        )
        # Mock find_matches so unique-site verification passes without real sequence content
        with patch.object(BridgeRNAReaction, "find_matches", return_value=[(1, 0, 7)]):
            reaction = BridgeRNAReaction("mock_path", [])
        # Assert all sites are set
        self.assertEqual(reaction.target_site, "CGCTGTA")
        self.assertEqual(reaction.lt, "CG")
        self.assertEqual(reaction.target_core, "CT")
        self.assertEqual(reaction.rt, "GTA")
        self.assertEqual(reaction.donor_site, "TACTTCG")
        self.assertEqual(reaction.ld, "TA")
        self.assertEqual(reaction.donor_core, "CT")
        self.assertEqual(reaction.rd, "TCG")

    # Add insertion test
    @patch("find_brna.find_target_donor_sites")
    @patch("Bio.SeqIO.parse")
    @patch("Bio.SeqIO.read")
    def test_detect_insertion(self, mock_read, mock_parse, mock_find):
        # Mock bRNA target/donor sites
        mock_find.return_value = (
            ("CGCTGTA", "CG", "CT", "GTA"), 
            ("TACTTCG", "TA", "CT", "TCG")
        )
        mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
        # Provide two additional sequence records
        mock_parse.return_value = [
            SeqRecord(Seq("TTTTTT"), id="target_contig"),
            SeqRecord(Seq("GGGGGG"), id="donor_contig"),
        ]

        # We have N=3 seqs total (brna + 2 others)
        # target matches only in index 1; donor matches only in index 2
        target_phase = [[], [(1, 100, 107)], []]
        donor_phase = [[], [], [(1, 200, 207)]]
        counters = {"target": 0, "donor": 0}

        def side_effect(site, seq):
            if site == "CGCTGTA":
                idx = counters["target"]
                counters["target"] += 1
                return target_phase[idx] if idx < len(target_phase) else []
            else:
                idx = counters["donor"]
                counters["donor"] += 1
                return donor_phase[idx] if idx < len(donor_phase) else []

        with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
            reaction = BridgeRNAReaction("brna.gbk", ["target.gbk", "donor.gbk"]) 
            self.assertEqual(reaction.mode, "insertion")

    @patch("find_brna.find_target_donor_sites")
    @patch("Bio.SeqIO.parse")
    @patch("Bio.SeqIO.read")
    def test_detect_excision(self, mock_read, mock_parse, mock_find):
        mock_find.return_value = (
            ("CGCTGTA", "CG", "CT", "GTA"), 
            ("TACTTCG", "TA", "CT", "TCG")
        )
        mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
        mock_parse.return_value = [
            SeqRecord(Seq("TTTTTT"), id="same_contig"),
        ]

        # N=2 (brna + one other). Both target and donor in same index 1, strands equal
        target_phase = [[], [(1, 50, 57)]]
        donor_phase = [[], [(1, 150, 157)]]
        counters = {"target": 0, "donor": 0}

        def side_effect(site, seq):
            if site == "CGCTGTA":
                idx = counters["target"]
                counters["target"] += 1
                return target_phase[idx] if idx < len(target_phase) else []
            else:
                idx = counters["donor"]
                counters["donor"] += 1
                return donor_phase[idx] if idx < len(donor_phase) else []

        with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
            reaction = BridgeRNAReaction("brna.gbk", ["same.gbk"]) 
            self.assertEqual(reaction.mode, "excision")

    @patch("find_brna.find_target_donor_sites")
    @patch("Bio.SeqIO.parse")
    @patch("Bio.SeqIO.read")
    def test_detect_inversion(self, mock_read, mock_parse, mock_find):
        mock_find.return_value = (
            ("CGCTGTA", "CG", "CT", "GTA"), 
            ("TACTTCG", "TA", "CT", "TCG")
        )
        mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
        mock_parse.return_value = [
            SeqRecord(Seq("TTTTTT"), id="same_contig"),
        ]

        # N=2. Both at index 1, but opposite strands to trigger inversion
        target_phase = [[], [(1, 50, 57)]]
        donor_phase = [[], [(-1, 150, 157)]]
        counters = {"target": 0, "donor": 0}

        def side_effect(site, seq):
            if site == "CGCTGTA":
                idx = counters["target"]
                counters["target"] += 1
                return target_phase[idx] if idx < len(target_phase) else []
            else:
                idx = counters["donor"]
                counters["donor"] += 1
                return donor_phase[idx] if idx < len(donor_phase) else []

        with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
            reaction = BridgeRNAReaction("brna.gbk", ["same.gbk"]) 
            self.assertEqual(reaction.mode, "inversion")

    # # PRODUCT TESTS
    # @patch("find_brna.find_target_donor_sites")
    # @patch("Bio.SeqIO.parse")
    # @patch("Bio.SeqIO.read")
    # def test_insertion_product_same_strand(self, mock_read, mock_parse, mock_find):
    #     # Use RT length 3 so split site = end-3
    #     mock_find.return_value = (
    #         ("TTTTTT", "L", "CC", "RTT"),
    #         ("GGGGGG", "L", "CC", "RTT")
    #     )
    #     mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
    #     target_rec = SeqRecord(Seq("TTTTAAAACCCC"), id="target_contig")
    #     donor_rec = SeqRecord(Seq("GGGGCCCCAAAA"), id="donor_contig")
    #     mock_parse.return_value = [target_rec, donor_rec]

    #     # Matches: target in index 1 at 4..10; donor in index 2 at 3..9
    #     target_phase = [[], [(1, 4, 10)], []]
    #     donor_phase = [[], [], [(1, 3, 9)]]
    #     counters = {"target": 0, "donor": 0}

    #     def side_effect(site, seq):
    #         if site == "TTTTTT":
    #             idx = counters["target"]; counters["target"] += 1
    #             return target_phase[idx] if idx < len(target_phase) else []
    #         else:
    #             idx = counters["donor"]; counters["donor"] += 1
    #             return donor_phase[idx] if idx < len(donor_phase) else []

    #     with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
    #         reaction = BridgeRNAReaction("brna.gbk", ["target.gbk", "donor.gbk"]) 
    #         self.assertEqual(reaction.mode, "insertion")
    #         products = reaction.calculate_product()
    #         self.assertEqual(len(products), 1)
    #         # Expected: target[:7] + donor[6:] + donor[:6] + target[7:]
    #         expected = "TTTTAAA" + "CCCAAAA" + "GGGGCC" + "ACCCC"
    #         self.assertEqual(str(products[0].seq), expected)

    # @patch("find_brna.find_target_donor_sites")
    # @patch("Bio.SeqIO.parse")
    # @patch("Bio.SeqIO.read")
    # def test_insertion_product_opposite_strand(self, mock_read, mock_parse, mock_find):
    #     # RT length 3
    #     mock_find.return_value = (
    #         ("TTTTTT", "L", "CC", "RTT"),
    #         ("GGGGGG", "L", "CC", "RTT")
    #     )
    #     mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
    #     target_rec = SeqRecord(Seq("TTTTAAAACCCC"), id="target_contig")
    #     donor_rec = SeqRecord(Seq("GGGGCCCCAAAA"), id="donor_contig")
    #     mock_parse.return_value = [target_rec, donor_rec]

    #     # Target match on + strand; donor match on - strand
    #     target_phase = [[], [(1, 4, 10)], []]
    #     donor_phase = [[], [], [(-1, 3, 9)]]
    #     counters = {"target": 0, "donor": 0}

    #     def side_effect(site, seq):
    #         if site == "TTTTTT":
    #             idx = counters["target"]; counters["target"] += 1
    #             return target_phase[idx] if idx < len(target_phase) else []
    #         else:
    #             idx = counters["donor"]; counters["donor"] += 1
    #             return donor_phase[idx] if idx < len(donor_phase) else []

    #     with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
    #         reaction = BridgeRNAReaction("brna.gbk", ["target.gbk", "donor.gbk"]) 
    #         self.assertEqual(reaction.mode, "insertion")
    #         products = reaction.calculate_product()
    #         self.assertEqual(len(products), 1)
    #         # After RC donor, donor becomes TTTTGGGGCCCC; split still 6
    #         # Expected: target[:7] + donor[6:] + donor[:6] + target[7:]
    #         expected = "TTTTAAA" + "GGGCCCC" + "TTTTGG" + "ACCCC"
    #         self.assertEqual(str(products[0].seq), expected)

    # @patch("find_brna.find_target_donor_sites")
    # @patch("Bio.SeqIO.parse")
    # @patch("Bio.SeqIO.read")
    # def test_excision_product_sequences(self, mock_read, mock_parse, mock_find):
    #     # RT length 3
    #     mock_find.return_value = (
    #         ("TTTTTT", "L", "CC", "RTT"),
    #         ("GGGGGG", "L", "CC", "RTT")
    #     )
    #     mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
    #     target_rec = SeqRecord(Seq("TTTTAAAACCCC"), id="same_contig")
    #     mock_parse.return_value = [target_rec]

    #     target_phase = [[], [(1, 4, 10)]]
    #     donor_phase = [[], [(1, 1, 5)]]
    #     counters = {"target": 0, "donor": 0}

    #     def side_effect(site, seq):
    #         if site == "TTTTTT":
    #             idx = counters["target"]; counters["target"] += 1
    #             return target_phase[idx] if idx < len(target_phase) else []
    #         else:
    #             idx = counters["donor"]; counters["donor"] += 1
    #             return donor_phase[idx] if idx < len(donor_phase) else []

    #     with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
    #         reaction = BridgeRNAReaction("brna.gbk", ["same.gbk"]) 
    #         self.assertEqual(reaction.mode, "excision")
    #         products = reaction.calculate_product()
    #         self.assertEqual(len(products), 2)
    #         leftovers, excised = products
    #         # Cuts at 2 and 7
    #         self.assertEqual(str(leftovers.seq), "TT" + "ACCCC")
    #         self.assertEqual(str(excised.seq), "TTAAA")

    # @patch("find_brna.find_target_donor_sites")
    # @patch("Bio.SeqIO.parse")
    # @patch("Bio.SeqIO.read")
    # def test_inversion_product_sequence(self, mock_read, mock_parse, mock_find):
    #     # RT length 3
    #     mock_find.return_value = (
    #         ("TTTTTT", "L", "CC", "RTT"),
    #         ("GGGGGG", "L", "CC", "RTT")
    #     )
    #     mock_read.return_value = SeqRecord(Seq("AAAAAA"), id="brna")
    #     target_rec = SeqRecord(Seq("TTTTAAAACCCC"), id="same_contig")
    #     mock_parse.return_value = [target_rec]

    #     target_phase = [[], [(1, 4, 10)]]
    #     donor_phase = [[], [(-1, 1, 5)]]
    #     counters = {"target": 0, "donor": 0}

    #     def side_effect(site, seq):
    #         if site == "TTTTTT":
    #             idx = counters["target"]; counters["target"] += 1
    #             return target_phase[idx] if idx < len(target_phase) else []
    #         else:
    #             idx = counters["donor"]; counters["donor"] += 1
    #             return donor_phase[idx] if idx < len(donor_phase) else []

    #     with patch.object(BridgeRNAReaction, "find_matches", side_effect=side_effect):
    #         reaction = BridgeRNAReaction("brna.gbk", ["same.gbk"]) 
    #         self.assertEqual(reaction.mode, "inversion")
    #         products = reaction.calculate_product()
    #         self.assertEqual(len(products), 1)
    #         expected = "TT" + "TTTAA" + "ACCCC"
    #         self.assertEqual(str(products[0].seq), expected)


if __name__ == "__main__":
    unittest.main()
