import argparse
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
# from find_brna import find_target_donor_sites
import find_brna
import genbank_utils
import copy
import os

class bRNAReaction():
    def __init__(self, brna_path, other_sequence_paths, out_dir):
        self.brna_path = brna_path
        self.other_sequence_paths = other_sequence_paths
        self.out_dir = out_dir
        self._process_brna()
        # self._read_dna_files()
        self._identify_sequences()

    def _process_brna(self):
        # Determines target and donor sites from bRNA in bRNA plasmid
        self.brna_plasmid_record = SeqIO.read(self.brna_path, "genbank")
        target, donor = find_brna.find_target_donor_sites(self.brna_plasmid_record)
        self.target_site, self.lt, self.target_core, self.rt = target
        self.donor_site, self.ld, self.donor_core, self.rd = donor

        print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")

    # TODO: more descriptive error messages
    def _verify_unique_sites(self, target_matches, donor_matches):
        unique_donor_plasmid = len([i for i in donor_matches if i]) == 1
        unique_donor_site = sum([len(i) for i in donor_matches]) == 1
        unique_target_plasmid = len([i for i in target_matches if i]) == 1
        unique_target_site = sum([len(i) for i in target_matches]) == 1
        if not unique_donor_plasmid:
            raise ValueError(f"Expected exactly one donor plasmid match, but found {len(donor_matches)}.")

        if not unique_donor_site:
            total_donor_sites = sum(len(i) for i in donor_matches)
            raise ValueError(f"Expected exactly one donor site match, but found {total_donor_sites} total.")

        if not unique_target_plasmid:
            raise ValueError(f"Expected exactly one target plasmid match, but found {len(target_matches)}.")

        if not unique_target_site:
            total_target_sites = sum(len(i) for i in target_matches)
            raise ValueError(f"Expected exactly one target site match, but found {total_target_sites} total.")
        
    def _identify_sequences(self):
        # Read other sequence records
        sequence_paths = [self.brna_path] + self.other_sequence_paths
        other_sequence_records = [
            SeqIO.read(file, "genbank") for file in self.other_sequence_paths
            ]
        sequence_records = [self.brna_plasmid_record] + other_sequence_records
        seqs = [str(record.seq) for record in sequence_records]
        # Search all plasmid records for target and donor sites
        target_matches = [
            self.find_matches(self.target_site, i) for i in seqs]
        donor_matches = [
            self.find_matches(self.donor_site, i) for i in seqs]
        # Make sure only one of each site is found
        self._verify_unique_sites(target_matches, donor_matches)
        # Determine rearrangement mode
        target_index = np.where([bool(i) for i in target_matches])[0].item()
        donor_index = np.where([bool(i) for i in donor_matches])[0].item()
        self.target_strand, self.target_start, self.target_end = target_matches[target_index][0]
        self.donor_strand, self.donor_start, self.donor_end = donor_matches[donor_index][0]
        # Store target contig, flip to + strand if necessary
        print(f"Found target plasmid: {sequence_paths[target_index]}\n\tstart: {self.target_start}\n\tend: {self.target_end}")
        target_contig = sequence_records[target_index]
        # if self.target_strand == -1:
        #     target_contig = genbank_utils.reverse_complement_record(target_contig)
        # Add annotations to target record
        self.annotate_site(target_contig, self.target_start, self.target_site, 
                           self.lt, self.rt, self.target_strand, mode="Target", inplace=True)
        self.target_contig = target_contig
        # Store split sites
        if self.target_strand == 1:
            self.target_split_site = self.target_end - len(self.rt)
        else:
            self.target_split_site = self.target_start + len(self.rt)
        if self.donor_strand == 1:
            self.donor_split_site = self.donor_end - len(self.rt)
        else:
            self.donor_split_site = self.donor_start + len(self.rt)
        print(f"Split sites: DONOR: {self.donor_split_site} TARGET: {self.target_split_site}")
        # self.donor_split_site = self.donor_start + 2 + len(self.lt)
        # self.donor_split_site = self.donor_start + 2 + len(self.ld)
        # Mode-specific preprocessing
        if target_index == donor_index:
            if self.target_strand == self.donor_strand:
                self.mode = "excision"
            else:
                self.mode = "inversion"
                # re-calculate split site for whichever site is on the (-) strand
                # L = len(self.target_contig.seq)
                # if self.target_strand == -1:
                #     self.target_split_site = L - self.target_split_site
                # if self.donor_strand == -1:
                #     self.donor_split_site = L - self.donor_split_site
        else:
            self.mode = "insertion"
            print(f"Found donor plasmid: {sequence_paths[donor_index]}")
            # Store separate donor record for insertions
            donor_contig = sequence_records[donor_index]
            # Add annotations to donor record
            # if self.donor_strand == -1:
            #     donor_contig = genbank_utils.reverse_complement_record(donor_contig)
            self.annotate_site(
                donor_contig, self.donor_start, self.donor_site, self.ld, 
                self.rd, self.donor_strand, mode="Donor", inplace=True)
            self.donor_contig = donor_contig

    def calculate_product(self):
        print(f"Calculating product in {self.mode} mode...")
        if self.mode == "insertion":
            # donor_split_site = self.donor_start + 2 + len(self.ld)
            # target_split_site = self.target_start + 2 + len(self.lt)
            result = self.calculate_insertion()
            return result
        elif self.mode == "excision":
            return self.calculate_excision()
        elif self.mode == "inversion":
            return self.calculate_inversion()
        else:
            raise ValueError("Invalid mode, no product calculated!")
    
    def calculate_insertion(self):
        """
        Calculates the product of an insertion reaction.
        """
        insertion_size = len(self.donor_contig.seq)
        insertion_record = genbank_utils.orient_record_around_pos(self.donor_contig, self.donor_split_site)
        # If the donor is on a different strand in its native plasmid than the 
        # target strand in the target plasmid, reverse complement the donor 
        # plasmid so we can add its sequence directly
        if self.donor_strand != self.target_strand:
            insertion_record = genbank_utils.reverse_complement_record(insertion_record)
        # Figure out the new sequence
        new_seq = (
            self.target_contig.seq[:self.target_split_site]
            + insertion_record.seq
            + self.target_contig.seq[self.target_split_site:]
            )
        # Transfer features
        features = []
        for feature in insertion_record.features:
            shifted_feature = genbank_utils.shift_feature(feature, self.target_split_site)
            features.append(shifted_feature)
        
        for feat in self.target_contig.features:
            feature = copy.deepcopy(feat)
            ends_before_insertion = feature.location.end <= self.target_split_site
            starts_before_insertion = feature.location.start < self.target_split_site
            before_insertion = starts_before_insertion & ends_before_insertion
            # If feature overlaps insertion, just drag it across the insertion. 
            if starts_before_insertion and not ends_before_insertion:
                feature.location = FeatureLocation(
                    start=int(feature.location.start),
                    end=int(feature.location.end) + insertion_size,
                    strand=feature.location.strand
                )
            # Shift the feature position up if it's to the right of the insertion
            elif not before_insertion:
                feature = genbank_utils.shift_feature(feature, insertion_size)
            features.append(feature)
        # Create the final record
        result_id = "recombined product" # I'll fill this in later.
        result_name = "recombined_product" # I'll fill this in later.
        result_description = f"Predicted recombination product using target site = {self.target_site} and donor site = {self.donor_site}"
        recombinant_record = SeqRecord(new_seq, id=result_id, name=result_name, 
                                    description=result_description)
        # Mark the record as circular
        print("Assuming circular topology for insertion")
        recombinant_record.annotations["topology"] = "circular"
        recombinant_record.annotations["molecule_type"] = "DNA"
        # Apply features to genbank record
        recombinant_record.features = features

        # # DEBUGGING
        # self.insert = insertion_record
        return [recombinant_record]

    def calculate_excision(self):
        cut_points = [self.target_split_site, self.donor_split_site]
        upstream_cut = min(cut_points)
        downstream_cut = max(cut_points)
        excision_length = downstream_cut - upstream_cut


        leftovers = copy.deepcopy(self.target_contig)
        excised = copy.deepcopy(self.target_contig)
        # Determine excision sequences
        leftovers_seq = (self.target_contig.seq[:upstream_cut] 
                         + self.target_contig.seq[downstream_cut:])
        excised_seq = self.target_contig.seq[upstream_cut:downstream_cut]

        leftovers_features, excised_features = [], []
        for feature in self.target_contig.features:
            starts_before_le = feature.location.start < upstream_cut
            ends_before_le = feature.location.end <= upstream_cut

            starts_before_re = feature.location.start < downstream_cut
            ends_before_re = feature.location.end <= downstream_cut

            if starts_before_le & ends_before_le:
                leftovers_features.append(feature)
            elif (not starts_before_re) & (not ends_before_re):
                # Shift features upstream by length of excision
                shifted_feat = genbank_utils.shift_feature(feature, -excision_length)
                leftovers_features.append(shifted_feat)
            elif starts_before_le & (not ends_before_le):
                # Break upstream overlap
                # print("Breaking upstream: ", feature)
                pass
            elif starts_before_re & (not ends_before_re):
                # Break downstream overlap
                # print("Breaking downstream: ", feature)
                pass
            else:
                # These are the features that end up in the excised product
                shifted_feat = genbank_utils.shift_feature(feature, -upstream_cut)
                excised_features.append(shifted_feat)
        leftovers.features = leftovers_features
        excised.features = excised_features
        # Add sequences
        leftovers.seq = leftovers_seq
        excised.seq = excised_seq
        return [leftovers, excised]

    def calculate_inversion(self):
        cut_points = [self.target_split_site, self.donor_split_site]
        # print(cut_points, [self.target_start, self.donor_start])
        upstream_cut = min(cut_points)
        downstream_cut = max(cut_points)
        result = copy.deepcopy(self.target_contig)
        original_features = self.target_contig.features
        features = []
        # Calculate new sequence
        result_seq = (
            self.target_contig.seq[:upstream_cut]
            + reverse_complement(self.target_contig.seq[upstream_cut:downstream_cut])
            + self.target_contig.seq[downstream_cut:]
        )
        print("Junctions before inversion")
        print(str(self.target_contig.seq[upstream_cut-10:upstream_cut+10]))
        print(str(self.target_contig.seq[downstream_cut-10:downstream_cut+10]))
        print("Junctions after inversion")
        print(result_seq[upstream_cut-10:upstream_cut+10])
        print(result_seq[downstream_cut-10:downstream_cut+10])

        # print(result_seq == self.target_contig.seq)
        # print(reverse_complement(self.target_contig.seq[upstream_cut:downstream_cut]))
        result.seq = result_seq
        # Store indices to use in genbank_utils.reverse_complement_feature
        invert_indices = np.concatenate([
            np.arange(upstream_cut),
            np.arange(upstream_cut, downstream_cut)[::-1],
            np.arange(downstream_cut, len(self.target_contig.seq)),
        ])
        # print(invert_indices[upstream_cut-10:upstream_cut+10])
        # print(invert_indices[downstream_cut-10:downstream_cut+10])
        # for feature in original_features:
        for feature in original_features:
            starts_before_le = feature.location.start < upstream_cut
            ends_before_le = feature.location.end <= upstream_cut

            starts_before_re = feature.location.start < downstream_cut
            ends_before_re = feature.location.end <= downstream_cut

            if (starts_before_le & ends_before_le) or ((not starts_before_re) & (not ends_before_re)):
                # print(feature)
                features.append(feature)
                # print("Appended! ", len(features))
            elif starts_before_le & (not ends_before_le):
                # Break upstream overlap
                # print("Breaking upstream: ", feature)
                pass
            elif starts_before_re & (not ends_before_re):
                # Break downstream overlap
                # print("Breaking downstream: ", feature)
                pass
            else:
                # These are the features that end up inverted
                inverted = genbank_utils.reverse_complement_feature(feature, invert_indices)
                features.append(inverted)
            result.features = features
            # print(len(features), len(self.target_contig.features))
        return [result]

    @classmethod    
    def find_matches(cls, site:str, seq:str):
        """
        Find matches for site in seq. 
        """
        # Find matches on positive and negative strands
        matches_pos_strand = [match for match in re.finditer(site, seq)]
        rc_contig = reverse_complement(seq)
        matches_neg_strand = [match for match in re.finditer(site, rc_contig)]
        matches = matches_pos_strand + matches_neg_strand
        match_strands = np.concatenate([
            np.ones_like(matches_pos_strand), np.ones_like(matches_neg_strand) * -1
            ])
        match_locations = []
        L = len(seq)
        for idx, match in enumerate(matches):
            strand = match_strands[idx]
            if strand == 1:
                site_start, site_end = match.start(), match.end()
            else:
                site_start, site_end = sorted([L - match.start(), L - match.end()])
            match_locations.append((strand, site_start, site_end))
        return match_locations

    @classmethod
    def annotate_site(
            cls, record, site_start, site_seq, le_seq, re_seq, strand=1, 
            mode="Target", inplace=False):
        """Adds left end, right end, and core annotations of either target or 
        donor site to a sequence record.

        Args:
            record (Bio.SeqRecord.SeqRecord): Record to annotate
            site_start (int): Start position of LE of site to annotate
            site_seq (str): Sequence of site to annotate
            le_seq (str): LE sequence
            re_seq (str): RE sequence
            strand (int): strand site is found on. Defaults to +1.
            mode (str, optional): Whether to annotate as Target or Donor site. 
                Defaults to "Target".
            inplace (bool, optional): Modify record inplace. Defaults to False.

        Returns:
            Bio.SeqRecord.SeqRecord: Input record with either Target or Donor
            site now annotated.
        """        
        allowed_modes = ["Target", "Donor"]
        assert mode in allowed_modes, f"mode must be one of: {allowed_modes}"

        # print(site_start, site_seq, le_seq, core, re_seq)
        colors = {
            "Target": "#00aeff",
            "Donor": "#ff8400",
            "Core": "#24951a",
        }
        # Extract indices of LT, RT relative to their respective contigs
        le_start, le_end = cls.find_partial_site(le_seq, site_seq, site_start, strand)
        re_start, re_end = cls.find_partial_site(re_seq, site_seq, site_start, strand)
        if strand == 1:
            core_start, core_end = sorted([le_end, re_start])
        else: 
            core_start, core_end = sorted([re_end, le_start])
        print(mode, " LE: ", le_start, le_end, " RE: ", re_start, re_end, " Strand:", strand, "core:", core_end - core_start)
        # Annotate the partial sites
        # total_feat = SeqFeature(
        #     location=FeatureLocation(le_start, re_end, strand), 
        #     type="misc_feature", 
        #     qualifiers = {
        #         "label": f"{mode}",
        #         "ApEinfo_fwdcolor": colors[mode],
        #         "ApEinfo_revcolor": colors[mode],
        #         "note": f"color: {colors[mode]}"
        #     })
        direction = "RIGHT" if strand == 1 else "LEFT"
        le_feat = SeqFeature(
            location=FeatureLocation(le_start, le_end, strand), 
            type="gene", 
            qualifiers = {
                "label": f"{mode} LE",
                "ApEinfo_fwdcolor": colors[mode],
                "ApEinfo_revcolor": colors[mode],
                "note": f"color: {colors[mode]}; direction: {direction}"
            })
        re_feat = SeqFeature(
            location=FeatureLocation(re_start, re_end, strand), 
            type="gene", 
            qualifiers = {
                "label": f"{mode} RE",
                "ApEinfo_fwdcolor": colors[mode],
                "ApEinfo_revcolor": colors[mode],
                "note": f"color: {colors[mode]}; direction: {direction}"
            })
        core_feat = SeqFeature(
            location=FeatureLocation(core_start, core_end, strand), 
            type="gene", 
            qualifiers = {
                "label": f"{mode} Core",
                "ApEinfo_fwdcolor": colors["Core"],
                "ApEinfo_revcolor": colors["Core"],
                "note": f"color: {colors['Core']}; direction: {direction}"
            })
        # total_feat = SeqFeature(
        #     location=CompoundLocation([
        #         le_feat.location,
        #         core_feat.location,
        #         re_feat.location
        #     ]),
        #     type="misc_feature",
        #     qualifiers={
        #         "label": f"{mode}",
        #         "ApEinfo_fwdcolor": colors[mode],
        #         "ApEinfo_revcolor": colors[mode],
        #         "note": f"color: {colors[mode]}; direction: {direction}"
        #     },
        #     # sub_features=[le_feat, core_feat, re_feat]
        # )
        # new_feats = [le_feat, core_feat, re_feat, total_feat]
        new_feats = [le_feat, core_feat, re_feat,]
        if inplace:
            record.features += new_feats
        else:
            new_record = copy.deepcopy(record)
            new_record.features += new_feats
            return new_record

    @classmethod
    def find_partial_site(cls, part, site, contig_pos, strand=+1):
        """Finds locations of parts of a given site in relation to the location of
        the overall site in its native contig

        Args:
            part (str): Part of site. e.g. LT, RT, LD, RD
            site (str): Full sequence of target or donor site
            contig_pos (str): Start index of site in contig
            strand (int): Strand in which the site is found on the contig. 
                Must be either 1 or -1. Defaults to 1.
        Returns:
            part_start, part_end: positions wrt the positive strand
        """    
        if strand == -1:
            part = reverse_complement(part)
            site = reverse_complement(site)
        part_start = site.index(part)
        part_end = part_start + len(part)
        print(part_start, part_end, contig_pos, part, site)
        # part_start, part_end = find_matches(part, site, complement_search=False, 
        #                                     return_strand=False)
        part_start, part_end = contig_pos + part_start, contig_pos + part_end
        part_start, part_end = sorted([part_start, part_end])
        return part_start, part_end

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = ("Script to predict IS110 recombination product given target "
        "and donor contigs, left and right ends of target and donor sites, and "
        "cores of each. ")
    )
    # TODO: make compatible with reading bRNA plasmid with regex
    # NOTE: handshake guide of last 2 bases of 3' of RTG and RDG is NOT targeting
    # NOTE 2: CA core at 3' of LDG?, CT core at 3' of LTG
    # NOTE: invariant AA/TT in the RDG
    # NOTE: the cores will be found as the 2 cores in the target guide
    # NOTE: it's fair to expect *FOR NOW* that cores are the same for a given recombination
    parser.add_argument('-s', '--other_sequence_records', nargs='+', 
                        help='Specify one or more additional sequence_records/genomes to ' \
                        'search for target and donor sites.')

    parser.add_argument(
        "-b",
        "--brna_path", 
        help="Path to bRNA plasmid",
        )

    parser.add_argument("-o", "--out_dir", required=True, help="Output directory",)

    # # # Test case
    # args = parser.parse_args(
    #     [
    #         # "--target_contig", "/groups/cress/projects/is110/sequences/SLMQRT/Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk",
    #         # # "--donor_contig", "/groups/cress/projects/is110/sequences/SLMQRT/monomeric_BAC_for_plate_scrape_13kb.gbk",
    #         # "--donor_contig", "/groups/cress/projects/efactors/jack/insertion_mapping/misc/W9G3HH_30kbp_BAC_gSES0133.gbk",
    #         "--brna_path", "../data/sequence_records/P7YMY9_pSES0099.gbk",
    #         "-o", "../misc/insertion_test", 
    #         "--other_sequence_records", 
    # #        "/groups/cress/projects/is110/sequences/SLMQRT/Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk", 
    # #        "/groups/cress/projects/efactors/jack/insertion_mapping/misc/W9G3HH_30kbp_BAC_gSES0133.gbk",
    #         "../data/genomes/Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk", 
    #         "../data/sequence_records/W9G3HH_30kbp_BAC_gSES0133.gbk",
    #         ])
    args = parser.parse_args()
    def write_genbank(record, path):
        with open(path, "w") as output_handle:
            SeqIO.write(record, output_handle, "genbank")

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    Experiment = bRNAReaction(
        args.brna_path, args.other_sequence_records, args.out_dir
    )
    results = Experiment.calculate_product()
    print(f"Saving results to {args.out_dir}")
    if Experiment.mode == "insertion":
        write_genbank(Experiment.target_contig, "../misc/insertion_test/genome.gbk")
        write_genbank(Experiment.donor_contig, "../misc/insertion_test/insert.gbk")
        # write_genbank(Experiment.insert, "../misc/insertion_test/insert_oriented.gbk")
        with open(os.path.join(args.out_dir, "insertion.gbk"), "w") as output_handle:
            SeqIO.write(results[0], output_handle, "genbank")
    elif Experiment.mode == "excision":
        with open(os.path.join(args.out_dir, "leftovers.gbk"), "w") as output_handle:
            SeqIO.write(results[0], output_handle, "genbank")
        with open(os.path.join(args.out_dir, "excision.gbk"), "w") as output_handle:
            SeqIO.write(results[1], output_handle, "genbank")
    elif Experiment.mode == "inversion":
        with open(os.path.join(args.out_dir, "inversion.gbk"), "w") as output_handle:
            SeqIO.write(results[0], output_handle, "genbank")
