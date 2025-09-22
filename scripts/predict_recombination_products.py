import argparse
import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import find_brna
import copy
import os

class BridgeRNAReaction():
    def __init__(self, brna_path, other_sequence_paths, 
                 annotate_split_features=True, verify_unique_target=True):
        self.brna_path = brna_path
        self.other_sequence_paths = other_sequence_paths
        self.annotate_split_features = annotate_split_features
        self.verify_unique_target = verify_unique_target
        self._process_brna()
        self._identify_sequences()

    def _process_brna(self):
        # Determines target and donor sites from bRNA in bRNA plasmid
        # NOTE: currently assumes there's only one sequence record in brna_plasmid
        self.brna_plasmid_record = SeqIO.read(self.brna_path, "genbank")
        target, donor = find_brna.find_target_donor_sites(self.brna_plasmid_record)
        self.target_site, self.lt, self.target_core, self.rt = target
        self.donor_site, self.ld, self.donor_core, self.rd = donor
        print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")

    def _verify_unique_site(self, matches, mode):
        unique_plasmid = len([i for i in matches if i]) == 1
        unique_site = sum([len(i) for i in matches]) == 1
        if not unique_plasmid:
            raise ValueError(f"Expected exactly one {mode} plasmid match, but found {len(matches)}.")

        if not unique_site:
            total_sites = sum(len(i) for i in matches)
            raise ValueError(f"Expected exactly one {mode} site match, but found {total_sites} total.")
        
    def _identify_sequences(self):
        # Read other sequence records
        sequence_paths = [self.brna_path]
        # Store all sequence records in additional sequence files
        other_sequence_records = []
        for filepath in self.other_sequence_paths:
            records = SeqIO.parse(filepath, "genbank")
            for record in records:
                    other_sequence_records.append(record)
                    sequence_paths.append(filepath)
        # other_sequence_records = [
        #     SeqIO.read(file, "genbank") for file in self.other_sequence_paths
        #     ]
        sequence_records = [self.brna_plasmid_record] + other_sequence_records
        seqs = [str(record.seq) for record in sequence_records]
        # Search all sequence records for target and donor sites
        target_matches = [
            self.find_matches(self.target_site, i) for i in seqs]
        donor_matches = [
            self.find_matches(self.donor_site, i) for i in seqs]
        if self.verify_unique_target:
            self._verify_unique_site(target_matches, "target")
        self._verify_unique_site(donor_matches, "donor")
        target_index = np.where([bool(i) for i in target_matches])[0].item()
        donor_index = np.where([bool(i) for i in donor_matches])[0].item()
        self.target_strand, self.target_start, self.target_end = target_matches[target_index][0]
        self.donor_strand, self.donor_start, self.donor_end = donor_matches[donor_index][0]
        # Store target contig
        self.target_file = os.path.splitext(os.path.basename(sequence_paths[target_index]))[0]
        print(f"Found target sequence: {self.target_file}\n\tstart: {self.target_start}\n\tend: {self.target_end}\n\tstrand: {self.target_strand}")
        self._find_split_sites()
        self.target_contig = sequence_records[target_index]
        self.annotate_site(
            record=self.target_contig, 
            site_start=self.target_start, 
            site_seq=self.target_site, 
            le_seq=self.lt, 
            re_seq=self.rt, 
            strand=self.target_strand, 
            mode="Target", 
            inplace=True)
        # Determine rearrangement mode
        if target_index == donor_index:
            self.annotate_site(
                record=self.target_contig, 
                site_start=self.donor_start, 
                site_seq=self.donor_site, 
                le_seq=self.ld, 
                re_seq=self.rd, 
                strand=self.donor_strand, 
                mode="Donor", 
                inplace=True
                )
            if self.target_strand == self.donor_strand:
                self.mode = "excision"
            else:
                self.mode = "inversion"
        else:
            self.mode = "insertion"
            # Store separate donor record for insertions
            self.donor_file = os.path.splitext(os.path.basename(sequence_paths[donor_index]))[0]
            print(f"Found donor sequence: {self.donor_file}\n\tstrand: {self.donor_strand}")
            self.donor_contig = sequence_records[donor_index]
            self.annotate_site(
                record=self.donor_contig, 
                site_start=self.donor_start, 
                site_seq=self.donor_site, 
                le_seq=self.ld, 
                re_seq=self.rd, 
                strand=self.donor_strand, 
                mode="Donor", 
                inplace=True
                )
            # Make sure donor and target are in the same orientation
            if self.donor_strand != self.target_strand:
                print("Reverse complementing donor contig...")
                self.donor_contig = reverse_complement_record(self.donor_contig)
                donor_length = len(self.donor_contig.seq)
                self.donor_split_site = donor_length - self.donor_split_site
                self.donor_start = donor_length - self.donor_start
                self.donor_end = donor_length - self.donor_end
                self.donor_strand = -self.donor_strand

    def _find_split_sites(self):
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
    
    def calculate_product(self):
        print(f"Calculating product in {self.mode} mode...")
        if self.mode == "insertion":
            result = self.calculate_insertion()
            return result
        elif self.mode == "excision":
            return self.calculate_excision()
        elif self.mode == "inversion":
            return self.calculate_inversion()
        else:
            raise ValueError("Invalid mode, no product calculated!")
    
    def calculate_insertion(self):
        if self.annotate_split_features:
            self.add_partial_sites(self.target_contig, self.target_split_site)
            self.add_partial_sites(self.target_contig, self.target_split_site)
            self.add_partial_sites(self.donor_contig, self.donor_split_site)
            self.add_partial_sites(self.donor_contig, self.donor_split_site)

        result = (
            self.target_contig[:self.target_split_site]
            + self.donor_contig[self.donor_split_site:]
            + self.donor_contig[:self.donor_split_site]
            + self.target_contig[self.target_split_site:]
        )
        prefix = (f"Predicted recombination product using target site = "
                  f"{self.target_site} and donor site = {self.donor_site}.\n"
                  f"Donor sequence info:\n{self.donor_contig.description}\n"
                  "Target sequence info:")
        result = self.transfer_metadata(
            old_record=self.target_contig,
            new_record=result,
            new_description_prefix=prefix
        )
        return [result]

    def calculate_excision(self):
        cut_points = [self.target_split_site, self.donor_split_site]
        upstream_cut = min(cut_points)
        downstream_cut = max(cut_points)
        if self.annotate_split_features:
            self.add_partial_sites(self.target_contig, upstream_cut)
            self.add_partial_sites(self.target_contig, downstream_cut)
        leftovers = (self.target_contig[:upstream_cut] 
                     + self.target_contig[downstream_cut:])
        excised = self.target_contig[upstream_cut:downstream_cut]
        # Transfer metadata
        common_prefix = (
            f"excision between {self.target_split_site+1} and "
            f"{self.donor_split_site+1} with target site = {self.target_site}, "
            f"donor site = {self.donor_site}\nDerived from:"
                      )
        leftovers = self.transfer_metadata(
            old_record=self.target_contig,
            new_record=leftovers,
            new_description_prefix=f"Backbone sequence from {common_prefix}"
        )
        excised = self.transfer_metadata(
            old_record=self.target_contig,
            new_record=excised,
            new_description_prefix=f"Excised sequence from {common_prefix}"
        )
        return [leftovers, excised]

    def calculate_inversion(self):
        cut_points = [self.target_split_site, self.donor_split_site]
        upstream_cut = min(cut_points)
        downstream_cut = max(cut_points)
        if self.annotate_split_features:
            self.add_partial_sites(self.target_contig, upstream_cut)
            self.add_partial_sites(self.target_contig, downstream_cut)
        inverted_record = (
            self.target_contig[:upstream_cut]
            + self.target_contig[upstream_cut:downstream_cut].reverse_complement()
            + self.target_contig[downstream_cut:]
        )
        # Transfer metadata
        inversion_prefix = (
            f"Inversion between {self.target_split_site+1} and "
            f"{self.donor_split_site+1} with target site = {self.target_site}, "
            f"donor site = {self.donor_site}\nDerived from:"
                      )
        inverted_record = self.transfer_metadata(
            old_record=self.target_contig, 
            new_record=inverted_record, 
            new_description_prefix=inversion_prefix)
        return [inverted_record]

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
        # print(mode, " LE: ", le_start, le_end, " RE: ", re_start, re_end, " Strand:", strand, "core:", core_end - core_start)
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
    def find_partial_site(cls, part:str, site:str, contig_pos:int, strand=+1):
        """Finds locations of parts of a given site in relation to the location of
        the overall site in its native contig

        Args:
            part (str): Part of site. e.g. LT, RT, LD, RD
            site (str): Full sequence of target or donor site
            contig_pos (int): Start index of site in contig
            strand (int): Strand in which the site is found on the contig. 
                Must be either 1 or -1. Defaults to 1.
        Returns:
            (part_start, part_end): positions wrt the positive strand
        """    
        if strand == -1:
            part = reverse_complement(part)
            site = reverse_complement(site)
        part_start = site.index(part)
        part_end = part_start + len(part)
        # print(part_start, part_end, contig_pos, part, site)
        # part_start, part_end = find_matches(part, site, complement_search=False, 
        #                                     return_strand=False)
        part_start, part_end = contig_pos + part_start, contig_pos + part_end
        part_start, part_end = sorted([part_start, part_end])
        return part_start, part_end

    @classmethod
    def transfer_metadata(cls, old_record:SeqRecord, new_record:SeqRecord, 
                          new_description_prefix:str=""):
        """
        Transfer metadata from one SeqRecord to another, with an optional prefix
        added to the new record. 
        """
        new_description = f"{new_description_prefix}\n{old_record.description}"
        new_record.id = old_record.id
        new_record.name = old_record.name
        new_record.description = new_description
        new_record.annotations = old_record.annotations.copy()
        return new_record

    @classmethod
    def split_feature(cls, feature:SeqFeature, split_site:int):
        """Splits a sequence feature at a provided split site into an upstream
        and downstream component. 

        Args:
            feature (SeqFeature): Feature to split
            split_site (int): site to split at

        Returns:
            (upstream_feat, downstream_feat): New partial features
        """        
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        original_label = (
            feature.qualifiers.get("label", [feature.type])[0]
            if feature.qualifiers else feature.type
        )
        # Upstream part
        upstream_loc = FeatureLocation(start, split_site, strand=strand)
        upstream_feat = SeqFeature(upstream_loc, type="misc_feature", qualifiers={
            "label": [f"UPSTREAM SPLIT: {original_label}"],
            "partial": ["true"],
            "ApEinfo_fwdcolor": ["#C0C0C0"], 
            "ApEinfo_revcolor": ["#C0C0C0"],
            # "ApEinfo_graphicformat": ["arrow_data {{}}"], # Hide by default
        })
        # Downstream part
        downstream_loc = FeatureLocation(split_site, end, strand=strand)
        downstream_feat = SeqFeature(downstream_loc, type="misc_feature", qualifiers={
            "label": [f"DOWNSTREAM SPLIT: {original_label}"],
            "partial": ["true"],
            "ApEinfo_fwdcolor": ["#C0C0C0"],
            "ApEinfo_revcolor": ["#C0C0C0"],
            # "ApEinfo_graphicformat": ["arrow_data {{}}"], # Hide by default
        })
        return (upstream_feat, downstream_feat)
    
    # TODO: add logic to properly split compound location features that overlap origin of circular record
    @classmethod
    def add_partial_sites(cls, record:SeqRecord, cut_site:int, 
                          inplace:bool=True):
        """Splits all features in a SeqRecord that pass through a certain spot
        by adding new features corresponding to the upstream and downstream
        portions of that feature.

        Args:
            record (SeqRecord): Record to add sites to.
            cut_site (int): Site at which to split
            inplace (bool, optional): Option to modify record inplace. Defaults to True.

        Returns:
            record (SeqRecord): Modified record, only returned if inplace=False
        """        
        if not inplace:
            record = copy.deepcopy(record)
        split_features = []
        for feature in record.features:
            if feature.type != "source":
                starts_before_cut_site = feature.location.start < cut_site
                ends_before_cut_site = feature.location.end <= cut_site
                if starts_before_cut_site & (not ends_before_cut_site):
                    split_feats = cls.split_feature(feature, cut_site)
                    split_features.extend(split_feats)
        record.features.extend(split_features)
        if not inplace:
            return record


def write_genbank(record:SeqRecord, path:str):
    with open(path, "w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")
        print(f"Saved to {path}")


def reverse_complement_record(record:SeqRecord):
    # Reverse complement the sequence and features
    rev_record = record.reverse_complement()
    # Manually transfer metadata
    rev_record.id = record.id
    rev_record.name = record.name
    rev_record.description = record.description
    rev_record.dbxrefs = record.dbxrefs
    rev_record.annotations = record.annotations
    return rev_record


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
    parser.add_argument('-s', '--other_sequence_records', nargs='+', default=[],
                        help='Specify one or more additional sequence_records/genomes to ' \
                        'search for target and donor sites.')
    parser.add_argument("-b", "--brna_path", help="Path to bRNA plasmid",)
    parser.add_argument('--allow_linear', action='store_true',
                    help='Treat the output product as circular unless this is set.')
    parser.add_argument('--no_split_features', action='store_true',
                    help='Add new  partial features for any features that overlap the target/donor sites unless this is set.')
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory",)
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
    split_features = not args.no_split_features
    Experiment = BridgeRNAReaction(
        args.brna_path, args.other_sequence_records, split_features
    )
    results = Experiment.calculate_product()
    if not args.allow_linear:
        for record in results:
            record.annotations["topology"] = "circular"
    print(f"Saving results to {args.out_dir}")
    if Experiment.mode == "insertion":
        # write_genbank(Experiment.target_contig, "../misc/insertion_test/genome.gbk")
        write_genbank(Experiment.donor_contig, os.path.join(args.out_dir, "insert.gbk"))
        # write_genbank(Experiment.insert, "../misc/insertion_test/insert_oriented.gbk")
        insertion_handle = f"insertion_{Experiment.donor_file}_into_{Experiment.target_file}.gbk"
        # print(f"SAVE TO: {outfile}")
        write_genbank(results[0], os.path.join(args.out_dir, insertion_handle))
    elif Experiment.mode == "excision":
        # Experiment.target_split_site, Experiment.donor_split_site
        leftovers_handle = f"{Experiment.target_file}_{Experiment.target_split_site}-{Experiment.donor_split_site}_leftovers.gbk"
        excision_handle = f"{Experiment.target_file}_{Experiment.target_split_site}-{Experiment.donor_split_site}_excision.gbk"
        write_genbank(results[0], os.path.join(args.out_dir, leftovers_handle))
        write_genbank(results[1], os.path.join(args.out_dir, excision_handle))
    elif Experiment.mode == "inversion":
        inversion_handle = f"inversion_{Experiment.target_file}_{Experiment.target_split_site}-{Experiment.donor_split_site}.gbk"
        write_genbank(results[0], os.path.join(args.out_dir, inversion_handle))
