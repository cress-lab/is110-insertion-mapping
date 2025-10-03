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
    """
    A class for predicting Bridge RNA (bRNA) recombination products.
    
    This class analyzes Bridge RNA plasmids and target/donor sequences to predict
    the outcome of IS110-mediated recombination events. It supports three types
    of rearrangements: insertion, excision, and inversion.
    
    The class uses a modular approach with separate methods for:
    - Loading and parsing sequence files
    - Finding target and donor sites
    - Determining rearrangement mode
    - Annotating sites on sequences
    - Calculating recombination products
    
    Attributes:
        brna_path (str): Path to the Bridge RNA plasmid file
        other_sequence_paths (list): List of paths to additional sequence files
        annotate_split_features (bool): Whether to add split feature annotations
        verify_unique_target (bool): Whether to verify target site uniqueness
        mode (str): Type of rearrangement ("insertion", "excision", or "inversion")
        target_contig (SeqRecord): Target sequence record
        donor_contig (SeqRecord): Donor sequence record (for insertions)
        target_site (str): Target site sequence
        donor_site (str): Donor site sequence
        
    Example:
        >>> reaction = BridgeRNAReaction(
        ...     brna_path="brna_plasmid.gbk",
        ...     other_sequence_paths=["target.gbk", "donor.gbk"]
        ... )
        >>> products = reaction.calculate_product()
        >>> print(f"Recombination mode: {reaction.mode}")
    """
    
    def __init__(self, brna_path, other_sequence_paths, 
                 annotate_split_features=True, verify_unique_target=True):
        """
        Initialize a Bridge RNA recombination reaction.
        
        Args:
            brna_path (str): Path to the Bridge RNA plasmid file (GenBank format)
            other_sequence_paths (list): List of paths to additional sequence files
                containing target and/or donor sequences
            annotate_split_features (bool, optional): Whether to add annotations
                for features that are split by recombination sites. Defaults to True.
            verify_unique_target (bool, optional): Whether to verify that the target
                site appears only once across all sequences. Defaults to True.
                
        Raises:
            ValueError: If target or donor sites are not found or are not unique
            FileNotFoundError: If any of the input files cannot be found
        """
        self.brna_path = brna_path
        self.other_sequence_paths = other_sequence_paths
        self.annotate_split_features = annotate_split_features
        self.verify_unique_target = verify_unique_target
        self._process_brna()
        self._identify_sequences()

    def _process_brna(self):
        """
        Process the Bridge RNA plasmid to extract target and donor site information.
        
        This method reads the bRNA plasmid file and uses the find_brna module to
        identify the target and donor sites, along with their left and right ends
        and core sequences.
        
        Note:
            Currently assumes there's only one sequence record in the bRNA plasmid.
            
        Sets:
            brna_plasmid_record (SeqRecord): The bRNA plasmid sequence record
            target_site (str): Full target site sequence
            donor_site (str): Full donor site sequence
            lt, rt (str): Left and right ends of target site
            ld, rd (str): Left and right ends of donor site
            target_core, donor_core (str): Core sequences of target and donor sites
        """
        self.brna_plasmid_record = SeqIO.read(self.brna_path, "genbank")
        target, donor = find_brna.find_target_donor_sites(self.brna_plasmid_record)
        self.target_site, self.lt, self.target_core, self.rt = target
        self.donor_site, self.ld, self.donor_core, self.rd = donor
        print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")

    def _verify_unique_site(self, matches, mode):
        """
        Verify that a site appears exactly once across all sequences.
        
        Args:
            matches (list): List of match results for each sequence
            mode (str): Type of site being verified ("target" or "donor")
            
        Raises:
            ValueError: If the site is not found exactly once
        """
        unique_plasmid = len([i for i in matches if i]) == 1
        unique_site = sum([len(i) for i in matches]) == 1
        if not unique_plasmid:
            raise ValueError(f"Expected exactly one {mode} plasmid match, but found {len(matches)}.")

        if not unique_site:
            total_sites = sum(len(i) for i in matches)
            raise ValueError(f"Expected exactly one {mode} site match, but found {total_sites} total.")
        
    def _identify_sequences(self):
        """
        Main method to identify target and donor sequences and set up the reaction.
        
        This method orchestrates the entire sequence identification process:
        1. Loads all sequence records from files
        2. Finds target and donor sites in the sequences
        3. Sets up target and donor properties
        4. Calculates split sites
        5. Determines the type of rearrangement
        6. Annotates sites on sequences
        7. Handles donor orientation for insertions
        
        This is the core method that sets up all the necessary attributes
        for calculating recombination products.
        """
        sequence_records, sequence_paths = self._load_sequence_records()
        target_index, donor_index = self._find_target_donor_sites(sequence_records)
        self._set_target_donor_properties(sequence_records, sequence_paths, target_index, donor_index)
        self._find_split_sites()
        self._determine_rearrangement_mode(target_index, donor_index)
        self._annotate_target_donor_sites()
        if self.mode == "insertion":
            self._handle_donor_orientation()

    def _load_sequence_records(self):
        """
        Load all sequence records from bRNA and other sequence files.
        
        Returns:
            tuple: (sequence_records, sequence_paths) where:
                - sequence_records: List of SeqRecord objects
                - sequence_paths: List of corresponding file paths
        """
        sequence_paths = [self.brna_path]
        other_sequence_records = []
        
        for filepath in self.other_sequence_paths:
            records = SeqIO.parse(filepath, "genbank")
            for record in records:
                other_sequence_records.append(record)
                sequence_paths.append(filepath)
        
        sequence_records = [self.brna_plasmid_record] + other_sequence_records
        return sequence_records, sequence_paths

    def _find_target_donor_sites(self, sequence_records):
        """
        Find target and donor sites in the sequence records.
        
        Args:
            sequence_records (list): List of SeqRecord objects to search
            
        Returns:
            tuple: (target_index, donor_index) - indices of sequences containing 
                sites in sequence_records
            
        Sets:
            target_strand, target_start, target_end: Target site coordinates
            donor_strand, donor_start, donor_end: Donor site coordinates
        """
        seqs = [str(record.seq) for record in sequence_records]
        
        # Search all sequence records for target and donor sites
        target_matches = [self.find_matches(self.target_site, seq) for seq in seqs]
        donor_matches = [self.find_matches(self.donor_site, seq) for seq in seqs]
        
        # Verify uniqueness
        if self.verify_unique_target:
            self._verify_unique_site(target_matches, "target")
        self._verify_unique_site(donor_matches, "donor")
        
        # Get indices of sequences containing the sites
        target_index = np.where([bool(matches) for matches in target_matches])[0].item()
        donor_index = np.where([bool(matches) for matches in donor_matches])[0].item()
        
        # Extract site coordinates
        self.target_strand, self.target_start, self.target_end = target_matches[target_index][0]
        self.donor_strand, self.donor_start, self.donor_end = donor_matches[donor_index][0]
        
        return target_index, donor_index

    def _set_target_donor_properties(self, sequence_records, sequence_paths, target_index, donor_index):
        """
        Set target and donor contig properties.
        
        Args:
            sequence_records (list): List of SeqRecord objects
            sequence_paths (list): List of corresponding file paths
            target_index (int): Index of sequence containing target site
            donor_index (int): Index of sequence containing donor site
            
        Sets:
            target_file, target_contig: Target sequence file name and record
            donor_file, donor_contig: Donor sequence file name and record (for insertions)
        """
        # Set target properties
        self.target_file = os.path.splitext(os.path.basename(sequence_paths[target_index]))[0]
        self.target_contig = sequence_records[target_index]
        print(f"Found target sequence: {self.target_file}\n\tstart: {self.target_start}\n\tend: {self.target_end}\n\tstrand: {self.target_strand}")
        
        # Set donor properties (for insertion mode)
        if target_index != donor_index:
            self.donor_file = os.path.splitext(os.path.basename(sequence_paths[donor_index]))[0]
            self.donor_contig = sequence_records[donor_index]
            print(f"Found donor sequence: {self.donor_file}\n\tstrand: {self.donor_strand}")

    def _determine_rearrangement_mode(self, target_index, donor_index):
        """
        Determine the type of rearrangement based on target and donor locations.
        
        Args:
            target_index (int): Index of sequence containing target site
            donor_index (int): Index of sequence containing donor site
            
        Sets:
            mode (str): Type of rearrangement:
                - "insertion": Target and donor on different sequences
                - "excision": Target and donor on same sequence, same strand
                - "inversion": Target and donor on same sequence, opposite strands
        """
        if target_index == donor_index:
            # Same sequence - either excision or inversion
            if self.target_strand == self.donor_strand:
                self.mode = "excision"
            else:
                self.mode = "inversion"
        else:
            # Different sequences - insertion
            self.mode = "insertion"

    def _annotate_target_donor_sites(self):
        """
        Annotate target and donor sites on their respective contigs.
        
        This method adds GenBank feature annotations to mark the target and donor
        sites, including left ends, right ends, and core sequences. The annotations
        are added in-place to the sequence records.
        """
        # Always annotate target site
        self.annotate_site(
            record=self.target_contig, 
            site_start=self.target_start, 
            site_seq=self.target_site, 
            le_seq=self.lt, 
            re_seq=self.rt, 
            strand=self.target_strand, 
            mode="Target", 
            inplace=True
        )
        
        # Annotate donor site
        if self.mode == "insertion":
            # Donor is on separate contig
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
        else:
            # Donor is on same contig as target
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

    def _handle_donor_orientation(self):
        """
        Handle donor orientation for insertion mode to ensure same orientation as target.
        
        For insertion mode, the donor sequence needs to be in the same orientation
        as the target sequence. If they are on opposite strands, this method
        reverse complements the donor sequence and updates all coordinates.
        """
        if self.donor_strand != self.target_strand:
            print("Reverse complementing donor contig...")
            self.donor_contig = reverse_complement_record(self.donor_contig)
            donor_length = len(self.donor_contig.seq)
            
            # Update coordinates after reverse complement
            self.donor_split_site = donor_length - self.donor_split_site
            self.donor_start = donor_length - self.donor_start
            self.donor_end = donor_length - self.donor_end
            self.donor_strand = -self.donor_strand

    def _find_split_sites(self):
        """
        Calculate the split sites for target and donor sequences.
        
        Split sites are the positions where the recombination will occur,
        determined by the position of the right end within the target or donor 
        site (RT or RD, not to be confused with the RE of the overall 
        insertion).
        
        Sets:
            target_split_site (int): Position where target sequence will be split
            donor_split_site (int): Position where donor sequence will be split
        """
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
        """
        Calculate the recombination product based on the determined mode.
        
        Returns:
            list: List of SeqRecord objects representing the recombination products
            
        Raises:
            ValueError: If the mode is invalid
        """
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
