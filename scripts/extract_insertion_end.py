"""
Script to create a fasta file with the RE + up to 200 nt of the donor sequence
for insertion mapping

NOTE: might be smart to wrap in the prediction of products here, OR add a non-
verbose option to not print all the target sites etc. I feel like it's ugly to
invoke BridgeRNAReaction twice with the same inputs if indeed I wanna predict 
products in the pipeline
"""
from Bio import SeqIO
from Bio.Seq import reverse_complement, Seq
from predict_recombination_products import BridgeRNAReaction
from os.path import join
import argparse
import numpy as np
import os
import pandas as pd

import find_brna
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description = ("""
        Script to extract the end of an IS110-mediated insertion reaction.
        """)
    )
    parser.add_argument('-s', '--other_sequence_records', nargs='+', default=[],
                        help='Specify one or more additional sequence_records/genomes to ' \
                        'search for target and donor sites.')
    parser.add_argument("-b", "--brna_path", nargs='+', help="Path to bRNA plasmid(s)",)
    parser.add_argument("-n", "--extract_length", 
                        type=int, default=200,
                        help="Number of nucleotides to extract from RE. Defaults to 200",)
    parser.add_argument("-o", "--out_dir", default=".", help="Output directory",)
    parser.add_argument("-m", "--mode", default="R", choices=["R", "L"],help="Which end to extract",)
    parser.add_argument("--save_targets", action="store_true", 
                        help="Option to save a file with all target sites")
    args = parser.parse_args()
    return args 


# TODO: finish refactoring to BridgeRNAPlasmid. 
# The goal is to have methods that allow you to take in a given sequence record 
# and process that for downstream stuff -- basically just use this thing as an 
# operator to find target/donor sites for downstream applications. 


# class BridgeRNAPlasmid():
#     def __init__(self, brna_path):
#         self.brna_paths = brna_path
#         self._process_brna()

#     def _process_brna(self):
#         """
#         Process the Bridge RNA plasmid to extract target and donor site information.
        
#         This method reads the bRNA plasmid file and uses the find_brna module to
#         identify the target and donor sites, along with their left and right ends
#         and core sequences.
        
#         Note:
#             Currently assumes there's only one sequence record in the bRNA plasmid.
            
#         Sets:
#             brna_plasmid_record (SeqRecord): The bRNA plasmid sequence record
#             target_site (str): Full target site sequence
#             donor_site (str): Full donor site sequence
#             lt, rt (str): Left and right ends of target site
#             ld, rd (str): Left and right ends of donor site
#             target_core, donor_core (str): Core sequences of target and donor sites
#         """
#         self.brna_plasmid_record = SeqIO.read(self.brna_paths, "genbank")
#         target, donor = find_brna.find_target_donor_sites(self.brna_plasmid_record)
#         self.target_site, self.lt, self.target_core, self.rt = target
#         self.donor_site, self.ld, self.donor_core, self.rd = donor
#         print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")
    
    
#     @classmethod    
#     def find_matches(cls, site:str, seq:str):
#         """
#         Find matches for site in seq. 
#         """
#         # Find matches on positive and negative strands
#         matches_pos_strand = [match for match in re.finditer(site, seq)]
#         rc_contig = reverse_complement(seq)
#         matches_neg_strand = [match for match in re.finditer(site, rc_contig)]
#         matches = matches_pos_strand + matches_neg_strand
#         match_strands = np.concatenate([
#             np.ones_like(matches_pos_strand), np.ones_like(matches_neg_strand) * -1
#             ])
#         match_locations = []
#         L = len(seq)
#         for idx, match in enumerate(matches):
#             strand = match_strands[idx]
#             if strand == 1:
#                 site_start, site_end = match.start(), match.end()
#             else:
#                 site_start, site_end = sorted([L - match.start(), L - match.end()])
#             match_locations.append((strand, site_start, site_end))
#         return match_locations


# class BridgeRNA():
#     """
#     A class to represent a Bridge RNA plasmid and extract its target and donor 
#     sites. Must provide a filepath to a bRNA plasmid in GenBank format.
#     """

#     def __init__(self, brna_path:str):
#         """
        
#         """

#         self.brna_path = brna_path
#         self._process_brna()


#     def _process_brna(self):
#         """
#         Process the Bridge RNA plasmid to extract target and donor site information.
        
#         This method reads the bRNA plasmid file and uses the find_brna module to
#         identify the target and donor sites, along with their left and right ends
#         and core sequences. If multiple bRNA plasmids are provided, it ensures that
#         all bRNAs have the same target and donor sites.
        
#         Note:
#             Currently assumes there's only one sequence record bRNA plasmid.
            
#         Sets:
#             brna_plasmid_record (SeqRecord): Sequence record for bRNA plasmid
#             target_site (str): Full target site sequence
#             donor_site (str): Full donor site sequence
#             lt, rt (str): Left and right ends of target site
#             ld, rd (str): Left and right ends of donor site
#             target_core, donor_core (str): Core sequences of target and donor sites
#         """
#         brna_plasmid_records = []
#         targets, donors = [], []
#         brna_plasmid_record = SeqIO.read(self.brna_path, "genbank")
#         target, donor = find_brna.find_target_donor_sites(brna_plasmid_record)
#         targets.append(target)
#         donors.append(donor)
#         brna_plasmid_records.append(brna_plasmid_record)
#         self.brna_plasmid_records = brna_plasmid_records
#         # Assert that all targets and donors are the same
#         assert all(
#             target == targets[0] for target in targets
#         ), "Provided bRNA plasmids have different target sites!"
#         assert all(
#             donor == donors[0] for donor in donors
#         ), "Provided bRNA plasmids have different donor sites!"
#         target, donor = targets[0], donors[0]
#         self.target_site, self.lt, self.target_core, self.rt = target
#         self.donor_site, self.ld, self.donor_core, self.rd = donor
#         print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")


class Insertion():
    """
    A workaround for the BridgeRNAReaction class to allow detection of insertion
    sequence ends without providing the sequence to insert into.
    """
    def __init__(self, brna_paths:list, other_sequence_paths:list, 
                 annotate_split_features=True, verify_unique_target=True):
        self.brna_paths = brna_paths
        self.other_sequence_paths = other_sequence_paths
        self.annotate_split_features = annotate_split_features
        self.verify_unique_target = verify_unique_target
        self._process_brnas()
        # self._identify_sites()

    def _process_brnas(self):
        """
        Process the Bridge RNA plasmid(s) to extract target and donor site information.
        
        This method reads the bRNA plasmid file(s) and uses the find_brna module to
        identify the target and donor sites, along with their left and right ends
        and core sequences. If multiple bRNA plasmids are provided, it ensures that
        all bRNAs have the same target and donor sites.
        
        Note:
            Currently assumes there's only one sequence record in each bRNA plasmid.
            
        Sets:
            brna_plasmid_records (SeqRecord): Sequence records for bRNA plasmid(s)
            target_site (str): Full target site sequence
            donor_site (str): Full donor site sequence
            lt, rt (str): Left and right ends of target site
            ld, rd (str): Left and right ends of donor site
            target_core, donor_core (str): Core sequences of target and donor sites
        """
        brna_plasmid_records = []
        targets, donors = [], []
        for brna_path in self.brna_paths:
            brna_plasmid_record = SeqIO.read(brna_path, "genbank")
            target, donor = find_brna.find_target_donor_sites(brna_plasmid_record)
            targets.append(target)
            donors.append(donor)
            brna_plasmid_records.append(brna_plasmid_record)
        self.brna_plasmid_records = brna_plasmid_records
        # Assert that all targets and donors are the same
        assert all(
            target == targets[0] for target in targets
        ), "Provided bRNA plasmids have different target sites!"
        assert all(
            donor == donors[0] for donor in donors
        ), "Provided bRNA plasmids have different donor sites!"
        target, donor = targets[0], donors[0]
        self.target_site, self.lt, self.target_core, self.rt = target
        self.donor_site, self.ld, self.donor_core, self.rd = donor
        print(f"Using:\tTarget site: {self.target_site}\n\tDonor site: {self.donor_site}")


    def extract_split_positions(self, match_list:list, site_right_end:str):
        sites = []
        for match in match_list:
            if match[0] == 1:
                split_site = match[2] - len(site_right_end)
            else:
                split_site = match[1] + len(site_right_end)
            sites.append(split_site)
        return sites

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
            

    def find_and_orient_donor(self, seq):
        """Searches for donor sequence in provided sequence. If donor sequence
        is found, returns that sequence oriented to how it would appear relative
        to the target site in the insertion reaction. Otherwise, returns None.

        Args:
            seq (str): Sequence to search for donor site

        Raises:
            ValueError: This function currently assumes that there's only one 
            donor site in the provided sequence. 

        """
        donor_matches = self.find_matches(self.donor_site, seq)
        donor_split_sites = self.extract_split_positions(donor_matches, self.rd)
        if len(donor_matches) == 1:
            oriented_donor_seq = (
                seq[donor_split_sites[0]:] 
                + seq[:donor_split_sites[0]]
            )
            strand = donor_matches[0][0]
            if strand == -1:
                oriented_donor_seq = reverse_complement(oriented_donor_seq)
            return oriented_donor_seq
        elif len(donor_split_sites) == 0:
            return None
        else:
            raise ValueError("Multiple donor sites found in provided sequence!")

    
    def find_target_sites(self, seq):
        target_matches = self.find_matches(self.target_site, seq)
        sites = []
        for match in target_matches:
            if match[0] == 1:
                split_site = match[2] - len(self.rt)
            else:
                split_site = match[1] + len(self.rt)
            sites.append(split_site)
        return sites


def main():
    args = parse_args()

    Experiment = Insertion(
            args.brna_path, args.other_sequence_records, verify_unique_target=False
        )

    insertion_ends, donor_names, experiment_target_sites = [], [], []


    left_ends, right_ends = [], []

    
    # Search all provided sequence records for donor and target sites
    for seqfile in args.other_sequence_records + args.brna_path:
        for record in SeqIO.parse(seqfile, "genbank"):
            oriented_donor_seq = Experiment.find_and_orient_donor(str(record.seq))
            target_sites = Experiment.find_target_sites(str(record.seq))
            if oriented_donor_seq is not None:
                assert len(target_sites) == 0, "Insertion reaction cannot contain both donor and target ends!"
                le = oriented_donor_seq[:args.extract_length]
                re = oriented_donor_seq[-args.extract_length:]



                left_ends.append(le)
                right_ends.append(re)



                if args.mode == "R":
                    end = re
                elif args.mode == "L":
                    end = le
                insertion_ends.append(end)
                donor_names.append(record.name)
            elif len(target_sites) > 0:
                site_for_df = list(zip([record.name]*len(target_sites), target_sites))
                experiment_target_sites += site_for_df
    # If multiple donors found, ensure they're all the same
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO? Save both left and right ends? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if len(insertion_ends) > 1:
        print(f"Multiple donors found for {args.mode}E extraction...", end=" ")
        multiple_donors_error = f"Multiple different insertion ends found with N = {args.extract_length}!"
        assert all(i == insertion_ends[0] for i in insertion_ends), multiple_donors_error
        print(f" all donors have same end")
    end = insertion_ends[0]
    # Save insertion end to fasta
    print(f"found {args.mode}E: {end}")
    donor_info = ','.join(donor_names)
    if len(donor_names) > 1:
        donor_info += " (conserved)"
    header = f"{args.mode}E from {donor_info}"
    record = SeqIO.SeqRecord(id=header, description="", seq=Seq(end))
    with open(join(args.out_dir, f"{args.mode}E.fasta"), "w") as file:
        SeqIO.write(record, file, "fasta")
    # OPTIONAL: Save target indices as a tsv
    if args.save_targets:
        target_df = pd.DataFrame(experiment_target_sites, columns=["contig", "target_site"])
        target_df.to_csv(join(args.out_dir, "target_sites.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
