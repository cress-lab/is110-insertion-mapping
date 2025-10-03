"""
Script to create a fasta file with the RE + up to 200 nt of the donor sequence
for insertion mapping

NOTE: might be smart to wrap in the prediction of products here, OR add a non-
verbose option to not print all the target sites etc. I feel like it's ugly to
invoke BridgeRNAReaction twice with the same inputs if indeed I wanna predict 
products in the pipeline
"""
from Bio import SeqIO
from Bio.Seq import reverse_complement
from predict_recombination_products import BridgeRNAReaction
from os.path import join
import argparse
import numpy as np
import os


class Insertion(BridgeRNAReaction):
    """
    A workaround for the BridgeRNAReaction class to allow detection of insertion
    sequence ends without providing the sequence to insert into.
    """
    def __init__(self, brna_path, other_sequence_paths, 
                 annotate_split_features=True, verify_unique_target=True):
        self.brna_path = brna_path
        self.other_sequence_paths = other_sequence_paths
        self.annotate_split_features = annotate_split_features
        self.verify_unique_target = verify_unique_target
        self._process_brna()
        self._identify_sites()

    def _identify_sites(self):
        # Read other sequence records
        sequence_paths = [self.brna_path]
        # Store all sequence records in additional sequence files
        other_sequence_records = []
        for filepath in self.other_sequence_paths:
            records = SeqIO.parse(filepath, "genbank")
            for record in records:
                    other_sequence_records.append(record)
                    sequence_paths.append(filepath)
        sequence_records = [self.brna_plasmid_record] + other_sequence_records
        seqs = [str(record.seq) for record in sequence_records]
        # Search all sequence records for target and donor sites
        donor_index = self._identify_donor(seqs, sequence_paths, sequence_records)
        target_index = self._identify_target(seqs, sequence_paths, sequence_records)
        self.target_matches = self.target_matches[target_index] if target_index is not None else self.target_matches
        self.donor_matches = self.donor_matches[target_index] if target_index is not None else self.donor_matches
        print(f"\n\n\n\t\tFOUND TARGET MATCHES: {self.target_matches}")
        assert (target_index != donor_index), "Target and Donor found on same contig, so not an insertion!"
        self.mode = "insertion"
        self._find_donor_site()
        

    def _identify_donor(self, seqs, sequence_paths, sequence_records):
        """
        Identifies and annotates donor site in order to extract insertion 
        ends. 

        Args:
            seqs (list): List of sequences to search
            sequence_paths (list): List of paths to sequences to search
            sequence_records: List of sequence records

        Returns:
            donor_index: Index of donor sequence in seqs
        Sets:
            donor_matches: List of tuples of (strand, start, end) for donor sites
            donor_contig (SeqRecord): Donor sequence record
        """        
        self.donor_matches = [
            self.find_matches(self.donor_site, i) for i in seqs]
        donor_index = np.where([bool(i) for i in self.donor_matches])[0].item()
        self.donor_strand, self.donor_start, self.donor_end = self.donor_matches[donor_index][0]
        self.donor_file = os.path.splitext(os.path.basename(sequence_paths[donor_index]))[0]
        self.donor_contig = sequence_records[donor_index]
        print(f"Found donor sequence: {self.donor_file}\n\tstrand: {self.donor_strand}")
        return donor_index
    
    def _identify_target(self, seqs, sequence_paths, sequence_records):
        """
        Identifies and annotates target site(s). Since this class isn't
        calculating the full reaction product, we only need to know which 
        sequences have target sites and where those are. 

        Args:
            seqs (list): List of sequences to search
            sequence_paths (_type_): List of paths to sequences to search
            sequence_records: List of sequence records

        Returns:
            target_index: Index of target sequence in seqs
        Sets:
            target_matches: List of tuples of (strand, start, end) for target sites
            target_contig (SeqRecord): Target sequence record (if target found)
        """        
        self.target_matches = [
            self.find_matches(self.target_site, i) for i in seqs]
        target_match_contigs = [bool(i) for i in self.target_matches]
        target_index = (
            np.where(target_match_contigs)[0].item() 
            if any(target_match_contigs) 
            else None
            )
        if target_index is not None:
            self.target_contig = sequence_records[target_index]
        return target_index


    def _find_donor_site(self):
        """
        Store donor site as class attribute
        """
        if self.donor_strand == 1:
            self.donor_split_site = self.donor_end - len(self.rt)
        else:
            self.donor_split_site = self.donor_start + len(self.rt)
    
    def extract_target_positions(self):
        sites = []
        for match in self.target_matches:
            if match[0] == 1:
                split_site = match[2] - len(self.rt)
            else:
                split_site = match[1] + len(self.rt)
            sites.append(split_site)
        return sites

            

def main():
    parser = argparse.ArgumentParser(
        description = ("""
        Script to extract the end of an IS110-mediated insertion reaction.
        """)
    )
    parser.add_argument('-s', '--other_sequence_records', nargs='+', default=[],
                        help='Specify one or more additional sequence_records/genomes to ' \
                        'search for target and donor sites.')
    parser.add_argument("-b", "--brna_path", help="Path to bRNA plasmid",)
    parser.add_argument("-n", "--extract_length", 
                        type=int, default=200,
                        help="Number of nucleotides to extract from RE. Defaults to 200",)
    parser.add_argument("-o", "--out_dir", default=".", help="Output directory",)
    parser.add_argument("-m", "--mode", default="R", choices=["R", "L"],help="Which end to extract",)
    parser.add_argument("--save_targets", action="store_true", 
                        help="Option to save a file with all target sites")
    args = parser.parse_args()

    Experiment = Insertion(
            args.brna_path, args.other_sequence_records, verify_unique_target=False
        )
    NotInsertionMessage = f"Provided set of plasmids was expected to be for an insertion reaction, got {Experiment.mode} instead"
    assert Experiment.mode == "insertion", NotInsertionMessage

    oriented_donor_seq = (
        Experiment.donor_contig.seq[Experiment.donor_split_site:] 
        + Experiment.donor_contig.seq[:Experiment.donor_split_site]
    )

    if Experiment.donor_strand == 1:
        le = oriented_donor_seq[:args.extract_length]
        re = oriented_donor_seq[-args.extract_length:]
    else:
        le = reverse_complement(oriented_donor_seq)[:args.extract_length]
        re = reverse_complement(oriented_donor_seq)[-args.extract_length:]
    if args.mode == "R":
        end = re
    elif args.mode == "L":
        end = le

    print(f"found {args.mode}E: {end}")
    header = f"{args.mode}E from {Experiment.donor_file}"
    record = SeqIO.SeqRecord(id=header, description="", seq=end)
    with open(join(args.out_dir, f"{args.mode}E.fasta"), "w") as file:
        SeqIO.write(record, file, "fasta")
    # OPTIONAL: Save target indices
    if args.save_targets:
        sites = Experiment.extract_target_positions()
        print(Experiment.target_contig.name)
        with open(join(args.out_dir, "target_sites.txt"), "w") as f:
            f.writelines(f"{Experiment.target_contig.name}\n")
            for i in sites:
                f.writelines(f"{i}\n")


if __name__ == "__main__":
    main()
