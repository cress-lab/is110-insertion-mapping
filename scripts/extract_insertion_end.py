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
        self._identify_donor()

    def _identify_donor(self):
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

        target_match_contigs = [bool(i) for i in target_matches]
        target_index = np.where(target_match_contigs)[0].item() if any(target_match_contigs) else None
        donor_index = np.where([bool(i) for i in donor_matches])[0].item()
        assert (target_index != donor_index), "Target and Donor found on same contig, so not an insertion!"
        self.donor_strand, self.donor_start, self.donor_end = donor_matches[donor_index][0]
        self.mode = "insertion"
        # Store separate donor record for insertions
        self.donor_file = os.path.splitext(os.path.basename(sequence_paths[donor_index]))[0]
        print(f"Found donor sequence: {self.donor_file}\n\tstrand: {self.donor_strand}")
        self.donor_contig = sequence_records[donor_index]
        self._find_donor_site()

    def _find_donor_site(self):
        # Store donor site
        if self.donor_strand == 1:
            self.donor_split_site = self.donor_end - len(self.rt)
        else:
            self.donor_split_site = self.donor_start + len(self.rt)
            

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

    print(Experiment.donor_strand)
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

if __name__ == "__main__":
    main()
