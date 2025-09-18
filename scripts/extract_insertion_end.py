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

parser = argparse.ArgumentParser(
    description = ("Script to predict IS110 recombination product given target "
    "and donor contigs, left and right ends of target and donor sites, and "
    "cores of each. ")
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

Experiment = BridgeRNAReaction(
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
