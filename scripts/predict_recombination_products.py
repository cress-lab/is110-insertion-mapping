import argparse
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
# from find_brna import find_target_donor_sites
import find_brna
import genbank_utils
import copy

# class bRNAReaction():
#     def __init__(self, brna_path, other_plasmid_paths, out_dir):
#         self.brna_path = brna_path
#         self.other_plasmid_paths = other_plasmid_paths
#         self.out_dir = out_dir
#         self._process_brna()
#         self._read_dna_files()
#         pass

#     def _process_brna(self):
#         brna_plasmid_record = SeqIO.read(self.brna_path, "genbank")
#         target, donor = find_brna.find_target_donor_sites(brna_plasmid_record)
#         self.target_site, self.lt, self.target_core, self.rt = target
#         self.donor_site, self.ld, self.donor_core, self.rd = donor

#         print(f"Using:\tTarget site: {target_site}\n\tDonor site: {donor_site}")

#     def _read_dna_files(self):
#         self.other_plasmid_records = [
#             SeqIO.read(file, "genbank") for file in self.other_plasmid_paths
#             ]


def determine_mode(*args):
    return "insertion"


def count_sites(site, sequence):
    matches_pos_strand = [match for match in re.finditer(site, sequence)]
    matches_neg_strand = [match for match in re.finditer(
                            site, reverse_complement(sequence))]
    pos_strand_hits = len(matches_pos_strand)
    neg_strand_hits = len(matches_neg_strand)
    return (pos_strand_hits, neg_strand_hits)


def determine_excision_or_inversion():
    # If target and donor are on the same plasmid, figure out if it's gonna be 
    # an excision or inversion
    pass


def annotate_plasmid_as_target_or_donor_or_both(record, target_site, donor_site):
    target_cts = sum(count_sites(record, target_site))
    donor_cts = sum(count_sites(record, donor_site))

    unique_donor = donor_cts == 1
    unique_target = target_cts == 1
    is_donor, is_target = False, False
    if unique_donor:
        is_donor = True
    elif donor_cts > 1:
        raise ValueError(f"Multiple donor sites found in {record.id}!")
    if unique_target:
        is_target = True
    elif target_cts > 1:
        raise ValueError(f"Multiple target sites found in {record.id}!")
    return is_donor, is_target

    # both = False

    # if both:
    #     determine_relative_orientation()


def identify_plasmids(
        brna_plasmid_record, other_plasmids, target_site, donor_site):
    plasmids = [brna_plasmid_record] + other_plasmids
    identities = [
        annotate_plasmid_as_target_or_donor_or_both(i, target_site, donor_site)
        for i in plasmids
    ]
    donor_cts = [i[0] for i in identities]
    target_cts = [i[1] for i in identities]

    one_donor_site = donor_cts.count(True) == 1
    one_target_site = target_cts.count(True) == 1

    if one_donor_site and one_target_site:
        donor_idx = np.where(donor_cts)[0].item()
        target_idx = np.where(target_cts)[0].item()
        if donor_idx == target_idx:
            target_plasmid = plasmids[target_idx]
            mode = determine_excision_or_inversion(target_plasmid)
            return mode, (target_plasmid)
        else:
            target_plasmid = plasmids[target_idx]
            donor_plasmid = plasmids[donor_idx]
            mode = "insertion"
            return mode, (target_plasmid, donor_plasmid)

    # brna_donor_hits = sum(count_sites(brna_plasmid_record, donor_site))
    # brna_target_hits = sum(count_sites(brna_plasmid_record, target_site))


def find_matches(site, contig):
    """
    Find matches for site in contig. 
    """
    # Find matches on positive and negative strands
    matches_pos_strand = [match for match in re.finditer(site, contig)]
    rc_contig = reverse_complement(contig)
    matches_neg_strand = [match for match in re.finditer(site, rc_contig)]
    matches = matches_pos_strand + matches_neg_strand
    match_strands = np.concatenate([
        np.ones_like(matches_pos_strand), np.ones_like(matches_neg_strand) * -1
        ])
    match_locations = []
    for idx, match in enumerate(matches):
        site_start, site_end = match.start(), match.end()
        strand = match_strands[idx]
        match_locations.append((strand, site_start, site_end))
    return match_locations


def find_unique_match(site, contig):
    match_locations = find_matches(site, contig)
    n_matches = len(match_locations)
    if n_matches == 1:
        strand, site_start, site_end = match_locations[0]
        return strand, site_start, site_end
    elif n_matches == 0:
        raise ValueError("Found 0 sites in contig!")
    else:
        raise ValueError("Found multiple sites in contig!")


def find_partial_site(part, site, contig_pos):
    """Finds locations of parts of a given site in relation to the location of
    the overall site in its native contig

    Args:
        part (str): Part of site. e.g. LT, RT, LD, RD
        site (str): Full sequence of target or donor site
        contig_pos (str): Start index of site in contig
    """    
    part_start = site.index(part)
    part_end = part_start + len(part)
    # part_start, part_end = find_matches(part, site, complement_search=False, 
    #                                     return_strand=False)
    part_start, part_end = contig_pos + part_start, contig_pos + part_end
    return part_start, part_end


def annotate_site(record, site_start, site_seq, le_seq, re_seq, mode="Target", inplace=False):
    allowed_modes = ["Target", "Donor"]
    assert mode in allowed_modes, f"mode must be one of: {allowed_modes}"

    # print(site_start, site_seq, le_seq, core, re_seq)
    colors = {
        "Target": "#00aeff",
        "Donor": "#ff8400",
        "Core": "#24951a",
    }
    # Extract indices of LT, RT relative to their respective contigs
    le_start, le_end = find_partial_site(le_seq, site_seq, site_start)
    re_start, re_end = find_partial_site(re_seq, site_seq, site_start)
    print(mode, le_start, le_end, re_start, re_end, re_end - le_start)
    # Annotate the partial sites
    total_feat = SeqFeature(
        location=FeatureLocation(le_start, re_end), 
        type="misc_feature", 
        qualifiers = {
            "label": f"{mode}",
            "ApEinfo_fwdcolor": colors[mode],
            "ApEinfo_revcolor": colors[mode],
            "note": f"color: {colors[mode]}"
        })
    le_feat = SeqFeature(
        location=FeatureLocation(le_start, le_end), 
        type="misc_feature", 
        qualifiers = {
            "label": f"{mode} LE",
            "ApEinfo_fwdcolor": colors[mode],
            "ApEinfo_revcolor": colors[mode],
            "note": f"color: {colors[mode]}"
        })
    re_feat = SeqFeature(
        location=FeatureLocation(re_start, re_end), 
        type="misc_feature", 
        qualifiers = {
            "label": f"{mode} RE",
            "ApEinfo_fwdcolor": colors[mode],
            "ApEinfo_revcolor": colors[mode],
            "note": f"color: {colors[mode]}"
        })
    core_feat = SeqFeature(
        location=FeatureLocation(le_end, re_start), 
        type="misc_feature", 
        qualifiers = {
            "label": f"{mode} Core",
            "ApEinfo_fwdcolor": colors["Core"],
            "ApEinfo_revcolor": colors["Core"],
            "note": f"color: {colors['Core']}"
        })
    # new_feats = [le_feat, core_feat, re_feat, total_feat]
    new_feats = [le_feat, core_feat, re_feat,]
    if inplace:
        record.features += new_feats
    else:
        new_record = copy.deepcopy(record)
        new_record.features += new_feats
        return new_record
    

def shift_feature(feature, offset):
    shifted_feature = copy.deepcopy(feature)
    shifted_feature.location = FeatureLocation(
        start=int(feature.location.start) + offset,
        end=int(feature.location.end) + offset,
        strand=feature.location.strand
    )
    return shifted_feature


def calculate_insertion(donor_record, target_record, donor_split_pos, target_split_pos):
    """

    Args:
        donor_record (_type_): _description_
        target_record (_type_): _description_
        donor_split_pos (_type_): donor_start + len(ld) + 2
        target_split_pos (_type_): target_start + len(lt) + 2
    """
    insertion_size = len(donor_record.seq)
    donor_record = genbank_utils.orient_record_around_pos(donor_record, donor_split_pos)
    with open("../misc/shifted_donor_test.gbk", "w") as f:
        SeqIO.write(donor_record, f, "genbank")
    with open("../misc/target.gbk", "w") as f:
        SeqIO.write(target_record, f, "genbank")
    # Figure out the new sequence
    new_seq = (
        target_record.seq[:target_split_pos]
        + donor_record.seq
        + target_record.seq[target_split_pos:]
        )

    # Transfer features
    features = []
    for feature in donor_record.features:
        # if "Donor Core" in feature.qualifiers.get("label", ""):
        #     print(feature.qualifiers.get("label", ""))
        #     print("Original: ", feature)
        #     print("Shifted: ", shift_feature(feature, target_split_pos))
        shifted_feature = shift_feature(feature, target_split_pos)
        features.append(shifted_feature)
        # feature.location.start = feature.location.start + target_split_pos
        # feature.location.end += target_split_pos
    
    for feat in target_record.features:
        feature = copy.deepcopy(feat)
        ends_before_insertion = feature.location.end <= target_split_pos
        starts_before_insertion = feature.location.start < target_split_pos
        before_insertion = starts_before_insertion & ends_before_insertion
        # if ends_before_insertion:
        #     continue
        # If feature overlaps insertion, just drag it across the insertion. 
        if starts_before_insertion and not ends_before_insertion:
            feature.location = FeatureLocation(
                start=int(feature.location.start),
                end=int(feature.location.end) + insertion_size,
                strand=feature.location.strand
            )
        elif not before_insertion:
            feature = shift_feature(feature, insertion_size)
        features.append(feature)
    # features = donor_record.features + target_record.features
    # Create the final record
    result_id = "recombined product" # I'll fill this in later.
    result_name = "recombined_product" # I'll fill this in later.
    result_description = f"Predicted recombination product using target site = {target_site} and donor site = {donor_site}"

    recombinant_record = SeqRecord(new_seq, id=result_id, name=result_name, 
                                   description=result_description)
    # Mark the record as circular
    print("Assuming circular topology for insertion")
    recombinant_record.annotations["topology"] = "circular"
    recombinant_record.annotations["molecule_type"] = "DNA"
    # Apply features to genbank record and save file
    recombinant_record.features = features
    return recombinant_record



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
parser.add_argument("-t", "--target_contig", help="Path to the genome/contig with target site")
parser.add_argument("-d", "--donor_contig", help="Path to the genome/contig with donor site", required=False)
parser.add_argument('--other_sequences', nargs='+', 
                    help='Specify one or more additional plasmids/genomes to ' \
                    'search for target and donor sites.')

parser.add_argument(
    "--brna_plasmid", 
    help="Path to bRNA plasmid",
    )

parser.add_argument("-o", "--output", required=True, help="Output .gbk file path",)

# Test case
args = parser.parse_args(
    [
        "--target_contig", "/groups/cress/projects/is110/sequences/SLMQRT/Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk",
        # "--donor_contig", "/groups/cress/projects/is110/sequences/SLMQRT/monomeric_BAC_for_plate_scrape_13kb.gbk",
        "--donor_contig", "/groups/cress/projects/efactors/jack/insertion_mapping/misc/W9G3HH_30kbp BAC_gSES0133.gbk",
        "--brna_plasmid", "/groups/cress/projects/is110/sequences/SLMQRT/P7YMY9_pSES0099.gbk",
        "-o", "/groups/cress/projects/efactors/jack/insertion_mapping/misc/TEST.gbk", 
        ])
# args = parser.parse_args()

# Read target, donor sites etc from bRNA plasmid
brna_plasmid_record = SeqIO.read(args.brna_plasmid, "genbank")
target, donor = find_target_donor_sites(brna_plasmid_record)
target_site, lt, target_core, rt = target
donor_site, ld, donor_core, rd = donor

print(f"Using:\tTarget site: {target_site}\n\tDonor site: {donor_site}")

mode = determine_mode()
print(f"\tTarget contig: {args.target_contig}")
print(f"\tDonor contig: {args.donor_contig}")

print("Reading target contig...", end=" ")
target_contig_record = SeqIO.read(args.target_contig, "genbank")
target_contig = str(target_contig_record.seq)

print("Reading donor contig...", end=" ")
donor_contig_record = SeqIO.read(args.donor_contig, "genbank")
donor_contig = str(donor_contig_record.seq)

print("Finding target site in target contig")
target_strand, target_start, target_stop = find_unique_match(target_site, target_contig)
if target_strand == -1:
    target_contig_record = genbank_utils.reverse_complement_record(target_contig_record)
# Add target site features
target_contig_record = annotate_site(
    target_contig_record, target_start, target_site, lt, rt, 
    mode="Target")


print("Finding donor site in donor contig")
donor_strand, donor_start, donor_stop = find_unique_match(donor_site, donor_contig) # NOTE: this is wrt - strand if applicable
if donor_strand == -1:
    donor_contig_record = genbank_utils.reverse_complement_record(donor_contig_record)
# Add donor site features
donor_contig_record = annotate_site(
    donor_contig_record, donor_start, donor_site, ld, rd, 
    mode="Donor")

if mode == "insertion":
    print("Calculating results in insertion mode...")
    donor_split_site = donor_start + 2 + len(ld)
    target_split_site = target_start + 2 + len(lt)
    recombinant_record = calculate_insertion(
        donor_contig_record, target_contig_record,
        donor_split_site, target_split_site
        )
    with open(args.output, "w") as output_handle:
        SeqIO.write(recombinant_record, output_handle, "genbank")
    print(f"Saved results to {args.output}")