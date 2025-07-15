import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq

# 5' ----- LT - target core - RT ----- 3' 
# 5' ----- LD - donor core - RD ----- 3' 


parser = argparse.ArgumentParser(
    description = "Script to predict IS110 recombination product given target "
    "and donor contigs, left and right ends of target and donor sites, and " \
    "cores of each. "
)
parser.add_argument("-tg", "--target_contig", help="Path to the genome/contig with target site")
parser.add_argument("-dg", "--donor_contig", help="Path to the genome/contig with donor site", required=False)


parser.add_argument("-lt", "--left_target", help="",)
parser.add_argument("-rt", "--right_target", help="",)
parser.add_argument("-ld", "--left_donor", help="",)
parser.add_argument("-rd", "--right_donor", help="",)
parser.add_argument("-tc", "--target_core", default="CT",
                    help="",)
parser.add_argument("-dc", "--donor_core", required=False,
                    help="",)

# Test case
args = parser.parse_args(
    [
        "--target_contig", "/groups/cress/projects/is110/sequences/SLMQRT/Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk",
        "--donor_contig", "/groups/cress/projects/is110/sequences/SLMQRT/monomeric_BAC_for_plate_scrape_13kb.gbk",
        "-lt", "attcat",
        "-tc", "ct",
        "-rt", "gcggat",
        "-ld", "CCTATTATGTT",
        "-dc", "CT",
        "-rd", "ATGACAATA", # NOTE: not sure if this donor seq is correct
        ])


same_contig = False

if not args.donor_core:
    args.donor_core = args.target_core
    print(f"Using donor core = target core = {args.target_core}")
else:
    print(f"donor core = {args.donor_core}; target core = {args.target_core}")


if not args.donor_contig:
    print("No donor contig provided. Assuming donor and target on same contig...")
    # Need logic to look for donor on target contig
    same_contig = True
    ...
    # If target is on one strand and donor is on the other, will need to complement all elements of donor.




lt = args.left_target.upper()
rt = args.right_target.upper()
ld = args.left_donor.upper()
rd = args.right_donor.upper()
tc = args.target_core.upper()
dc = args.donor_core.upper()

target_site = lt + tc + rt
donor_site = ld + dc + rd
print(f"Using:\tTarget site: {target_site}\n\tDonor site: {donor_site}")
print(f"\tTarget contig: {args.target_contig}")
print(f"\tDonor contig: {args.donor_contig}")



# Read sequences of target and donor files
print("Reading target contig...", end=" ")
target_contig_records = list(SeqIO.parse(args.target_contig, "genbank"))
assert len(target_contig_records) == 1, "Found more than one record in target contig file!"

for record in target_contig_records:
    target_contig = str(Seq(record.seq))
    print("done")

print("Reading donor contig...", end=" ")
donor_contig_records = list(SeqIO.parse(args.donor_contig, "genbank"))
assert len(donor_contig_records) == 1, "Found more than one record in donor contig file!"

for record in donor_contig_records:
    donor_contig = str(Seq(record.seq))
    print("done")

# NOTE: cell 2

def complement(seq):
    return str(Seq(seq).complement())


def find_matches(site, contig, complement_search=True, return_contig=True):
    """
    Find matches for site in contig. If none found, look for it in the complement
    of the contig instead if complement_search=True.
    """
    matches_in_contig = [match for match in re.finditer(site, contig)]
    if (len(matches_in_contig) > 1):
        raise ValueError("Found multiple sites in contig!")
    elif (len(matches_in_contig) == 0):
        if complement_search:
            print("\tsearching for site in complement...", end=" ")
            contig = complement(contig)
            matches_in_contig = [match for match in re.finditer(site, contig)]
            if (len(matches_in_contig) > 1):
                raise ValueError("Found multiple sites in contig!")
            elif (len(matches_in_contig) == 0):
                raise ValueError("Found 0 sites in contig!")
            else:
                print("found!")
        else:
            raise ValueError("Found 0 sites in contig!")

    site_start, site_end = matches_in_contig[0].start(), matches_in_contig[0].end()
    # contig[site_start:site_end] == site
    # NOTE: safe to index like [site_end:] or [:site_start] to extract contig 
    if return_contig:
        return contig, site_start, site_end
    else:
        return site_start, site_end


print("Finding target site in target contig")
target_contig, target_start, target_stop = find_matches(target_site, target_contig)
print("Finding donor site in donor contig")
donor_contig, donor_start, donor_stop = find_matches(donor_site, donor_contig)
# Get sequences of donor + cargo besides target / donor site
target_external_seq = target_contig[target_stop:] + target_contig[:target_start]
donor_cargo_seq = donor_contig[donor_stop:] + donor_contig[:donor_start]


# Calculate insertion product
if not same_contig:
    print("Calculating results in insertion mode...")
    result = (complement(lt) + tc + rd 
            + complement(donor_cargo_seq) 
            + complement(ld) + dc + rt + target_external_seq)
    print(result)
else:
    # TODO: use seq logic to determine if rearrangement or excision







# REARRANGEMENT AND EXCISION LOGIC!!!!!!!!!!!!!!!!!!! NOT COMPLETE!
#
#
#
#
#
#
#
#
#
# NOTE: only useful for excisions / rearrangements ? 
def find_partial_site(part, site, contig_pos):
    """_summary_

    Args:
        part (str): Part of site. e.g. LT, RT, LD, RD
        site (str): Full sequence of target or donor site
        contig_pos (str): Start index of site in contig
    """    
    part_start, part_end = find_matches(part, site, complement_search=False, return_contig=False)
    part_start, part_end = contig_pos + part_start, contig_pos + part_end
    return part_start, part_end

# Extract indices of LT, RT, LD, RD
lt_start, lt_end = find_partial_site(lt, target_site, target_start)
rt_start, rt_end = find_partial_site(rt, target_site, target_start)
ld_start, ld_end = find_partial_site(ld, donor_site, donor_start)
rd_start, rd_end = find_partial_site(rd, donor_site, donor_start)
# if same_contig:
#     ...

