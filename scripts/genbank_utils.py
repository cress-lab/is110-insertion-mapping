"""
Script to write utils to modify genbank records
"""
import argparse
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from find_brna import find_target_donor_sites

import sys
from Bio.SeqFeature import FeatureLocation
import copy


# REVERSE COMPLEMENT UTILS
def flip_direction_note(qualifiers):
    """
    flip the direction from left to right or vice versa in a misc_feature
    in a snapgene genbank file
    """
    if "note" not in qualifiers:
        return qualifiers

    new_notes = []
    for note in qualifiers["note"]:
        if "direction: RIGHT" in note:
            new_note = note.replace("direction: RIGHT", "direction: LEFT")
        elif "direction: LEFT" in note:
            new_note = note.replace("direction: LEFT", "direction: RIGHT")
        else:
            new_note = note  # leave 'BOTH' or others unchanged
        new_notes.append(new_note)

    qualifiers["note"] = new_notes
    return qualifiers


def reverse_complement_feature(feat, rc_locations):
    """Calculates the reverse complemented feature for a given feature from 
    a snapgene genbank file

    Args:
        feat (Bio.SeqFeature.SeqFeature): Feature to edit
        rc_locations (np.array): numpy array of the equivalent minus strand
        positions for all positions in the original record

    Returns:
        Bio.SeqFeature.SeqFeature: Updated feature
    """    
    feature = copy.deepcopy(feat)
    start, end, strand = (
        feature.location.start, 
        feature.location.end, 
        feature.location.strand,
        )
    rc_end, rc_start = rc_locations[start], rc_locations[end]
    if feature.type != "source":
        strand = -1 * strand if strand is not None else None
    new_location = FeatureLocation(
        int(min(rc_start, rc_end)),
        int(max(rc_start, rc_end)),
        strand=int(strand) if strand is not None else None
    )
    feature.location = new_location
    # # Update direction for certain features
    # if feature.type in ["misc_feature", "gene", "ncRNA"]:
    feature.qualifiers = flip_direction_note(feature.qualifiers)
    return feature


def reverse_complement_record(record):
    """Given a genbank record, reverse complement the entire record + all 
    features
    """    
    new_record = copy.deepcopy(record)
    rc_locations = np.arange(0, len(new_record.seq)+1)[::-1]
    rc_feats = []
    for feature in new_record.features:
        rc_feature = reverse_complement_feature(feature, rc_locations)
        rc_feats.append(rc_feature)
    new_record.seq = reverse_complement(new_record.seq)
    new_record.features = rc_feats
    return new_record


# INSERTION UTILS


# OTHER
def orient_record_around_pos(record, new_start):
    seq_len = len(record.seq)
    new_seq = record.seq[new_start:] + record.seq[:new_start]

    new_features = []
    for feature in record.features:
        old_start = int(feature.location.start)
        old_end = int(feature.location.end)
        strand = feature.location.strand

        # Shift both coordinates
        new_start_pos = (old_start - new_start) % seq_len
        new_end_pos = (old_end - new_start) % seq_len
        # If new end position is 0, change it to end at the last position of the sequence
        if new_end_pos == 0:
            new_end_pos = seq_len
        # Handle wrap-around
        if new_start_pos > new_end_pos:
            # Feature crosses the origin now
            # new_location = FeatureLocation(new_start_pos, new_end_pos + seq_len, strand=strand)
            part1 = FeatureLocation(new_start_pos, seq_len, strand=strand)
            part2 = FeatureLocation(0, new_end_pos, strand=strand)
            # print("ORIENTING: ", feature.qualifiers.get("label", ""), part1, part2)
            new_location = CompoundLocation([part1, part2])
        else:
            new_location = FeatureLocation(new_start_pos, new_end_pos, strand=strand)

        # Build new feature
        new_feature = copy.deepcopy(feature)
        new_feature.location = new_location
        new_features.append(new_feature)

    # Rebuild record
    new_record = copy.deepcopy(record)
    new_record.seq = new_seq
    new_record.features = new_features
    new_record.annotations['comment'] = (
        f"Sequence rotated to start at position {new_start} (0-based offset)."
    )
    return new_record


target_path = "/groups/cress/projects/efactors/jack/insertion_mapping/misc/gbtest/NATIVE_Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk"
out_path = "/groups/cress/projects/efactors/jack/insertion_mapping/misc/gbtest/RC_Manually_Adjusted_Seq_of_HalfKan_Reporter_Strain_SLMQRT_4_dnaSES0085_no_BAC.gbk"
# donor_path = "/groups/cress/projects/efactors/jack/insertion_mapping/misc/W9G3HH_30kbp BAC_gSES0133.gbk"

record = SeqIO.read(target_path, "genbank")
record = reverse_complement_record(record)
# donor_record = SeqIO.read(donor_path, "genbank")

# print()
# print("LENGTH: ", len(record.seq))




with open(out_path, "w") as out:
    SeqIO.write(record, out, "genbank")


