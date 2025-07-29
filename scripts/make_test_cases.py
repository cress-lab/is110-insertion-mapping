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
import os
import random
from predict_recombination_products import bRNAReaction

def make_random_seq(n_bases):
    bases = ['A', 'T', 'C', 'G']
    random_sequence_list = [random.choice(bases) for _ in range(n_bases)]
    random_sequence_str = "".join(random_sequence_list)
    return random_sequence_str

def format_as_seqrecord(seq):
    record = SeqRecord(Seq(seq), id="inv_test_1", name="test1", )
    record.annotations["topology"] = "circular"
    record.annotations["molecule_type"] = "DNA"
    return record


def make_inversion_testcase(target_site, donor_site):
    left_of_target = make_random_seq(1000)
    right_of_target = make_random_seq(500)
    left_of_donor = make_random_seq(500)
    right_of_donor = make_random_seq(1000)

    test1 = (
        left_of_target
        + target_site
        + right_of_target
        + left_of_donor
        + reverse_complement(donor_site)
        + right_of_donor
    )
    test2 = (
        left_of_target
        + reverse_complement(target_site)
        + right_of_target
        + left_of_donor
        + donor_site
        + right_of_donor
    )
    test3 = (
        left_of_target
        + donor_site
        + right_of_target
        + left_of_donor
        + reverse_complement(target_site)
        + right_of_donor
    )
    test4 = (
        left_of_target
        + reverse_complement(donor_site)
        + right_of_target
        + left_of_donor
        + target_site
        + right_of_donor
    )
    return [format_as_seqrecord(i) for i in [test1, test2, test3, test4]]

def make_excision_testcase(target_site, donor_site):
    left_of_target = make_random_seq(1000)
    right_of_target = make_random_seq(500)
    left_of_donor = make_random_seq(500)
    right_of_donor = make_random_seq(1000)

    test1 = (
        left_of_target
        + target_site
        + right_of_target
        + left_of_donor
        + donor_site
        + right_of_donor
    )
    test2 = (
        left_of_target
        + reverse_complement(target_site)
        + right_of_target
        + left_of_donor
        + reverse_complement(donor_site)
        + right_of_donor
    )
    return [format_as_seqrecord(i) for i in [test1, test2,]]

"#66FF00"

def add_features(record):
    upstream_feature = SeqFeature(
            location=FeatureLocation(10, 300), 
            type="gene", 
            qualifiers = {
                "label": f"Upstream feature",
                "ApEinfo_fwdcolor": "#002DF7",
                "ApEinfo_revcolor": "#002DF7",
                "note": "color: " + "#002DF7"
            })
    span_feat = SeqFeature(
            location=FeatureLocation(950, 1050), 
            type="gene", 
            qualifiers = {
                "label": f"Feature to break",
                "ApEinfo_fwdcolor": "#F70000",
                "ApEinfo_revcolor": "#F70000",
                "note": "color: " + "#F70000"
            })
    inter_feat = SeqFeature(
            location=FeatureLocation(1100, 1200), 
            type="gene", 
            qualifiers = {
                "label": f"Feature to flip or excise",
                "ApEinfo_fwdcolor": "#66FF00",
                "ApEinfo_revcolor": "#66FF00",
                "note": "color: " + "#66FF00"
            })
    inter_feat2 = SeqFeature(
            location=FeatureLocation(1600, 1700, strand=-1), 
            type="gene", 
            qualifiers = {
                "label": f"Feature to flip or excise",
                "ApEinfo_fwdcolor": "#FF00EE",
                "ApEinfo_revcolor": "#FF00EE",
                "note": "color: " + "#FF00EE",
            })
    span_feat = SeqFeature(
            location=FeatureLocation(1915, 2115), 
            type="gene", 
            qualifiers = {
                "label": f"Feature to break downstream",
                "ApEinfo_fwdcolor": "#F70000",
                "ApEinfo_revcolor": "#F70000",
                "note": "color: " + "#F70000"
            })
    downstream_feature = SeqFeature(
            location=FeatureLocation(1000 + len(ts) + 1100, 1000 + len(ts) + 1300), 
            type="gene", 
            qualifiers = {
                "label": f"Downstream feature",
                "ApEinfo_fwdcolor": "#F6FF00",
                "ApEinfo_revcolor": "#F6FF00",
                "note": "color: " + "#F6FF00"
            })
    new_features = [
        upstream_feature,
        span_feat,
        inter_feat,
        inter_feat2,
        span_feat,
        downstream_feature
        ]
    record.features += new_features
    

ts = "TAGGCGTCTACTTA"
ds = "ACAGTATCTTGTAT"



tests = make_inversion_testcase(ts, ds)
for i, test in enumerate(tests):
    bRNAReaction.annotate_site(test, 1000, ts, 'TAGGCGT', 'ACTTA', inplace=True)
    bRNAReaction.annotate_site(test, 1000 + len(ts) + 1000, ds, 'ACAGTAT', 'TGTAT', mode="Donor", inplace=True)
    add_features(test)
    with open(f"../misc/invex/inv/test{i}.gbk", "w") as f:
        SeqIO.write(test, f, "genbank")

tests = make_excision_testcase(ts, ds)
for i, test in enumerate(tests):
    bRNAReaction.annotate_site(test, 1000, ts, 'TAGGCGT', 'ACTTA', inplace=True)
    bRNAReaction.annotate_site(test, 1000 + len(ts) + 1000, ds, 'ACAGTAT', 'TGTAT', mode="Donor", inplace=True)
    add_features(test)
    with open(f"../misc/invex/ex/test{i}.gbk", "w") as f:
        SeqIO.write(test, f, "genbank")
