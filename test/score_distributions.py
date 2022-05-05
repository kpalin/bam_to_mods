#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 

Created at Thursday 05 May 2022  11:21 by Kimmo Palin <kpalin@helsinki.fi>
"""

__version__ = "0.1"

from argparse import Namespace
from email.policy import default
from os import EX_OK
import sys
from typing import List, Optional


def parse_args(args: List[str]) -> Namespace:
    import argparse

    description = "\n".join(__doc__.splitlines()[:-1]).strip()

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i", "--input", help="Input file [default:%(default)s]", default="/dev/stdin"
    )

    version = f"%(prog)s {__version__}"
    parser.add_argument("--version", action="version", version=version)

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        action="store_true",
        help="Be more verbose with output",
    )

    args = parser.parse_args()

    import logging

    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s:%(funcName)s:%(levelname)s:%(message)s",
        )

    return args


def main(argv: List[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)

    import pysam
    from collections import Counter, defaultdict
    import logging as log

    samfile = pysam.AlignmentFile(args.input, "r")

    distribs = dict()
    for n_reads, r in enumerate(samfile):
        if r.is_supplementary or r.is_secondary:
            continue
        ref_name = r.reference_name
        if ref_name != "chr1":
            break
        if ref_name not in distribs:
            distribs[ref_name] = Counter()
        q_to_ref_pos = r.get_reference_positions(full_length=True)

        for k, mC in r.modified_bases.items():
            ref_key = ref_name, k
            if ref_key not in distribs:
                distribs[ref_key] = Counter()

            d = distribs[ref_key]
            for q_pos, mod_prob in mC:
                ref_pos = q_to_ref_pos[q_pos]
                if ref_pos is not None:
                    d[(ref_pos, mod_prob)] += 1

    return distribs


if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt

    ret = main()
    r = pd.Series(r["chr1", ("C", 0, "m")]).reset_index()
    r.columns = ["Pos", "Score", "Count"]
    D = (
        r.query("Count>2")
        .Score.value_counts()
        .sort_index()
        .reset_index()
        .rename(columns={"index": "ProbInt"})
    )
    D["Prob"] = (D.ProbInt + 0.5) / 256.0
    D.plot.scatter(x="Prob", y="Score")
    plt.yscale("log")
    plt.gca().locator_params(axis="x", nbins=18)
    plt.axvline(0.1)
    plt.axvline(0.8)
# %run /home/kpalin/software/bam_to_mods/test/score_distributions.py -i /home/kpalin/software/bam_to_mods/subset.cram
