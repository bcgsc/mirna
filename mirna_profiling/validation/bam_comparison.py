#!/usr/bin/env python

import pysam
import json
import os
import argparse


def main():
    """
    A script to compare tags in two bams
    """
    args = parse_args()
    compare_sam_tag(args.old_bam, args.new_bam)


def compare_sam_tag(inbam, outbam):
    """
    Compare read annotations (tags: XC, XI and XD) in two bams

    Args:
        inbam:
        outbam:

    Returns:

    """
    # Assume same number of colon-separated items in annotation_tags
    annotation_tags = ["XC", "XI", "XD"]
    same_count = 0
    diff_count = 0
    libname = os.path.basename(inbam)[0:6]
    samfile = pysam.AlignmentFile(inbam, "r")
    samfile2 = pysam.AlignmentFile(outbam, "r")
    with open("{}_diff_tags.txt".format(libname), "w") as of:
        for (read1, read2) in zip(
            samfile.fetch(until_eof=True), samfile2.fetch(until_eof=True)
        ):
            if any(
                any(not x.has_tag(y) for y in annotation_tags) for x in [read1, read2]
            ):
                continue

            x1s = [
                dict(zip(annotation_tags, x))
                for x in zip(*[read1.get_tag(y).split(";") for y in annotation_tags])
            ]
            x2s = [
                dict(zip(annotation_tags, x))
                for x in zip(*[read2.get_tag(y).split(";") for y in annotation_tags])
            ]
            for x in zip(x1s, x2s):
                if x[0].get("XC") == x[1].get("XC") and x[0].get("XD") == x[1].get(
                    "XD"
                ):
                    same_count += 1
                elif "UTR" in x[0].get("XC"):
                    same_count += 1
                elif (
                    x[0].get("XC") == "Coding Exon"
                    and x[0].get("XD") == "No CDS"
                    and x[1].get("XC") == "snoRNA"
                ):
                    same_count += 1
                elif x[0] == "Unknown" and x[1] == "rmsk Low complexity":
                    same_count += 1
                else:
                    diff_count += 1
                    of.write("{}\t{}\t{}\n".format(read1.query_name, x[0], x[1]))

    samfile2.close()
    samfile.close()
    print("Same : {} Diff: {}".format(same_count, diff_count))

    return diff_count == 0


def parse_args():
    parser = argparse.ArgumentParser(description="Compare Read Annotations in Bams")
    parser.add_argument(
        "--old_bam", help="Sam annotated with old module", required=True
    )
    parser.add_argument(
        "--new_bam", help="Bam/Sam annotated with new module", required=True
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
