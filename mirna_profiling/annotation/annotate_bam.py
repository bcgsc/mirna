#!/usr/bin/env python

import pysam
import argparse
import logging
import os
from pathlib import Path
import re

from mirna_profiling.annotation.ucsc import UCSC
from mirna_profiling.annotation.mirbase import Mirbase
from mirna_profiling.annotation.feature_db import FeatureDatabase


def main():
    args = parse_args()

    feature_db = UCSC(args.ucsc, config_file=args.config_file)
    mirbase = Mirbase(args.mirbase, args.species, config_file=args.config_file)
    feature_db.add_features(mirbase.get_features())
    if os.path.isdir(args.bam):
        bams = (str(x) for x in Path(args.bam).rglob("*.[bs]am"))
    else:
        bams = (args.bam,)

    for bam in bams:
        if args.out_dir:
            out_bam = os.path.join(
                args.out_dir, "{}.annot".format(os.path.basename(bam))
            )
        else:
            out_bam = "{}.annot".format(bam)
        annotate_sam(feature_db, bam, out_bam, config_file=args.config_file)


def annotate_sam(feature_db, inbam, outbam, config_file=None):
    """
    Annotate bam file using feature info retrieved from UCSC and miRBase

    Args:
        feature_dbs (FeatureDatabase) : an instance of FeatureDatabase class
        inbam (str): input bam path
        outbam (str): annotated bam path

    """
    logging.info("Start annotating bam: {}".format(inbam))
    samfile = pysam.AlignmentFile(inbam, "rb")
    pairedreads = pysam.AlignmentFile(outbam, "w", template=samfile)
    progress = 0
    for read in samfile.fetch(until_eof=True):
        annotate_read(feature_db, read)
        pairedreads.write(read)
        progress += 1
        if progress % 100000 == 0:
            logging.info(f"Annotated {progress:,} reads")

    pairedreads.close()
    samfile.close()
    logging.info("Finished Annotation. Output bam: {}".format(outbam))


def is_annotated(inbam):
    """
    Check if a bam file has been annotated using this package, i.e., XC tags are present.

    Args:
        inbam (string): input bam

    Returns:
        annot_stat (bool): True is bam is annotated, False otherwise
    """

    annot_stat = False
    samfile = pysam.AlignmentFile(inbam, "rb")
    for read in samfile.fetch(until_eof=True):
        if read.cigarstring is None or read.reference_name == "*":
            continue
        if read.has_tag("XC"):
            annot_stat = True
        break

    return annot_stat


def annotate_read(feature_db, read):
    """
    Annotate an aligned read, i.e., save aligned feature info in XC, XI and XD tags

    Args:
        feature_db (FeatureDatabase):
        read (pysam.AlignedSegment):

    """
    if read.cigarstring is None or read.reference_name == "*":
        return

    strand = "+" if (read.flag & 0x0010) == 0 else "-"
    # pysam coordinate is 0-based, different from the fields in the alignment file
    ref_start = read.reference_start + 1
    (xc, xi, xd) = feature_db.annotate_region(
        read.reference_name,
        ref_start,
        ref_start + parse_cigar_string(read.cigarstring) - 1,
        strand,
    )
    xc += ";"
    xi += ";"
    xd += ";"
    if read.has_tag("XA"):
        (tag, tag_type) = read.get_tag("XA", with_value_type=True)
        if tag_type == "Z":
            for x in tag.rstrip().rstrip(";").split(";"):
                (xchr, strand_start, xcigar) = x.split(",")[0:3]
                if not xcigar:
                    continue
                xstart = int(strand_start[1:])
                (xxc, xxi, xxd) = feature_db.annotate_region(
                    xchr,
                    xstart,
                    xstart + parse_cigar_string(xcigar) - 1,
                    strand_start[0],
                )
                xc += "{};".format(xxc)
                xi += "{};".format(xxi)
                xd += "{};".format(xxd)

    read.set_tag("XC", xc)
    read.set_tag("XI", xi)
    read.set_tag("XD", xd)


def parse_cigar_string(cigar):
    """
    Parse the number of aligned bases (MNDP) from cigar string

    Args:
        cigar (str):

    Returns:
        bases (int): number of aligned bases

    """
    bases = 0
    if cigar is not None:
        matches = re.findall("(\\d+)[MNDP]", cigar)
        bases = sum(int(x) for x in matches)

    return bases


def parse_args():
    parser = argparse.ArgumentParser(description="micro RNA profiling analyis")
    parser.add_argument(
        "-u", "--ucsc", help="UCSC database name, e,g. hg38, hg19, mm10", required=True
    )
    parser.add_argument(
        "-m",
        "--mirbase",
        help="miRBase MySQL database name, e.g. mirna_21, mirna_16",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--species",
        help="miRBase species code, e.g., hsa for human, mmu for mouse",
        required=True,
    )
    group_annot = parser.add_mutually_exclusive_group(required=True)
    group_annot.add_argument(
        "-p",
        "--proj_dir",
        dest="bam",
        help="Project directory path with annotated b/sam files",
    )
    group_annot.add_argument(
        "-b",
        "--bam",
        dest="bam",
        help="bam or sam file to annotate",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        help="Output directory",
    )
    parser.add_argument(
        "-c",
        "--config_file",
        help="Custom annotation config file in yaml format",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
