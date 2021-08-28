import yaml
import pysam
import logging
import os

import pandas


class FeatureDatabase(object):
    """
    A class for loading, indexing and querying features of interest

    """

    feature_attrs = ["name", "feature_type", "chr", "strand", "start", "end"]
    feature_indices = ["chr", "strand", "feature_type"]

    def __init__(self, config_file=None):
        """
        Load config file and initialize object

        # TODO: Validate annotation.config
        """
        if not config_file:
            config_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                "configuration",
                "annotation_config.yaml",
            )

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        self.features = []
        self.feature_df = None
        self.feature_dict = None
        self.feature_priority = self.config["Feature Annotation"]["feature_priority"]
        self.feature_type_display_label = self.config["Feature Annotation"].get(
            "feature_type_display_label"
        )
        self.mito_chrs = self.config["Feature Annotation"]["mito_chrs"]
        self.bin_size = self.config["Feature Annotation"]["bin_size"]
        self.minimum_overlap = self.config["Feature Annotation"]["minimum_overlap"]

    def get_features(self):
        return self.features

    def get_feature_df(self, refresh=False):
        if self.feature_df is None or refresh:
            self.feature_df = pandas.DataFrame(self.get_features())
            self.feature_df.set_index(
                FeatureDatabase.feature_indices, drop=False, inplace=True
            )

        return self.feature_df

    def get_feature_dict(self, refresh=False):

        if self.feature_dict is None or refresh:
            self.feature_dict = self.build_feature_hierarchy()

        return self.feature_dict

    def add_features(self, features):

        mfeatures = self.get_features()
        mfeatures.extend(filter(FeatureDatabase.validate_feature, features))
        self.features = mfeatures

    def get_overlapping_features(
        self, chr, chr_start, chr_end, strand, filter_func=None
    ):
        """
        Get feature(s) overlapping the specified chromosome region
        TODO: Find an efficient way to slice df using feature_type, chr, strand and a range of bins

        Args:
            chr (str):
            chr_start (int):
            chr_end (int):
            strand (str):
            filter_func (func): function used to filter features
        Returns:
            (list): feature(s) overlapping the specified region
        """

        df = self.get_feature_df().xs(chr + strand, drop_level=False)
        features = df[(chr_start <= df["end"]) & (chr_end >= df["start"])].to_dict(
            "records"
        )
        if filter_func is not None:
            features = list(filter(filter_func, features))
        return features

    def annotate_region(self, chr, chr_start, chr_end, strand):
        """
        Get the highest priority feature overlapping the specified chromosome region

        Args:
            chr (str):
            chr_start (int):
            chr_end (int):
            strand (str):

        Returns:
            (dict): feature overlapping the specified region
        """

        chr_start += self.minimum_overlap
        chr_end -= self.minimum_overlap

        features = self.get_overlapping_features_from_hierarchy(
            chr, chr_start, chr_end, strand
        )

        annot = ("Unknown", "", "")
        if features:
            # features = sorted(
            #     features,
            #     key=lambda x: self.feature_priority.index(x.get("feature_type"))
            #     if x.get("feature_type") in self.feature_priority
            #     else 1000000,
            # )
            feature = features[0]
            feature_type = feature.get("feature_type")
            if feature_type in [
                "precursor",
                "stemloop",
                "mature",
                "star",
                "unannotated",
            ]:
                if feature_type in ["star", "mature"]:
                    (mirna_id, mature_acc) = features[0].get("name").split(",")
                    mirna_type = "{},{}".format(feature_type, mature_acc)
                else:
                    mirna_id = features[0].get("name")
                    mirna_type = feature_type
                annot = ("miRNA", mirna_id, mirna_type)
            else:
                annot = (feature_type, features[0].get("name"), "")

        return annot

    def get_overlapping_features_from_hierarchy(
        self, chr, chr_start, chr_end, strand, filter_func=None, find_all=False
    ):
        """
        Get feature(s) overlapping the specified chromosome region

        Args:
            chr (str):
            chr_start (int):
            chr_end (int):
            strand (str):
            filter_func (func): function used to filter features
        Returns:
            (list): feature(s) overlapping the specified region
        """

        all_feature_dict = self.get_feature_dict()

        overlapping_features = []
        start_bin = int(chr_start / self.bin_size)
        end_bin = int(chr_end / self.bin_size)
        chr_strand = chr + strand
        for feature_type in self.feature_priority:
            feature_dict = all_feature_dict.get(feature_type)
            if feature_dict and feature_dict.get(chr_strand):
                feature_chr_dict = feature_dict.get(chr_strand)
                for bin_idx in range(start_bin, end_bin + 1):
                    features = feature_chr_dict.get(bin_idx)
                    if features:
                        for key, feat in features.items():
                            if (
                                chr_start <= feat.get("end")
                                and chr_end >= feat.get("start")
                                and (not filter_func or filter_func(feat))
                            ):
                                overlapping_features.append(feat)
                                if not find_all:
                                    return overlapping_features

        overlapping_features = unique_dicts(overlapping_features)

        return overlapping_features

    def get_alternate_chrs(self, chr):
        alt_chrs = []
        if chr in self.mito_chrs:
            alt_chrs.extend(self.mito_chrs)
            alt_chrs.remove(chr)
        elif chr.startswith("chr"):
            alt_chrs.append(chr[3:])
        else:
            alt_chrs.append("chr{}".format(chr))

        return alt_chrs

    def conform_features_chr_to_bam(self, inbam):
        samfile = pysam.AlignmentFile(inbam, "rb")
        chrs = samfile.references
        samfile.close()

        with_prefix = True if any(x.startswith("chr") for x in chrs) else False
        mito = "chrMT" if with_prefix else "MT"
        use_mt = True if mito in chrs else False
        for feature in self.get_features():
            self.reformat_chr(feature, with_prefix, use_mt)

    def reformat_chr(self, feature, with_prefix=True, use_mt=False):
        """
        Reformat chr name to confirm with UCSC

        Args:
            with_prefix (bool): chromosome name starts with "chr"
            use_mt (bool): MT for the mitochondrial

        """

        if feature["chr"] in self.mito_chrs:
            feature["chr"] = ("chr" if with_prefix else "") + ("MT" if use_mt else "M")
        elif with_prefix and not feature["chr"].startswith("chr"):
            feature["chr"] = "chr" + feature["chr"]
        elif not with_prefix and feature["chr"].startswith("chr"):
            feature["chr"] = feature["chr"][3:]

    def build_feature_hierarchy(self, features=None):
        """
        Build a feature hash ordered by chromosome, strand and bin for fast access

        """
        all_feature_dict = {}
        feature_ids = {}
        idx = 1
        if features is None:
            features = self.get_features()
        for feat in features:
            feature_type = feat.get("feature_type")
            if not all_feature_dict.get(feature_type):
                all_feature_dict[feature_type] = {}
            feature_dict = all_feature_dict.get(feature_type)
            chr_strand = feat.get("chr") + feat.get("strand")
            if not feature_dict.get(chr_strand):
                feature_dict[chr_strand] = {}
                for alt_chr in self.get_alternate_chrs(feat.get("chr")):
                    feature_dict[alt_chr + feat.get("strand")] = feature_dict[
                        chr_strand
                    ]

            start_bin = int(feat.get("start") / self.bin_size)
            end_bin = int(feat.get("end") / self.bin_size)
            refname = feat.get("name")
            if feature_ids.get(refname):
                refname += "_{}".format(idx)
            else:
                feature_ids[refname] = 1
            for i in range(start_bin, end_bin + 1):
                if not feature_dict[chr_strand].get(i):
                    feature_dict[chr_strand][i] = {}
                feature_dict[chr_strand][i][refname] = feat
            idx += 1

        return all_feature_dict

    def update_display_label(self, features):
        for x in features:
            if (
                self.feature_type_display_label
                and x.get("feature_type") in self.feature_type_display_label
            ):
                x["feature_type"] = self.feature_type_display_label.get(
                    x["feature_type"]
                )

    def load_features_from_csv(self, csv):
        """
        Load features from (custom) csv file

        Args:
            csv (str):

        Returns:
            (list): list of dictionaries
        """
        self.features = load_csv(csv, filter_func=FeatureDatabase.validate_feature)
        return self.features

    def save_features(self):
        if hasattr(self, "feature_csv"):
            self.dump_features(self.feature_csv)

    def dump_features(self, csv, sep="\t"):

        pandas.DataFrame(self.get_features()).to_csv(csv, sep=sep, index=False)

    @staticmethod
    def validate_feature(feat):
        """
        Validate a feature
        Args:
            feat (dict):

        Returns:
            is_valid (bool):
        """
        is_valid = True
        if any(x not in feat for x in FeatureDatabase.feature_attrs) or feat.get(
            "strand"
        ) not in ["+", "-"]:
            is_valid = False

        return is_valid


def load_csv(csv, required_attrs=None, filter_func=None):
    """
    Load csv file with header, and return a list of records

    Args:
        csv (str):
        required_attrs (list:
        filter_func (func):

    Returns:
        features(list): list of dictionaries
    """
    with open(csv, "r") as f:
        first_line = f.readline().strip()
    sep = "\t" if "\t" in first_line else ","
    df = pandas.read_csv(csv, header=0, sep=sep)
    if required_attrs:
        missing_attrs = set(required_attrs) - set(df.columns)
        if missing_attrs:
            logging.error(
                "Invalid feature CSV file: {}. Missing required fields: {}".format(
                    csv, missing_attrs
                )
            )
            return []

    features = df.to_dict("records")
    if filter_func is not None:
        features = list(filter(filter_func, features))

    return features


def unique_dicts(mlist):

    return [dict(t) for t in {tuple(sorted(d.items())) for d in mlist}]
