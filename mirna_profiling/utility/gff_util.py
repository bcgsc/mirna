import gzip
import logging

from mirna_profiling.utility.url_util import url_to_list

GFF_COLUMNS = [
    "chr",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def gff_to_records(gff, gff_columns=None):
    """
    Parses the content of gff into a list of features

    Args:
        gff (str): file path or URL

    Returns:
        records (list): list of feature dictionaries, with keys in GFF_COLUMNS
    """
    if gff_columns is None:
        gff_columns = GFF_COLUMNS

    records = []
    if "://" in gff:
        # get GFF entries from URL
        content = url_to_list(gff)
    elif gff.endswith(".gz"):
        with gzip.open(gff, "rb") as gf:
            content = gf.read().decode("utf-8").splitlines()
    else:
        with open(gff, "r") as gf:
            content = gf.read().splitlines()

    if not content:
        logging.warning("Failure reading GFF: {}".format(gff))
        return records

    for line in content:
        if line.startswith("#"):
            continue
        tokens = line.strip().split("\t")
        entry = {gff_columns[idx]: val for idx, val in enumerate(tokens[:8])}
        if len(tokens) > 8:
            for x in tokens[8].strip().rstrip(";").split(";"):
                if "=" in x:
                    [key, val] = x.split("=", maxsplit=1)
                else:
                    [key, val] = x.split(maxsplit=1)
                entry[key] = val.strip('"')
        records.append(entry)

    return records
