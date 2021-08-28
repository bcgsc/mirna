import json
import os
import pytest

from mirna_profiling.annotation.ucsc import UCSC


def test_databae_setup():
    with pytest.raises(Exception):
        UCSC("hg38", user="dummy", password="dummy", host="ucsc")


def test_annotate_gene_features():
    gene_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", "gene.json"
    )
    with open(gene_file, "r") as jf:
        gene = json.load(jf)
    features = UCSC.annotate_gene_features(gene)
    assert features[0].get("feature_type") == "5 UTR"
    assert features[-1].get("feature_type") == "3 UTR"
    assert (
        len(list(filter(lambda x: x.get("feature_type") == "Coding Exon", features)))
        == 6
    )
    assert len(list(filter(lambda x: x.get("feature_type") == "Intron", features))) == 5

    gene_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mock_input",
        "gene_w_utr_intron.json",
    )
    with open(gene_file, "r") as jf:
        gene = json.load(jf)
    features = UCSC.annotate_gene_features(gene)

    assert features[0].get("feature_type") == "5 UTR"
    assert features[-1].get("feature_type") == "3 UTR"
    assert (
        len(list(filter(lambda x: x.get("feature_type") == "Coding Exon", features)))
        == 2
    )
    assert len(list(filter(lambda x: x.get("feature_type") == "Intron", features))) == 5
