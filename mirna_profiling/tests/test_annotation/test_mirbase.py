import os
import json
import pytest

from mirna_profiling.annotation.mirbase import Mirbase

test_dir = os.path.dirname(os.path.realpath(__file__))
out_dir = os.path.join(test_dir, "mock_output")


def test_databae_setup():
    with pytest.raises(Exception):
        Mirbase("mirna_21", "hsa", user="dummy", password="dummy", host="mirbase.org")


def test_load_mirna_from_gff(mirbase_21_db):
    os.makedirs(out_dir, 0o775, exist_ok=True)
    mirnas = mirbase_21_db.load_mirnas_from_gff()
    mature_accs = mirbase_21_db.get_mature_mirna_acc()
    assert len(mature_accs) == 2588
    assert len(mirnas) == 1881
    features = mirbase_21_db.get_features()
    mirbase_21_db.dump_features(
        os.path.join(out_dir, "mirbase_mirna_21_hsa_features.csv")
    )
    assert len(features) == 9879


def test_annotate_mirna_features():
    mirna_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", "mirna_from_db.json"
    )
    with open(mirna_file, "r") as jf:
        mirna = json.load(jf)
    features = Mirbase.annotate_mirna_feature(Mirbase.transform_mirna_db_rec(mirna))

    assert [x.get("feature_type") for x in features] == [
        "precursor",
        "star",
        "unannotated",
        "stemloop",
        "mature",
        "precursor",
    ]
    assert mirna.get("start") == features[0].get("start")
    assert mirna.get("end") == features[-1].get("end")

    mirna_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mock_input",
        "mirna_mature_from_db.json",
    )
    with open(mirna_file, "r") as jf:
        mirna = json.load(jf)
    features = Mirbase.annotate_mirna_feature(Mirbase.transform_mirna_db_rec(mirna))

    assert [x.get("feature_type") for x in features] == [
        "precursor",
        "mature",
        "stemloop",
        "unannotated",
    ]
    assert mirna.get("start") == features[0].get("start")
    assert mirna.get("end") == features[-1].get("end")

    mirna_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", "mirna.json"
    )
    with open(mirna_file, "r") as jf:
        mirna = json.load(jf)
    features = Mirbase.annotate_mirna_feature(Mirbase.transform_mirbase_gff_rec(mirna))
    assert [x.get("feature_type") for x in features] == [
        "unannotated",
        "stemloop",
        "mature",
    ]


def test_get_mature_mirnas_by_id(mirbase_21_db):

    mature_id = "MIMAT0000227"
    mature_name = "hsa-miR-197-3p"
    mature_mirnas = mirbase_21_db.get_mature_mirnas_by_id(mature_id)

    assert (
        mature_mirnas
        and mature_mirnas[0].get("mature_acc") == mature_id
        and mature_mirnas[0].get("mature_name") == mature_name
    )
