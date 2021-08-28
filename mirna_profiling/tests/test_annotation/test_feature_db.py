import os
import json
from mirna_profiling.annotation.feature_db import FeatureDatabase, load_csv


def test_get_overlapping_features():
    ucsc = FeatureDatabase()
    feature_csv = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mock_input",
        "ucsc_hg38_chr1_features.csv",
    )
    ucsc.load_features_from_csv(feature_csv)

    feature = {
        "chr": "chr1",
        "strand": "-",
        "end": 1115021,
        "feature_type": "Intron",
        "name": "C1orf159",
        "start": 1092104,
    }

    features = ucsc.get_overlapping_features_from_hierarchy(
        "chr1", 1111655, 1111670, "-"
    )
    assert len(features) == 1
    assert features[0].items() >= feature.items()

    features = ucsc.get_overlapping_features_from_hierarchy(
        "chr1", 1111655, 1111670, "-", find_all=True
    )

    assert len(features) == 4
