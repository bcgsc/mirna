import os
import json

from mirna_profiling.utility.gff_util import gff_to_records

test_dir = os.path.dirname(os.path.realpath(__file__))


def test_gff_to_records():
    ucsc_genome = "hg38"
    mirbase_version = "mirna_21"
    species_code = "hsa"

    gff = os.path.join(
        test_dir, "mock_input", "{}_{}.gff3".format(mirbase_version, species_code)
    )
    recs_gtf = gff_to_records(gff)

    assert len(recs_gtf) == 4694

    gff = os.path.join(
        test_dir, "mock_input", "{}_knownGene.gtf.gz".format(ucsc_genome)
    )
    recs_gtf = gff_to_records(gff)

    assert len(recs_gtf) == 1000
