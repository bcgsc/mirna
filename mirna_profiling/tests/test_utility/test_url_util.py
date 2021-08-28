import os
from mirna_profiling.utility.url_util import (
    url_exists,
    url_to_list,
    url_to_str,
    get_url_content,
    download_url,
)

gzip_url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz"
ftp_url = "ftp://mirbase.org/pub/mirbase/21/genomes/hsa.gff3"
dummy_url = "http://dummy"

test_dir = os.path.dirname(os.path.realpath(__file__))
out_dir = os.path.join(test_dir, "mock_output")

VALDI_URL_SIZE = 1000


def test_url_exists():

    assert not url_exists(dummy_url)

    assert url_exists(ftp_url)


def test_get_url_content():

    # TEST valid URL
    content = url_to_list(ftp_url)
    assert len(content) > VALDI_URL_SIZE

    # Test invalid URL
    assert not get_url_content(dummy_url)


def test_get_gzip_url_content():

    # TEST gzip URL
    content = url_to_str(gzip_url)
    assert len(content) > VALDI_URL_SIZE


def test_download_url():
    os.makedirs(out_dir, 0o775, exist_ok=True)
    outfile = download_url(gzip_url, out_dir)
    assert outfile and os.path.exists(outfile)

    outfile = download_url(dummy_url, out_dir)
    assert not outfile
