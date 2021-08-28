# requests library does not work with FTP schema, so use urllib.request instead
import urllib.request

import logging
import gzip
import os


def url_exists(url):
    website_exists = False
    try:
        with urllib.request.urlopen(url) as response:
            website_exists = True
    except Exception as e:
        logging.error("Error Retrieving URL {}: {}".format(url, e))

    return website_exists


def get_url_content(url, max_retry=3):

    data = None
    retries = 0
    status = False
    while not status and retries <= max_retry:
        try:
            with urllib.request.urlopen(url) as response:
                if url.endswith(".gz"):
                    data = gzip.decompress(response.read())
                else:
                    data = response.read()
                status = True
        except Exception as e:
            retries += 1
            logging.error(
                "Error Retrieving URL {}: {}. Retry {}".format(url, e, retries)
            )

    return data


def url_to_str(url):
    data = get_url_content(url)
    if data:
        data = data.decode("utf-8")

    return data


def url_to_list(url):
    data = get_url_content(url)
    if data:
        data = data.decode("utf-8").splitlines()

    return data


def download_url(url, outdir):
    if not url_exists(url):
        return None

    outfile = os.path.join(outdir, url.split("/")[-1])
    with open(outfile, "wb") as f:
        with urllib.request.urlopen(url) as response:
            f.write(response.read())

    return outfile
