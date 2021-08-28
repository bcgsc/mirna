#!/usr/bin/env python

import math

ADAPTER_CATEGORIES = [0, 14, 25, 35, math.inf]

ADAPTERS = dict(
    zip(
        ADAPTER_CATEGORIES,
        ["Adapter dimers"]
        + [
            "Adapter at {}-{}bp".format(
                ADAPTER_CATEGORIES[i - 1] + 1, ADAPTER_CATEGORIES[i]
            )
            for i in range(1, len(ADAPTER_CATEGORIES) - 1)
        ]
        + ["Adapter after {}bp".format(ADAPTER_CATEGORIES[-2])],
    )
)


def get_adapter_stat(adapter_file):
    """Get read count from adapter report file"""
    adapter_stat = {"% Adapter dimers": 0}
    adapter_stat.update({ADAPTERS[x]: 0 for x in ADAPTERS.keys()})
    total_reads = 0
    with open(adapter_file, "r") as f:
        for line in f:
            (pos, read_count) = [int(x) for x in line.rstrip().split(" ")]
            bin = get_bin(pos, ADAPTER_CATEGORIES)
            adapter_stat[ADAPTERS.get(bin)] += read_count
            total_reads += read_count

    if total_reads > 0:
        adapter_stat["% Adapter dimers"] = "{:0.2f}%".format(
            adapter_stat["Adapter dimers"] * 100 / total_reads
        )
    return adapter_stat


def get_bin(pos, bins):
    """Given a sorted list of bins/numbers, return the bin that pos belongs to"""
    for i in range(0, len(bins) - 1):
        if pos == bins[i]:
            return bins[i]
        if pos > bins[i] and pos <= bins[i + 1]:
            return bins[i + 1]

    return bins[-1]
