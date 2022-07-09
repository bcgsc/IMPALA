#!/usr/bin/env python
"""sampletsvtoyaml.py: Converts tab-separated sample information to configuration file"""

import yaml
import argparse
import os

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Converts a tab-separated file of sample information to a samples.yaml file")

parser.add_argument("-t", "--tsv_file", help = "Path to input tsv file", default = "./samples.tsv")
parser.add_argument("-o", "--output_path", help = "Path to config directory. Defaults to ./config/", default = "./config/")
parser = parser.parse_args()

if not os.path.exists(parser.output_path):
    os.mkdir(parser.output_path)

with open(parser.tsv_file) as f:
    f = [l.strip().split("\t") for l in f.readlines()]

def extract(lst, n):
    return [item[n] for item in lst]

samples = list(set(extract(f, 0)))
out_dict = {"samples":{}}

for s in samples:
    must_have = [s]
    must_have_set = set(must_have)
    sub = [x for x in f if set(x) & must_have_set == must_have_set]
    lst = []
    for n in sub:
        identity = n[1]
        path = n[2]
        element = {identity:path}
        lst.append(element)
        out_dict["samples"][s] = lst

with open(os.path.join(parser.output_path, "samples.yaml"), "w") as f:
    yaml.dump(out_dict, f, default_flow_style=False)