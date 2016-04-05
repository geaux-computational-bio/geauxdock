#!/usr/bin/env python

import os
import argparse


def split(sdf_path, dirname=""):
    """split the molecules in one sdf into individual files
    Keyword Arguments:
    sdf_path -- file path
    dirname  -- output directory, the same directory as the input file by default
    """
    section = []
    if dirname == "":
        dirname = os.path.dirname(sdf_path)
    ofn = ""
    with open(sdf_path, 'r') as ifs:
        lines = ifs.readlines()
        for idx, line in enumerate(lines):
            section.append(line)
            if '>  <MOLID>' in line:
                molecule_name = lines[idx + 1].rstrip()
                ofn = os.path.join(dirname, molecule_name + '.sdf')
            if '$$$$' in line:
                with open(ofn, 'w') as ofs:
                    ofs.write("".join(section))
                section = []


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="split large sdf file into individual files")
    parser.add_argument(
        "-d",
        "--directory",
        default="",
        type=str,
        help="output directory, directory of the sdf file by default")
    parser.add_argument("-s", "--sdf", type=str, help="ligand sdf file")

    args = parser.parse_args()

    work_dir = args.directory
    lig_sdf = args.sdf
    split(lig_sdf, dirname=work_dir)
