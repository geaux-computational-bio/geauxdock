#!/usr/bin/env python2

import glob
import os
import sys


def func1 (lig_dir):
    ffs = sorted (glob.glob (lig_dir + '/*.ff'))
    sdfs = sorted (glob.glob (lig_dir + '/*.sdf'))
    sdf_ffs = zip (sdfs, ffs)

    for sdf_ff in sdf_ffs:
        sdf_path = os.path.abspath (sdf_ff[0])
        ff_path = os.path.abspath (sdf_ff[1])
        lig_id = os.path.basename (sdf_path).split ('.')[0]
        str = lig_id + ' ' + sdf_path + " " + ff_path
        print str


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: ./mkliglist.py LIG_DIR_PATH'
    else:
        lig_dir = sys.argv[1]
        func1 (lig_dir)

