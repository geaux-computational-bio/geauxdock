#!/usr/bin/env python2

import glob
import os
import sys


def func1 (prt_dir):
    prts = sorted (glob.glob (prt_dir + '/*.pdb'))

    for prt in prts:
        prt_path = os.path.abspath (prt)
        prt_id = os.path.basename (prt_path).split ('.')[0]
        str = prt_id + ' ' + prt_path
        print str


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: ./mkprtlist.py PRT_DIR_PATH'
    else:
        prt_dir = sys.argv[1]
        func1 (prt_dir)


