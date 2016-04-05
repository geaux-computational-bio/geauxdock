#!/usr/bin/env python2

import re
import os
import glob
import pandas as pd


def readReport (fn):
#    print 'reading', fn
    lines = open (fn, 'r').readlines ()
    if len (lines) > 0:
        ligfile = [l.split()[-1] for l in lines if re.match('ligand file', l)][0]
        prtfile = [l.split()[-1] for l in lines if re.match('protein file', l)][0]
        n_lig = [int(l.split()[-1]) for l in lines if re.match('ligand conformations', l)][0]
        n_prt = [int(l.split()[-1]) for l in lines if re.match('prt conformations', l)][0]
        n_tmp = [int(l.split()[-1]) for l in lines if re.match('temperatures', l)][0]
        n_rep = [int(l.split()[-1]) for l in lines if re.match('replica ensembles', l)][0]
        size_lig = [int(l.split()[-1]) for l in lines if re.match('size_lig', l)][0]
        size_prt = [int(l.split()[-1]) for l in lines if re.match('size_prt', l)][0]
        size_pnk = [int(l.split()[-1]) for l in lines if re.match('size_pnk', l)][0]
        size_mcs = [int(l.split()[-1]) for l in lines if re.match('size_mcs', l)][0]
        kernel_time = [float(l.split()[-1]) for l in lines if re.match('monte carlo', l)][0]
        kernel_sweep= [float(l.split()[-1]) for l in lines if re.match('MC sweeps per second', l)][0]

        lig_id = os.path.basename(ligfile)
        prt_id = os.path.basename(prtfile)
        complex_id = lig_id + '-' + prt_id
#       print complex_id

        return complex_id, n_lig, n_prt, n_tmp, n_rep, size_lig, size_prt, size_pnk, size_mcs, kernel_time, kernel_sweep




def readReportDir(dir):
    columns = ["n_lig", "n_prt", "n_tmp", "n_rep", "size_lig", "size_prt", "size_pnk", "size_mcs", "kernel_time", "kernel_sweep"]
    fns = glob.glob (dir + "/*.txt")
    data = {}
    for fn in fns:
        result = readReport (fn)
        if result is not None:
            data[result[0]] = result[1:]

    dset = pd.DataFrame(data).T
    dset.columns = columns
    csv_file = './tmp/' + os.path.basename (dir) + '.csv'
    dset.to_csv (csv_file)





if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print 'Usage: ./readReport.py DIR'
    else:
        readReportDir (sys.argv[1])

