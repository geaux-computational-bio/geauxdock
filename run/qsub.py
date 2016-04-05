#!/usr/bin/env python2

import os


def func1 (bin, nt, n):
    x = \
"""
#!/bin/sh
#PBS -q workq
#PBS -l nodes=1:ppn=16
#PBS -l walltime=3:00:00
#PBS -j oe

date
/home/yfang11/work/geauxdock_multibackend/bin/batch.py ~APP~ ../../astex/proteins/prt.txt ../../astex/ligands/ligs.txt ~NT~
date

exit 0
"""
    x = x.replace ('~APP~', bin)
    x = x.replace ('~NT~', nt)

    file = './tmp/' + bin + '_' + nt + '.pbs'
    f = open (file, 'w')
    f.write (x)
    f.close ()

    cmd = 'qsub' + ' ' + (file)
    print (cmd)
    #os.system (command)
    return



if __name__ == '__main__':
    nt_vars = [1, 2, 4, 8]
    nt_vars = map (str, nt_vars)          # to string
    dock_vars = ['dock-gpu-full', 'dock-gpu-kde', 'dock-gpu-mcs', 'dock-gpu-prt']

    n = 0
    for dock_var in dock_vars:
       for nt_var in nt_vars:
          n += 1
          func1 (dock_var, nt_var, n)


