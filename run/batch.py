#!/usr/bin/env python3

import os
import sys
import socket
import re
#import subprocess32    # a more robust replacement of  os.system ()

BIN = ''
OUTPUT_DIR = ''
OUTPUT_DIR_BASE = './out'
NT_LOW = 1
NT_HI = 1
NT_STEP = 1
ITER_LOW = 0
ITER_HI = 0


def Run2 (cmd, dir):
  for i in range (ITER_LOW, ITER_HI):
      ii = str (i).zfill (2)
      file =  dir + '/' + ii
      cmd2 = cmd
      cmd2 += ' | tee ' + file + '.txt'
      cmd3 = 'grep "\<Benchmark\>" ' + file + '.txt > ' + file + '_bench.csv'
      cmd4 = 'grep "\<Benchmark_papi\>" ' + file + '.txt > ' + file + '_benchpapi.csv'

      #print(cmd2)
      #print(cmd3)
      #print(cmd4)
      os.system (cmd2)
      os.system (cmd3)
      os.system (cmd4)



def Run1 (bin, output_dir2):
    for nt in range (NT_LOW, NT_HI + 1, NT_STEP):
        cmd = bin
        cmd = cmd + ' -nt ' + str (nt)

        cmd = cmd + ' -ll /worka/work/yfang11/geauxdock_cs_v2/data/astex/ligands/list.txt'
#        cmd = cmd + ' -ll /worka/work/yfang11/geauxdock_cs_v2/data/astex/ligands/lig_5cppA1.txt'

        cmd = cmd + ' -lp /worka/work/yfang11/geauxdock_cs_v2/data/astex/proteins/prt1.txt'
#       cmd = cmd + ' -lp /worka/work/yfang11/geauxdock_cs_v2/data/astex/proteins/prt11.txt'

        dir = output_dir2 + '_nt' + str (nt)
        #print(dir)
        os.mkdir (dir)

        Run2 (cmd, dir)



if __name__ == '__main__':
    bin = sys.argv[1]
    NT_LOW = sys.argv[2]
    NT_HI = sys.argv[3]
    NT_STEP = sys.argv[4]
    ITER_LOW = sys.argv[5]
    ITER_HI = sys.argv[6]
    NT_LOW = int (NT_LOW)
    NT_HI = int (NT_HI)
    NT_STEP = int (NT_STEP)
    ITER_LOW = int (ITER_LOW)
    ITER_HI = int (ITER_HI)



    hostname = socket.gethostname ()
    hostname = re.sub ("\d+$", "", hostname)
    binname = os.path.basename (bin)
    #output_dir2 = os.path.join (OUTPUT_DIR_BASE, hostname + '_' + binname)
    output_dir2 = OUTPUT_DIR_BASE + '_' + hostname + '_' + binname

    Run1 (bin, output_dir2)


