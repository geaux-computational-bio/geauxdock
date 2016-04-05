./mkliglist.py ../data/astex/ligands/ > ../data/astex/ligands/list.txt
grep "5cppA1.sdf" ../data/astex/ligands/list.txt > ../data/astex/ligands/lig_5cppA1.txt
head -n 1 ../data/astex/ligands/list.txt > ../data/astex/ligands/ligs.txt
head -n 1 ../data/astex/ligands/list.txt >> ../data/astex/ligands/ligs.txt
head -n 1 ../data/astex/ligands/list.txt >> ../data/astex/ligands/ligs.txt



./mkprtlist.py ../data/astex/proteins/ > ../data/astex/proteins/list.txt
grep "1a07C.pdb" ../data/astex/proteins/list.txt > ../data/astex/proteins/prt.txt
grep "1a07C-1.pdb" ../data/astex/proteins/list.txt > ../data/astex/proteins/prt1.txt
grep "1a07C-3.pdb" ../data/astex/proteins/list.txt > ../data/astex/proteins/prt3.txt
grep "1a07C.pdb" ../data/astex/proteins/list.txt > ../data/astex/proteins/prt11.txt

