#!/bin/bash
L1=$(grep -n 'Atoms' large_0.lmp | tr -dc '0-9')
L1=$((${L1} + 2))
L0=$(grep 'atoms' large_0.lmp | tr -dc '0-9')
L2=$((${L1}+${L0}-1))
echo $L0 $L1 $L2

awk -v l1=${L1} 'NR < l1 {print $0}' large_0.lmp > 0in.crystal
awk -v l1=${L1} -v l2=${L2} 'NR>=l1 && NR<=l2 {print $1, int(($1+3)/4), $3, $4, $5, $6}' large_0.lmp > 1in.crystal
#paste 1in.crystal newlattice.xyz > temp
#awk '{print $1, $2, $3, $8, $9, $10}' temp > temp2
awk  -v l2=${L2} 'NR>l2 {print $0}' large_0.lmp > 2in.crystal
cat 0in.crystal 1in.crystal 2in.crystal > large-s.lmp
rm -f 0in.crystal 1in.crystal 2in.crystal
