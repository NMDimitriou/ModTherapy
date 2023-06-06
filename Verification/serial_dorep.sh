#! /bin/bash

pwd
mkdir -p repdir.$1;
cd model; cp -r * ../repdir.$1
cd ../repdir.$1

awk -v mv=$1 '{if(NR==mv) print $0}' ../../../Verification/Parameter_Replicates/rnd_$2_$3.txt > params.txt
./main params.txt num_tp.txt 
mv output*.txt ../output_rep_$2_vs_$3/output_rep_$1.txt
cd ../
