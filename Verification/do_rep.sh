#! /bin/bash 

# number of parameter replicates are passed as $1
# time points are passed as $2
# data are passed as $3
# model is passed as $4

#cd ../src/parent
#./main_parent start ../../../Data $3.txt $2
##########################################################

cp serial_dorep.sh ../src/$4/
cd ../src/$4

mkdir -p output_rep_$3_vs_$4

parallel ./serial_dorep.sh ::: $(seq 1 $1) ::: $3 ::: $4 

#for ((i=1; i<=$1; i++))
#do
#		mkdir -p repdir.$i; 
#		cd model; cp -r * ../repdir.$i
#		cd ../repdir.$i	   	

#		awk -v mv=$i '{if(NR==mv) print $0}' ../../../Verification/Parameter_Replicates/rnd_$3_$4.txt > params.txt 
#		./main params.txt num_tp.txt  #> verification_$3_vs_$4_out.txt
#		mv output*.txt ../output_rep_$3_vs_$4/output_rep_$i.txt
#	   cd ../	
#done
rm -rf repdir.* params.txt





#######################################################33
#cd ../parent/
#./main_parent stop ../../../Data $3.txt $2
#cd ../

#mkdir output_rep_$3_vs_$4
#mv output_rep*.txt output_rep_$3_vs_$4/.

#echo "-> Output is written in: src/$4/output_rep_$3_vs_$4"
