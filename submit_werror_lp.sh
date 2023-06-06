#!/bin/bash

#UNCOMMENT THE LINE FOR THE MODEL AND CONDITION TO BE EXAMINED

#SLURM_JOB_NAME=DSR00_treatment_expanded_inv_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_expanded_inv_RUN_1_C33A_x100_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_expanded_sym_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_expanded_sym_RUN_1_C33A_x100_Rall_obs_sigma

#SLURM_JOB_NAME=ID_DSR00_treatment_expanded_sym_no_gen_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=ID_DSR00_treatment_expanded_sym_no_gen_RUN_1_C33A_x100_Rall_obs_sigma

#SLURM_JOB_NAME=DSRF00_RUN_1_C33A_controls_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_RUN_1_C33A_x0_01_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_RUN_1_C33A_x0_1_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_RUN_1_C33A_x1_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_RUN_1_C33A_x100_Rall_obs_sigma

#SLURM_JOB_NAME=DSRF00_treatment_expanded_inv_RUN_1_C33A_x100_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_expanded_inv_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_expanded_sym_RUN_1_C33A_x100_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF00_treatment_expanded_sym_RUN_1_C33A_x10_Rall_obs_sigma

#SLURM_JOB_NAME=DSRF01_treatment_expanded_inv_RUN_1_C33A_x100_Rall_obs_sigma
#SLURM_JOB_NAME=DSRF01_treatment_expanded_RUN_1_C33A_x100_Rall_obs_sigma

SLURM_JOB_NAME=DSR00_RUN_1_C33A_controls_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_RUN_1_C33A_x0_01_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_RUN_1_C33A_x0_1_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_RUN_1_C33A_x1_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSR00_treatment_RUN_1_C33A_x100_Rall_obs_sigma

#SLURM_JOB_NAME=DSRG00_RUN_1_C33A_controls_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_RUN_1_C33A_controls_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_treatment_RUN_1_C33A_x0_01_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_treatment_RUN_1_C33A_x0_1_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_treatment_RUN_1_C33A_x1_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_treatment_RUN_1_C33A_x10_Rall_obs_sigma
#SLURM_JOB_NAME=DSR02_treatment_RUN_1_C33A_x100_Rall_obs_sigma

dnam=C33A_control
#dnam=C33A_x0_01
#dnam=C33A_x0_1
#dnam=C33A_x1
#dnam=C33A_x10
#dnam=C33A_x100

data=(C33A_control_R1 C33A_control_R2 C33A_control_R3)
#data=(C33A_x0_01_R1 C33A_x0_01_R2 C33A_x0_01_R3)
#data=(C33A_x0_1_R1 C33A_x0_1_R2 C33A_x0_1_R3)
#data=(C33A_x1_R1 C33A_x1_R2 C33A_x1_R3)
#data=(C33A_x10_R1 C33A_x10_R2 C33A_x10_R3)
#data=(C33A_x100_R1 C33A_x100_R2 C33A_x100_R3)

suff=_werr

model=00
#model=G00
#model=00_treatment
#model=00_treatment_expanded_sym
#model=00_treatment_expanded_sym_no_gen
#model=00_treatment_expanded_inv
#model=F00
#model=F00_treatment
#model=F00_treatment_expanded_inv
#model=F00_treatment_expanded_sym
#model=F01_treatment_expanded_inv
#model=F01_treatment_expanded
#model=01
#model=02
#model=02_treatment
dim=2
attempt=1

vf=Verification

tp_arr_1=(1 2 3 4 5 6 7 8 9 10)
tp_arr_2=9
#tp_pred=10

tp=('tp_arr_1' 'tp_arr_2')
study=(Calibration Validation)

echo "------------------------------------------------------"
echo "        DSR${model}_RUN_${attempt}_${data}            "
echo "------------------------------------------------------"
echo " "

#for st in {1} #..2
#do
st=1

	echo "=========================="
	echo " Beginning ${study[st-1]} "
	echo "=========================="
	echo " "
	echo " -- Specifying time-points -- "
	prf=${study[st-1],}

	tmp=tp_arr_$st[@]
	tp_arr=(`echo "${!tmp[0]}"`)
	time_points=$(echo ${#tp_arr[@]})
	num_sets=$(echo ${#data[@]})
	echo "  - # time-points = ${time_points}"
	echo "  - # sets = ${num_sets}"
	echo "  - time-point indices = ${tp_arr[@]}"

	rm -rf Data/${prf}_*.txt
	echo ${time_points} >  src/${model}/model/num_tp.txt
	echo ${num_sets}    >> src/${model}/model/num_tp.txt 
	echo ${dnam}	    >> src/${model}/model/num_tp.txt

	for u in ${data[@]}
	do
		for i in ${tp_arr[@]}; 
		do  
			awk -v line=${i} 'NR==line' Data/${u}${suff}.txt >> Data/${prf}_${u}${suff}.txt 
		done
	done

	echo " -- ${study[st-1]} -- "
	cp src/tmcmc/sample src/${model}/.
	cd src/parent/
	./main_parent start ../../Data ${prf}_${data[0]}${suff}.txt ${prf}_${data[1]}${suff}.txt ${prf}_${data[2]}${suff}.txt ${time_points}
	cd ../${model}

	sleep 5
	mpirun -np 4 ./sample 
	sleep 5

	mkdir -p ../../${study[st-1]}
	mv *.txt *.out ../../${study[st-1]}
	echo "-- ${study[st-1]} Fitting finshed --"
	echo "-- Results saved in ${study[st-1]}"

	echo "-- Plot posteriors --"
	cd ../../Visualization/
	Rscript plot_samples.R "../${study[st-1]}/final.txt" "${model}"

	echo "-- Verification - Draw parameters from posterior distributions --"
	cd ../${study[st-1]}
	#Last three columns are useless for now, so delete them (two columns if sigma is from experimental data)
	awk '{NF-=2}1' final.txt > final_${prf}_${dnam}${suff}_${model}.txt
	cp final_${prf}_${dnam}${suff}_${model}.txt ../$vf/Parameter_Replicates/


	cd ../$vf/Parameter_Replicates/
	# Draw the replicates from the posteriors
	./hist ${dim} ${prf}_${dnam}${suff}_${model} ${model}

	cd ../
	./do_rep.sh 5000 ${time_points} ${prf}_${dnam}${suff} ${model}
	cd ../src/parent
    ./main_parent stop ../../Data ${prf}_${data[0]}${suff}.txt ${prf}_${data[1]}${suff}.txt ${prf}_${data[2]}${suff}.txt ${time_points}
    rm -rf ../inExp*.txt
	cd ../../${vf}

	echo "-- Parameter Optmization step --"
	mv ../src/${model}/output_rep_${prf}_${dnam}${suff}_vs_${model} .
	cd ../Optimization/
	julia --sysimage sys_complete.so optimize.jl "../Data/${prf}_${dnam}*.txt" "../${vf}/output_rep*/*.txt" "../${study[st-1]}/final_tmcmc_params.txt" "${dim}" "${model}"

	# Put the posterior predictive p-value here
	# for each rep:
	# for each time point: (sim_data - exp_data > 0)/reps
	#echo "-- Estimating Posterior Predictive p-values for 10 time-points--"
	#cd ../src/parent
    #./main_parent start ../../Data ${data}${suff}.txt 10
    #cd ../../${vf}/Post_p_value
	#./ppp 5000 10 ${prf} ${data}${suff} ${model}
	#cd ../../src/parent
	#./main_parent stop ../../Data ${data}${suff}.txt 10
	#rm -rf ../inExp*.txt


	echo "-- Finishing with Summary folder --"
	cd ../
	mkdir -p Summary/$SLURM_JOB_NAME
    cd Summary/$SLURM_JOB_NAME	
	mv ../../${study[st-1]}/ .
	cd ${study[st-1]}
	mv ../../../Optimization/*.png .
	mv ../../../Optimization/*.txt .
	#mv ../../../${vf}/Post_p_value/ppp_${prf}_${data}${suff}_${model}.txt .
	cd ../../../
	echo "-- Finished with ${study[st-1]} dataset --"

#done

#echo " "
#echo "==============================="
#echo " Statistical comparison of QOI "
#echo "==============================="
#echo " "

#This a 2-sample Kolmogorov-Smirnov test
# for values of cell count at t=72
# between calibration and validation 
# datasets. 

#cd QOI_comparisons/
#matlab -nodisplay -r "qoi_comp('${data}_werr','../${vf}','output_rep_calibration_${data}${suff}_vs_${model}','output_rep_validation_${data}${suff}_vs_${model}','output_rep_',5000); exit"
#mv plot_qoi_comp_output_rep_validation_${data}${suff}_vs_${model}.png ../Summary/$SLURM_JOB_NAME
#mv ks_qoi_output_rep_validation_${data}${suff}_vs_${model}.txt ../Summary/$SLURM_JOB_NAME
#cd ../
echo "-- Finished --"
#make clear
