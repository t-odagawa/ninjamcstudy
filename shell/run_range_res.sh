for i in `seq 1 999`;do
inputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/momch/momch_ecc5_${i}_pid.momch
outputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/mom_root_file/mom_file_${i}.root
bsub -q s ../build/src/Background/MomentumConsistency ${inputfile} ${outputfile}
done
