for i in `seq 1 999`;do
inputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_anu_20220714/momch/momch_ecc5_${i}_pid.momch
outputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_anu_20220714/output/bg_dist_${i}.root
bsub -q s ../build/src/Background/AnuIntBG ${inputfile} ${outputfile}
done
