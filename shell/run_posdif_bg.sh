for i in `seq 1 999`;do
ntbmfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/trackmatch/ninja_mc_h2o_ninjamatch_${i}.root
inputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/momch/momch_ecc5_${i}_pid.momch
outputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/bg_posdif_dist_${i}.root
bsub -q s ../build/src/Background/MuonMisMatchCheck ${ntbmfile} ${inputfile} ${outputfile}
done
