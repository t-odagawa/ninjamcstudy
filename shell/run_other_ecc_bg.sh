for i in `seq 1 999`;do
#inputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/momch/momch_ecc5_${i}_pid.momch
#outputfile=/hsm/nu/ninja/pra_tmp/mc_tmp_20220505/output/bg_dist_${i}.root
inputfile=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/ecc/momch/momch_ecc5_${i}_sig.momch
outputfile=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/ecc/output/bg_dist_${i}.root
bsub -q s ../build/src/Background/OtherEccBG ${inputfile} ${outputfile}
done
