for i in `seq 30001 39999`;do
inputfile=/hsm/nu/ninja/pra_tmp/wall_mc_20220616/momch/momch_ecc5_${i}_pid.momch
outputfile=/hsm/nu/ninja/pra_tmp/wall_mc_20220616/output/bg_dist_${i}.root
bsub -q s ../build/src/Background/WallBG ${inputfile} ${outputfile}
done
