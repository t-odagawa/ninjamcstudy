#for i in `seq 30001 39999`;do
for i in `seq 1 49999`;do
#inputfile=/hsm/nu/ninja/pra_tmp/wall_mc_20220616/momch/momch_ecc5_${i}_pid.momch
#outputfile=/hsm/nu/ninja/pra_tmp/wall_mc_20220616/output/bg_dist_${i}.root
inputfile=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/wall/momch/momch_ecc5_${i}_sig.momch
outputfile=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/wall/output/bg_dist_${i}.root
bsub -q s ../build/src/Background/WallBG ${inputfile} ${outputfile} 
done
