#prefix=/hsm/nu/ninja/pra_tmp/mc_tmp_20220505
prefix=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620

for i in `seq 0 1000`;do
#i=1
bsub -q s ../build/src/ReconDistribution/MomchDraw ${prefix} ${i} 5 2
done 
