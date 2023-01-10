for i in `seq 1 49999`;do
#for i in `seq 1 2`;do
b2file=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/wall/ninja_mc_wall_${i}.root
momchfile=/group/nu/ninja/work/odagawa/20220930-bkg-mc-new-matching/wall/momch/momch_ecc5_${i}_pid.momch
ofile=/group/nu/ninja/work/odagawa/20221027-shingakujutsu-neutrino-preliminary/wall-2d-hist/output_${i}.root
bsub -q s ../build/src/ReconDistribution/WallMomchDraw ${b2file} ${momchfile} ${ofile}
#../build/src/ReconDistribution/WallMomchDraw ${b2file} ${momchfile} ${ofile}
done
