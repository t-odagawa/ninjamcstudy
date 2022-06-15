for i in `seq 0 1000`;do
inputfile=/home/t2k/odagawa/data/mc_data/wall/ninja_mc_wall_${i}.root
outputfile=/hsm/nu/ninja/pra_tmp/wall_angle_corr/angle_corr_${i}.root
bsub -q s ../build/src/ReconDistribution/AngleCorrection $inputfile $outputfile
done
