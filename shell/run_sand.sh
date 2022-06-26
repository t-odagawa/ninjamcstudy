for i in `seq 30001 40000`;do
inputfile=/hsm/nu/ninja/pra_tmp/wall_mc_20220616/track/ninja_mc_wall_track_${i}.root
outputfile=/hsm/nu/ninja/pra_tmp/wall_angle_corr/angle_corr_${i}.root
bsub -q s ../build/src/ReconDistribution/AngleCorrection $inputfile $outputfile
done
