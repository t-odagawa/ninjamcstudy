for i in `seq 0 1000`;do
#i=1
dirname=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620
b2filename=${dirname}/track/ninja_mc_h2o_track_${i}.root
momchfilename=${dirname}/momch/momch_ecc5_${i}.momch
momchbmfilename=${dirname}/momch/momch_ecc5_${i}_pid.momch
outputfilename=${dirname}/eff_br/efficiency_breakdown_${i}.root
bsub -q s ../build/src/Efficiency/EventEfficiencyBreakdown ${b2filename} ${momchfilename} ${momchbmfilename} ${outputfilename}
#../build/src/Efficiency/EventEfficiencyBreakdown ${b2filename} ${momchfilename} ${outputfilename}
done
