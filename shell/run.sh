inputdir=${HOME}/data/mc_data/pra1
outputdir=${HOME}/data/plots/muon_true_distribution

for i in `seq 0 400`;do
inputfile=${inputdir}/ninja_mc_h2o_${i}.root
outputfile=${outputdir}/result_muon_${i}.root

bsub -q s "../build/src/TrueDistribution/MuonTrueDistribution ${inputfile} ${outputfile} > /dev/null 2>&1"
done
