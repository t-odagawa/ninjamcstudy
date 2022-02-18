#inputdir=${HOME}/data/mc_data/pra1
inputdir=${HOME}/data/mc_data/pra1/trackmatch
#outputdir=${HOME}/data/plots/muon_true_distribution
outputdir=${HOME}/data/mc_data/pra1/momch

for i in `seq 0 1000`;do
inputfile=${inputdir}/ninja_mc_h2o_ninjamatch_b2_${i}.root
#outputfile=${outputdir}/result_muon_${i}.root
outputfile=${outputdir}/kink_interaction_${i}.momch

#bsub -q s "../build/src/TrueDistribution/MuonTrueDistribution ${inputfile} ${outputfile} > /dev/null 2>&1"
bsub -q l "../build/src/Kink/KinkInteraction ${inputfile} ${outputfile}"
done
