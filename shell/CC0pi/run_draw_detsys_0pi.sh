prefix=$1

systematics=/detsys/$2

variation=/$3

for kinematics in 0 2 3;do
#for kinematics in 6 7;do
    #mcdir=${prefix}${systematics}${variation}
    mcdir=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620 # Hit Threshold
    #mcdir=/hsm/nu/ninja/pra_tmp/detsys/Nominal
    momchdir=${prefix}${systematics}${variation}/momch

    if [ ! -d /hsm/nu/ninja/pra_tmp/CC0pi_20221213${systematics}${variation}/output ]; then
	mkdir -p /hsm/nu/ninja/pra_tmp/CC0pi_20221213${systematics}${variation}/output
    fi

    outputdir=/hsm/nu/ninja/pra_tmp/CC0pi_20221213${systematics}${variation}/output
    for i in `seq 1 1000`;do
	b2filename=${mcdir}/ninja_mc_h2o_${i}.root
	momchfilename=${momchdir}/momch_ecc5_${i}_sig.momch
	outputfilename=${outputdir}/output_mode${kinematics}_${i}.root
	if [ -f ${momchfilename} ];then
	    bsub -q s ../../build/src/ReconDistribution/MomchDraw ${b2filename} ${momchfilename} ${outputfilename} ${kinematics}
	fi
    done
done
