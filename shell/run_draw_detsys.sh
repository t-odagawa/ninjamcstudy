prefix=$1

systematics=/detsys/$2

variation=/$3

for kinematics in 0 2 3;do
    #mcdir=${prefix}${systematics}${variation}
    mcdir=/hsm/nu/ninja/pra_tmp/detsys/Nominal
    momchdir=${prefix}${systematics}${variation}/momch
    outputdir=${prefix}${systematics}${variation}/output
    for i in `seq 1 1000`;do
	b2filename=${mcdir}/ninja_mc_h2o_${i}.root
	momchfilename=${momchdir}/momch_ecc5_${i}_sig.momch
	outputfilename=${outputdir}/output_mode${kinematics}_${i}.root
	if [ -f ${momchfilename} ];then
	    bsub -q s ../build/src/ReconDistribution/MomchDraw ${b2filename} ${momchfilename} ${outputfilename} ${kinematics}
	fi
    done
done
