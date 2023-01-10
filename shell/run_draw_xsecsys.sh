mcprefix=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620

momchprefix=$1

systematics=/xsecsys/$2

variation=/$3

if [ ! -d ${momchprefix}${systematics}${variation}/output ]; then
    mkdir ${momchprefix}${systematics}${variation}/output
fi

for kinematics in 0 2 3;do  # 0 multi 2 momentum 3 angle
#for kinematics in 0 2 3 4;do # 4 PID
    momchdir=${momchprefix}${systematics}${variation}/momch
    outputdir=${momchprefix}${systematics}${variation}/output
    for i in `seq 1 1000`;do
#    for i in `seq 1 1`;do
	b2filename=${mcprefix}/ninja_mc_h2o_${i}.root
	momchfilename=${momchdir}/momch_ecc5_${i}_sig.momch
	outputfilename=${outputdir}/output_mode${kinematics}_${i}.root
	bsub -q s ../build/src/ReconDistribution/MomchDraw ${b2filename} ${momchfilename} ${outputfilename} ${kinematics}
#	echo ${b2filename} ${momchfilename} ${outputfilename} ${kinematics}
    done 
done
