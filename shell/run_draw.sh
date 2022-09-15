#prefix=/hsm/nu/ninja/pra_tmp/mc_tmp_20220505
prefix=/hsm/nu/ninja/pra_tmp/mc_tmp_20220620

kinematics=0 #0 multi 1 momentum 2 angle

#systematics=0
#systematics=mcs_scale_syst
systematics=angres_syst
#variation=nominal
variation=minus
#variation=plus

for i in `seq 0 1000`;do
#i=1
bsub -q s ../build/src/ReconDistribution/MomchDraw ${prefix} ${i} 5 ${kinematics} ${systematics} ${variation}
done 
