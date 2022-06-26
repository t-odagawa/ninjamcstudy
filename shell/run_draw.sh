for i in `seq 0 1000`;do
#i=1
bsub -q s ../build/src/ReconDistribution/MomchDraw ${i} 5 5
done
