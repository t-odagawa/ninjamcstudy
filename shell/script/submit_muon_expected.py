#/usr/bin/env python3

import subprocess

filepath = '../run_muon_expected.sh'
nrun = 1000
for i in range(1, nrun+1) :
    subprocess.call(["bsub -q s " + filepath + " " + str(i)], shell=True)
