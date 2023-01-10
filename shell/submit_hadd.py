#/usr/bin/env python3

import subprocess
from submit_draw import get_dial_from_id

filepath="./run_hadd.sh"

for isyst in range(0,20):
    dial, syst = get_dial_from_id(isyst)
    for imode in [0, 2, 3]:
        subprocess.call([filepath + " " + str(syst) + " " + str(imode)], shell=True)
