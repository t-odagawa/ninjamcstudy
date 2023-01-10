#/usr/bin/env python3

import subprocess

filepath = './run_draw_detsys_0pi.sh'

def main():
    #for isyst in [1, 4, 11]:
    for isyst in [2]:
    #for isyst in [3]:
    #for isyst in [5, 6, 8, 9]:
    #for isyst in [12]:
        dial = get_dial_from_id(isyst)
        if dial == "Nominal" :
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial)], shell=True)
        elif dial == "MPPCNoise" :
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " minus"], shell=True)
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " plus"], shell=True)
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " nominal"], shell=True)
        else : 
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " minus"], shell=True)
            subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " plus"], shell=True)

def get_dial_from_id(i):
    if i == 0 :
        dial = "BabyMind"
    elif i == 1 :
        dial = "NinjaBabyMindDistance"
    elif i == 2 :
        dial = "HitThreshold"
    elif i == 3 :
        dial = "MPPCNoise"
    elif i == 4 :
        dial = "Alignment"
    elif i == 5 :
        dial = "AngularResolution"
    elif i == 6 : 
        dial = "DetectionEfficiency"
    elif i == 7 :
        dial = "McsScaling"
    elif i == 8 :
        dial = "VphMean"
    elif i == 9 :
        dial = "VphSigma"
    elif i == 10 :
        dial = "MaterialThickness"
    elif i == 11 :
        dial = "PhysicsList"
    elif i == 12 :
        dial = "Nominal"

    return dial

if __name__ == "__main__":
    main()
