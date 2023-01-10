#/usr/bin/env python3

import subprocess

filepath = './run_draw_xsecsys.sh'

def main():
    for isyst in range(0, 21):
        dial = get_dial_from_id(isyst)
        subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " minus"], shell=True)
        subprocess.call([filepath + " /hsm/nu/ninja/pra_tmp " + str(dial) + " plus"], shell=True)
    #   subprocess.call(["echo " + str(i)], shell=True)  

def get_dial_from_id(i):
    if i == 0 :
        dial = "MaCCQE"
    elif i == 1 :
        dial = "QETwk_HighQ2Weight_1"
    elif i == 2 :
        dial = "QETwk_HighQ2Weight_2"
    elif i == 3 :
        dial = "QETwk_HighQ2Weight_3"
    elif i == 4 :
        dial = "SF_OptPotTwkDial_O16"
    elif i == 5 :
        dial = "MECTwkDial_Norm_O16"
    elif i == 6 :
        dial = "MECTwkDial_PDDWeight_O16_NN"
    elif i == 7 :
        dial = "MECTwkDial_PDDWeight_O16_np"
    elif i == 8 :
        dial = "MECTwkDial_PNNN_Shape"
    elif i == 9 :
        dial = "RES_Eb_O_numu"
    elif i == 10 :
        dial = "BgSclRES"
    elif i == 11 :
        dial = "CA5RES"
    elif i == 12 :
        dial = "MaRES"
    elif i == 13 :
        dial = "PionFSI_AbsProb"
    elif i == 14 :
        dial = "PionFSI_CExHighMomProb"
    elif i == 15 :
        dial = "PionFSI_CExLowMomProb"
    elif i == 16 :
        dial = "PionFSI_InelProb"
    elif i == 17 :
        dial = "PionFSI_QEHighMomProb"
    elif i == 18 :
        dial = "PionFSI_QELowMomProb"
    elif i == 19 :
        dial = "TwkDial_FateNucleonFSI"
    elif i == 20 :
        dial = "CC_DIS_MultiPi_Norm_Nu"

    return dial

if __name__ == "__main__":
    main()
