#/usr/bin/env python3

import subprocess

filepath = './run_draw.sh'

def main():
    for isyst in range(1, 20):
        dial, syst = get_dial_from_id(isyst)
        #subprocess.call(["echo " + str(syst)], shell=True)
        subprocess.call(["bsub -q s " + filepath + " " + str(syst)], shell=True)

def get_dial_from_id(i):
    if i == 0 :
        dial = "MaCCQE"
        syst = "maccqe"
    elif i == 1 :
        dial = "QETwk_HighQ2Weight_1"
        syst = "highq2_1"
    elif i == 2 :
        dial = "QETwk_HighQ2Weight_2"
        syst = "highq2_2"
    elif i == 3 :
        dial = "QETwk_HighQ2Weight_3"
        syst = "highq2_3"
    elif i == 4 :
        dial = "SF_OptPotTwkDial_O16"
        syst = "sf_opt_pot"
    elif i == 5 :
        dial = "MECTwkDial_Norm_O16"
        syst = "mec_norm"
    elif i == 6 :
        dial = "MECTwkDial_PDDWeight_O16_NN"
        syst = "mec_pdd_nn"
    elif i == 7 :
        dial = "MECTwKDial_PDDWeight_O16_np"
        syst = "mec_pdd_np"
    elif i == 8 :
        dial = "MECTwKDial_PNNN_Shape"
        syst = "mec_pnnn_shape"
    elif i == 9 :
        dial = "RES_Eb_O_numu"
        syst = "res_eb"
    elif i == 10 :
        dial = "BgScIRES"
        syst = "iso_nonres"
    elif i == 11 :
        dial = "CA5RES"
        syst = "ca5_res"
    elif i == 12 :
        dial = "MaRES"
        syst = "mares"
    elif i == 13 :
        dial = "PionFSI_AbsProb"
        syst = "pi_fsi_abs"
    elif i == 14 :
        dial = "PionFSI_CExHighMomProb"
        syst = "pi_fsi_cex_high"
    elif i == 15 :
        dial = "PionFSI_CExLowMomProb"
        syst = "pi_fsi_cex_low"
    elif i == 16 :
        dial = "PionFSI_InelProb"
        syst = "pi_fsi_inel"
    elif i == 17 :
        dial = "PionFSI_QEHighMomProb"
        syst = "pi_fsi_qe_high"
    elif i == 18 :
        dial = "PionFSI_QELowMomProb"
        syst = "pi_fsi_qe_low"
    elif i == 19 :
        dial = "T2kDial_FateNucleonFSI"
        syst = "nucl_fsi"

    return dial, syst


if __name__ == "__main__":
    main()

