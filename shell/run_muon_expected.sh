#!/bin/sh

cd ${HOME}/NinjaMCStudy/build

subrun=pra1
material=h2o

b2datadir=${HOME}/data/mc_data/${subrun}/track
b2datafile=${b2datadir}/ninja_mc_${material}_track_$1.root

matchdir=${HOME}/data/mc_data/${subrun}/trackmatch
matchfile=${matchdir}/ninja_mc_${material}_ninjamatch_$1.root

outputdir=${HOME}/data/plots/${subrun}/muon_expected
outputfile=${outputdir}/muon_expected_distribution_$1.root

./src/MatchAcceptance/MuonExpectedDistribution ${b2datafile} ${matchfile} ${outputfile}
