#!/bin/bash

./bin/fitmyy data/ws2.root data/mc_8TeV.root fit_8TeV_unweighted.root
./bin/fitmyy data/ws2.root data/mc_8TeV.root fit_8TeV_weighted.root -w
./bin/fitmyy data/ws2.root data/mc_13TeV_h009.root fit_13TeV_unweighted.root
./bin/fitmyy data/ws2.root data/mc_13TeV_h009.root fit_13TeV_weighted.root -w

./bin/fix_stats fit_8TeV_weighted.root fit_8TeV_unweighted.root
./bin/fix_stats fit_13TeV_weighted.root fit_13TeV_unweighted.root
