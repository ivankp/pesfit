#!/bin/bash

mkdir -p nov16

# 8 TeV

rm -v .build/pesfit.o .build/confnoteplot1.o
make -j4 Old8TeVFile=1

./bin/pesfit data/PAPERAC_Final_mc12_skimmed_new_weight.root -w data/ws2.root -o nov16/fit_8TeV_nosumw2_norm.root -f cb --norm
./bin/confnoteplot1 nov16/fit_8TeV_nosumw2_norm.root nov16/fit_8TeV_nosumw2_norm.pdf minsig

./bin/pesfit data/PAPERAC_Final_mc12_skimmed_new_weight.root -w data/ws2.root -o nov16/fit_8TeV_norm.root -f cb --norm --sumw2
./bin/confnoteplot1 nov16/fit_8TeV_norm.root nov16/fit_8TeV_norm.pdf minsig

./bin/pesfit data/PAPERAC_Final_mc12_skimmed_new_weight.root -w data/ws2.root -o nov16/fit_8TeV_nosumw2.root -f cb
./bin/confnoteplot1 nov16/fit_8TeV_nosumw2.root nov16/fit_8TeV_nosumw2.pdf minsig

./bin/pesfit data/PAPERAC_Final_mc12_skimmed_new_weight.root -w data/ws2.root -o nov16/fit_8TeV.root -f cb --sumw2
./bin/confnoteplot1 nov16/fit_8TeV.root nov16/fit_8TeV.pdf minsig

# 13 TeV

rm -v .build/pesfit.o .build/confnoteplot1.o
make -j4 Old8TeVFile=0

./bin/pesfit data/h009/* -w data/ws2.root -o nov16/fit_13TeV_nosumw2_norm.root -f cb --norm
./bin/confnoteplot1 nov16/fit_13TeV_nosumw2_norm.root nov16/fit_13TeV_nosumw2_norm.pdf minsig

./bin/pesfit data/h009/* -w data/ws2.root -o nov16/fit_13TeV_norm.root -f cb --norm --sumw2
./bin/confnoteplot1 nov16/fit_13TeV_norm.root nov16/fit_13TeV_norm.pdf minsig

./bin/pesfit data/h009/* -w data/ws2.root -o nov16/fit_13TeV_nosumw2.root -f cb
./bin/confnoteplot1 nov16/fit_13TeV_nosumw2.root nov16/fit_13TeV_nosumw2.pdf minsig

./bin/pesfit data/h009/* -w data/ws2.root -o nov16/fit_13TeV.root -f cb --sumw2
./bin/confnoteplot1 nov16/fit_13TeV.root nov16/fit_13TeV.pdf minsig
