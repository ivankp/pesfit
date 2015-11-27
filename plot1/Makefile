mc_8TeV  := ../data/PAPERAC_Final_mc12_skimmed_new_weight.root
mc_13TeV := $(shell find ../data/h009 -type f -name '*.root')

m_yy_range := 105:160

.SECONDEXPANSION:

.SECONDARY: \
	mc_8TeV.root mc_13TeV.root \
	fit_8TeV_unweighted.root fit_8TeV_weighted.root \
	fit_13TeV_unweighted.root fit_13TeV_weighted.root \
	fit_8TeV.root fit_13TeV.root

all: plot1.pdf

mc_%.root: $$(mc_%)
	@echo $@
	@../bin/data_select $@ $^ -x $(m_yy_range)
	@chmod 444 $@

fit_%_unweighted.root: mc_%.root
	@echo $@
	@../bin/fitmyy ../data/ws2.root $^ $@

fit_%_weighted.root: mc_%.root
	@echo $@
	@../bin/fitmyy ../data/ws2.root $^ $@ -w

fit_%.root: fit_%_weighted.root fit_%_unweighted.root
	@echo $@
	@../bin/fix_stats $^

plot1.pdf: fit_8TeV.root fit_13TeV.root
	@echo $@
	@../bin/plot1 $^ $@
