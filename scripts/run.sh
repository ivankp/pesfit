#!/bin/bash

./bin/pesfit data/h008/P* -l -f cb --ws-setRange='mean_offset_bin0:-1.5:1.' --prec=4 \
-o float_alpha.pdf
./bin/pesfit data/h008/P* -l -f cb --ws-setRange='mean_offset_bin0:-1.5:1.' --prec=4 \
--fix-alpha -o fix_alpha.pdf
