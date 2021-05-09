#!/bin/bash

./process-all-inputs-in-dir.sh inputs/g16/t1 data/g16-t1
perl plot.pl inputs/plotparams-16-t1.cfg 
perl plot-cuts.pl inputs/plotcutparams-16-t1.cfg 
perl find-thinning-velocity-from-sheet-profile.pl inputs/findvelparams-g16-beta0.2.cfg
