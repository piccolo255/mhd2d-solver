#!/bin/bash

###############################################################################################################################################
# g32 tau2 beta0.2
mkdir -p "data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition"
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0004-00.1250000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0004-00.1250000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0008-00.2500000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0008-00.2500000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0016-00.5000000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0016-00.5000000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0024-00.7500000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0024-00.7500000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0032-01.0000000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0032-01.0000000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0064-02.0000000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0064-02.0000000.wdat
perl scripts/find-wave-decomposition.pl -i data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/shlf-32-t2-beta0.2-u1.0-0128-04.0000000.dat -o data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decomposition/waves-shlf-32-t2-beta0.2-u1.0-0128-04.0000000.wdat
