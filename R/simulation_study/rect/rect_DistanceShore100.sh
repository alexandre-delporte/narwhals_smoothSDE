#!/bin/bash
#OAR -n rect_DistanceShore100
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore100.out
#OAR --stderr rect_DistanceShore100.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore100.R
