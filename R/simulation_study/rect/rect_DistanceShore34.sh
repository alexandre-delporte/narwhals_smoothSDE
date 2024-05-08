#!/bin/bash
#OAR -n rect_DistanceShore34
#OAR -l /core=16,walltime=01:00:00
#OAR --stdout rect_DistanceShore34.out
#OAR --stderr rect_DistanceShore34.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore34.R
