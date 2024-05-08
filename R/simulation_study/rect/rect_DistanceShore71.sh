#!/bin/bash
#OAR -n rect_DistanceShore71
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore71.out
#OAR --stderr rect_DistanceShore71.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore71.R
