#!/bin/bash
#OAR -n rect_DistanceShore52
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore52.out
#OAR --stderr rect_DistanceShore52.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore52.R
