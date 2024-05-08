#!/bin/bash
#OAR -n rect_DistanceShore25
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore25.out
#OAR --stderr rect_DistanceShore25.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore25.R
