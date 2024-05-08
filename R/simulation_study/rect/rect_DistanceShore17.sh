#!/bin/bash
#OAR -n rect_DistanceShore17
#OAR -l /core=16,walltime=01:00:00
#OAR --stdout rect_DistanceShore17.out
#OAR --stderr rect_DistanceShore17.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore17.R
