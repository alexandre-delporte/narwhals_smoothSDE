#!/bin/bash
#OAR -n rect_DistanceShore72
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore72.out
#OAR --stderr rect_DistanceShore72.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore72.R
