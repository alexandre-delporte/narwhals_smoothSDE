#!/bin/bash
#OAR -n rect_DistanceShore85
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore85.out
#OAR --stderr rect_DistanceShore85.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore85.R
