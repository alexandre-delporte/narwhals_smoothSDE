#!/bin/bash
#OAR -n rect_DistanceShore87
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore87.out
#OAR --stderr rect_DistanceShore87.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore87.R
