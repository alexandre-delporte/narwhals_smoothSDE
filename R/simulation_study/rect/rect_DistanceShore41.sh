#!/bin/bash
#OAR -n rect_DistanceShore41
#OAR -l /core=16,walltime=01:00:00
#OAR --stdout rect_DistanceShore41.out
#OAR --stderr rect_DistanceShore41.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore41.R
