#!/bin/bash
#OAR -n rect_DistanceShore6
#OAR -l /core=16,walltime=01:00:00
#OAR --stdout rect_DistanceShore6.out
#OAR --stderr rect_DistanceShore6.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore6.R