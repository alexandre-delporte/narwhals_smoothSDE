#!/bin/bash
#OAR -n rect_DistanceShore60
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore60.out
#OAR --stderr rect_DistanceShore60.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore60.R
