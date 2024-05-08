#!/bin/bash
#OAR -n rect_DistanceShore88
#OAR -l /core=16,walltime=01:30:00
#OAR --stdout rect_DistanceShore88.out
#OAR --stderr rect_DistanceShore88.err
#OAR --project pr-formation-ced-calcul
source /applis/site/guix-start.sh
Rscript rect_DistanceShore88.R
