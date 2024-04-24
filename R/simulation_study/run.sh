#!/bin/bash
# -- Un nom pour le job --
#OAR -n simulation_study_rect
# -- Description des ressources souhaitées -
#OAR -l /core=352,walltime=03:00:00
# -- Et bien sur, le projet Perseus --
#OAR --project pr-formation-ced-calcul

# -- Puis la liste des commandes à exécuter --
source /applis/site/guix-start.sh
R -e "rmarkdown::render('simulation_study_rect.Rmd')"
