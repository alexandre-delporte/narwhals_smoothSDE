echo "Preprocessing..."
cd preprocessing
Rscript process_boundaries.R
Rscript -e "rmarkdown::render('data_preprocessing.Rmd')"

echo "Simulation study..."
cd ..
cd simulation_study
Rscript -e "rmarkdown::render('simulation_study_rect.Rmd')"
Rscript -e "rmarkdown::render('simulation_study_fjords.Rmd')"



echo "Narwhals data analysis..."
cd ..
cd mixed_effects
cd no_constraints
echo "Unconstrained SDEs..."



echo "Constrained SDEs..."



