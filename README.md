# The structure of the hematopoietic system can explain chronic myeloid leukemia progression
This repository provides the data and the code necessary to reproduce the results published in the article:

**Include link to the article**

## Computer simulation of the CML progression model

An single simulation of CML progression can be produced by running the code cml_simulator.py, with the proper
arguments. The code will produce two files a bone and blood txt files, containing a sample of the status
of the bone marrow and bloodstream at different timepoints. 

## Simulated data used in the manuscript figures

The data simulated using the cml_simulator.py is presented in the folder simulated_data.
This folder contains two files with the aggregated results from the optimization
of the CML model parameters. This files cannot be reproduced directly with the
cml_simulator.py, the entire optimization process has a high resource consumption
and will not be easy to reproduce. For more detailed instructions of how to 
reproduce the entire optimization data please contact the corresponding author.

## Figures in the manuscript

Five figures are currently included in the manuscript which can be found in the manuscript_figures folder.
- mutational_dynamics.svg 
  - This figure can be reproduced using the figure_creator.py code. This is Figure 6 of the manuscript,
  showing the mutational dynamics at different timepoints.
- schematic.eps
  - This figure cannot be reproduced from any code.
  In the figure we represent qualitatively CML progression.
- median_durations.eps
  - This figure can be reproduced using the figure_creator.py code. This is Figure 4 in the manuscript,
  showing the boxplots for the CML progression phases durations.
- hierarchical_model.eps
  - This figure cannot be reproduced from any code.
  In the figure we sketch the hierarchical structure representing the hematopoiesis.
  We also present the interaction between the bone marrow and bloodstream,
  and the effect of the driver mutations on the cell rates.
- bcr_sc_relation.svg
  - This figure can be reproduced using the figure_creator.py code. This is Figure 5 in the manuscript,
  showing the filtered evolutionary parameter combinations for the computer simulations.
- individual_time_graph.eps
  - This figure can be reproduced using the figure_creator.py code. This is Figure 3 in the manuscript,
  the three panels show the result of an individual simulation of the model, where the background shadow 
  indicates the different phases of progression of CML

## Reproducing figures from the manuscript

The figures in the manuscript can be reproduced running the file create_manuscript_figures.py
This code will create a new folder in the working directory with 4 images for the 
figures [3-6] in the manuscript. The figures are created with the simulated data 
available in the simulated_data folder.