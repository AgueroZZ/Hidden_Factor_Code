# Hidden_Factor_Code
The online web-page of the manuscript "The hidden factor: accounting for covariate effects in power and sample size computation for a binary trait"

## Data

- **Section 4.1-4.2**: simulated datasets are generated in the R scripts.

- **Section 5**: datasets are obtained from UK Biobank.

## Scripts

The simulation/examples in the main paper are replicated in the following scripts. Specific points to pay attention to:
1. In each simulation/example, please make sure the working directory is correctly specified.
2. All the required R packages will be loaded in the R script, please make sure they are installed in the system before proceed.

The results in the paper are reproduced as follows (please make sure each cpp file is compiled before running the R script in the corrsponding folder):

- **Section 4.1**:  
     - *script.R* generates the dataset and runs the simulation.
     - *script.R* will also produce the figures used in the manuscript, and save them in the current working directory.

- **Section 4.2**:
	 - *script.R* runs the simulation and then *plot.R* does the plotting works in the main paper.
