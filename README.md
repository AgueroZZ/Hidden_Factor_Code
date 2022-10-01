# Hidden_Factor_Code
The online web-page of the manuscript "The hidden factor: accounting for covariate effects in power and sample size computation for a binary trait"

## Data

- **Section 4.1-4.2**: simulated datasets are generated in the R scripts.

- **Section 5**: datasets are obtained from UK Biobank.

## Scripts

The simulation/examples in the main paper are replicated in the following scripts. Specific points to pay attention to:
1. In each simulation/example, please make sure the working directory is correctly specified.
2. All the required R packages will be loaded in the R script, please make sure they are installed in the system before proceed.
3. For the codes of section 5, please run them in the following orders "01_QC.R" to "02_GWAS.R" to "03_Computation.R".

The results in the paper are reproduced as follows (please make sure each cpp file is compiled before running the R script in the corrsponding folder):

- **Section 4.1-4.2**:  
     - *script.R* generates the dataset and runs the simulation. It will also produce the figures used in the manuscript, and save them in the current working directory.

- **Section 5**:
	 - *01_QC.R* does the quality control works that are preliminary to the GWAS analysis in section 5, and produces the corresponding PC plots.
	 - *02_GWAS.R* uses the output from *01_QC.R* to produce the GWAS analysis and figures in the section 5.
	 - *03_Computation.R* does the power/sample size computation with the GWAS summary statistics as input, using the R package SPCompute.
