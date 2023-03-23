# MetabolicModeling
## Bayesian models for predicting phenotypes based on metabolic networks
The scripts and files here are required to reproduce the results of the following paper:

	A Bayesian model for genomic prediction using metabolic networks
	Akio Onogi
	bioRxiv 2023.03.12.532311; doi: https://doi.org/10.1101/2023.03.12.532311

To run the scripts, files provided at https://github.com/Hao-Tong/netGS/tree/master/netGS are required.

The directory "/netGS" in the scripts need to include all files provided at the link.

The descriptions of files are as follows:

- EditData.R  
	Create input files. This script should be run first.
- MetaboliteOrder.txt  
	The order of metabolites in the stoichiometry matrix of Tong et al. (2020) in that of Arnold et al. (2014).
- TongH2020SupplementaryData2.csv  
	Supplementary Data2 of Tong et al. (2020). Abbreviations of reactions were added by Onogi. The supplementary data of Tong et al. (2020) is distributed under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
- CreateSimulation.R  
	Simulate data. CSV files are output.
- Simulation1.R  
	Analyzed the data simulated by CreateSimulation.R. Scripts for rrBLUP, quadratic programming, MegaLMM, and the proposed Bayesian models are included. The number (1) in the file name indicates the serial number of simulations. In the paper, 10 simulations (number 1-10) were conducted.
- CrossValidation1.R  
	Analyzed the real data of Arabidopsis provided at https://github.com/Hao-Tong/netGS/tree/master/netGS. Scripts for rrBLUP, quadratic programming, MegaLMM, and the proposed Bayesian models are included. The number (1) in the file name indicates the serial number of cross validations. In the paper, 20 cross validations (number 1-20) were conducted.
- MetabolicModeling10rho.stan  
	Stan scripts for the proposed Bayesian model.

The versions of R and R packages used are:  

- R: ver. 4.1.3 on Windows and 4.2.2 on Ubuntu 18.04  
- rstan: ver. 2.21.5 or 2.21.7
- rrBLUP: ver. 4.6.1
- snow: ver. 0.4-4
- CVXR: ver. 1.0-10
- MegaLMM: ver. 0.1.0

rstan may have troubles on R ver. 4.2 and later on Windows.
