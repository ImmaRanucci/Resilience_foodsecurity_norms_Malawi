# Weather Shocks and Resilience to Food Insecurity: Exploring the Role of Gender and Kinship Norms - Replication Code
This README describes the directory structure & Stata code necessary to replicate the main analysis presented in the paper "Weather Shocks and Resilience to Food Insecurity: Exploring the Role of Gender and Kinship Norms", authored by Immacolata Ranucci, Donato Romano, and Luca Tiberti, to be published in World Development. 

The dataset available in this replication package contains cleaned data, including constructed indicators. MIHS data is merged with SPEI data. This dataset includes all variables necessary to produce the main results and also all additional variables used in the robustness checks. We provide the STATA code to reproduce the main results presented in the paper (code for robustness checks is available upon request). 

[![DOI](https://zenodo.org/badge/895882376.svg)](https://doi.org/10.5281/zenodo.14242898)


## Index
 - [Contributors](#contributors)
 - [Data](#data-sources)
 - [STATA Code](#stata-code)

## Contributors
* Immacolata Ranucci
* Donato Romano
* Luca Tiberti

## Data Sources
The work relies on the Malawi Integrated Household Panel Survey (MIHPS) 2010-2013-2016, part of the Living Standards Measurement Study (LSMS) project of the World Bank, publicly accessible in the World Bank Microdata Library. Data is spatially merged with the Standardised Precipitation-Evapotranspiration Index (SPEI), publicly available in the database of the SPEI Global Drought Monitor.

The World Bank LSMS survey data is publicly available for each survey round together with questionnaires, basic information documents, interview manuals, and metadata. Panel survey data and additional material can be accessed at the following link:
- Integrated Household Panel Survey 2010-2013-2016: https://microdata.worldbank.org/index.php/catalog/2939

The Standardised Precipitation-Evapotranspiration Index (SPEI) time series can be downloaded as a single grid cell index for a selected geographical area at the following link:
- Global Drought Monitor: https://spei.csic.es/map/maps.html#months=1#month=9#year=2024

The cleaning process of household-level survey data, preparation of SPEI index, and merging of the two datasets is explained in the paper. We provide here the final version of the cleaned dataset (Data.dta), which contains all relevant variables used in the main analysis and robustness checks. 

## STATA code
The Stata do file (Code.do) available in this replication package allows the reproduction of the tables and figures containing the main results presented in the paper. 

### Instructions for Replication
Open the Code.do file and update the global "PATH" with your username in the Filepath section. 
The Data.dta file should be saved in the folder that you indicate as PATH. 

The code includes commands to generate sub-folders for storing results within the folder PATH.

The following user-written STATA programs are necessary to run the analysis:
- exbsample: https://www.stata.com/meeting/uk22/slides/UK22_Van_Kerm.pdf
- estout: https://repec.sowi.unibe.ch/stata/estout/
- grc1leg2: https://ideas.repec.org/c/boc/bocode/s459360.html

The .do file includes lines of code to install these packages.

Run the Code.do file. Output tables and graphs will be saved into the Tables and Figures sub-folders within the PATH folder. 




