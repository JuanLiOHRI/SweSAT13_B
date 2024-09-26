# Information Theory for Tests

Code for analysis in the manuscript entitled "*An Information Manifold Perspective for Analyzing Test Data*".

doi: [link to be inserted]

## Data files

Data files are inside the `data` folder.

1.  `Quant_13B_sample.csv`: A subset of the SweSAT Quant 13B test data (N = 5379) that maintains the sum score distribution as the full data.
2.  `Quant_13B_key.txt`: Answer keys of the test.

Intermediate results will be saved to the `data/saved results` folder. See details below **Note:** This folder is empty in the Github version.

## Code files

__Important:__ This project used TestGardener version 3.2.6. TestGardener has had a major update recently (version 3.3.3) which includes various changes that are incompatible with code here. We will update the repo to accommodate the changes, but until that completes, please use TestGardener version 3.2.6 to run this code.

1.  `1. APM 2024 Quant_13b TG.Rmd`
    1.  Code for the TestGardener analysis.
    2.  Need to read
        -   `data/Quant_13B_sample.csv`
        -   `data/Quant_13B_key.txt`
    3.  Generate files:
        -   `data/saved results/Quant_13B_dataList.rds`: Results for setting up the data
        -   `data/saved results/WfdResult0.rds`: Preliminary displays and plots before `Analyze` show initial mean fitting function and initial surprisal curves.
        -   `data/saved results/AnalyzeResult_20.rds"`: Results after the 20 cycles.
        -   `data/saved results/Quant_13B_infoList.rds`: Results after transforming score indices to arc length values.
2.  `2. APM 2024 Quant_13b mirt.Rmd`
    1.  Code for the nominal (mirt) analysis.
    2.  Need to read
        -   `data/Quant_13B_sample.csv`
        -   `data/Quant_13B_key.txt`
    3.  Generate files:
        -   `data/saved results/mirt_mod_13.rds`: Results of the nominal model.
        -   `data/saved results/Quant_13B_dataList_nom.rds`: Results for setting up the data
        -   `data/saved results/Quant_13B_parList_nom.rds"`: Results for setting up for TG style plotting.
        -   `data/saved results/Quant_13B_infoList_nom.rds`: Results after transforming theta values to arc length values.
3.  `3 APM 2024 Quant_13b Results.Rmd`
    1.  Code for generating results and figures.
    2.  Need to read all the saved results from the above two files.

## Supportive function

Supportive functions are placed in the `R` folder.

1.  `nominal_surpmat.R`: calculate the surprisal matrix for the nominal model.
2.  `ICC.plot.R`: Modified `TestGardener` functions to generate to make the background in white.
