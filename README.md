# Methylation_Analysis
DNA methylation workflow, from beta-values to signatures

################################################

2017-12-01

This workflow contains 4 modules.
1. Normalization.R
2. DMP_DMR.R
3. RFBoot.R
4. Lasso.R

Instructions
* Requirement of each module can be seen at the beginning of the scripts.
* Rscripts and data can be found at "/lustre/rdi/user/hantc/tools/Methylation", just copy them to your own directory and test whether it can be run.
* It is to be noted that packages required by Normalization and DMP_DMR cannot be easily installed on the cluster, so I recommend you to run these 2 modules on your PC.
* RFBoot takes a long time. Running it on a Cluster is a better option. 
* If the sample size is small, Lasso may return an error.
