# Methylation Analysis
DNA methylation workflow, from beta-values to signatures  

This workflow contains 4 modules.
1. Normalization
    - Normalization.R
2. DMP & DMR
    - DMP_DMR.R
    - 2013Method.R & 2013Method_cv.R  
    
　Test the methods in *Exploring genome-wide DNA methylation profiles altered in hepatocellular carcinoma using Infinium HumanMethylation 450 BeadChips. Epigentics 2013*.  
    
　By using the methods, the hyper/hypo-methylated sites can be distinguished from DMPs. The dmrcate function was re-defined and the number of hyper/hypo-methylated sites covered by each DMR are also a part of the output now.  
3. Random Forests
    - RFBoot.R
4. Lasso
    - Lasso.R
 
# Instructions
* Requirement of each module can be seen at the beginning of the scripts.
* Rscripts and data can be found at "/lustre/rdi/user/hantc/tools/Methylation", just copy them to your own directory and test whether it can be run.
* It is to be noted that packages required by Normalization and DMP_DMR cannot be easily installed on the cluster, so I recommend you to run these 2 modules on your PC.
* RFBoot takes a long time. Running it on a Cluster is a better option. 
* If the sample size is small, Lasso may return an error.


# Change Logs
################################################

2017-12-11

Uploaded new scripts for testing the methods in *Exploring genome-wide DNA methylation profiles altered in hepatocellular carcinoma using Infinium HumanMethylation 450 BeadChips. Epigentics 2013*.
* By using the methods, the hyper/hypo-methylated sites can be distinguished from DMPs. The dmrcate function was re-defined and the number of hyper/hypo-methylated sites covered by each DMR are also a part of the output now.


################################################

2017-12-01

Upload 4 modules.
