# RAD_R_Functions

This project will be a collection of useful R code that our group at the Bedford Institute of Oceanography will be using for
genetic analysis using RAD-seq and SNP data. 

# Discliamer 

This code is use at own risk. We are starting to use GitHub for collaborating and starting a uniform RAD-Seq data pipeline. This code will evolve significantly over the next few months. 

## GenePopConvert.R
GenePop convert is a function which will prune a GenePop file for outliers or neutral markers.
This function permits a quick division of SNP data in the GenePop format .txt file. The output can be converted to 
other file types (STRUCTURE, Arelquin, Genodiver etc) using conventional conversion programs such as PGDSpider. 

  #Example GenePop Subset.R
  This file will document how the code works. We use an example dataset named "Sal.txt" which is availble for download:
  https://www.dropbox.com/s/jz79ed39o6tynpm/SAL.txt?dl=0
  
# 2Plink.R 
This function will help to create the command line interface code to check Stacks output for Hardy-Weinberg and MAF data checks. 
For now the output from this code can be copied and pasted into PLINK for data filtering. 
