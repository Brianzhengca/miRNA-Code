# miRNA-Code
This project is the source code for identifying three prognostic and diagnostic biomarkers for human colorectal cancer. Certain miRNA expression files referenced in the code are too large to be uploaded to GitHub and will be available upon request (Contact at briansiyuanzheng@gmail.com)

**differentialexpressionanalysis.R**: The R code that was used in Limma to conduct differential expression analysis, requires a pre-built design and  expression matrix, both of which can be obtained by contacting the author. 

**get_uniquely_expressed_miRNAs.py**: Finds miRNAs that are expressed at unique levels between colorectal cancer and all other cancers.

**linreg.py**: Linear Regression to Find miRNAs that are expressed at different levels across stages. Requires two files: "colorectaltoalldesign" and "colorectaltoallexp", which are design matrix and expression matrix for colorectal cancer against all other cancers plus healthy controls. Can be obtained by contacting the author. 

**logisticregressioncode.py**: Code for Logistic Regression Classifier. Requires two files: "colorectaltoalldesign" and "colorectaltoallexp", which are design matrix and expression matrix for colorectal cancer against all other cancers plus healthy controls. Can be obtained by contacting the author. This code is also used for pairwise classification tasks involving colorectal cancer and each of the other cancers plus healthy controls. The design matrix and expression matrix for those pairs can be obtained by contacting the author. 
