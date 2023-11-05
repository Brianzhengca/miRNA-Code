# miRNA-Code
This project is the source code for identifying three prognostic and diagnostic biomarkers for human colorectal cancer. Certain miRNA expression files referenced in the code are too large to be uploaded to GitHub and will be available upon request (Contact at briansiyuanzheng@gmail.com)

**differentialexpressionanalysis.R**: Code that was used in Limma to conduct differential expression analysis (moderated t-tests), requires a pre-built design and  expression matrix. As the design and expression matrices are too large to be included in this repository, they are available upon contact. 

**get_differentially_expressed_miRNAs.py**: Generates violin plots for miRNAs that are expressed at different levels between stages.

**linreg.py**: Linear regression t-test to Find miRNAs that are expressed at different levels across stages. 

**logisticregressioncode.py**: Code for Logistic Regression Classifier, the mirna_feature file contacts the name of the miRNA that is used as the sole feature for the classifier. This file also graphs the ROC curves for the classifier. It requires a design and expression matrix, both of which are available upon request. 
