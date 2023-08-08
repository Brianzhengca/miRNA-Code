  library(limma)
  raw <- read.csv("C://Users//xDxDxD//Desktop//MiRNA-Colorectal//cancer1exp", sep="\t")
  colorectalexp <- raw[,-1]
  rownames(colorectalexp) <- raw[,1]
  raw <- read.csv("C://Users//xDxDxD//Desktop//MiRNA-Colorectal//cancer1design", sep="\t")
  colorectaldesign <- raw[,-1]
  rownames(colorectaldesign) <- raw[,1]
  colorectal_cont <- makeContrasts(C_v_H=Y-N, levels=colorectaldesign)
  fit <- lmFit(colorectalexp, colorectaldesign)
  fit_contrast <- contrasts.fit(fit, colorectal_cont)
  fit_contrast <- eBayes(fit_contrast)
  volcanoplot(fit_contrast)
  fit_contrast.treat <- treat(fit_contrast,lfc=3)
  #top_genes <- topTable(fit_contrast, number = 100, adjust = "BH")
  #top_genes
  top_genes_treated <- topTreat(fit_contrast.treat,coef=1,sort.by="p", number=20000)
  res.treat <- decideTests(fit_contrast.treat)
  #summary(res.treat)
  #top_genes_treated
  sig_genes <- subset(top_genes_treated, top_genes_treated$adj.P.Val < 0.01)
  sig_genes_untreated <- subset(top_genes, top_genes$adj.P.Val < 0.01)
  write.csv(sig_genes_untreated, "C://Users//xDxDxD//Desktop//MiRNA-Colorectal//output_file.csv")
  volcanoplot(fit_contrast.treat)
  