library(ggpubr)    
library(patchwork)

data <- read.csv("GSE117613.data.csv") 

#boxplot
ANGPT2 <- ggboxplot( data,  x="Group",  y="ANGPT2",  fill = "Group",   xlab = FALSE,
			  ylab = "Normalized Expression",  bxp.errorbar = TRUE,   bxp.errorbar.width = 0.4,
			  short.panel.labs = TRUE,  linetype = "solid",   width = 0.7,  notch = FALSE,  outlier.shape = 19,
			  add = "jitter",  error.plot = "pointrange",  repel = TRUE,  label.rectangle = FALSE) + theme_pubr(legend = "none")     
	
	
CA4 <- ggboxplot(data,  x="Group",  y="CA4",  fill = "Group",   xlab = FALSE,
			  ylab = "Normalized Expression",  bxp.errorbar = TRUE,   bxp.errorbar.width = 0.4,
			  short.panel.labs = TRUE,  linetype = "solid",   width = 0.7,  notch = FALSE,  outlier.shape = 19,
			  add = "jitter",  error.plot = "pointrange",  repel = TRUE,  label.rectangle = FALSE)+ theme_pubr(legend = "none") 
	

PADI4 <- ggboxplot( data,  x="Group",  y="PADI4",  fill = "Group",   xlab = FALSE,
			  ylab = "Normalized Expression",  bxp.errorbar = TRUE,   bxp.errorbar.width = 0.4,
			  short.panel.labs = TRUE,  linetype = "solid",   width = 0.7,  notch = FALSE,  outlier.shape = 19,
			  add = "jitter",  error.plot = "pointrange",  repel = TRUE,  label.rectangle = FALSE)+ theme_pubr(legend = "none") 
	
	
SLC2A3 <- ggboxplot( data,  x="Group",  y="SLC2A3",  fill = "Group",   xlab = FALSE,
			  ylab = "Normalized Expression",  bxp.errorbar = TRUE,   bxp.errorbar.width = 0.4,
			  short.panel.labs = TRUE,  linetype = "solid",   width = 0.7,  notch = FALSE,  outlier.shape = 19,
			  add = "jitter",  error.plot = "pointrange",  repel = TRUE,  label.rectangle = FALSE)+ theme_pubr(legend = "none") 
		
(ANGPT2 | CA4)  / (PADI4 | SLC2A3)
