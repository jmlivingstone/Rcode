# target file has the names of the .cel files and the model classification
load("Target.rda")
names(target) < c("filename","Cy3","Cy5","mouse.model")
data <- read.maimages(target$filename ,source="agilent")

# normalize data using limma package
data.bg <- backgroundCorrect(data,method="rma")
data.w  <- normalizeWithinArrays(data.bg,method="loess",bc.method="none")
MA <- normalizeBetweenArrays(data.w,method="quantile")

design <- modelMatrix(target,ref="Ref")

contrast.matrix <- makeContrasts((G1+G2-F1-F2)/2,
        (N1+N2-NE1-NE.2)/2,
        (NDL1+NDL2-NDLE1-NDLE.2)/2,
        levels=design)

fit <- lmFit(MA, design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

top1 <- topTable(fit2,coef=1,number=nrow(fit2))
top2 <- topTable(fit2,coef=2,number=nrow(fit2))
top3 <- topTable(fit2,coef=3,number=nrow(fit2))

#plot(top1$logFC)
  #want to see a peak and then equal
#hist(top1$P.Value,100)
  #$M is the fold change - will be opposite for dye swaps
#heatmap(julie.w$M[as.numeric(rownames(top1)[top1$adj.P.Val <= 0.005]),])
#heatmap(fit$coefficients[as.numeric(rownames(top1)[top1$adj.P.Val <= 0.005]),])

# because there are three different groups being compared - index only applicable patients
sean.heatmap(fit$coefficients[as.numeric(rownames(top1)[top1$adj.P.Val <= 0.05][1:200]),c(1,2,7,8)],
        midcolour = "black")

sean.heatmap(fit$coefficients[as.numeric(rownames(top2)[top2$adj.P.Val <= 0.05][1:200]),c(7,8,9,10)],
        bottom.margin = 12, midcolour = "black")

sean.heatmap(fit$coefficients[as.numeric(rownames(top3)[top3$adj.P.Val <= 0.05][1:127]),c(3,4,5,6)],
        bottom.margin = 12, midcolour = "black")
