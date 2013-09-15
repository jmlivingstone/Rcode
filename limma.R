#########
# cut patients into two groups based on expression ###
########

# vector of gene indices to use for clustering
m1 <- match(genes, GeneName)

#create expression matrix
e <- exprs[m1,]

#get the clustering indices
output <- sean.heatmap(e,rowInd=1:nrow(e),showplot=FALSE)

#create vector with the indices of the split
if(mean(e[,cutree(output$hcc,k=2)==1],na.rm=T) > mean(e[,cutree(output$hcc,k=2)==2],na.rm=T)){
        split <- ifelse(cutree(output$hcc,k=2)==1,"red","blue")
} else {
        split <- ifelse(cutree(output$hcc,k=2)==2,"red","blue")
}

#####

mostVar <- function(data) {
        m <- match(unique(data$GeneName), data$GeneName)
        e <- data$exprs[m,]
        ord <- order(apply(e, 1, var, na.rm=T), decreasing = T)
        var <- e[ord,]

#plot the top 500 most variable genes
        sean.heatmap(var[1:500,], midcolour = "black", collim = c(-2,2), main = "500.most.variable", clinical = data$heatmap.clinical)
}

###

limma <- function(exprs, variable) {
# targets = groups
        targets <- rep(NA, ncol(data$exprs))
        targets[which(variable == 1)] <- "GROUP1"
        targets[which(variable == 0)] <- "GROUP2"
        gtargets <- unique(targets)
        m <- match(targets,gtargets)
        design <- model.matrix(~ -1 + factor(m))
        colnames(design) <- gtargets
        contrast.matrix <- makeContrasts(GROUP1-GROUP2, levels = design)
        fit <- lmFit(exprs, design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit3 <- eBayes(fit2)
        top1 <- topTable(fit3,coef=1,number=nrow(fit3))
        m1 <- as.numeric(rownames(top1)[(which(top1$adj.P.Val < limit))])

#plot the top 500 most differential expressed genes              
	sean.heatmap(exprs[m1[1:500],], midcolour = "black", collim = c(-2,2))
        return(top1)
}
