# naive bayes predictor #
samples <- colnames(exprs)
# binary vector of recurrence/survival
recurrence.vec <- variable
names(recurrence.vec) <- colnames(exprs)
classvector<-as.factor(recurrence.vec.wang[samples])
trainarray <- t(exprs[best.pred,names(recurrence.vec[samples])])

predictor <-naiveBayes(trainarray, classvector)
res <- predict(predictor, trainarray,type = 'class' )
prob <- predict(predictor, trainarray,type = 'raw' )

PREDICTOR <- ifelse(res1 == 1, "red","blue")
o <- order(prob[,1])


generate.test.predictor <- function(uniq.genes.inds # inds of ordered unique univariate predictors
                                ,tau # number of genes to select the random sets from
                                ,pred.size #predictor size
                                ,num.errors #number of allowed errors
                                ,cpc.max # number of predictors generated
                                ,samples # samples to select from (e.g only ER postivie samples)
                                ,recurrence.vec= recurrence.vec # recurrence information
                                ,exprdata=data.set$exprdata # exprs data
                                ,pos.var="Yes",neg.var="No"
                                ,subtype="all subtypes"
                                        )
{

        tmp1=names(recurrence.vec[samples])
        tmp2=colnames(exprdata[, samples])
        if(sum(tmp1==tmp2)!= length(recurrence.vec[samples])) {
                stop("ordering of the patients in the recurrence vector and expression data is not the same")
        }
	if(length(which(is.na(recurrence.vec[samples])))>0) {
                stop("fix NA values in recurrence vector")
        }

	classvector.all<-as.factor(recurrence.vec[samples])

        prob = matrix(data=0, nrow=length(samples), ncol=2)
        colnames(prob) = c(pos.var,neg.var)
        rownames(prob) = samples
        res.list = list()
        res=c()
        num.runs=0
        cpc.counter=1

        library(e1071)
# maximum number of runs set to 100000 which seems to be what people report in publications
#        while( (cpc.counter<= cpc.max)&&(num.runs<= cpc.max*5)) {
        while(cpc.counter<= cpc.max) {
                num.runs=num.runs+1
                cat(paste("\n cpc.counter:",cpc.counter," num.runs:",num.runs,"in subtype",subtype))
                # sample from the top tau genes, for the predictor of size
                currentpredictor <- sample(uniq.genes.inds[1:tau],pred.size, replace=FALSE)
                trainarray <- t(exprdata[currentpredictor,names(recurrence.vec[samples])])
                #as.matrix(exprdata[currentpredictor,names(recurrence.vec)],nrow=length(currentpredictor),
                ncol=length(recurrence.vec)
                trainarray = data.frame(trainarray)
                FP= 0
                FN= 0
                TP=0
                TN=0
                totalerrors = 0
                current.patient.error= matrix(0,nrow=1,ncol=length(recurrence.vec[samples]))
                colnames(current.patient.error)=samples
                for (ii in 1:length(recurrence.vec[samples])) {
                        looguy = ii
                        predictor1 <-naiveBayes(trainarray[-looguy,], as.factor(classvector.all[-looguy]) )
                        res1 <- predict(predictor1, trainarray[looguy,],type = 'class' )
                        prob[ii,] <- predict(predictor1, trainarray[looguy,],type = 'raw' )
                        if ((res1==pos.var) & (classvector.all[looguy]==neg.var)) {
                                FP <- FP + 1
                                current.patient.error[1,looguy]=1
                        }
                        if ((res1==neg.var) & (classvector.all[looguy]==pos.var)) {
                                FN <- FN + 1
                                current.patient.error[1,looguy]=1
                        }
                        if ((res1==pos.var) & (classvector.all[looguy]==pos.var)) {
                                TP <- TP + 1
                        }
                        if ((res1==neg.var) & (classvector.all[looguy]==neg.var)) {
                                TN <- TN + 1
                        }
                } # end of ii (loosteps)
                totalerrors <- FN + FP
                if (totalerrors <= num.errors){
                        test.tbl = matrix(nrow=2,ncol=2,data=0)
                        test.tbl[1,1] = TP
                        test.tbl[1,2]=  FP
                        test.tbl[2,1]=  FN
                        test.tbl[2,2] = TN
                        tmp = fisher.test(test.tbl, simulate.p.value=FALSE)
                        fishers.pvalue = tmp$p.value
                        res$currentpredictor = currentpredictor
                        res$FP= FP
                        res$FN= FN
                        res$TN= TN
                        res$TP= TP
                        res$totalerrors = totalerrors
                        res$current.patient.error = current.patient.error
                        res$fishers.pvalue = fishers.pvalue
                        res$prob = prob
                        res.list[[cpc.counter]]=res
                        cpc.counter=cpc.counter+1
                } #if (totalerrors<= num.errors){

                }

        return(list(result=res.list,num.runs=num.runs))

}

###################################
### example of how to run function ###

meta.pred = generate.test.predictor(uniq.genes.inds=uniq.ind.ord,
tau=100,pred.size=25,num.errors=50,cpc.max=100,
samples=samples,recurrence.vec= recurrence.vec,
exprdata=data$exprs,pos.var=1,neg.var=0,subtype="all subtypes")


getOutcome <- function(data, recurrence.vec){
        for(i in 1:length(data$result)) {
                b <- which(data$result[[i]]$current.patient.error == 1)
                outcome <- recurrence.vec
                for(j in 1:length(b)) {
                        if(outcome[b][j] == 1) {
                                outcome[b][j] <- 0
                        }
                        else {
                              	outcome[b][j] <- 1
                        }
                }
                data$result[[i]]$outcome.colour <- data.frame(ifelse(outcome == 1, "red", "blue"), stringsAsFactors=F)
                names(data$result[[i]]$outcome.colour) <- "PREDICTOR"
                data$result[[i]]$outcome <- outcome
        }
return(data)
}

####### TEST in another dataset ###
test.predictor <- function(currentpredictor = currentpredictor,samples = samples,
                        recurrence.vec= recurrence.vec,
                        exprdata=exprdata,
                        pos.var="yes",
                        neg.var="no")
        {
        tmp1=names(recurrence.vec[samples])
        tmp2=colnames(exprdata[, samples])

        classvector.all<-as.factor(recurrence.vec[samples])

        #samples=c(1:length(recurrence.vec))
        prob = matrix(data=0, nrow=length(samples), ncol=2)

        rownames(prob) = samples
        res = c()
        library(e1071)

        tmp.length=length(currentpredictor)
        if (tmp.length>1){
                trainarray <- t(exprdata[currentpredictor,samples])
                #as.matrix(exprdata[currentpredictor,names(recurrence.vec)],nrow=length(currentpredictor), ncol=length(recurrence.vec))
                trainarray = data.frame(trainarray)
                FP <- 0
                FN <- 0
                TP<-0
                TN<-0
                totalerrors <- 0
                current.patient.error<- matrix(0,nrow=1,ncol=length(recurrence.vec[samples]))
                for (ii in 1:length(recurrence.vec[samples])) {
                        looguy <- ii
                        #predictor1 <-naiveBayes(trainarray, as.factor(classvector.all) )
                        predictor1 <-naiveBayes(trainarray[-looguy,], as.factor(classvector.all[-looguy]) )
                        res1 <- predict(predictor1, trainarray[looguy,],type = 'class' )
                        tmp = predict(predictor1, trainarray[1,],type = 'raw' )
                        colnames(prob) = colnames(tmp)

                        prob[ii,] <- predict(predictor1, trainarray[looguy,],type = 'raw' )
                        if ((res1==pos.var) & (classvector.all[looguy]==neg.var)) {
                                FP <- FP + 1
                                current.patient.error[1,looguy]=1
                        }
                        if ((res1==neg.var) & (classvector.all[looguy]==pos.var)) {
                                FN <- FN + 1
                                current.patient.error[1,looguy]=1
                        }
                        if ((res1==pos.var) & (classvector.all[looguy]==pos.var)) {
                                TP <- TP + 1
                        }
                        if ((res1==neg.var) & (classvector.all[looguy]==neg.var)) {
                                TN <- TN + 1
                        }
                } # end of ii (loosteps)
                totalerrors <- FN + FP
                test.tbl = matrix(nrow=2,ncol=2,data=0)
                test.tbl[1,1] = TP
                test.tbl[1,2]=  FP
                test.tbl[2,1]=  FN
                test.tbl[2,2] = TN
                tmp = fisher.test(test.tbl, simulate.p.value=FALSE)
                fishers.pvalue = tmp$p.value
                res$FP=FP
                res$FN= FN
                res$TN= TN
                res$TP= TP
                res$totalerrors = totalerrors
                res$current.patient.error = current.patient.error
                res$fishers.pvalue = fishers.pvalue
                res$prob = prob
                res$current.predictor = currentpredictor
        }
return(res)

}

### fishers exact test
for(i in 1:5000) {
        test.tbl = matrix(nrow=2,ncol=2,data=0)
        test.tbl[1,1] = results[[i]]$TP + results2[[i]]$TP
        test.tbl[1,2]=  results[[i]]$FP + results2[[i]]$FP
        test.tbl[2,1]=  results[[i]]$FN + results2[[i]]$FN
        test.tbl[2,2] = results[[i]]$TN + results2[[i]]$TN
        tmp = fisher.test(test.tbl, simulate.p.value=FALSE)

        fishers.pvalue.luminal[i] = tmp$p.value
}
