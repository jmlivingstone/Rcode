concordance<-function(x,y){2 * cov(x,y,use="complete") /
(var(x,na.rm=TRUE) + var(y,na.rm=TRUE) +
(mean(x,na.rm=T)-mean(y,na.rm=T))^2 )}

concordanceM <- function(x){3
    n <- ncol(x)
    cm <- matrix(0,ncol=n,nrow=n)
    for(i in 1:n){
        print(i)
        for(j in i:n){
            cm[i,j] <- concordance(x[,i],x[,j])
            cm[j,i] <- cm[i,j]
        }
    }
    return(cm)
}

# example
values <- concordanceM(exprs)

