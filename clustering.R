sean.heatmap(exprs[m,], collim=c(-2,2), midcolour = "black",
distfun=function(d) { dist(d, method="manhattan")}, hclustfun = function(d) { hclust(d, method = "complete")})

#library(survival)

#returns the survival coefficient
survivalAnalysis_Coeff <- function(e,time,event){
  s <- NA
  try( s <- coxph(Surv(time,event) ~ e))
  if(is.na(s)){
    return(NA)
  }
  else {
    return(s$coefficients)
  }
}

#returns the logtest pvalue
survivalAnalysis <- function(e,time,event){
  s <- NA
  try( s <- coxph(Surv(time,event) ~ e))
  if(is.na(s)){
    return(NA)
  }
  else {
    return(summary(s)$logtest[3])
  }
}


# plot the kaplan meier survival curve
orderedKaplanMeierPlot <- function(time,event,factors,main=""){
        p <- coxph(Surv(time,event)~factors)
        n <- summary(p)$n
        p2 <- summary(p)$logtest[3]
        hr <- summary(p)$coef[2]
        fit <- survfit(Surv(time,event)~factors,type="kaplan-meier" )
        plot(fit,
                col=c("red","blue"),
                xlab="TIME (MONTHS)",
                ylab="OUTCOME",
                cex.lab=1,
                main=main,
        )
	red = sum(factors=="red")
        blue = sum(factors=="blue")
        legend("bottomleft",
                c(paste("HR =",format(hr,digits=2)),
                paste("p =",format(p2,digits=2)),
                paste("BAD OUTCOME (n=",red,")"),
                paste("GOOD OUTCOME (n=",blue ,")")),
                col=c("white","white","red","blue"),
                lty=1,
                #lwd=3,
                bty="n")
return(p2)
}
