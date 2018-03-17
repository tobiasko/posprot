#code block to compare results for LS and WLS
#candidate 1 ####
#LS
poi <- names(model.list[464])
poi.subset <- (droplevels(CasTC_MasterPools_PeptideQuantDefault[CasTC_MasterPools_PeptideQuantDefault$EG.ModifiedSequence == poi,]))
str(poi.subset)
plot.linear.model(model.list[[464]], poi.subset)
#WLS
plot.linear.model(w.model.list[[464]], poi.subset)
#candidate 2 ####
poi <- names(model.list[538])
poi.subset <- (droplevels(CasTC_MasterPools_PeptideQuantDefault[CasTC_MasterPools_PeptideQuantDefault$EG.ModifiedSequence == poi,]))
str(poi.subset)
plot.linear.model(model.list[[538]], poi.subset)
#WLS
plot.linear.model(w.model.list[[538]], poi.subset)

#global picture for WLS models
#basic try to catch exception using class argument
get.model.p <- function(linear.model){
  f.stat <- summary(linear.model)$fstatistic
  if(class(f.stat) == "numeric"){
    p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail=FALSE)
    attributes(p) <- NULL
    return(p)
  }
  else return(NA)
}
w.model.quality <- ldply(w.model.list, function(x){c(rsqare = get.Rsquared(x), p.value=get.model.p(x))}, .progress = "text")
q.value <- p.adjust(w.model.quality$p.value, method = "BH")
w.model.quality <- transform(w.model.quality, q.value=q.value)
histogram(~p.value+q.value, data=w.model.quality, main="WLS")
good.w.model.index <- w.model.quality$q.value <= 0.05
summary(good.w.model.index)

w.model.2ndCoef <- ldply(w.model.list, failwith(NA, function(x){get.coef.1(x, 2)}), .progress = "text")
add.q <- function(df){
  names <- colnames(df)
  names[5] <- "p.value"
  colnames(df) <- names
  q.value <- p.adjust(df$p.value, method = "BH")
  return(cbind(df,q.value))
}
w.model.2ndCoef <- add.q(w.model.2ndCoef)
xyplot(-log(q.value)~Estimate, data = w.model.2ndCoef, main="WLS regression - slope estimation")

#select good models with positive slope, here abs(slope)>1
wls.candidates.pos <- na.omit(w.model.2ndCoef[good.w.model.index & (w.model.2ndCoef$Estimate > 0.75 & w.model.2ndCoef$q.value < 0.05),])
wls.candidates.neg <- na.omit(w.model.2ndCoef[good.w.model.index & (w.model.2ndCoef$Estimate < -0.75 & w.model.2ndCoef$q.value < 0.05),])

