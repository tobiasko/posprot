#testing robust regression
library("MASS", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
#read in peptide std. report from Spectronaut 9
#create subset for a specific peptide of ineterst
poi <- "_[+28]VSK[+28]PDLTAALR_"
poi.subset <- droplevels(all_in_edition[all_in_edition$EG.ModifiedSequence == poi,])

# fit lm, weighted lm and rlm
fmla <- as.formula("log2(FG.Quantity)~R.Condition + EG.PrecursorId")
poi.w.lm <- lm(fmla, weights = FG.Quantity, data = poi.subset)
summary(poi.w.lm)
poi.lm <- lm(fmla, data = poi.subset)
summary(poi.lm)
poi.r.lm <- rlm(fmla, data = poi.subset, method = "MM")
summary(poi.r.lm)

xyplot(log2(FG.Quantity)~R.Condition|EG.PrecursorId, panel = function(x,y){
  panel.xyplot(x, y)
  panel.abline(coef = coef(poi.lm), col = "blue")
  panel.abline(coef = coef(poi.w.lm), col = "red")
  panel.abline(coef = coef(poi.r.lm), col = "green")
}, data=poi.subset)

xyplot(log2(FG.Quantity)~R.Condition, groups = EG.PrecursorId, data = poi.subset)
xyplot(FG.Quantity~R.Condition, groups = EG.PrecursorId, data = poi.subset)


#Maybe generalized additive models would be better
library("mgcv")
poi.subset.3 <- droplevels(poi.subset[poi.subset$EG.PrecursorId == "_[+28]VSK[+28]PDLTAALR_.3",])
str(poi.subset.3)
poi.gam <- gam(log2(FG.Quantity)~R.Condition, data = poi.subset.3)
summary(poi.gam)
xv <- 0:5
yv <- predict(poi.gam, newdata = list(R.Condition=xv))
plot(log2(FG.Quantity)~R.Condition, data = (poi.subset.3))
lines(xv,yv)
#same-same
