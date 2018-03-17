#libs ####
library("plyr")
library("lattice")
library("tidyr")
library("qvalue")
library("MASS")

# exploratory plotting functions for raw data ####
plot.raw <- function(df, ...){
  xyplot(FG.Quantity~R.Condition | EG.PrecursorId, main="Raw scale measurements", panel = function(x, y, ...){
    panel.xyplot(x, y, ...)
    panel.loess(x, y, ...)
    panel.lmline(x, y, ...)
  }, data = df, auto.key = TRUE)
}
plot.ylog <- function(df, ...){
  xyplot(log2(FG.Quantity)~R.Condition | EG.PrecursorId, panel = function(x, y, ...){
    panel.xyplot(x, y, ...)
    panel.loess(x, y, ...)
    panel.lmline(x, y, ...)
  }, data = df, auto.key = TRUE, main="log-transformed measurements")
}
plot.xylog <- function(df, v = 2, ...){
  if(v==1){
    xyplot(log2(FG.Quantity)~log2(R.Condition), groups = EG.PrecursorId, panel = function(x, y, ...){
      panel.xyplot(x, y, ...)
      panel.loess(x, y, ...)
      #panel.lmline(x, y, ...)
    }, data = df, auto.key = TRUE, main="log-log transformed data (grouped)")
  }
  else{
    xyplot(log2(FG.Quantity)~log2(R.Condition) | EG.PrecursorId, panel = function(x, y, ...){
      panel.xyplot(x, y, ...)
      panel.loess(x, y, ...)
      #panel.lmline(x,y, ...)
    }, data = df, auto.key = TRUE, main="log-log transformed data (conditioned)")
  }
}

# co-ploting functions for fitted models and data
plot.linear.model <- function(lm, df){
  #expects linear model object and df containing observations used for model fitting
  #limitations: uses only coef 1 and 2 for drawing model line
  xyplot(log2(FG.Quantity)~R.Condition|EG.PrecursorId, panel = function(x, y, ...){
    panel.xyplot(x, y)
    panel.abline(coef = coef(lm))
    }, data=df)
}
plot.loglog.linear.model <- function(lm, df, c =1){
  #expects linear model object and df containing observations used for model fitting
  #limitations: uses only coef 1 and 2 for drawing model line
  xyplot(log2(FG.Quantity)~log2(R.Condition+c)|EG.PrecursorId, panel = function(x, y, ...){
    panel.xyplot(x, y)
    panel.abline(coef = coef(lm))
  }, data=df)
}

#ploting functions for model parameter estimates
vplot <- function(df, fc = 1){
  #assumes df to contain parameter estimates and q-values and grouping var
  #red line is 5% FDR
  #green line is 1% FDR
  xyplot(-log10(q.value)~Estimate | groups, panel = function(x,y){
  panel.xyplot(x, y)
  panel.abline(v=c(-fc, fc))
  panel.abline(h=-log10(0.05), col="red")
  panel.abline(h=-log10(0.01), col="green")},
  data = df, main="Vulcano plot")
}
zplot <- function(x){
  #x is ecpected to be a ecdf instance
  plot(x, pch=1)
  abline(h=0.01, col="red")
  abline(h=0.99, col="blue")
  abline(h=0.5, col="black")
  abline(v=0, col="black")
}

#########################
# prepare input data ####
#########################

#read in Peptide Quant Default export from Spectronaut 9 ####
# pd workflow (3 pass SEQUESTHT search + Percolator validation) was executed using all DDA runs 
# Spectronaut 9 was used to query DIA maps and resulst were exported using default scheme
PepQuantDef <- read.csv("~/Documents/RStudio/p2095_CasTC/CasTC_all_in_edition_PeptideQuantDefault.csv", sep=";")
str(PepQuantDef)
#diagnostic plots
bwplot(log10(FG.Quantity)~R.FileName, data = PepQuantDef, main="Signal distribution per run")
histogram(~EG.Qvalue, data = PepQuantDef, main="Q-value distribution per runs")

# load pd 2.1 results export for Peptide Group level####
# spectronaut assays were build using pd results file and fasta DB used for SEQUESTHT searches
PeptideGroups <- read.delim("~/Documents/RStudio/p2095_CasTC/CasTC_3pass_analysis_PeptideGroups.txt")
str(PeptideGroups)

# load degrabase 1.0 ####
# was downloaded from degrabase website and saved as csv
degrabase <- read.csv("~/Documents/RStudio/p2095_CasTC/degrabase_071712.csv")
str(degrabase)
get.stripped.peptides <- function(x){
  y <- sub("^.*-", "", as.character(x), perl = TRUE)
  z <- sub("\\(.*\\)", "", y, perl=TRUE)
  return(z)
}
degrabase <- mutate(degrabase, stripped.peptide = get.stripped.peptides(degrabase$peptide))


####################
# model fitting ####
####################

# LS - ordinary linear regression using least squares estimation ####
# Version 1 : incl. charge as independent variable
fit.lm <- function(df){
  df <- droplevels(df)
  if (nlevels(df$EG.PrecursorId)>1) {
    lm(log2(FG.Quantity)~R.Condition + EG.PrecursorId, data = df)
  }
  else {
    lm(log2(FG.Quantity)~R.Condition, data = df)
  }
}
model.list <- dlply(.variables = "EG.ModifiedSequence", .fun = failwith(NA, fit.lm), .data = PepQuantDef, .parallel=FALSE, .progress = "text")

# Version 2 : fit a model for every precursor using split-apply-combine strategy
fit.lm <- function(df){
  df <- droplevels(df)
  return(lm(log2(FG.Quantity)~R.Condition, data = df))
}
model.list <- dlply(.variables = "EG.PrecursorId", .fun = failwith(NA, fit.lm), .data = PepQuantDef, .parallel=FALSE, .progress = "text")

head(model.list)
save(model.list, file = "lm_models.RData")

# WLS - estimation using weighted least squares according to signal intensity ####
fit.w.lm <- function(df){
  df <- droplevels(df)
  if (nlevels(df$EG.PrecursorId)>1) {
    lm(log2(FG.Quantity)~R.Condition + EG.PrecursorId, weights = FG.Quantity, data = df)
  }
  else {
    lm(log2(FG.Quantity)~R.Condition, weights = FG.Quantity, data = df)
  }
}
w.model.list <- dlply(.variables = "EG.ModifiedSequence", .fun = failwith(NA, fit.w.lm), .data = PepQuantDef, .parallel=FALSE, .progress = "text")
head(w.model.list)
save(w.model.list, file = "wlm_models.RData")

# RL - robust linear regression using M estimators ####
# not there yet

#log-log model ####
#log2 transformation of response and explanatory variable
#normaly used as a transformation for power functions : y = axˆb
# parameter a : slope
# parameter b: shape
#since log(0) is not def. -> log(0+c)

fit.loglog.lm <- function(df, c = 1){
  df <- droplevels(df)
  if (nlevels(df$EG.PrecursorId)>1) {
    lm(log2(FG.Quantity)~log2(I(R.Condition+c)) + EG.PrecursorId, data = df)
  }
  else {
    lm(log2(FG.Quantity)~log2(I(R.Condition+c)), data = df)
  }
}
loglog.model.list <- dlply(.variables = "EG.ModifiedSequence", .fun = failwith(NA, fit.loglog.lm), .data = PepQuantDef, .parallel=FALSE, .progress = "text")
save(loglog.model.list, file = "logloglm_models.RData")


###################
# test section ####
###################
iRT.protein <- c("LGGNEQVTRYILAGVENSKGTFIIDPGGVIRGTFIIDPAAVIRGAGSSEPVTGLDAKTPVISGGPYEYRVEATFGVDESNAKTPVITGAPYEYRDGLDAASYYAPVRADVTPADFSEWSKLFLQFGAQGSPFLK")
iRT <- gsub(x = iRT.protein, pattern = "R", replacement = "R,")
iRT <- gsub(x = iRT, pattern = "K", replacement = "K,")
iRT <- unlist(strsplit(iRT, split = ","))

# easy to expand to protein digest function
#digest.protein <- function(x, spec="RK"){
#  
#}

poi <- "_[+28]VSK[+28]PDLTAALR_" #two charge states
poi <- "_[+28]GQVITIGNER_" #one charge state

# select subset according to modified sequence
poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.ModifiedSequence == poi,])
summary(poi.subset)
# select subset according to precursor ID
poi <- "_[+28]GQVITIGNER_.2"
poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.PrecursorId == poi,])
# select according to stripped sequence
poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.StrippedSequence %in% iRT,])

plot.raw(poi.subset)
plot.xylog(poi.subset)
lm1 <- fit.lm(poi.subset)
lm2 <- fit.w.lm(poi.subset)
lm3 <- fit.loglog.lm(poi.subset)
summary(lm1)
summary(lm2)
summary(lm3)
#is the shape parameter of the log-log model sig different from 1?
t.value <- abs(0.38284-1)/0.07296
dt(t.value, df=1) #not sure about correct df

#compare the difference between log-log model and linear model in raw space
plot(FG.Quantity~R.Condition, data=poi.subset[poi.subset$EG.PrecursorId=="_[+28]VSK[+28]PDLTAALR_.2",])
plot(FG.Quantity~R.Condition, data=poi.subset[poi.subset$EG.PrecursorId=="_[+28]VSK[+28]PDLTAALR_.3",])
abline(lm(FG.Quantity~R.Condition, data=poi.subset[poi.subset$EG.PrecursorId=="_[+28]VSK[+28]PDLTAALR_.2",]))
abline(lm(FG.Quantity~R.Condition, data=poi.subset[poi.subset$EG.PrecursorId=="_[+28]VSK[+28]PDLTAALR_.3",]))
xv <- seq(0,5,0.1)
a <- coef(lm3)[1]
a <- 2^a  
b <- coef(lm3)[2]
yv <- a*xv^b 
lines(xv, yv, col="red")

##########################
# Examine the results ####
##########################

# helper functions
get.R2 <- function(linear.model){summary(linear.model)$r.squared}
get.model.p <- function(linear.model){
  f.stat <- summary(linear.model)$fstatistic
  if(class(f.stat) == "numeric"){
    p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail=FALSE)
    attributes(p) <- NULL
    return(p)
  }
  else return(NA)
}
get.coef <- function(linear.model, n){
  if(class(linear.model) != "lm"){
    warning("First argument is not an instance of class lm!")
    return(c(NA, NA, NA, NA))
  }
  else {
    s <- summary(linear.model)$coefficients
    n.coef <- c(NA, NA, NA, NA)
    tryCatch(expr = n.coef <- s[n,], finally = return(n.coef))
  }
}

# compile model quality table with model R2 and model p and adjusted p (here q for simplicity) ####
model.quality <- ldply(model.list, function(x){c(R2 = get.R2(x), p.value=get.model.p(x))}, .progress = "text")
xyplot(R2~p.value, data = model.quality, main="Model p-value vs. R2")
# choose multiple testing adjustment method by # usage
#model.quality <- transform(model.quality, q.value=qvalue(model.quality$p.value)$qvalues)
model.quality <- transform(model.quality, q.value=p.adjust(model.quality$p.value, method = "BH"))
histogram(~p.value + q.value, data=model.quality, main="Model p-value distribution before and after multiple testing adjustment", breaks = 50)
summary(model.quality)
good.model.index <- model.quality$q.value <= 0.05 #5% FDR

# compile table with 2nd coef (slope) plus q-value ####
model.2ndCoef <- ldply(model.list, failwith(NA, function(x){get.coef(x, 2)}), .progress = "text")
str(model.2ndCoef)
#choose multiple testing adjustment
add.q <- function(df){
  names <- colnames(df)
  names[5] <- "p.value"
  colnames(df) <- names
  return(transform(df, q.value = p.adjust(df$p.value, method = "BH")))
}
add.q <- function(df){
  names <- colnames(df)
  names[5] <- "p.value"
  colnames(df) <- names
  q <- qvalue(df$p.value)
  return(cbind(df, q.value=q$qvalues))
}
model.2ndCoef <- add.q(model.2ndCoef)
summary(model.2ndCoef)

# Disgnostic plots ####
# multiple testing adjustment
histogram(~q.value + p.value, data = model.2ndCoef, main="Slope p & q-value distribution", breaks = 50)
#vulcano plot
#vulcano plots for good and bad models and fc cutoff as used for candidate list
vplot(transform(model.2ndCoef, groups = as.character(good.model.index)), fc = 0.8)
# marginal distributions
histogram(~q.value, data = model.2ndCoef, main="slope q-value distribution", breaks = 50)
histogram(~Estimate, data = model.2ndCoef, main="slope effect size distribution", breaks = 50)

#############################
# compile candidate list ####
#############################
#example: select good models with abs(slope)>0.8 at 1% sig. level (look at vulcano plot)
candidates.pos <- na.omit(model.2ndCoef[good.model.index & (model.2ndCoef$Estimate > 0.8 & model.2ndCoef$q.value < 0.01),])
candidates.neg <- na.omit(model.2ndCoef[good.model.index & (model.2ndCoef$Estimate < -0.8 & model.2ndCoef$q.value < 0.01),])

#########################################################
# estimation of cutoff from the data using a 0-model ####
#########################################################
# workflow
# 1: select a set of peptides detected by DDA that should not show a dependence on time
# 2: examine the distribution of estimated coefs

#instead of using a fixed cutoff one can determine a 0 distribution
# set1 : natural N-term. peptides
# Problem : split variable has changed -> name of feature is now "EG.PrecursorId"
# solution could be: rename column to generic name: model.Id or use column number for reference
model.quality <- rename(model.quality, c("EG.PrecursorId" = "ModelId"))
Met.removed <- grepl("_[-131]", model.quality$ModelId, fixed=TRUE)
Nterm.acet <- grepl("_[+42]", model.quality$ModelId, fixed=TRUE)
Met.removed.acet <- grepl("_[-89]", model.quality$ModelId, fixed=TRUE)
natural <- Met.removed | Nterm.acet | Met.removed.acet

#set2: unmodified, semitryptic peptides without internal D
semi.tryptic <- grepl("R_", model.quality$ModelId)
free.nterm <- !grepl("_\\[", model.quality$ModelId)
noD <- !grepl("D", model.quality$ModelId)
trp <- semi.tryptic & free.nterm & noD

#diagnostic plots ####
#group TRUE is semi-tryptic, free peptides without internal Ds 
vplot(transform(model.2ndCoef, groups=as.character(trp)))
#marginal x distribution
densityplot(~Estimate, groups = groups, data = transform(model.2ndCoef, groups = as.character(trp)), auto.key = TRUE)
#compute empirical density function of zero distribution of slope effect size
zero.dist.ecdf <- ecdf(model.2ndCoef$Estimate[trp])
zplot(zero.dist.ecdf)

#compute empirical density function of zero distribution of slope p-value
zero.dist.ecdf.pvalue <- ecdf(model.2ndCoef$p.value[trp])
plot(zero.dist.ecdf.pvalue, pch = 1, main="slope p.value ECDF")
abline(a=0, b=1) #for zero we expect uniform distribution of p-value = diagonal
(quantile(zero.dist.ecdf.pvalue, 0.01))

#compute empirical density function of zero distribution of slope q-value
zero.dist.ecdf.qvalue <- ecdf(model.2ndCoef$q.value[trp])
plot(zero.dist.ecdf.qvalue, pch = 1, main="slope q-value ECDF")
abline(a=0, b=1) #for zero we expect uniform distribution of q-value
(quantile(zero.dist.ecdf.qvalue, 0.01))

# Determine the critical slope value from the data distribution ####
n <- length(model.list)
(n*0.01) #expected number of false positives at 1% FDR
(cutoff <- quantile(zero.dist.ecdf, 0.99)) #critical value for slope at 1% FDR
vplot(transform(model.2ndCoef, groups = as.character(good.model.index)), fc = cutoff)
#override with new candidate list using 0-model
candidates.pos <- na.omit(model.2ndCoef[good.model.index & (model.2ndCoef$Estimate > cutoff & model.2ndCoef$q.value < 0.01),])
candidates.pos.nterm <- na.omit(model.2ndCoef[good.model.index & (model.2ndCoef$Estimate > cutoff & model.2ndCoef$q.value < 0.01) & grepl("^_\\[\\+28]", model.2ndCoef$EG.ModifiedSequence, perl = TRUE),])
#write candidate lists to hard drive
write.csv(candidates.pos, file = "candidate_pos.csv", row.names = FALSE, quote = FALSE)


##########################################
# 1-Model - true cases from degrabase ####
##########################################
# workflow:
# 1. select peptides iden. in DDA runs that start after D and are dimethylated at the N-term
# 2. retrieve corresponding protein acc and position of peptide within protein
# 3. query degrabase based on protein acc and position (implemented by SQL left join)
# 4. query list of models and raw data - How do known Cas-dependent peptides behave? 

# helper functions
get.acc <- function(x){
  acc.pattern <- "^.*\\s"
  grep(pattern = acc.pattern, x=x, perl=TRUE)
  m <- regexpr(pattern = acc.pattern, x)
  n <- sub(pattern = " ", "", regmatches(x=x, m))
  return(n)
}
get.start <- function(x){
  start.pattern <- "\\[.*-"
  m <- regexpr(pattern = start.pattern, x)
  n <- sub("\\[","", sub(pattern = "-", "", regmatches(x=x, m)))
  return(as.numeric(n))
}
get.interval <- function(x){
  #expects x to be a character, returns numeric vector x1=start, x2=stop
  interval.pattern <- "\\[\\d*-\\d*\\]"
  m <- regexpr(pattern = interval.pattern, x)
  n <- strsplit(regmatches(x=x, m),split = "-")
  return(as.numeric(sub(pattern = "\\[|\\]", replacement = "", unlist(n))))
}

# select all peptide groups indentified in DDA starting after D and having N-term. dimethyl mark
# queries pd 2.1 peptide groups export
after.D.pattern <- "^\\[D\\]\\."
nterm.dim.pattern <- "Dimethyl\\s\\[N-Term"
subset <- PeptideGroups[grepl(after.D.pattern, PeptideGroups$Annotated.Sequence, perl = TRUE) & grepl(nterm.dim.pattern, PeptideGroups$Modifications, perl = TRUE),]
subset <- separate_rows(data = subset[,c(3:4,11,13)], Positions.in.Master.Proteins, sep = ";\\s") #ungroups
subset <- mutate(subset, acc = get.acc(subset$Positions.in.Master.Proteins), start = get.start(subset$Positions.in.Master.Proteins)) #extracts acc and start position of peptide
#SCL-type left join
joined.subset <- droplevels(join(x = subset, y = degrabase, by=c("acc","start"), type="left"))
summary(joined.subset)

# Diagnostic plots ####
histogram(~P4 + P3 + P2 + P1 + P1. + P2. + P3. + P4., joined.subset)
summary(is.na(joined.subset)) #check for missing values

# filter joined data for relevant experimental context ####
# CasTC exp has been performed using Jurkat cells and Stau. treatment
annotated <- !is.na(joined.subset$uniprot.ID) #select all features (queries) that are known to degrabase
jurkat <- joined.subset$cell.type  %in% c("Jurkat", "MM1S, DB, Jurkat") #select all fatures that are known in the context of Jurkat cells
staurosporine <- joined.subset$perturbagen %in% c("Staurosporine (2uM) or TRAIL (2nM), extracts combined", "Staurosporine, Bortezomib, Doxorubicin") #select all features detected in Staurosporine treatent
known <- as.character(unique(joined.subset[annotated & jurkat & staurosporine, "Annotated.Sequence"]))

#helper functions
# intended use : each model is named according to split factor (e.g. EG.PrecursorId), 
# version 1 : EG.ModifiedSequence gives pattern: _XX[mod]XX_
# version 2 : EG.PrecursorId gives pattern: _XX[mod]XX_.z
isolate.peptide <- function(x){
  y <- sub(x, pattern = "^\\[.\\].", replacement = "", perl = TRUE)
  z <- sub(y, pattern = "\\.\\[.\\]$", replacement = "", perl = TRUE)
  return(z)
}
modify.peptide <- function(x){
  y <- sub(x, pattern = "M", replacement = "M[+16]", perl = TRUE)
  z <- sub(y, pattern = "K", replacement = "M[+28]", perl = TRUE)
  return(z)
}
strip.peptide <- function(x, aa.only =FALSE){
  #strips off modifications
  if(aa.only){
    z <- gsub(x, pattern = "(!?[ABCDEFGHIJKLMNOPQRSTUVWYZX])", replacement = "", perl = TRUE)
  }
  else{
    y <- gsub(x, pattern = "\\[(\\+||\\-)\\d+\\]", replacement = "", perl = TRUE)
    z <- gsub(y, pattern = "_||_\\.\\d", replacement = "", perl = TRUE)
  }
  return(z)
}
#test statement
grep(pattern = isolate.peptide(known[1]), x = strip.peptide(names(model.list)), value = TRUE)

# See how known Cas-dependent peptides behave in CasTC exp ####
length(known)
i <- 30
known[i] #Peptide known to degrabase in the context of Jurkat cells treated with Stau.
l <- model.list[grep(pattern = isolate.peptide(known[i]), x = strip.peptide(names(model.list)))] #corresponding linear model, based on stripped sequence
# Model IDs
names(l)
# How does the raw data over these models look like?
#poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.ModifiedSequence %in% names(l), ])
poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.PrecursorId %in% names(l), ]) #substrings also match!
plot.raw(poi.subset)
plot.ylog(poi.subset)
plot.xylog(poi.subset, v = 2)

# Look at raw data for a specific model
j <- 14
#poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.ModifiedSequence == names(l)[j], ])
poi.subset <- droplevels(PepQuantDef[PepQuantDef$EG.PrecursorId == names(l)[j], ])
plot.raw(poi.subset)
plot.ylog(poi.subset)

#next step
#generate QC plots over a list of features to determine distribution of estimated coef

# grep protein group information for candidate ####
(pep <- as.character(candidates.pos.nterm$EG.ModifiedSequence[7]))
p <- strip.ends(strip.mods(pep))
grep(pattern = p, x = CasTC_3pass_analysis_PeptideGroups$Annotated.Sequence, fixed = TRUE, value = TRUE) #gets also substring matches :-)
CasTC_3pass_analysis_PeptideGroups[grep(pattern = p, x = CasTC_3pass_analysis_PeptideGroups$Annotated.Sequence, fixed = TRUE),]


