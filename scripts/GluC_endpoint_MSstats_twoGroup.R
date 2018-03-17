###########
# info ####
###########
# this script generates a peptide-centric analysis of the GluC-endpoint experiment using MSstats


######################
# install package ####
######################
# MSstats package install
install.packages(c("gplots","lme4","ggplot2","ggrepel","reshape","reshape2","data.table","Rcpp","survival","minpack.lm"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))
install.packages("~/Downloads/MSstats_3.7.1.tar.gz", repos = NULL, type = "source")
library("MSstats")
vignette("MSstats")

################
# load libs ####
################
# MSstats
library("MSstats")
# plyr
library("plyr")
# tidy data
library("tidyr")
# lattice
library("lattice")
# dplyr
library("dplyr")

##################
# data import ####
##################
# this script was written using MSstats 3.6
# data export was done using custom scheme in Spectronaut 10
MSstats.report <- read.csv("GluC_endpoint_with_protein_inference_MSstatsReport.csv", sep=";") #read csv file exported from Spectronaut 9
names <- colnames(SRMRawData)

#QC
str(MSstats.report)
head(MSstats.report)
tail(MSstats.report)

#reshape report according to MSstats conventions
FragmentIon <- paste(MSstats.report$F.FrgIon, MSstats.report$F.FrgLossType, sep = "") #problem: FragIon is not uniqued without FrgLossType
df <- ldply(strsplit(as.character(MSstats.report$EG.PrecursorId),".", fixed = TRUE)) #problem: MSstats wants precuorsor charge in a sep column
IsotopeLabelType <- "L" #label-free DIA experiment, no SIL peptide spike-in
#default: run.id us used as BioReplicate id
input.df <- data.frame(ProteinName = MSstats.report$PG.ProteinNames, PeptideSequence = df$V1, PrecursorCharge = df$V2, FragmentIon, ProductCharge = MSstats.report$F.FrgZ, IsotopeLabelType, Condition = MSstats.report$Condition, BioReplicate = MSstats.report$R.FileName, Run = MSstats.report$R.FileName, Intensity = MSstats.report$F.PeakArea)

# Problem: "." is also used for mod description XXX[+28.0]XXX, so either on replaces .d+ by " or "." can not be used as split character

#####################
# data import II ####
#####################
# new custom report in Spectronaut 10
# does not provide all the information in separate columns
# import report
MSstats.report <- read.csv("GluC_endpoint_inclProtInf_MSstatsReportII.csv", sep = ";")
str(MSstats.report)
names <- colnames(SRMRawData) # column names for 10 column format
# issue 1 - Condition and BioReplicate 
# Condition and bioreplicate are taken from fields in the file name, "_" is split character, gives 4 fields (V1 to V4)
MSstats.report <- separate(data = MSstats.report, col = R.FileName, into = c("V1", "V2", "V3", "V4"), remove = FALSE)
# issue 2 - .zero
MSstats.report$FG.Id <- gsub(pattern = ".0", replacement = "", MSstats.report$FG.Id)
# issue 3 - mod. peptide sequence is missing for some features, create from FG.Id
mod.pep <- gsub(pattern = "\\.\\d$|_", replacement = "", MSstats.report$FG.Id, perl = TRUE)

# next lines generates the 10 column format for peptide-level modeling
pepitdeLevel.input <- data.frame(
  ProteinName = mod.pep,
  PeptideSequence = mod.pep,
  PrecursorCharge = as.factor(MSstats.report$PrecursorCharge),
  FragmentIon=paste(MSstats.report$F.FrgIon, MSstats.report$F.FrgLossType, sep = ""),
  ProductCharge = as.factor(MSstats.report$F.FrgZ),
  IsotopeLabelType = "L",
  Condition = substr(MSstats.report$V4, 1, 1), #condition is 1th char in V4
  BioReplicate = MSstats.report$V4, #replicate number
  Run = MSstats.report$R.FileName,
  Intensity = MSstats.report$F.PeakArea
)
str(pepitdeLevel.input)
summary(pepitdeLevel.input)

######################
# data import III ####
######################
#crazy shit - finally the Vitec lab started talking to the Biognosys people. A new export scheme emerged!
MSstats.report <- read.csv("GluC_endpoint_inclProtInf_ReportIII.csv", sep = ";")
str(MSstats.report)
summary(MSstats.report)
histogram(~EG.Qvalue, data = MSstats.report, breaks = 50)
# control is not the reference level in R.Condition, need to relevel
MSstats.report$R.Condition <- relevel(MSstats.report$R.Condition, ref = "control")
str(MSstats.report)

# pre-process input
input <- SpectronauttoMSstatsFormat(MSstats.report, intensity = "PeakArea", qvalue_cutoff = 0.01)
summary(input)
quant <- dataProcess(raw = input, censoredInt = "0")

#QC
dataProcessPlots(data = quant, type = "QCPlot", which.Protein = c("Q9QZF2","E9Q616"), address = FALSE)
# protein summary and normalization looks pretty bad!
dataProcessPlots(data = quant, type = "ProfilePlot", which.Protein = c("Q9QZF2","E9Q616"), address = FALSE)
dataProcessPlots(data = quant, type = "ConditionPlot", which.Protein = "E9Q616", address = FALSE)

# two-group comparision
levels(quant$ProcessedData$GROUP_ORIGINAL)
c1 <- matrix(c(-1,1), nrow = 1)
c <- c1
row.names(c) <- c("control vs. GluC") 
GluC.endpoint <- groupComparison(contrast.matrix = c, data = quant)

######################
# data import IV ####
######################
# imports all available peptide features, not only proteotypic
MSstats.report <- read.csv("GluC_endpoint_inclProtInf_ReportIV.csv", sep = ";")
str(MSstats.report)
summary(MSstats.report)
# control is not the reference level in R.Condition, need to relevel
MSstats.report$R.Condition <- relevel(MSstats.report$R.Condition, ref = "control")

#check signal distribution in report
bwplot(log2(F.PeakHeight)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Signal distribution per run")
bwplot(log2(F.NormalizedPeakHeight)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Signal distribution per run")
bwplot(log2(F.NormalizedPeakHeight)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Signal distribution per run")
bwplot(log2(F.PeakArea)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Signal distribution per run")
bwplot(log2(F.NormalizedPeakArea)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Signal distribution per run")
histogram(~log2(F.NormalizedPeakArea)|R.FileName, data = MSstats.report, breaks = 25)

# check signal quality
bwplot(-log10(EG.Qvalue)~R.FileName, data = MSstats.report, scales=list(x=list(rot=45)), main="Q-value distribution")
(-log10(0.01))
histogram(~EG.Qvalue|R.FileName, data = MSstats.report, main="Q-value distribution per run")

#How many Q-values are above 0.01?
dotplot(MSstats.report[MSstats.report$EG.Qvalue > 0.01, "R.FileName"], main="#Q-values above 0.01")

#######################
# pre-process IV.I ####
#######################
# first try
input <- SpectronauttoMSstatsFormat(MSstats.report.proteotypic, intensity = "PeakArea", qvalue_cutoff = 0.01)
summary(input)
input$ProteinName <- input$PeptideSequence #peptide level summary
# check signal distribution per run before summary
bwplot(log2(Intensity)~Run | Condition, data = input)
densityplot(~log2(Intensity)|Run, groups = Condition, data = input)


########################
# pre-process IV.II ####
########################
# 2nd try - use normalization computed by SC10?
input <- MSstats::SpectronauttoMSstatsFormat(input = MSstats.report, intensity = "PeakArea", filter_with_Qvalue = 0.01, useUniquePeptide = FALSE, fewMeasurements = FALSE)
norm.input <- SpectronauttoMSstatsFormat(input = MSstats.report, intensity = "NormalizedPeakArea", filter_with_Qvalue = 0.01, useUniquePeptide = FALSE, fewMeasurements = FALSE)
summary(input)
summary(norm.input)
# Intensity distribution in pre-processed data
bwplot(log2(Intensity)~Run, data = input, scales=list(x=list(rot=45)), main="Signal distribution per run - raw")
bwplot(log2(Intensity)~Run, data = norm.input, scales=list(x=list(rot=45)), main="Signal distribution per run - SC10 norm.")

#####################
# summarize data ####
#####################
#Protein-level summary
protein.quant <- dataProcess(raw = norm.input, censoredInt = "0", normalization = FALSE)
save(protein.quant, file = "protein_quant.RData")

# Peptide-level summary
norm.input$ProteinName <- norm.input$PeptideSequence
peptide.quant <- dataProcess(raw = norm.input, censoredInt = "0", normalization = FALSE, n_top_feature = 5, remove_proteins_with_interference = FALSE)
save(peptide.quant, file = "peptide_quant.RData")

# QC protein-level
#check for the first 10 proteins
poi <- as.character(levels(protein.quant$ProcessedData$PROTEIN))[1:10]
dataProcessPlots(data = protein.quant, type = "QCPlot", which.Protein = poi, address = FALSE)
dataProcessPlots(data = protein.quant, type = "ProfilePlot", which.Protein = poi, address = FALSE)

# QC peptide-level
#check for the first 10 peptides
poi <- as.character(levels(peptide.quant$ProcessedData$PROTEIN))[1:10]
dataProcessPlots(data = peptide.quant, type = "QCPlot", which.Protein = poi, address = FALSE)
dataProcessPlots(data = peptide.quant, type = "ProfilePlot", which.Protein = poi, address = FALSE)


############################
# two-group comparision ####
############################
#contrast matrix
levels(quant$ProcessedData$GROUP_ORIGINAL)
c1 <- matrix(c(-1,1), nrow = 1)
c <- c1
row.names(c) <- c("control vs. GluC") 

#Protein-level
GluC.endpoint.protein <- groupComparison(contrast.matrix = c, data = protein.quant)
groupComparisonPlots(data = GluC.endpoint.protein$ComparisonResult, type = "VolcanoPlot", ProteinName = FALSE, address = "Protein-level_")

#Peptide-level
GluC.endpoint.peptide <- groupComparison(contrast.matrix = c, data = peptide.quant)
save(GluC.endpoint.peptide, file = "GluC_endpoint_peptide_comp.RData")

##############################
# check the comp. results ####
##############################
# output from MSstats group comparision
load("GluC_endpoint_peptide_comp.RData")

#test section ####
test.case.1 <- "_A[+28.0]AAATAGK[+28.0]MMTVR_" #nterm & internal label
test.case.2 <- "_AAAATAGK[+28.0]MMTVR_" # internal label
test.case.3 <- "_RAAATAGK[+28.0]MMTVR_" # ArgC
strip.ends(test.case.1)
strip.mods(test.case.2)
grep(pattern = nterm.label, x = strip.ends(test.case.1))
grep(pattern = nterm.label, x = strip.ends(test.case.2))
grep(pattern = cterm.R, x = strip.ends(test.case.1))
grep(pattern = nterm.R, x = strip.ends(test.case.2))

# make groups ####
# this should finally go into a grouping variable
index.1 <- grepl(pattern = nterm.label, x = strip.ends(GluC.endpoint.peptide$ComparisonResult$Protein))
index.2 <- grepl(pattern = cterm.R, x = strip.ends(GluC.endpoint.peptide$ComparisonResult$Protein))
index <- index.1 & index.2

#GluC.endpoint.peptide$ComparisionResult needs to be annotated
# get upstream and downstream sequence from protein DB

# add grouping variable
# need to match all peptides versus ref proteome

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
biocLite("BiocGenerics")
library("BiocGenerics")
library("plyr")


## sample from quant results table
## n <- 50
## s <- sample_n(GluC.endpoint.peptide$ComparisonResult, n)

## full data
## MSstats ComparisionResult is expected to be from a peptite-level analysis !!!
## Protein contains mod. peptide sequence
s <- GluC.endpoint.peptide$ComparisonResult #extract MSstats comparision results
simple.mods <- function(x){
  y <- gsub(pattern = "\\.0+", x = x, replacement = "")
  return(y)
}
strip.ends <- function(x){
  x <- as.character(x)
  y <- gsub(pattern = "^_|_$", x = x, replacement = "")
  return(y)
}
strip.mods <- function(x){
  #any.mod <- "\\[(\\+|\\-)\\d+\\.\\d+]"
  #any.mod <- "\\[(\\+|-)\\d+\\]"
  any.mod <- "\\[(\\+|-)\\d+\\.?\\d+\\]"
  y <- gsub(pattern = any.mod, x = x, replacement = "")
  return(y)
}
s <- mutate(s, s.pep = strip.mods(strip.ends(Protein)), m.pep = strip.ends(Protein))

## The Magic
## use AhoCorasick to search for perfect matches
library("AhoCorasickTrie")

## laod reference proteome
ref.proteome <- readAAStringSet(filepath = "uniprot-proteome_UP000000589.fasta")
proteins <- BiocGenerics::sapply(ref.proteome, toString) # full reference proteome as string vector
system.time(l <- AhoCorasickSearch(keywords = s$s.pep, text = proteins, alphabet = "aminoacid", groupByKeyword = TRUE))
#system.time(m <- AhoCorasickSearch(keywords = s$s.pep, text = proteins, alphabet = "aminoacid", groupByKeyword = FALSE))


# -2- use pattern matching functions from Biostrings class
# foo <- function(x){vmatchPattern(pattern = x, subject = ref.proteome)} #one-to-many matching
# system.time(m <- BiocGenerics::lapply(peptides, foo))
# works, but now it is pretty unclear to me, how to use the M index to get upstream and downstream sequences

## helper functions to transform list returned by ACS to data frame
## as.df(...) transforms the list-of-lists to a dataframe
sec.split <- function(x){
  y <- ldply(.data = x, .fun = function(x){data.frame(prot = x$Text, offset=x$Offset)})
  return(y)
}
as.df <- function(x){
  y <- ldply(.data = x, .fun = sec.split, .id = "s.pep")
  return(y)
}
df <- as.df(l)

anno <- mutate(.data = df,
             prot.seq = BiocGenerics::sapply(ref.proteome[prot], toString),
             first.aa = substr(prot.seq, start = offset, stop = offset),
             last.aa = substr(prot.seq, start = offset+nchar(as.character(s.pep))-1, stop = offset+nchar(as.character(s.pep))-1),
             up.aa = substr(prot.seq, start = offset-1, stop = offset-1),
             down.aa = substr(prot.seq, start = offset+nchar(as.character(s.pep)), stop = offset+nchar(as.character(s.pep))),
             nterm.spec = up.aa %in% c("R", ""),
             cterm.spec = last.aa == "R")

# patterns ####
label <- "\\[\\+28\\.?0*\\]" # dimethly
nterm.label <- paste("^.", label, sep = "") # nterm dimethyl
nterm.acetylated <- "^.\\[\\+42\\.?0*\\]" #nterm acetylation
cterm.R <- "R$"
nterm.R <- "^R"

results <- mutate(.data = s,
                  nterm.label = case_when(
                    grepl(pattern = nterm.label, x = "m.pep") == TRUE ~ "label",
                    grepl(pattern = nterm.acetylated, x = "m.pep") == TRUE ~ "acetylated",
                    TRUE ~ "free"
                  )
            )

## join annotations and comparision results
annotated.results <- dplyr::left_join(x = results, y = anno, by = "s.pep")

## need one function that takes MSstats comparision results and transforms into annotated comparision results

annotate.groupComparisionResults <- function(fasta.file, MSstats.GC){
  
  ## fasta.file : name of the fasta file, expected in working dir
  ## MSstats.GC : R object returned by grouopComparision() function
  ## dependencies : Biostrings, BiocGenerics, AhoCorasickTrie, dplyr
  
  # read fasta file from disc ------------------------------------------------------------------
  
  if (fasta.file %in% list.files()) {
    message(paste("Reading protein sequences from:", fasta.file))
    proteome <- readAAStringSet(filepath = fasta.file)
    proteins <- BiocGenerics::sapply(proteome, toString)
    message(paste("Found", length(proteins), "protein entries in fasta file."))
  } else {
    stop("Fasta file does not exist in current working directory.")
  }
  
  # integrity check for comparision result -----------------------------------------------------
  
  if (class(MSstats.GC) != "list") {
    stop("Object is not of class list.")
  } else { 
    if ("ComparisonResult" %in% names(MSstats.GC)) {
      message("Integrity check: ok")
    } else {
      stop("Slot ComparisionResult was not found.")
      }
    }
  
  # helper functions ---------------------------------------------------------------------------
  
  # slot extraction of group comp. result
  get.CR <- function(x){
    x$ComparisonResult
  }
  
  # better
  get.gcr <- function(x, label){
    stopifnot(is.MSstats_gcr(x), class(label) == "character")
    x$ComparisonResult[x$ComparisonResult$Label == label,]
  }
  
  
  # strip modified aa sequences
  simple.mods <- function(x){
    y <- gsub(pattern = "\\.0+", x = x, replacement = "")
    return(y)
  }
  strip.ends <- function(x){
    x <- as.character(x)
    y <- gsub(pattern = "^_|_$", x = x, replacement = "")
    return(y)
  }
  strip.mods <- function(x){
    #any.mod <- "\\[(\\+|\\-)\\d+\\.\\d+]"
    #any.mod <- "\\[(\\+|-)\\d+\\]"
    any.mod <- "\\[(\\+|-)\\d+\\.?\\d+\\]"
    y <- gsub(pattern = any.mod, x = x, replacement = "")
    return(y)
  }
  
  # list-of-list conversion to df
  sec.split <- function(x){
    y <- ldply(.data = x, .fun = function(x){data.frame(prot = x$Text, offset=x$Offset)})
    return(y)
  }
  as.df <- function(x){
    y <- ldply(.data = x, .fun = sec.split, .id = "s.pep")
    return(y)
  }
  
  # modification patterns ----------------------------------------------------------------------
  
  label <- "\\[\\+28\\.?0*\\]" # dimethly
  nterm.label <- paste("^.", label, sep = "") # nterm dimethyl
  nterm.acetylated <- "^.\\[\\+42\\.?0*\\]" #nterm acetylation
  
  # function body ------------------------------------------------------------------------------
  
  gCR <- get.CR(MSstats.GC)
  message(paste("Annotating", nrow(gCR), "features found in MSstats group comparision:", levels(gCR$Label)))
  
  message("step 1 - Stripping peptide sequences.")
  gCR <- mutate(gCR, s.pep = strip.mods(strip.ends(Protein)), m.pep = strip.ends(Protein))
  
  message("step 2 - Locating stripped peptides using AhoCorasickSearch.")
  ll <- AhoCorasickSearch(keywords = gCR$s.pep, text = proteins, alphabet = "aminoacid", groupByKeyword = TRUE)
  
  message("step 3 - Converting search results to data frame.")
  pos.anno <- as.df(ll)
  
  message("step 4 - Adding variables to features.")
  pos.anno <- mutate(pos.anno,
              prot.seq = BiocGenerics::sapply(proteome[prot], toString),
              first.aa = substr(prot.seq, start = offset, stop = offset),
              last.aa = substr(prot.seq, start = offset+nchar(as.character(s.pep))-1, stop = offset+nchar(as.character(s.pep))-1),
              up.aa = substr(prot.seq, start = offset-1, stop = offset-1),
              down.aa = substr(prot.seq, start = offset+nchar(as.character(s.pep)), stop = offset+nchar(as.character(s.pep))))
  
  gCR <- mutate(gCR,
         nterm.label = case_when(
                       grepl(pattern = nterm.label, x = "m.pep") == TRUE ~ "label",
                       grepl(pattern = nterm.acetylated, x = "m.pep") == TRUE ~ "acetylated",
                       TRUE ~ "free"))
  
  message("step 5 - Joining positional information with comparision result.")
  pos.anno$s.pep <- as.character(pos.anno$s.pep)
  annotated.results <- dplyr::left_join(x = gCR, y = pos.anno, by = "s.pep")
  
  message("Done!")
  message(paste("Returning", nrow(annotated.results), "positionally annotated peptides."))
  return(annotated.results)
}

## test
annotated.results <- annotate.groupComparisionResults("uniprot-proteome_UP000000589.fasta", GluC.endpoint.peptide)


###################
# vulcano plot ####
###################
# take care of Inf cases
i.neg.inf <- x$log2FC == -Inf
i.neg.inf[is.na(i.neg.inf)] <- FALSE
i.pos.inf <- x$log2FC == Inf
i.pos.inf[is.na(i.pos.inf)] <- FALSE
x[i.pos.inf, c("log2FC")] <- 5
x[i.neg.inf, c("log2FC")] <- -5
x[i.pos.inf|i.neg.inf, c("adj.pvalue")] <- 0.00001

# add grouping variable for vplot or grouped summary
annotated.results <- mutate(annotated.results, pep.group = as.factor(up.aa %in% c("E","D")):as.factor(nterm.label):as.factor(cterm.spec))
annotated.results <- mutate(annotated.results, pep.group = as.factor(up.aa %in% c("E","D")):as.factor(nterm.label))
# add numeric ID column
annotated.results <- tibble::rownames_to_column(annotated.results, var="ID")

vplot <- function(x, alpha = 0.01, fc.cutoff = 2, version = ""){
  # filter for def. cases
  x <- filter(x, !is.na(adj.pvalue), !is.na(log2FC))
  # replace infinite cases
  x[is.infinite(x$log2FC), "adj.pvalue"] <- 0.00001
  x[x$log2FC == -Inf, "log2FC"] <- -5
  x[x$log2FC == Inf, "log2FC"] <- 5
  
  message(paste("n =", dim(x)[1]))
  # plot
  
  if(version == ""){
    message("Drawing basic vulcano plot")
    xyplot(-log10(adj.pvalue) ~ log2FC, data = x, main="Basic Vulcano Plot")
  }
  
  if(version == "c"){
    message("Drawing conditioned vulcano plot.")
    xyplot(-log10(adj.pvalue) ~ log2FC,
         groups = cut(log2FC, breaks= c(-Inf,-fc.cutoff, fc.cutoff, Inf)),
         panel = function(...){
           panel.xyplot(..., cex = 1, col.symbol = c("red", "grey", "green"))
           panel.abline(h=-log10(alpha), col="grey", lty = "dotted")
           panel.abline(v=c(-fc.cutoff, fc.cutoff), col="grey", lty = "dotted")
         },
         main="Conditioned Vulcano plot", data = x
    )
  }
  
  if(version == "g"){
    message("Drawing grouped vulcano plot.")
    xyplot(-log10(adj.pvalue) ~ log2FC,
           groups = pep.group,
           panel = function(...){
             panel.xyplot(..., cex = 1, jitter.y = FALSE, jitter.x = FALSE, amount = 0.5)
             panel.abline(h=-log10(alpha), col="grey", lty = "dotted")
             panel.abline(v=c(-fc.cutoff, fc.cutoff), col="grey", lty = "dotted")
           },
           main="Grouped Vulcano plot", auto.key = TRUE, data = x
    )
  }
}
vplot(annotated.results, version = "g")


#######################
# grouped analysis ####
#######################

make.pep.groups <- function(x){
  # creates peptides groups
  # expects data frame
  message(paste(c(" ### Grouping ", nrow(x), " features. ### ")))
  #y <- dplyr::summarise(group_by(x, Protein), mean.log2FC = mean(log2FC), mean.adj.pvalue = mean(adj.pvalue), up.aa = toString(levels(as.factor(up.aa))), down.aa = toString(levels(as.factor(down.aa))), nterm.label = toString(levels(as.factor(nterm.label))))
  y <- dplyr::summarise(group_by(x, Protein),
                        mean.log2FC = mean(log2FC, na.rm = TRUE),
                        mean.adj.pvalue = mean(adj.pvalue, na.rm = TRUE),
                        up.aa = paste(unique(up.aa), collapse = ","),
                        down.aa = paste(unique(down.aa), collapse = ","),
                        nterm.label = paste(unique(nterm.label), collapse = ","),
                        IDs = paste(unique(ID), collapse = ","),
                        prots = paste(unique(sub(pattern = "\\s.*", replacement = "", prot)), collapse = ","))
  return(y)
}
pep.g.annotated.results <- make.pep.groups(annotated.results)
pep.g.annotated.results <- mutate(pep.g.annotated.results, s.pep = strip.mods(strip.ends(Protein)))

grouped.vplot <- function(x){
  # !!! expects a grouped input !!!
  c1 <- 0.001
  c2 <- 5
  # filter for def. cases
  message(paste(c("### Input contains ", nrow(x), " cases. ###")))
  y <- filter(x, !is.na(mean.adj.pvalue), !is.na(mean.log2FC))
  message(paste(c("### Filtering for ", nrow(y), " complete cases. ###")))
  # replace infinite cases
  y[is.infinite(y$mean.log2FC), "mean.adj.pvalue"] <- c1
  y[y$mean.log2FC == -Inf, "mean.log2FC"] <- -c2
  y[y$mean.log2FC == Inf, "mean.log2FC"] <- c2
  message(paste("Inf anf -Inf cases are plotted at: ", c2,",", c1))
  # vulcano plot
  xyplot(-log10(mean.adj.pvalue)~mean.log2FC, groups = nterm.label, data = y, auto.key = TRUE, main="Vulcano plot for grouped peptides")
  xyplot(-log10(mean.adj.pvalue)~mean.log2FC | nterm.label, data = y, auto.key = TRUE, main="Vulcano plot for prouped peptides")
}

# neoN ####
#neoN <- filter(pep.g.annotated.results, mean.log2FC > 2, mean.adj.pvalue < 0.01, nterm.label == "label", substr(as.character(Protein), nchar(as.character(Protein)), nchar(as.character(Protein))) == "R")
filter.neoN <- function(x){
  y <- filter(x,
          mean.log2FC > 2,
          mean.adj.pvalue < 0.01,
          nterm.label == "label",
          substr(s.pep, nchar(s.pep), nchar(s.pep)) == "R"
        )
  return(y)
}
neoN <- filter.neoN(pep.g.annotated.results)
histogram(~as.factor(up.aa), data = neoN, subset = !grepl(",", up.aa), main="upstream aa of neo-N peptide")
xtabs(~up.aa, data = neoN)

# neoC ####
neoC <- filter(pep.g.annotated.results, mean.log2FC > 2, mean.adj.pvalue < 0.01, up.aa == "R")
filter.neoC <- function(x){
  y <- filter(x,
          mean.log2FC > 2,
          mean.adj.pvalue < 0.01,
          up.aa == "R"
        )
  return(y)
}
neoC <- filter.neoC(pep.g.annotated.results)
histogram(~as.factor(substr(s.pep, nchar(s.pep), nchar(s.pep))), data = neoC, main="last aa of neo-C peptide")

# spanning ####
filter.spanning <- function(x){
  y <- filter(x,
              up.aa == "R",
              substr(s.pep, nchar(s.pep), nchar(s.pep)) == "R",
              nterm.label == "free")
  return(y)
}
spanning <- filter.spanning(pep.g.annotated.results)

# neoN as substring in spanning? ####
sp <- spanning$s.pep
names(sp) <- spanning$s.pep
neoN.l <- AhoCorasickSearch(keywords = neoN$s.pep, text = sp, alphabet = "aminoacid", groupByKeyword = TRUE)

# neoC as substring in spanning? ####
neoC.l <- AhoCorasickSearch(keywords = neoC$s.pep, text = sp, alphabet = "aminoacid", groupByKeyword = TRUE)

#################
# use of tbl ####
#################

# dplyr grouping of tbl
prot.g.annotated.results <- group_by(annotated.results, prot)
hist(tally(prot.g.annotated.results)$n)
# plyr split+apply+combine on df to list
test <- dlply(.data = annotated.results, .variables = "prot", .fun = nrow)

#############
# ranges ####
#############

## AIM: creates a ranged representation of peptides mapped to proteins using
## IRanges and Biostrings
## list structure, each entry corresponds to a protein, ranges correspond to peptides
## 

## initial idea was to use views, but that is not as versatile as ranges
.make.views <- function(x){
  aa <- AAString(x$prot.seq[1])
  v <- Views(subject = aa, start = x$offset, width = nchar(x$s.pep), names = x$ID)
  return(v)
}

## central function to confert tidy df into a list of ranges
.make.ranges <- function(x){
  ## helper function for split+apply+combine
  ## transforms df into ranges object with metadata
  r <- IRanges(start = x$offset, width = nchar(x$s.pep), names = x$ID)
  elementMetadata(r) <- x[c(4,9,13:15,19:22)]
  return(r)
}

## above version is hard coded...better select columns based on names

.make.ranges <- function(x){
  ## helper function for split+apply+combine
  ## transforms df into ranges object with metadata
  r <- IRanges(start = x$offset, width = nchar(x$s.pep), names = x$ID)
  elementMetadata(r) <- x[c("log2FC", "adj.pvalue", "s.pep", "m.pep", "nterm.label", "first.aa", "last.aa", "up.aa", "down.aa")]
  return(r)
}

## unit test

poi <- annotated.results$prot[1]
df <- annotated.results[!is.na(annotated.results$prot),] #remove cases with missing protein
test <- df[df$prot == poi,] #single protein case
rl <- dlply(.data = test, .variables = "prot", .fun = failwith(NA, make.ranges), .progress = "text")

## split + apply + combine solution
## rl is a primitive list of IRanges
rl <- dlply(.data = df, .variables = "prot", .fun = failwith(NA, make.ranges), .progress = "text")

## coerse primitive list to RangedList object
real.rl <- as(rl, "RangesList")
class(real.rl)
head(start(real.rl))
universe(real.rl) <- "UP000000589"

## should be transformed into a function - and here it is:

df2IRangesList <- function(x, u = "proteome"){
  x <- x[!is.na(x$prot),] #remove cases with missing protein annotation, can not split here!!!
  y <- dlply(.data = x, .variables = "prot", .fun = failwith(NA, .make.ranges), .progress = "text")
  y <- as(y, "RangesList")
  universe(y) <- u
  return(y)
}

rl <- df2IRangesList(annotated.results, "UP000000589")
universe(rl)


## test applying ranges to proteins to make Views
n <- 10
poi <- names(rl[n])
ref.proteome[[poi]]
Views(ref.proteome[[poi]], rl[[n]])
(rl[[n]])
# works!

#######################
# finding patterns ####
#######################

## now that we have a list of ranges, we use this structure to find patterns!

## case 1 - finding adjacent ranges
adjacent <- function(x){
  # helper function
  # x should be instance of IRanges
  # function returns a Hits object
  s <- match(end(x)+1, start(x))
  q <- 1:length(x)
  y <- Hits(from = q[!is.na(s)], to = s[!is.na(s)], nLnode = length(x), nRnode = length(x))
  return(y)
}

## case 2 - finding paired ranges
## extension of adjacent, taking metadata into account
paired.ranges <- function(x, A = "free", B = "label"){
  # x is expected to be an instance of IRanges with metadata column nterm.label
  # functions returns ranges of type B that are directly precedet by a range of type A
  # A and B must be factor levels of nterm.label
  # ----AAAABBBB-----
  i <- (start(x)-1) %in% end(x) & elementMetadata(x)$nterm.label == B
  j <- elementMetadata(x[follow(x[i], x)])$nterm.label == A
  return(x[i&j])
}

## extension of case 2
# we need a global vars to check

test.prot.spec <- c("E", "D")
work.prot.spec <- c("R")

paired.ranges <- function(x, A = "free", B = "label"){
  # x is expected to be an instance of IRanges with ext. metadata
  # functions returns ranges of type B that are directly precedet by a range of type A
  # looks at working protease specificity defined in global variable work.prot.spec -> []
  # A and B must be factor levels of nterm.label
  # ---[-]AAAABBB[B]-----
  i <- (start(x)-1) %in% end(x) & elementMetadata(x)$nterm.label == B & elementMetadata(x)$last.aa %in% work.prot.spec
  j <- elementMetadata(x[follow(x[i], x)])$nterm.label == A & elementMetadata(x[follow(x[i], x)])$up.aa %in% work.prot.spec
  return(x[i&j])
}

## 2nd extension of case 2
paired.ranges <- function(x, A = "free", B = "label"){
  # x is expected to be an instance of IRanges with ext. metadata 
  # functions returns ranges of type B that are directly precedet by a range of type A
  # looks at working protease specificity defined in global variable work.prot.spec -> []
  # looks at test protease specificity -> {}
  # A and B must be factor levels of nterm.label
  # ---[R]XXX{E|D}XXX[R]---
  i <- (start(x)-1) %in% end(x) & elementMetadata(x)$nterm.label == B & elementMetadata(x)$last.aa %in% work.prot.spec & elementMetadata(x)$up.aa %in% test.prot.spec
  j <- elementMetadata(x[follow(x[i], x)])$nterm.label == A & elementMetadata(x[follow(x[i], x)])$up.aa %in% work.prot.spec 
  return(x[i&j])
}

## unit tests
# ---[R]XXX{E|D}XXX[R]---
# ---[R]XXXXXXXXXXX[R]
test <- IRanges(start = c(1, 6, 1), end = c(5, 10, 10))
elementMetadata(test) <- data.frame(nterm.label = c("free", "label", "free"),
                                    first.aa = c("X"),
                                    last.aa = c("E", "R", "R"),
                                    up.aa = c("R","E","R"),
                                    down.aa = c("X")
                                    )
plotRanges(test)
expect_is(paired.ranges(test), "IRanges")
expect_equal(length(paired.ranges(test)), 1)

# split + apply + combine solution for finding paired ranges over all proteins
pl <- llply(.data = rl, .fun = failwith(default = NA, f = paired.ranges), .progress = "text")

## examine the results I
## counting events
n.pairs <- ldply(.data = pl, .fun = failwith(NA, length), .progress = "text")
head(n.pairs[n.pairs$V1 != 0,], n = 20)
dim(n.pairs[n.pairs$V1 != 0,]) #n of proteins
sum(n.pairs$V1) #number of ranges
histogram(~V1, data = n.pairs[n.pairs$V1 != 0,])

## examine results for a specific poi
n <- which(n.pairs$V1 != 0)
i <- 16
pl[[n[i]]]
sort(rl[[n[i]]])
plotRanges(rl[[n[i]]], main=names(rl[n[i]]))
abline(v=start(pl[[n[i]]]), col="blue") #draws marker for potential cleavage site

head(as.data.frame(as(pl, "RangesList")))

###################################################
# find neoN-neoC pairs                         ####
###################################################

.neoCN.pairs <- function(x){
   # functions returns pairs of ranges of type label that are directly precedet by a range of type free
   # looks at working protease specificity defined in global variable work.prot.spec -> []
   # looks at test protease specificity -> {}
   # ---[R]XXX{E|D}XXX[R]---
   i <- (start(x)-1) %in% end(x) & elementMetadata(x)$nterm.label == "label" & elementMetadata(x)$last.aa %in% work.prot.spec & elementMetadata(x)$up.aa %in% test.prot.spec
   j <- elementMetadata(x[follow(x[i], x)])$nterm.label == "free" & elementMetadata(x[follow(x[i], x)])$up.aa %in% work.prot.spec 
   #return(x[i&j])
   return(c(x[i&j], x[follow(x[i&j], x)]))
}

test_that(".neoCN.pair is found in positive example", {
  ## test instance
  ## ---[R]XXX{E|D}XXX[R]---
  ## ---[R]XXXXXXXXXXX[R]
  test <- IRanges(start = c(1, 6, 1), end = c(5, 10, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label", "free"),
                                    first.aa = c("X"),
                                    last.aa = c("E", "R", "R"),
                                    up.aa = c("R","E","R"),
                                    down.aa = c("X")
  )
  #plotRanges(test)
  expect_is(.neoCN.pairs(test), "IRanges")
  expect_equal(length(.neoCN.pairs(test)), 2)
  #expect_equal(.neoCN.pairs(test), test[1:2])
  #order matters for testthat
})

## split + apply + combine
## implements filtering over lists of IRanges
filter.neoCN.pairs <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .neoCN.pairs), .progress = "text")
  return(y)
}
expect_is(filter.neoCN.pairs(rl[1]),"list")


## works, but does not generate a representation that keeps pairs together
## need better version that returns pairs or hits object
## keyfunction : findOverlapPairs()

## we need global vars to check specificity
## helper function
## retrieves enzyme specificity from a table and returns aa that are cleaved

lookup.specificity <- function(protease = "Trypsin"){
  t <- data.frame(name = c("GluC", "Trypsin", "Trypsin*", "unknown", "ArgC"), aa = c("ED", "RK", "R", NA, "R"))
  message(paste("Retrieving "), protease, " specificity.")
  i <- match(protease, t$name)
  return(unlist(strsplit(x = as.character(t[i, "aa"]), split = "")))
}

## unit test
lookup.specificity()
lookup.specificity("GluC")

test.prot.spec <- lookup.specificity("GluC")
work.prot.spec <- lookup.specificity("Trypsin*")

test.prot.spec <- lookup.specificity("unknown")


## suggestion by M. Lawrence

.experimental <- function(x){
  m <- elementMetadata(x)
  a <- x[m$nterm.label == "label" & m$last.aa %in% work.prot.spec & m$up.aa %in% test.prot.spec]
  b <- x[m$nterm.label == "free" & m$last.aa %in% test.prot.spec & m$up.aa %in% work.prot.spec]
  p <- findOverlapPairs(a, b, maxgap=1L)
  #p <- punion(subset(p, start(first) == end(second) + 1L | end(first) == start(second) - 1L))
  p <- punion(subset(p, start(first) == end(second) + 1L))
  elementMetadata(p) <- data.frame(nterm.label = "neoNC") 
  return(p)
}

test_that(".neoCN.pair is found in positive example", {
  ## test instance
  ## ---[R]XXX{E|D}XXX[R]---
  ## ---[R]XXXXXXXXXXX[R]
  test <- IRanges(start = c(1, 6, 1), end = c(5, 10, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label", "free"),
                                      first.aa = c("X"),
                                      last.aa = c("E", "R", "R"),
                                      up.aa = c("R","E","R"),
                                      down.aa = c("X")
  )
  expect_is(.experimental(test), "IRanges")
  expect_equal(length(.experimental(test)), 1)
})


## split + apply + combine
## implements finding pairs over lists of IRanges

filter.neoCN.pairs <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .experimental), .progress = "text")
  return(y)
}

## unit test
test_that("3 neoCN pairs are found in list of IRanges", {
  test.list <- list(a = test, b = test, c = test)
  #test.list <- as(test.list, "RangesList")
  expect_is(filter.neoCN.pairs(test.list), "list")
  expect_length(filter.neoCN.pairs(test.list), 3)
})

## okidoki - now see what we get on the complete data

neoCN <- filter.neoCN.pairs(rl)
sum(count.ranges(neoCN))
plotRanges(rl[[145]])
neoCN[[145]][1]

poi <- names(rl[145])
ref.proteome[[poi]]
Views(ref.proteome[[poi]], neoCN[[145]])
(rl[[n]])

############################################
### extension to missing test prot spec  ###
############################################

.experimental <- function(x){
  m <- elementMetadata(x)
  message(paste("Evaluating using working protease specificity:", paste(work.prot.spec, collapse = "")))
  message(paste("Evaluating using test protease specificity:", paste(test.prot.spec, collapse = "")))
  if (!is.na(test.prot.spec)){
    a <- x[m$nterm.label == "label" & m$last.aa %in% work.prot.spec & m$up.aa %in% test.prot.spec]
    b <- x[m$nterm.label == "free" & m$last.aa %in% test.prot.spec & m$up.aa %in% work.prot.spec]
  } else {
    a <- x[m$nterm.label == "label" & m$last.aa %in% work.prot.spec]
    b <- x[m$nterm.label == "free" & m$up.aa %in% work.prot.spec]
  }
  p <- findOverlapPairs(a, b, maxgap=1L)
  p <- punion(subset(p, start(first) == end(second) + 1L))
  elementMetadata(p) <- data.frame(nterm.label = "neoNC") 
  return(p)
}

test_that(".neoCN.pair is found in positive example", {
  ## test instance
  ## ---[R]XXX{E|D}XXX[R]---
  ## ---[R]XXXXXXXXXXX[R]
  test <- IRanges(start = c(1, 6, 1), end = c(5, 10, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label", "free"),
                                      first.aa = c("X"),
                                      last.aa = c("E", "R", "R"),
                                      up.aa = c("R","E","R"),
                                      down.aa = c("X")
  )
  test.prot.spec <- lookup.specificity("ArgC")
  work.prot.spec <- lookup.specificity("Trypsin")
  #plotRanges(test)
  #View(test)
  expect_is(.experimental(test), "IRanges")
  expect_equal(length(.experimental(test)), 1)
})


###################################################
# find pairs of neoN-spanning or neoC-spanning ####
###################################################

## case 1 - homo pairs
## same start, could also be a homo B, C, D, ...-type pair
## ----AAAAAAAAA----
## ----AAAAA--------

# returns ranges forming homo pair
homo.pairs <- function(x){
  y <- findOverlaps(x, type = "start", drop.self = TRUE, drop.redundant = TRUE)
  same <- elementMetadata(x)$nterm.label[from(y)] == elementMetadata(x)$nterm.label[to(y)]
  #return(y[same]) #hits object
  i <- c(from(y[same]), to(y[same]))
  return(x[i]) #returns all ranges contributing to hits
}

## unit test
test1 <- IRanges(start = c(1,1,4), end = c(5,3,5))
elementMetadata(test1) <- data.frame(nterm.label = c("free", "free", "label"))
plotRanges(test1)
expect_is(homo.pairs(test1), "IRanges")
expect_equal(length(homo.pairs(test1)), 2)

## case 2 - hetero pairs
## same start, but different types, need to be def., can not be any
## ----AAAAAAAAA----
## ----BBBBB--------

# returns ranges forming hetero pairs
hetero.pairs <- function(x, a = "free", b = "label"){
  x.a <- x[elementMetadata(x)$nterm.label == a]
  x.b <- x[elementMetadata(x)$nterm.label == b]
  h <- findOverlaps(query = x.a, subject = x.b, type = "start")
  #return(h) #for diagnostics
  return(c(x.a[from(h)], x.b[to(h)])) #returns all ranges contributing to hits
}

## unit test
test2 <- IRanges(start = c(1,1,4), end = c(5,3,5))
elementMetadata(test2) <- data.frame(nterm.label = c("free", "label", "free"))
plotRanges(test2)
expect_is(hetero.pairs(test2), "IRanges")
expect_equal(length(hetero.pairs(test2)), 2)

## case 3 - neoC-S pairs
## same start, same N-term label (free)
## ---[R].XXXXXXX[R].---- free
## ---[R].XXX{E|D}.------ free

.neoCS.pairs <- function(x){
   s.subset <- x[elementMetadata(x)$nterm.label == "free" & elementMetadata(x)$up.aa %in% work.prot.spec & elementMetadata(x)$last %in% work.prot.spec]
   c.subset <- x[elementMetadata(x)$nterm.label == "free" & elementMetadata(x)$up.aa %in% work.prot.spec & elementMetadata(x)$last %in% test.prot.spec]
   #message(length(s.subset))
   #message(length(c.subset))
   h <- findOverlaps(query = s.subset, subject = c.subset, type = "start")
   #return(h) #for diagnostics
   return(c(s.subset[from(h)], c.subset[to(h)])) #returns all ranges contributing to hits
}

## unit test 
test3 <- IRanges(start = c(1,1), end = c(5,3))
elementMetadata(test3) <- data.frame(nterm.label = c("free", "free"),
                                     first = c("X", "X"),
                                     last = c("R", "E"),
                                     up.aa = c("R", "R"),
                                     down.aa = c("X", "X"))
plotRanges(test3)
.neoCS.pairs(test3)
expect_is(.neoCS.pairs(test3), "IRanges")
expect_equal(length(.neoCS.pairs(test3)), 2)
#okidoki!

## split + apply + combine
filter.neoCS.pairs <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .neoCS.pairs), .progress = "text")
  return(y)
}
neoCS <- filter.neoCS.pairs(rl)
neoCS <- as(neoCS, "RangesList")
length(neoCS)
## examine the results
sum(count.ranges(neoCS))


## case 4 - neoN overlaps spanning
## same end
## ---[R]XXXXXXXX[R]----
## --------{E}XXX[R]----

.neoNS.pairs <- function(x){
   s.subset <- x[elementMetadata(x)$nterm.label == "free" & elementMetadata(x)$up.aa %in% work.prot.spec & elementMetadata(x)$last %in% work.prot.spec]
   n.subset <- x[elementMetadata(x)$nterm.label == "label" & elementMetadata(x)$up.aa %in% test.prot.spec & elementMetadata(x)$last %in% work.prot.spec]
   #message(length(s.subset))
   #message(length(c.subset))
   h <- findOverlaps(query = s.subset, subject = n.subset, type = "end")
   #return(h) #for diagnostics
   return(c(s.subset[from(h)], n.subset[to(h)])) #returns all ranges contributing to hits
}

## unit test 
test4 <- IRanges(start = c(1,3), end = c(5,5))
elementMetadata(test4) <- data.frame(nterm.label = c("free", "label"),
                                     first = c("X", "X"),
                                     last = c("R", "R"),
                                     up.aa = c("R", "E"),
                                     down.aa = c("X", "X"))
plotRanges(test4)
expect_that(.neoNS.pairs(test4), is_a(class(test4)))
expect_equal(length(.neoNS.pairs(test4)), 2)

## split + apply + combine
filter.neoNS.pairs <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .neoNS.pairs), .progress = "text")
  return(y)
}
neoNS <- filter.neoNS.pairs(rl)
neoNS <- as(neoNS, "RangesList")
length(neoNS)
## examine the results
sum(count.ranges(neoNS))

###################################################
# find ragged                                  ####
###################################################

## case 5 - ragged neoN series
## same end, all labeled
## ----{E}XXXXX[R]----
## --------XXXX[R]----
## ---------XXX[R]----
## ----------XX[R]----

.raggedN <- function(x){
  n.subset <- x[elementMetadata(x)$nterm.label == "label" & elementMetadata(x)$last %in% work.prot.spec]
  h <- findOverlaps(query = n.subset, type = "end", drop.self = TRUE, drop.redundant = TRUE)
  #View(h) #for diagnostics
  i <- unique(c(from(h), to(h)))
  return(n.subset[i]) #returns all ranges contributing to hits
}

## unit testing 
test5 <- IRanges(start = c(1, 2, 3, 4), end = c(5, 5, 5, 5))
elementMetadata(test5) <- data.frame(nterm.label = c("label", "free", "label", "label"),
                                     first = c("X"),
                                     last = c("R"),
                                     up.aa = c("E", "X", "X", "X"),
                                     down.aa = c("X")
                                     )
plotRanges(test5)
expect_is(.raggedN(test5), "IRanges")
expect_equal(length(.raggedN(test5)), 3)
plotRanges(.raggedN(test5))

## split + apply + combine
filter.raggedN <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .raggedN), .progress = "text")
  return(y)
}
ragged <- filter.raggedN(rl)
ragged <- as(ragged, "RangesList")
length(ragged)
## examine the results
count.ranges <- function(x){
  d <- ldply(.data = x, .fun = failwith(NA, length), .progress = "text")
  return(d$V1)
}
count.ranges(ragged)
plotRanges(rl[[380]])
plotRanges(ragged[[380]])
rl[380]

# there is an actin cluster, showing a lot of ragged sequences.

# find overlapping pairs
findOverlaps(query = rl[[n]], maxgap = 1, drop.self=TRUE, drop.redundant=TRUE)
findOverlaps(query = rl[[n]], maxgap = 0, type = "start", drop.self=TRUE, drop.redundant=TRUE) #same start
findOverlaps(query = rl[[n]], maxgap = 0, type = "end", drop.self=TRUE, drop.redundant=TRUE) #same end

make.views(test)
ranges(make.views(test))
make.ranges(test)

###################################################
# find neoC-neoN + spanning triplet            ####
###################################################

## strategy:
## 1. find neoC-neoN pairs
## 2. find overlab with spanning peptide

## unit tests
# ---[R]XXX{E|D}XXX[R]---
# ---[R]XXXXXXXXXXX[R]
test6 <- IRanges(start = c(1, 6, 1), end = c(5, 10, 10))
elementMetadata(test6) <- data.frame(nterm.label = c("free", "label", "free"),
                                     first.aa = c("X"),
                                     last.aa = c("E", "R", "R"),
                                     up.aa = c("R","E","R"),
                                     down.aa = c("X"))
plotRanges(test6)

## helper function to select peptides generated by working protease
## ---[R].XXXXXXX[R].---
## takes IRanges and returns only subset

.spanning <- function(x){
  e.meta <- elementMetadata(x)
  return(x[e.meta$up.aa == work.prot.spec & e.meta$last.aa == work.prot.spec & e.meta$nterm.label == "free"])
}

# unit tests
expect_is(.spanning(test6), "IRanges")
expect_equal(length(.spanning(test6)), 1)
expect_equal(.spanning(test6), test[3])

findOverlaps(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal") #returns hits object
subsetByOverlaps(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal") #returns IRanges query
findOverlapPairs(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal") #returns pairs of ranges

pintersect(findOverlapPairs(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal"))
expect_is(findOverlapPairs(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal"), "Pairs")
p <- Pairs(first = IRanges(1,10),  second = IRanges(1,10))
expect_equal(findOverlapPairs(reduce(.neoCN.pairs(test6)), .spanning(test6), type = "equal"), p)

## split + apply + combine
filter.spanning <- function(x){
  y <- llply(.data = x, .fun = failwith(default = NA, f = .spanning, .progress = "text"))
  return(y)
}
## unit test


filtered.rl <- as(filter.neoCN.pairs(rl), "RangesList")
count.ranges(filtered.rl)

#####################################
# sand box for playing with code ####
#####################################

# old, old, old... first ideas
# proteins as XStrings instance, AAString
poi <- ref.proteome[[1]]
alphabetFrequency(poi)
# peptides as Views, Ranges (like a clipping mask for XStrings)
v1 <- Views(poi, start = c(1,11,25), width = c(10,8,5), names = c("pepA", "pepB", "pepC"))
# subsetting
v1[1] #still a views instance
v1[[1]] #returns XString
v1[-3] #subsetting using neg i
# get slots
width(v1)
start(v1)
end(v1)
subject(v1)
ranges(v1)
#coercion
as(v1, "IRanges")

v1[1] == v1[2] # same AA in both ranges?

# compare views on AAString
hits <- findOverlaps(query = v1, maxgap = 1, drop.self=TRUE, drop.redundant=TRUE) #returns a Hits object
c <- countOverlaps(query = v1, maxgap = 1, drop.self=TRUE, drop.redundant=TRUE) #returns # hits per query

#compare ranges
pcompare(r[1], r[-1]) #returns code
subsetByOverlaps(query = r[1], subject = r[-1],  maxgap = 1)
poverlaps(query = r[1], subject = r[-1], maxgap = 1)
findOverlapPairs(query = r[1], r[-1], maxgap=1)
findOverlaps(query=ranges(v1), maxgap=1, drop.self=TRUE, drop.redundant=TRUE)


#####################
# filter results ####
#####################
# !!! these are the none grouped results !!!
# filter for n-term labeled peptides after D or E ending on R
test <- filter(annotated.results, up.aa %in% c("D","E"), nterm.label == "label", last.aa == "R")
# filter for sig, upregulated ones with label
test <- filter(annotated.results, log2FC > 2, adj.pvalue < 0.01, nterm.label == "label")
histogram(~as.factor(up.aa), data=filter(annotated.results, log2FC > 2, adj.pvalue < 0.01, nterm.label == "label"))
