###########
# info ####
###########
# this script uses dplyr to access the SQLite based spectral libraries build by BiblioSpec
# so far, it simply checks for iRT-related spectra
# https://skyline.ms/project/home/software/BiblioSpec/begin.view


######################
# install package ####
######################

# dplyr
source("https://bioconductor.org/biocLite.R")
biocLite("dplyr")
library("dplyr")

# create SQLlite connection
#my_db <- src_sqlite("F245177.blib", create=TRUE) #only one .dat file
my_db <- src_sqlite("joined_Gluc_control.blib", create = TRUE) # 2 .dat per .raw file
my_red.db <- src_sqlite("joined_Gluc_control.redundant.blib", create =TRUE) #2. dat per .raw file

my_db <- src_sqlite("dimethyl_GluC.blib", create = TRUE)
my_red.db <- src_sqlite("dimethyl_GluC.redundant.blib", create =TRUE)


# select specific tables
my_RefSpectra <- tbl(my_db, "RefSpectra")
my_Modifications <- tbl(my_db, "Modifications")
my_red.RefSpectra <- tbl(my_red.db, "RefSpectra")
my_red.SpectrumSourceFiles <- tbl(my_red.db, "SpectrumSourceFiles")
head(my_RefSpectra)
head(my_Modifications)

head(my_red.RefSpectra)
head(my_red.SpectrumSourceFiles)

# look for peptides from iRT set
# Generate iRT-C18 scale
Sequence <- c("LGGNEQVTR", "GAGSSEPVTGLDAK", "VEATFGVDESNAK", "YILAGVENSK", "TPVISGGPYEYR", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "LFLQFGAQGSPFLK")
Name <- c("iRT-pep a", "iRT-pep b", "iRT-pep c", "iRT-pep d", "iRT-pep e", "iRT-pep f", "iRT-pep g", "iRT-pep h", "iRT-pep i", "iRT-pep k", "iRT-pep l")
iRT <- c(-24.92, 0.00, 12.39, 19.79, 28.71, 33.38, 42.26, 54.62, 70.52, 87.23, 100.00)
iRT.C18 <- data.frame(Sequence, iRT, row.names = Name)

# retrieve table and plot RT vs iRT and regression line
test <- collect(filter(my_RefSpectra, peptideSeq %in% iRT.C18$Sequence))
index <- match(test$peptideSeq, iRT.C18$Sequence)
plot(test$retentionTime, iRT.C18[index, "iRT"])
abline(lm(iRT.C18[index, "iRT"]~test$retentionTime))

# retrieve all PSMs for iRT peptides and nth run
n <- 6
test <- collect(filter(my_red.RefSpectra, peptideSeq %in% iRT.C18$Sequence & fileID == n))
#summary(test)
#head(test)
#plot(test$retentionTime)
levels(as.factor(test$peptideSeq))
dotchart(test$retentionTime, test$peptideSeq)


# how are mod sequences coded
pep.mod.seq <- as.character(collect(select(my_table, peptideModSeq)))
