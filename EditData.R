#Read information from files provided by Tong et al. (2020) and Arnold et al. (2014)
#And output files required for the experiments of Onogi (2023)
#The directory "/netGS/" includes all files provided at https://github.com/Hao-Tong/netGS/tree/master/netGS

#Metabolite network information
Xml <- readLines("netGS/model.xml")

#Biomass reaction for each variety
Biomass <- read.csv("netGS/biomreaction.csv", header = FALSE)

#Number of metabolites
##Metabolites are called as species
Metabolite <- matrix(unlist(strsplit(Xml[grep("species id", Xml)], "\"")), nrow = 13)[2, ]
K <- length(Metabolite)

#Number of reactions
Reaction <- matrix(unlist(strsplit(Xml[grep("reaction id", Xml)], "\"")), nrow = 7)[2, ]
J <- length(Reaction)

#Stoichiometry matrix
Smatrix <- matrix(0, K, J)
rownames(Smatrix) <- Metabolite
colnames(Smatrix) <- Reaction

#Lower and upper bounds of fluxes
Bound <- matrix(0, 2, J)
rownames(Bound) <- c("Lower", "Upper")
colnames(Bound) <- Reaction

#Assign metabolite names to biomass reaction
rownames(Biomass) <- Metabolite

#Fill matrices
Reaction.start <- grep("reaction id", Xml)
Reaction.start <- c(Reaction.start, length(Xml) + 1)
names(Reaction.start) <- c(Reaction, "End")
for(j in 1:J){
  
  reaction <- Reaction[j]
  v <- Xml[Reaction.start[j]:(Reaction.start[j + 1] - 1)]
  
  w <- grep("listOfReactants", v)
  if(length(w) > 0){
    Reactant <- v[(w[1] + 1):(w[2] - 1)]
    for(i in 1:length(Reactant)){
      r <- unlist(strsplit(Reactant[i], "\""))[c(2, 4)]
      Smatrix[r[1], reaction] <- as.numeric(r[2]) * (-1)
    }
  }
  
  w <- grep("listOfProducts", v)
  if(length(w) > 0){
    Product <- v[(w[1] + 1):(w[2] - 1)]
    for(i in 1:length(Product)){
      r <- unlist(strsplit(Product[i], "\""))[c(2, 4)]
      Smatrix[r[1], reaction] <- as.numeric(r[2])
    }
  }
  
  w <- grep("BOUND", v)
  Bound["Lower", reaction] <- as.numeric(unlist(strsplit(v[w[1]], "\""))[4])
  Bound["Upper", reaction] <- as.numeric(unlist(strsplit(v[w[2]], "\""))[4])
}
any(colSums(Smatrix != 0) == 0)
FALSE
any(rowSums(Smatrix != 0) == 0)
FALSE

#Order the metabolites according to Arnold et al. (2014) SupplementaryData3.csv
MetaboliteOrder <- unlist(read.table("MetaboliteOrder.txt"))
#=>This text file includes the order in Arnold et al. (2014) and created by Onogi
Smatrix <- Smatrix[MetaboliteOrder, ]

#Remove reactions with zero fluxes in the reference genotype
DataS2 <- read.csv("TongH2020SupplementaryData2.csv", header = TRUE)
#=>This csv file was created from supplementary data 2 of Tong et al. (2020) with modifications
#=>and includes the fluxes of each genotype for the reactions that show non-zero fluxes in the reference
#=>The abbreviations of the reaction names were added by Onogi
w <- which(!is.element(paste("R", DataS2$Abbreviation, sep = "_"), colnames(Smatrix)))
#Some discrepancies in abbreviations of reaction names between Smatrix and DataS2
w
#4  21  22 174 184
DataS2$Abbreviation[w]
#"Fd-NADPR_h"  "StS_H2"      "StS_H3"      "Asp-SeADH_h" "5M-THFOR_c" 
colnames(Smatrix)[grep("NADPR", colnames(Smatrix))]
#"R_Fd_DASH_NADPR_h"
colnames(Smatrix)[grep("StS_h2",colnames(Smatrix))]
#"R_StS_h2"
colnames(Smatrix)[grep("StS_h3",colnames(Smatrix))]
"R_StS_h3"
colnames(Smatrix)[grep("ADH",colnames(Smatrix))]
#[5] "R_Asp_DASH_SeADH_h"
colnames(Smatrix)[grep("THFOR",colnames(Smatrix))]
#"R_5M_DASH_THFOR_c"

#Rename abbreviations of DataS2
DataS2$Abbreviation[c(4, 21, 22, 174, 184)] <- 
  c("Fd_DASH_NADPR_h", "StS_h2", "StS_h3", "Asp_DASH_SeADH_h", "5M_DASH_THFOR_c")
any(!is.element(paste("R", DataS2$Abbreviation, sep = "_"), colnames(Smatrix)))
FALSE

Smatrix.nonzero <- Smatrix[, is.element(colnames(Smatrix), paste("R", DataS2$Abbreviation, sep = "_"))]
dim(Smatrix.nonzero)
#407 336
Smatrix.nonzero <- Smatrix.nonzero[rowSums(abs(Smatrix.nonzero)) > 0, ]
dim(Smatrix.nonzero)
#350 336
write.csv(Smatrix.nonzero, "Smatrix.nonzero.csv")

#Bounds of fluxes
Bound.nonzero <- Bound[, is.element(colnames(Bound), paste("R", DataS2$Abbreviation, sep = "_"))]
identical(colnames(Bound.nonzero), colnames(Smatrix.nonzero))
TRUE
write.csv(Bound.nonzero, "Bound.nonzero.csv") 

#Biomass reaction unique to varieties
Biomass.nonzero <- Biomass[rownames(Smatrix.nonzero), ]
dim(Biomass.nonzero)
#350  67
write.csv(Biomass.nonzero, "Biomass.nonzero.csv")

#Create the relationship matrix
Geno<-read.csv("netGS/snp.csv", header = FALSE)
dim(Geno)
#67 1824

MAF <- colSums(Geno + 1)/(nrow(Geno) * 2)
MAF[MAF > 0.5] <- 1 - MAF[MAF > 0.5]
range(MAF)
#0.01492537 0.49253731
library(rrBLUP)
Gmatrix <- A.mat(Geno, shrink = TRUE, min.MAF = 0)#Use all SNPs
hist(diag(Gmatrix))
hist(Gmatrix[upper.tri(Gmatrix, diag = FALSE)])
write.csv(Gmatrix, "Gmatrix.csv", row.names = FALSE)
