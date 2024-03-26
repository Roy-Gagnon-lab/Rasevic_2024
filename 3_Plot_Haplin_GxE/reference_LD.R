library(LDlinkR)
library(readr)
library(stringr)

load("mat.Rdata")
load("p_origin.Rdata")

region_LD <- function(snps){
  x <- LDmatrix(snps = snps, pop="EUR", r2d = "r2", token = "6fd23f7590fe", genome_build = "grch38")
  return(x)
}

load("m_LD1.Rdata") #plink R2 matrix
m_LD2 <- m_LD1^2 #will replace the entries with the LDlinkR values. some SNPs are missing from the results so this will keep those as the plink values

for (i in unique(mat$region)){
  index = which(mat$region == i)
  snps = mat[index,]$rsid
  
  x <- region_LD(snps)
  index = which(x[,1] %in% mat$rsid)
  
  x <- as.matrix(x[,-1])
  m_LD2[index, index] <- x
  Sys.sleep(0.5)
}

save(m_LD2, file="m_LD2.Rdata")

load("p_LD1.Rdata")
p_LD2 <- p_LD1^2 #will replace the entries with the LDlinkR values. some SNPs are missing from the results so this will keep those as the plink values
for (i in unique(p_origin$region)){
  index = which(p_origin$region == i)
  snps = p_origin[index,]$rsid
  
  x <- region_LD(snps)
  index = which(x[,1] %in% p_origin$rsid)
  
  x <- as.matrix(x[,-1])
  p_LD2[index, index] <- x
  Sys.sleep(0.5)
}
save(p_LD2, file="p_LD2.Rdata")

#----------get the genetics R matrix  
library(genetics)
library(readr)
mat_ped <- read_table("maternal.ped", col_names=FALSE)
load("mat.Rdata")

m_genotypes <- list()
for (i in 1:nrow(mat)){
  var = mat$snp[i]
  m_genotypes[[var]] <- paste(mat_ped[[(2*i)+5]], mat_ped[[(2*i)+6]], sep="/")
}
m_genotypes <- as.data.frame(m_genotypes)
m_genotypes <- makeGenotypes(m_genotypes)
m_LD3 <- LD(m_genotypes)$r

