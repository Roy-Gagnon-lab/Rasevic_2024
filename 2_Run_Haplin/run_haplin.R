library(Haplin)

for (i in c("maternal", "parental")){
    for (j in c("OFC", "CPO", "CLPP")){
        for(k in c("folate", "smoke")){
            pedfile <- paste0(i, "_", j, "_", k, ".ped")
            covfile <- paste0(k, "_", j, ".cov")

            gendata <- genDataRead(file.in = pedfile, format = "ped", cov.file.in = covfile, na.strings = "0",
                                    overwrite = TRUE)
            gendata <- genDataPreprocess(data.in = gendata, design = "triad",
                                         overwrite=TRUE)
            poo     <- ifelse(i=="parental", TRUE, FALSE)
            output  <- haplinSlide(gendata, maternal = TRUE, poo = poo, strata = k,
                                design = "triad", use.missing = TRUE, response = "mult", reference = "ref.cat", winlength = 1 )          
            
            save(output, paste0("haplin_",i, "_", j, "_", k,".Rdata")) 
            if (i=="maternal"){
                pval <-  sapply(output, function(x) `[`(x,"gxe_maternal_pval")[[1]][1])
            }else{
                pval <-  sapply(output, function(x) `[`(x,"gxe_poo_pval")[[1]][1])
            }
            # print(c(i,j,k))
            # print(pval)
            save(pval, paste0("pval_",i, "_", j, "_", k,".Rdata"))             
        }
    }
}
