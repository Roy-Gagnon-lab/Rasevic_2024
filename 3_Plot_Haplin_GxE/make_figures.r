
library(readr)
library(dplyr)
library(ggplot2)
library(poolr)
# library(matrixcalc)
library(STAAR)
library(ggnewscale)
library(ggsci)

#-----------GGPLOT THEME------------------
my_theme <- function()
{
  res <- theme(
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey"),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom")
  res
}


#--------------- plotting function -------------------
LD_lead_SNP <- function(df, LD_matrix, prune_list, title){
  df <- na.omit(df)
  
  chipped_regions <- df %>% group_by(region) %>% summarise(any(Genotyped == "Genotyped"))
  chipped_regions <- as.data.frame(chipped_regions)
  chipped_regions <- as.vector(chipped_regions[,1][chipped_regions[,2]])
  
  df <- df[df$region %in% chipped_regions, ]
  df <- df[order(df$chr, df$position),]
  
  df$region <- df$region2
  df$region <- factor(df$region, levels=unique(df$region))
  
  pruned_results <- df[df$snp %in% prune_list,]
  
  #add a distance column to help plot the x-axis
  df <- df %>% group_by(region) %>% summarise(start_BP = min(position)) %>% left_join(df , . , by=c("region"="region"))
  df$distance <- df$position - df$start_BP
  increment <- (df %>% group_by(region) %>% summarise(increment = max(distance)))$increment
  increment <- c(0, head(increment, -1) + 5000)
  increment <- cumsum(increment)
  
  df$distance <- df$distance + rep(increment, as.vector(table(df$region)))
  
  #put the x axis labels in the right spot
  xaxis <- df %>% group_by(region) %>% 
    summarise(center=(max(distance) + min(distance))/2)
  
  #put vertical line between regions
  v_line <- df %>% group_by(region) %>% 
    summarise(x_int = max(distance) + 2500)
  
  # get pooled values
  extra_layer <- df %>% group_by(region) %>%
    summarise (Cauchy = CCT(pval),
               Empirical = fisher(pval, adjust="empirical", R = LD_matrix[unname(snp), unname(snp)])$p)
  
  # extra_layer <- pruned_results %>% group_by(region) %>%
  #   summarise(Fisher = fisher(pval)$p) %>%
  #   left_join(extra_layer, . , by="region")
  
  #turn dataframe long
  extra_layer <- extra_layer %>% tidyr::pivot_longer(c(!region), names_to="pooled", values_to="pval")
  extra_layer$pooled[extra_layer$pooled == "Empirical"] <- "Empirical Fisher"
  extra_layer <- xaxis %>% left_join(extra_layer, . , by="region")
  
  LDvec <- c()
  for (region in levels(df$region)){
    subdf   <- df[df$region == region,]
    snps    <- subdf$snp
    leadSNP <- subdf[which(subdf$pval == min(subdf$pval)), "snp"][1]
    LDvec   <- c(LDvec, unname(LD_matrix[leadSNP, snps]))
  }
  
  df$LD <- LDvec^2
  p <-ggplot(df, aes(distance, -log10(pval), color=LD)) +
    geom_point(shape=18) +
    # geom_point(data=df[df$region %in% unique(df$region)[c(T,F)],], mapping = aes(x=distance, y=-log10(pval)), shape=18) +
    # geom_point(data=df[df$region %in% unique(df$region)[c(F,T)],], mapping = aes(x=distance, y=-log10(pval)), shape=20) +
    geom_point(data=df[df$Genotyped == "Genotyped", ], mapping = aes(x=distance, y = -log10(pval), color= LD), shape=23, fill="white") +
    geom_point(data=extra_layer, aes(x=center, y=-log10(pval), shape=pooled), color="black") +
    scale_shape_manual(values=c(1,2)) + 
    scale_colour_stepsn(breaks = c(0.2,0.4,0.6,0.8), colours = pal_locuszoom("default")(7)[c(6,4,3,2,1)]) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#071E22") +
    geom_hline(yintercept=-log10(0.05/length(unique(df$region))), linetype="dashed", color = "#679289") +
    geom_hline(yintercept=-log10(0.05/nrow(df)), linetype="dashed", color = "#F4C095") +
    ylim(c(0,5)) + my_theme() +  labs( x="region", shape="") +
    scale_x_continuous(label=xaxis$region, breaks=xaxis$center)
  return(p)
}

#--------------- load pre-processed data ------------------------------
m_keep <- read_table("maternal.in", col_names = FALSE)[[1]]
p_keep <- read_table("p_origin.in", col_names=FALSE)[[1]]

load("mat.Rdata")
load("p_origin.Rdata")
#LD matrix
load("m_LD1.Rdata")
load("p_LD1.Rdata")

mat <- mat[mat$impute > 0.8, ]
p_origin <- p_origin[p_origin$impute > 0.8 , ]
p_origin <- p_origin[p_origin$region2 != "chr14\n88.3Mb", ]


#--------------- make figures. fig1 = maternal effects. fig2= parent of origin.
sub_pheno <- c("OFC", "CPO", "CLPP")
for (i in sub_pheno){
  in_file        <- paste0("../../NewAnalyses/New_Cleft_Haplin_Analysis/Input_and_Output_Files/maternal_smoke_", i, "_haplin.txt")
  haplin_results <- read_table(in_file)
  haplin_results <- merge(mat, haplin_results[,c("snp", "maternal_effectRR", "gxe_smoke_mat", "maf", "no_env_mat_RR", "env_mat_RR")], by="snp")
  names(haplin_results)[names(haplin_results) == 'gxe_smoke_mat'] <- 'pval'
  
  title  <- paste0(i, ": maternal effect smoking GxE")
  p_A    <- LD_lead_SNP(haplin_results, m_LD1, m_keep, "CPO: maternal effect smoking GxE")
  legend <- cowplot::get_legend(p_A + theme(legend.box.margin = margin(0, 0, 0, 0)))
  
  in_file        <- paste0("../../NewAnalyses/New_Cleft_Haplin_Analysis/Input_and_Output_Files/maternal_folate_", i , "_haplin.txt")
  haplin_results <- read_table(in_file)
  haplin_results <- merge(mat, haplin_results[,c("snp", "maternal_effectRR", "gxe_folate_mat", "maf", "no_env_mat_RR", "env_mat_RR")], by="snp")
  names(haplin_results)[names(haplin_results) == 'gxe_folate_mat'] <- 'pval'
  
  title  <- paste0(i, ": maternal effect folate GxE")
  p_B <- LD_lead_SNP(haplin_results, m_LD1, m_keep, title)
  
  plot <- cowplot::plot_grid(p_A + theme(legend.position="none"), p_B + theme(legend.position="none"), ncol=1, labels=c("A","B"))
  plot <- cowplot::plot_grid(plot, legend, ncol=1, rel_heights = c(1,0.1))
  out  <- paste0("./plots_and_tables/figure1_", i, ".pdf")
  cowplot::ggsave2(out, width=7.3, height=5, plot)
}

sub_pheno <- c("OFC", "CPO", "CLPP")
for (i in sub_pheno){
  in_file        <- paste0("../../NewAnalyses/New_Cleft_Haplin_Analysis/Input_and_Output_Files/parental_smoke_", i, "_haplin.txt")
  haplin_results <- read_table(in_file)
  haplin_results <- merge(p_origin, haplin_results[,c("snp", "gxe_smoke_poo", "env_poo_RR", "no_env_poo_RR")], by="snp")
  haplin_results <- rename(haplin_results, 'pval'='gxe_smoke_poo', 'no_env_RR'='no_env_poo_RR', 'env_RR'='env_poo_RR')
  p_A <- LD_lead_SNP(haplin_results, p_LD1, p_keep, "")
  legend <- cowplot::get_legend(p_A + theme(legend.box.margin = margin(0, 0, 0, 0)))
  
  in_file        <- paste0("../../NewAnalyses/New_Cleft_Haplin_Analysis/Input_and_Output_Files/parental_folate_", i, "_haplin.txt")
  haplin_results <- read_table(in_file)
  haplin_results <- merge(p_origin, haplin_results[,c("snp", "gxe_folate_poo", "env_poo_RR", "no_env_poo_RR")], by="snp")
  haplin_results <- rename(haplin_results, 'pval'='gxe_folate_poo', 'no_env_RR'='no_env_poo_RR', 'env_RR'='env_poo_RR')
  p_B <- LD_lead_SNP(haplin_results, p_LD1, p_keep, "")
  
  plot <- cowplot::plot_grid(p_A + theme(legend.position="none"), p_B + theme(legend.position="none"), ncol=1, labels=c("A","B"))
  plot <- cowplot::plot_grid(plot, legend, ncol=1, rel_heights = c(1,0.1))
  out  <- paste0("./plots_and_tables/figure2_", i, ".pdf")
  cowplot::ggsave2(out, width = 7.3, height=5, plot)
}


