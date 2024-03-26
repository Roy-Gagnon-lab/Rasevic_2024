library(ggplot2)
library(ggsci)

my_theme <- function()
{
  res <- theme(
    plot.title = element_text(hjust=0),
    plot.title.position = "plot",
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey"),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom")
  res
}

load("chr14_regdata.Rdata")

p1 <- ggplot(chr14_regdata, aes(pos, -log10(POO_folate))) +
    geom_point(color="grey") +
    geom_point(data=chr14_regdata[!is.na(chr14_regdata$regulomeDB),], aes(pos, -log10(POO_folate), color=regulomeDB)) +
    scale_colour_viridis_d() + ylim(c(0,5)) + xlim(76100000, 76140000)+ labs(y="-log10(pval)") +my_theme()

p2 <- ggplot(chr14_regdata, aes(pos, -log10(POO_folate))) +
  geom_point(color="grey") +
  geom_point(data=chr14_regdata[!is.na(chr14_regdata$forgeDB),], aes(pos, -log10(POO_folate), color=forgeDB)) +
  scale_colour_viridis_d() + ylim(c(0,5)) + xlim(76100000, 76140000)  + labs(y="-log10(pval)") +my_theme()

cowplot::plot_grid(p1,p2, ncol=1)
cowplot::ggsave2("chr14_reg_score.pdf", width=6, height=6)

#--------------------------
LDX <- read.table("./LDexpress_08R2_500kb.txt", header=TRUE, sep='\t', quote="")
LDX$Position <- sapply(strsplit(LDX$Position, ":"), "[[", 2)

LDX$Gene.Symbol <- as.factor(LDX$Gene.Symbol)
LDX$Tissue <- as.factor(LDX$Tissue)
#turn rsid into factor but keep current ordernig which is postion based
LDX$RS.ID <- factor(LDX$RS.ID, levels = unique(LDX$RS.ID)) 

p1 <- ggplot(LDX, aes(x=RS.ID, y=P.value, fill=Tissue)) + geom_col() + facet_wrap(vars(Gene.Symbol), ncol =1) +
    theme(axis.text.x =element_text(angle=90)) 
# cowplot::ggsave2("chr14_expression_13SNPs.pdf", width=8, height=5)

#------
chr14_regdata$rsid <- NA
chr14_regdata$rsid[chr14_regdata$pos %in% LDX$Position] <- unique(LDX$RS.ID)
p2 <- ggplot(chr14_regdata, aes(pos, -log10(POO_folate), label=rsid)) +
  geom_point(color="grey") +
  geom_point(data=chr14_regdata[!is.na(chr14_regdata$rsid),], aes(pos, -log10(POO_folate), color="red")) +
  ylim(c(0,5)) + xlim(76100000, 76140000)+ labs(y="-log10(pval)") +my_theme()

cowplot::plot_grid(p2, p1, ncol=1, rel_heights = c(1,4))
cowplot::ggsave2("chr14_13SNPs_LDX.pdf", width=8, height=8)

# load("../p_LD1.RData")
# p_LD1 <- p_LD1[grep("^chr14:",colnames(p_LD1)), grep("^chr14:",colnames(p_LD1))]
# colnames(p_LD1) <- sapply(strsplit(colnames(p_LD1), ":"), "[[", 2) 
# rownames(p_LD1) <- colnames(p_LD1)







