library(edgeR)
library(ggplot2)
library(statmod)
library(cowplot)
library(tidyverse)


feature_sum <- readRDS(file = "feature_total_uniq_counts.rds")
metaD <- read.csv(file = "./2629_ProjectSummary.csv")
metaD <- metaD[order(metaD$File.Name),]
metaD <- metaD[2:13,]
metaD$Treatment <- c("Latency_Input", "Lytic_Input", "Lytic_IgG", "Lytic_UPF1", "Latency_IgG","Latency_UPF1",
                     "Latency_Input", "Latency_IgG","Latency_UPF1",
                     "Lytic_Input", "Lytic_IgG", "Lytic_UPF1")


# differentially test------------------------------------------------------------------------------------------
y <- DGEList(counts = feature_sum$counts, 
             samples = metaD,
             group = metaD$Treatment)

# cpm = 1 is count of 36~37 counts in smallest sample
# keep <- rowSums(cpm(y)>1) >= 2 not run

#y <- y[keep, , keep.lib.sizes=FALSE] not run
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + Treatment, metaD)  
y <- estimateDisp(y, design, robust = TRUE)

fit <- glmQLFit(y, design = design, robust = TRUE)
#write.csv(fit$counts, file = "filtered_counts.csv")


con1 <- makeContrasts(TreatmentLatency_UPF1 - TreatmentLatency_Input, levels = design)
con2 <- makeContrasts(TreatmentLatency_UPF1 - TreatmentLatency_IgG, levels = design)
con3 <- makeContrasts(TreatmentLytic_UPF1 - TreatmentLytic_Input, levels = design)
con4 <- makeContrasts(TreatmentLytic_UPF1 - TreatmentLytic_IgG, levels = design)


con <- list(con1, con2, con3, con4)



qlf <- lapply(con, function(x, fit){glmQLFTest(fit, contrast = x)}, fit = fit)
qlf <- lapply(qlf, function(x){
  x$table$padj <- p.adjust(x$table$PValue, method = "BH")
  x$table$threshold <- as.factor(ifelse(x$table$padj < 0.05 & abs(x$table$logFC) > 1,
                                        ifelse(x$table$logFC > 1, 'Up', 'Down'), 'Not'))
  return(x)})



 export gene_type ------------------------------------------------------------------------------------------------
# only host genes
library(rtracklayer)
SAF <- import.gff2(con = "./gencode.v24.primary_assembly_GQ994935.annotation.gtf")
SAF <- data.frame(SAF)


qlf_anno <- lapply(qlf, function(x){
  #as_tibble(qlf[[1]]$table) %>% rownames_to_column()
  x_tb <- as_tibble(x$table) %>% rownames_to_column()
  colnames(x_tb) <- c("gene_id", "logFC","logCPM","F","PValue","padj","threshold")
  x_tb <- left_join(x_tb, SAF, by = "gene_id")
  x_tb <- x_tb[!duplicated(x_tb$gene_id),]
  write_csv(x_tb, path = paste0(x$comparison, "all.csv"))
  x_tb <- x_tb[x_tb$padj < 0.05 & x_tb$logFC > 1,]
  x$table <- x_tb
  return(x)
})


# NMD_target wilcox test -----------------------------------------------------
NMD_targets <- read.csv(file = "NMD_targets.csv", header = TRUE, skip = 1)

all_genes <- qlf_anno

NMD_tb <- lapply(all_genes, function(x){
  x$table <- x$table[x$table$gene_name %in% NMD_targets$Gene.ID,]
  return(x)
})

cum_plot <- function(qlf, NMD_tb){
  x <- qlf
  y <- NMD_tb
  pl <- data.frame(cs = c(cumsum(rep(100,17516))/17516, cumsum(rep(100,89))/89), 
                   fc = c(x$table$logFC, y$table$logFC),
                   set = c(rep("all", 17516), rep("NMD_targets", 89))
  )
  p <- ggplot(data = pl, aes(x=fc,y=cs, group=set))+
    geom_line(aes(color=set))+
    ylab("Cumulative fraction (%)")+
    xlab("log2 ratio of fold-change")+
    xlim(c(-5,5))
  ggsave(p, filename = paste0(x$comparison, "_cum_plot.pdf"))
  return(NULL)
}

cum_plot(qlf_anno[[1]], NMD_tb[[1]])
cum_plot(qlf_anno[[2]], NMD_tb[[2]])
cum_plot(qlf_anno[[3]], NMD_tb[[3]])
cum_plot(qlf_anno[[4]], NMD_tb[[4]])

qlf_anno[[1]]$comparison
wilcox.test(qlf_anno[[1]]$table$logFC, NMD_tb[[1]]$table$logFC, paired = F)
qlf_anno[[2]]$comparison
wilcox.test(qlf_anno[[2]]$table$logFC, NMD_tb[[2]]$table$logFC, paired = F)
qlf_anno[[3]]$comparison
wilcox.test(qlf_anno[[3]]$table$logFC, NMD_tb[[3]]$table$logFC, paired = F)
qlf_anno[[4]]$comparison
wilcox.test(qlf_anno[[4]]$table$logFC, NMD_tb[[4]]$table$logFC, paired = F)
