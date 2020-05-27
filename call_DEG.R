library(edgeR)
library(ggplot2)
library(statmod)
library(cowplot)
library(tidyverse)



feature_sum <- readRDS(file = "feature_total_uniq_counts.rds")
# read in metadata for samples
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



# volcano plot -------------------------------------------------------------------------------
library(ggplot2)
library(cowplot)

p <- lapply(qlf, function(x){
  p <- ggplot(data = x$table, aes(x = logFC, 
                                  y = -log10(padj),
                                  colour = threshold))+
    geom_point()+
    xlab("log2 fold change")+
    ylab("-log10 padj")+
    ggtitle(x$comparison)+
    scale_color_manual(values = c("blue", "black", "red"))+
    theme_cowplot()
  return(p)
})

if(TRUE){
  # volcano plot
plot_grid(plotlist = list(p[[1]], p[[2]], p[[3]], p[[4]]))
ggsave("volcanoPlot.pdf")
}

lapply(qlf, function(x){write.csv(x$table, file=paste0(x$comparison,".csv"))})

# export gene_type ------------------------------------------------------------------------------------------------
# only host genes
library(rtracklayer)
SAF <- import.gff2(con = "./gencode.v24.primary_assembly_GQ994935.annotation.gtf")
SAF <- data.frame(SAF)


# abundant gene_type
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

# export report
lapply(qlf_anno, function(x){
  #ty <- table(qlf_anno[[1]]$table$gene_type) %>% as_tibble()
  x$table %>% as_tibble() %>% write_csv(path = paste0(".",x$comparison,"up.csv"))
  return(NULL)
})
lapply(qlf_anno, function(x){
  #ty <- table(qlf_anno[[1]]$table$gene_type) %>% as_tibble()
  p <- table(x$table$gene_type) %>% as_tibble() %>%
    ggplot(aes(x="", y = n, fill = Var1)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
    ggtitle(x$comparison) + theme_minimal()
  ggsave(filename = paste0(x$comparison,"gene_type.pdf"),p)
  table(x$table$gene_type) %>% as_tibble() %>% write_csv(path = paste0(".",x$comparison,"gene_type.csv"))
  return(NULL)
})


