library(Rsubread)



# generate count matrix
files <- list.files(path = ".",
                  pattern = "marked.bam$",
                  recursive = TRUE,
                  include.dirs = TRUE)
SAF <- read.table(file = "gencode_v24_GQ994935.SAF", header = T)
feature_sum <- featureCounts(files = files,
                             isPairedEnd = TRUE,
                             #isGTFAnnotationFile = TRUE,
                             annot.ext = SAF,
                             #GTF.attrType = "gene_id",
                             useMetaFeatures = TRUE,
                             strandSpecific = 2,
                             requireBothEndsMapped = TRUE,
                             nthreads = 12,
                             autosort = TRUE,
                             countMultiMappingReads = FALSE,
                             countChimericFragments = FALSE
)

saveRDS(feature_sum, file = "feature_total_uniq_counts.rds")
