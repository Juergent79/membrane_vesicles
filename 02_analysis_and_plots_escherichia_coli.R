library(Biostrings) # DNA seq loading and manipulation
library(zoo) # sliding windows
library(tidyverse)
library(RColorBrewer)
library(edgeR)

#### Load Data ####
pileup_path <- paste(file.path("./omv_pileup/"), list.files("./omv_pileup/"), sep = "/")

pileup_df <- data.frame(Replicon = rep("NZ_CP009273.1", 4631469),
                        Position = seq(1:4631469),
                        stringsAsFactors = FALSE)

for (i in 1:length(pileup_path)) {
  pileup_df[,i+2] <- read.table(pileup_path[i], header = FALSE, stringsAsFactors = FALSE)[,4]
  
}

colnames(pileup_df)[-(1:2)] <- gsub("_all.pileup", "",list.files("./omv_pileup/"))
colnames(pileup_df) <- gsub("-", "_", colnames(pileup_df))

#### Calculate mean of sliding window ####
# size 1kb, step 0.5kb 

covmeans_rep <- as.data.frame(rollapply(pileup_df[,-1],
                                        FUN = mean, width = 1000,
                                        by = 500, align = "left", partial = 50))

#### Extract dif region ####
# single nucleotide resolution

covmeans_dif <- pileup_df[(1585052-3999):(1585052+4000),-1]

covmeans_rep[,-1] <- round(covmeans_rep[,-1], 0)

rm(pileup_df)


#### Sliding window full genome ####

covmeans_edgeR <- DGEList(counts=covmeans_rep[,-1],
                          group=colnames(covmeans_rep)[-1] %>% str_replace("_[0-9]", ""),
                          genes=data.frame(replicon = rep("NZ_CP009273.1", nrow(covmeans_rep)),
                                           position=covmeans_rep[,1]))

covmeans_edgeR <- calcNormFactors(covmeans_edgeR, method = "RLE")

# design: group sampling point
design <- model.matrix(~0+group, data = covmeans_edgeR$samples)
colnames(design) <- gsub("group", "", colnames(design))

# estimate dispersion
covmeans_edgeR <- estimateDisp(covmeans_edgeR)

plotBCV(covmeans_edgeR)

cov_contrast <- makeContrasts(
  xerC-WT, xerD-WT,
  levels=design)

covmeans_fit <- glmQLFit(covmeans_edgeR, design,
                         robust=TRUE,
                         dispersion = covmeans_edgeR$common.dispersion,
                         abundance.trend = TRUE)


lfc_matA <- matrix(NA, nrow(covmeans_fit$genes), ncol(cov_contrast))
fdr_matA <- matrix(NA, nrow(covmeans_fit$genes), ncol(cov_contrast))

colnames(lfc_matA) <- paste("lfc", colnames(cov_contrast), sep = "_")
colnames(fdr_matA) <- paste("fdr", colnames(cov_contrast), sep = "_")

for(i in 1:ncol(cov_contrast)) {
  treat_intermediate <- glmQLFTest(covmeans_fit, contrast = cov_contrast[,i])
  lfc_matA[,i] <- treat_intermediate$table$logFC
  fdr_matA[,i] <- p.adjust(treat_intermediate$table$PValue, method = "BH")
}

result_table <- cbind(covmeans_fit$genes,
                      WT = rowMeans(cpm(covmeans_fit)[,1:3]),
                      xerC = rowMeans(cpm(covmeans_fit)[,4:5]),
                      xerD = rowMeans(cpm(covmeans_fit)[,6:8]),
                      lfc_matA, fdr_matA)


write.table(result_table, file = "coverage_edgeR.txt", row.names = FALSE, sep = "\t")

covmeans_long <- gather(result_table, strain, coverage, WT:xerD, factor_key=TRUE)

ggplot(covmeans_long,
       aes(x = position/1000000, y = coverage, color = strain)) +
  geom_vline(xintercept = c(1.3, 1.585052, 1.9)) +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  theme_bw()
ggsave(filename="./plots_jt/ec_full_chr_cov.pdf",
       width = 8.5,
       height = 4.5,
       units = c("cm"))


covmeans_long %>% 
  filter(position > 1.3e+06 & position < 1.9e+06) %>%
  ggplot(aes(x = position/1000000, y = coverage, color = strain)) +
  geom_vline(xintercept = c(1.585052)) +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="./plots_jt/ec_zoom_chr_cov.pdf",
       width = 8.5,
       height = 4.5,
       units = c("cm"))  

covFC_long <- gather(result_table, strain, logFC, 6:7, factor_key=TRUE)

covFC_long %>% 
  filter(position > 1.3e+06 & position < 1.9e+06) %>%
  ggplot(aes(x = position/1000000, y = logFC, color = strain)) +
  geom_vline(xintercept = c(1.585052)) +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="./plots_jt/ec_zoom_chr_lfc.pdf",
       width = 8.5,
       height = 4.5,
       units = c("cm")) 

#### dif site zoom ####

dif_edgeR <- DGEList(counts=covmeans_dif[,-1],
                          group=colnames(covmeans_dif)[-1] %>% str_replace("_[0-9]", ""),
                          genes=data.frame(replicon = rep("NZ_CP009273.1", nrow(covmeans_dif)),
                                           position=covmeans_dif[,1]))

dif_edgeR <- calcNormFactors(dif_edgeR, method = "RLE")

# design: group sampling point
design <- model.matrix(~0+group, data = dif_edgeR$samples)
colnames(design) <- gsub("group", "", colnames(design))

# estimate dispersion
dif_edgeR <- estimateDisp(dif_edgeR)

plotBCV(dif_edgeR)

cov_contrast <- makeContrasts(
  xerC-WT, xerD-WT,
  levels=design)

dif_fit <- glmQLFit(dif_edgeR, design,
                    robust=TRUE,
                    dispersion = covmeans_edgeR$common.dispersion,
                    abundance.trend = TRUE)


lfc_matA <- matrix(NA, nrow(dif_fit$genes), ncol(cov_contrast))
fdr_matA <- matrix(NA, nrow(dif_fit$genes), ncol(cov_contrast))

colnames(lfc_matA) <- paste("lfc", colnames(cov_contrast), sep = "_")
colnames(fdr_matA) <- paste("fdr", colnames(cov_contrast), sep = "_")

for(i in 1:ncol(cov_contrast)) {
  treat_intermediate <- glmQLFTest(dif_fit, contrast = cov_contrast[,i])
  lfc_matA[,i] <- treat_intermediate$table$logFC
  fdr_matA[,i] <- p.adjust(treat_intermediate$table$PValue, method = "BH")
}

result_table_dif <- cbind(dif_fit$genes,
                          WT = rowMeans(cpm(dif_fit)[,1:3]),
                          xerC = rowMeans(cpm(dif_fit)[,4:5]),
                          xerD = rowMeans(cpm(dif_fit)[,6:8]),
                          lfc_matA, fdr_matA)

covdif_long <- gather(result_table_dif, strain, coverage, WT:xerD, factor_key=TRUE)


covdif_long %>% 
  filter(position > 1583700 & position < 1586200) %>%
  ggplot(aes(x = position/1000000, y = coverage, color = strain)) +
  geom_vline(xintercept = c(1.585005, 1.585052)) +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  theme_bw()
ggsave(filename="./plots_jt/ec_difsite_zoom.pdf",
       width = 4.5,
       height = 4,
       units = c("cm"))


covFC_dif_long <- gather(result_table_dif, strain, logFC, 6:7, factor_key=TRUE)

covFC_dif_long %>% 
  filter(position > 1583700 & position < 1586200) %>%
  ggplot(aes(x = position/1000000, y = logFC, color = strain)) +
  geom_vline(xintercept = c(1.585005, 1.585052)) +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="./plots_jt/ec_zoom_chr_lfc.pdf",
       width = 8.5,
       height = 4.5,
       units = c("cm")) 

covfdr_dif_long <- gather(result_table_dif, strain, fdr, 8:9, factor_key=TRUE)

covfdr_dif_long %>% 
  filter(position > 1583700 & position < 1586200) %>%
  ggplot(aes(x = position/1000000, y = fdr, color = strain)) +
  geom_vline(xintercept = c(1.585005, 1.585052)) +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="./plots_jt/ec_zoom_chr_lfc.pdf",
       width = 8.5,
       height = 4.5,
       units = c("cm")) 