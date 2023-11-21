#### Load libraries ####

library(zoo)
library(tidyverse)
library(RColorBrewer)
library(edgeR)

#### Load pileup files ####

pileup_files <- list.files(path = "/Juergen/Projects/vesicles_Ecoli/ncbi_sra/",
                           full.names = TRUE) %>% grep(".pileup", ., value = TRUE)

pileup_list <- list()

for(i in 1:length(pileup_files)) {
  pileup_list[[i]] <- read_tsv(pileup_files[i], 
                               col_names = c("replicon", "position", "coverage"))
}

names(pileup_list) <- basename(pileup_files) %>% gsub(".pileup", "", .)

rm(pileup_files)

pileup_list$Vibrio_cholerae_chr1 <- pileup_list$Vibrio_cholerae %>%
  filter(replicon == "NC_009457.1")

pileup_list$Vibrio_cholerae_chr2 <- pileup_list$Vibrio_cholerae %>%
  filter(replicon == "NC_009456.1")

pileup_list <- pileup_list[-4]

covmeans_rep <- list()

for(i in 1:length(pileup_list)) {
  covmeans_rep[[i]] <- as.data.frame(rollapply(pileup_list[[i]][,-1],
                                                    FUN = mean, width = 500, by = 250, align = "center", partial = 50))
  colnames(covmeans_rep[[i]]) <- colnames(covmeans_rep[[i]])
}

names(covmeans_rep) <- names(pileup_list)

strain_info <- data.frame(strain = names(covmeans_rep),
                          dif = c(829000, 829000, 2443067, 1129240, 567632))

covmeans <- covmeans_rep %>% bind_rows(.id = "strain")


covmeans %>%
  filter(strain == "Prochlorococcus_Med4_01") %>%
ggplot(aes(x = position/1000000, y = log2(coverage+2))) +
  geom_vline(xintercept = c(strain_info$dif[1]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  ylim(0,15) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Prochlorococcus.pdf",
       width = 9,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Pseudomonas_aeruginosa_01") %>%
  ggplot(aes(x = position/1000000, y = log2(coverage+8))) +
  geom_vline(xintercept = c((strain_info$dif[3]-99999)/1000000,
                            (strain_info$dif[3]+100000)/1000000,
                            strain_info$dif[3]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 1) +
  ylim(7,11) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Pseudomonas.pdf",
       width = 9,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Pseudomonas_aeruginosa_01") %>%
  filter(position > strain_info$dif[3]-99999 & position < strain_info$dif[3]+100000) %>%
  ggplot(aes(x = position/1000000, y = log2(coverage+2))) +
  geom_vline(xintercept = strain_info$dif[3]/1000000, color = "red") +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Pseudomonas_aeruginosa_zoom.pdf",
       width = 5,
       height = 4.5,
       units = c("cm"))  

covmeans %>%
  filter(strain == "Vibrio_cholerae_chr1") %>%
  ggplot(aes(x = position/1000000, y = log2(coverage+2))) +
  geom_vline(xintercept = c(strain_info$dif[4]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  ylim(3,10) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Vibrio_cholerae_chr1.pdf",
       width = 6.5,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Vibrio_cholerae_chr2") %>%
  ggplot(aes(x = position/1000000, y = log2(coverage+2))) +
  geom_vline(xintercept = c(strain_info$dif[5]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  ylim(3,10) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Vibrio_cholerae_chr2.pdf",
       width = 3.5,
       height = 4.5,
       units = c("cm"))


covmeans %>%
  filter(strain == "Prochlorococcus_Med4_01") %>%
  ggplot(aes(x = position/1000000, y = coverage/1000)) +
  geom_vline(xintercept = strain_info$dif[1]/1000000, color = "red") +
  geom_line(show.legend = FALSE, lwd = 1) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Prochlorococcus_nt.pdf",
       width = 9,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Pseudomonas_aeruginosa_01") %>%
  ggplot(aes(x = position/1000000, y = coverage/100)) +
  geom_vline(xintercept = c((strain_info$dif[3]-99999)/1000000,
                            (strain_info$dif[3]+100000)/1000000,
                            strain_info$dif[3]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 1) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Pseudomonas_nt.pdf",
       width = 9,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Pseudomonas_aeruginosa_01") %>%
  filter(position > strain_info$dif[3]-99999 & position < strain_info$dif[3]+100000) %>%
  ggplot(aes(x = position/1000000, y = coverage/100)) +
  geom_vline(xintercept = strain_info$dif[3]/1000000, color = "red") +
  geom_line(show.legend = TRUE) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Pseudomonas_aeruginosa_zoom_nt.pdf",
       width = 5,
       height = 4.5,
       units = c("cm"))  

covmeans %>%
  filter(strain == "Vibrio_cholerae_chr1") %>%
  ggplot(aes(x = position/1000000, y = coverage/1000)) +
  geom_vline(xintercept = c(strain_info$dif[4]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  ylim(0,1) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Vibrio_cholerae_chr1_nt.pdf",
       width = 7.0,
       height = 4.5,
       units = c("cm"))

covmeans %>%
  filter(strain == "Vibrio_cholerae_chr2") %>%
  ggplot(aes(x = position/1000000, y = coverage/1000)) +
  geom_vline(xintercept = c(strain_info$dif[5]/1000000), color = "red") +
  geom_line(show.legend = FALSE, lwd = 0.5) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  ylim(0,1) +
  theme_bw()
ggsave(filename="/Juergen/Projects/vesicles_Ecoli/ncbi_sra/Vibrio_cholerae_chr2_nt.pdf",
       width = 3.75,
       height = 4.5,
       units = c("cm"))



plot(covmeans_rep$Pseudomonas_aeruginosa_01$position,
     covmeans_rep$Pseudomonas_aeruginosa_01$coverage,
     xlim = c(2e+06, 3e+06))
abline(v = 2443067)

plot(covmeans_rep$Prochlorococcus_Med4_01$position,
     covmeans_rep$Prochlorococcus_Med4_01$coverage)
abline(v = 829000)

plot(covmeans_rep$Prochlorococcus_Med4_02$position,
     covmeans_rep$Prochlorococcus_Med4_02$coverage)
abline(v = 829000)


plot(covmeans_rep$Vibrio_cholerae_chr1$position,
     covmeans_rep$Vibrio_cholerae_chr1$coverage, ylim = c(0,1000))
abline(v = 1129240)

plot(covmeans_rep$Vibrio_cholerae_chr2$position,
     covmeans_rep$Vibrio_cholerae_chr2$coverage, ylim = c(0,500))
abline(v = 567632)
