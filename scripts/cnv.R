rm(list=ls(all=TRUE))

args <- commandArgs(trailingOnly = TRUE)

omimcnv <- args[1]
coverage <- args[2]
pohlavi <- args[3]

omimgeny <- read.table(omimcnv, header = TRUE, sep = "\t")
coverage <- read.table(coverage, header = TRUE, sep = "\t",check.names = FALSE)
pohlavi <- read.table(pohlavi, header = TRUE, sep = ";")

pohlavi <- data.frame(pohlavi, row.names = 1)
  
  
  # rozdeleni tabulky na M a Z
  for (i in rownames(pohlavi)){
    ID <- pohlavi[i,]
    colnames(coverage)[which(names(coverage) == i)] <- paste(ID, i, sep="_")
  }
  
  M <- colnames(coverage)[grepl("M_",colnames(coverage))]
  Z <- colnames(coverage)[grepl("Z_",colnames(coverage))]
  coverage_M <- coverage[,c(M)]
  coverage_Z <- coverage[,c(Z)]
  
  # vypocet medianu pro kazdeho pacienta + normalizace vydelenim pokryti medianem -> standardizovane pokryti (NP)
  coverage_M_norm <- coverage_M/vapply(coverage_M, median, numeric(1))[col(coverage_M)]
  coverage_Z_norm <- coverage_Z/vapply(coverage_Z, median, numeric(1))[col(coverage_Z)]
  
  # prumer NP pro jednotlive oblasti +  vypocet mean(NP) - NP
  #coverage_M_mean <- rowMeans(coverage_M_norm) - coverage_M_norm
  #coverage_Z_mean <- rowMeans(coverage_Z_norm) - coverage_Z_norm
  
  # prumer NP pro jednotlive oblasti +  vypocet NP - mean(NP)
  coverage_M_mean <- coverage_M_norm - rowMeans(coverage_M_norm)
  coverage_Z_mean <- coverage_Z_norm - rowMeans(coverage_Z_norm)
  
  # vytvoreni finalnich tabulek s hlavickou
  coverage_M_final <- cbind(coverage[,c("chr","start","stop","name")], seq.int(nrow(coverage_M_mean)), coverage_M_mean)
  coverage_Z_final <- cbind(coverage[,c("chr","start","stop","name")], seq.int(nrow(coverage_M_mean)), coverage_Z_mean)
  
  # ulozeni tabulek
  #write.table(coverage_M_final,"coverage_M_final.csv",col.names=TRUE, row.names=FALSE, sep=";")
  #write.table(coverage_Z_final,"coverage_Z_final.csv",col.names=TRUE, row.names=FALSE, sep=";")
  
  
  
  # vytvoreni tabulky se vsemi a celymi radky s extremnimi hodnotami (>0.25) + ulozeni
  
  greater_M <- coverage_M_final[rowSums(coverage_M_final[6:ncol(coverage_M_final)] > 0.25 | coverage_M_final[6:ncol(coverage_M_final)] < -0.25) > 0, ]
  greater_Z <- coverage_Z_final[rowSums(coverage_Z_final[6:ncol(coverage_Z_final)] > 0.25 | coverage_Z_final[6:ncol(coverage_Z_final)] < -0.25) > 0, ]
  
  # vlozi omim info
  names_M <- omimgeny$phenotyp[match(greater_M$name,omimgeny$gene)]
  names_Z <- omimgeny$phenotyp[match(greater_Z$name,omimgeny$gene)]
  
  #prida sloupec names jako sloupec "gene" k expresi
  greater_M$OMIM <- names_M
  greater_Z$OMIM <- names_Z
  
  
  
  write.table(greater_M,"CNV_M.csv",col.names=TRUE, row.names=FALSE, sep=";")
  write.table(greater_Z,"CNV_Z.csv",col.names=TRUE, row.names=FALSE, sep=";")
  
