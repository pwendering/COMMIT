# compare model features
library(plotrix, warn.conflicts = F, quietly = T)
library(scales)

rm(list = ls())
habitats = c("Soil", "Leaf", "Root")
topDir = "~/ComGapFill/figures/model-features/"
writeToFile = F

for (i in 1:length(habitats)) {
  if (!exists("rxns")) {
    rxns <- read.table(paste(topDir, "rxns-", habitats[i], ".txt", sep = ""),
                     header = T)
  } else {
    rxns <- rbind(rxns, read.table(paste(topDir, "rxns-", habitats[i], ".txt", sep = ""),
                                       header = T))
  }
  
  if (!exists("mets")) {
    mets <- read.table(paste(topDir, "mets-", habitats[i], ".txt", sep = ""),
                       header = T)
  } else {
    mets <- rbind(mets, read.table(paste(topDir, "mets-", habitats[i], ".txt", sep = ""),
                                   header = T))
  }
  
  if (!exists("genes")) {
    genes <- read.table(paste(topDir, "genes-", habitats[i], ".txt", sep = ""),
                       header = T)
  } else {
    genes <- rbind(genes, read.table(paste(topDir, "genes-", habitats[i], ".txt", sep = ""),
                                   header = T))
  }
}

  # Boxplot comparing consensus with sums of all models for each feature
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#009E73")
  
  box_data <- data.frame(rxns$consensus_rxns, rxns$sums_rxns,
                         mets$consensus_mets, mets$sums_mets,
                         genes$consensus_genes, genes$sums_genes)
  
  x_pos <- c(1,1.8,4,4.8,7,7.8)
  y_max <- 1.1*max(box_data)
  cex <- 1.2
  
  if (writeToFile) {
    png(paste(topDir, "Comp_features.png", sep = ""), height = 20, width = 20, units = "cm", res = 300)
  }
  par(bty = "n", family = "Arial", mai = rep(1,4))
  boxplot(box_data, at = x_pos, ylim = c(0,y_max),
          col = rep(c(cbp1[c(1,2)]),3),
          xaxt = "n", bty = "n",
          ylab = "Number", cex.axis = cex, cex.lab = cex)
  
  # Points
  alpha = 0.1
  pch = 4
  # reactions
  rxns <- c(rxns$CarveMe_rxns, rxns$KBase_rxns, rxns$RAVEN_rxns, rxns$AuReMe_rxns)
  points(x = c(rep(x_pos[2]-0.1, length(rxns)/4),
               rep(x_pos[2]-0.05, length(rxns)/4),
               rep(x_pos[2]+0.05, length(rxns)/4),
               rep(x_pos[2]+0.1, length(rxns)/4)
               ),
         y = rxns, col = alpha(c(rep(cbp1[3], length(rxns)/4),
                           rep(cbp1[4], length(rxns)/4),
                           rep(cbp1[5], length(rxns)/4),
                           rep(cbp1[6], length(rxns)/4)), alpha),
         pch = pch
         )
  
  # metabolites
  mets <- c(mets$CarveMe_mets, mets$KBase_mets, mets$RAVEN_mets, mets$AuReMe_mets)
  points(x = c(rep(x_pos[4]-0.15, length(mets)/4),
               rep(x_pos[4]-0.05, length(mets)/4),
               rep(x_pos[4]+0.05, length(mets)/4),
               rep(x_pos[4]+0.15, length(mets)/4)
  ),
  y = mets, col = alpha(c(rep(cbp1[3], length(mets)/4),
                    rep(cbp1[4], length(mets)/4),
                    rep(cbp1[5], length(mets)/4),
                    rep(cbp1[6], length(mets)/4)), alpha),
  pch = pch
  )
  
  # genes
  genes <- c(genes$CarveMe_genes, genes$KBase_genes, genes$RAVEN_genes, genes$AuReMe_genes)
  points(x = c(rep(x_pos[6]-0.1, length(genes)/4),
               rep(x_pos[6]-0.05, length(genes)/4),
               rep(x_pos[6]+0.05, length(genes)/4),
               rep(x_pos[6]+0.1, length(genes)/4)
  ),
  y = genes, col = alpha(c(rep(cbp1[3], length(genes)/4),
                    rep(cbp1[4], length(genes)/4),
                    rep(cbp1[5], length(genes)/4),
                    rep(cbp1[6], length(genes)/4)), alpha),
  pch = pch
  )
  
  axis(1, at = c(1.4, 4.4, 7.4), labels = c("Reactions", "Metabolites", "Genes"),
       font = 2, cex.axis = cex)
  legend("topright", legend = c("Sum", "Consensus", "CarveMe", "KBase", "RAVEN 2.0", "AuReMe"),
         fill = cbp1[c(2,1,3,4,5,6)], bty = "n", cex = cex)

  if (writeToFile) dev.off()
