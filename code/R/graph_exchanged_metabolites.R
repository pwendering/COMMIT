# plot a growth displaying the exchanged metabolites between all bacterial families
library(igraph)
library(scales)

plotExchangeGraph <- function(fileBaseName, classFile, experiment, outPath, writeToFile = F) {
  
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  fileName <- paste(fileBaseName, experiment, ".txt", sep = "")
  data <- read.table(fileName, header = T)
  
  # classification
  fileName <- paste(classFile, experiment, ".txt", sep = "")
  classification <- read.table(fileName, header = T, sep = "\t", na.strings = "")
  row.names(classification) = classification[,1]
  classification <- classification[,-1]
  classes <- unique(as.vector(as.matrix(classification)))
  

  e_colors <- c()
  for (i in 1:length(classes)) {
    if (is.na(classes[i])) {
      idx <- which(is.na(classification))
    } else {
    idx <- which(as.vector(as.matrix(classification))==classes[i])
    }
    e_colors[idx] <- cbp1[i]
  }
  
  # graph
  adj_mat <- as.matrix(data)
  adj_mat <- adj_mat / mean(adj_mat)
  G <- graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = T)
  classes <- unique(as.vector(classification[adj_mat>0]))
  classes[is.na(classes)] <- "Other"
  fontsize <- 1.2
  
  E(G)$width <- 2*E(G)$weight
  palette <- colorRampPalette(c("lightblue", "firebrick"))
  E(G)$color <- alpha(as.vector(t(matrix(e_colors, nrow = nrow(classification), ncol = ncol(classification))))[which(t(adj_mat)>0)],0.8)
  E(G)$arrow.width = 0.2*E(G)$weight
  E(G)$arrow.size = 1.5
  
  V(G)$color <- "white"
  V(G)$frame.color <- "white"
  V(G)$size = 25
  # V(G)$label.dist <- c(2, 5, 2, 4, 4, 4, 1.5, 5)
  # V(G)$label.degree <- c(length(V(G)):1)/length(V(G))*2*pi
  # V(G)$label.degree[c(2, 4, 6, 8)] <- c(2*pi, pi, pi, 2.1*pi)
  V(G)$label.color <- "black"
  V(G)$label.font <- 4
  V(G)$label.cex <- fontsize
  V(G)$label.family <- "Arial"
  
  outFileName <- paste(outPath, "graph_", experiment, ".png", sep = "")
  if (writeToFile) {
    png(filename = outFileName,
        units = "cm", width = 20, height = 20, res = 300)
  }
  par(family = "Arial", mar = c(0, 4, 1, 2) + 0.1)
  
  plot.igraph(G, layout = layout_in_circle(G), order = c(1:length(V(G))), edge.curved = 0.2, margin = 0.25)
  legend(x=0.7,y=1.5, legend = classes,
         fill = cbp1[unique(as.vector(as.matrix(classification))) %in% classes],
         border = "NA", bty = "o", bg = NA, box.lwd = 0.7, cex = fontsize)
  if (writeToFile) dev.off()
}

writeToFile = T
topDir = "~/ComGapFill/figures"

### Consensus
# Soil
# modelType <- "all"
# fileBaseName <- paste(topDir, "/exchanged_metabolites/graph/", modelType, "_exchanged_metabolites_", sep = "")
# outPath <- paste(topDir, "/exchanged_metabolites/graph/", modelType, "_", sep = "")
# classFile <- paste(topDir, "/exchanged_metabolites/graph/", modelType,
#                    "_brite_exchanged_", sep = "")
# 
# plotExchangeGraph(fileBaseName = fileBaseName, classFile = classFile, experiment = "Schlaeppi", outPath = outPath, writeToFile = writeToFile)
# plotExchangeGraph(fileBaseName = fileBaseName, classFile = classFile, experiment = "Bulgarelli", outPath = outPath, writeToFile = writeToFile)

# Root
modelType <- "all"
fileBaseName <- paste(topDir, "/exchanged_metabolites/graph/Root_", modelType, "_exchanged_metabolites_", sep = "")
outPath <- paste(topDir, "/exchanged_metabolites/graph/Root_", modelType, "_", sep = "")
classFile <- paste(topDir, "/exchanged_metabolites/graph/Root_", modelType,
                   "_brite_exchanged_", sep = "")

plotExchangeGraph(fileBaseName = fileBaseName, classFile = classFile, experiment = "Schlaeppi", outPath = outPath, writeToFile = writeToFile)
plotExchangeGraph(fileBaseName = fileBaseName, classFile = classFile, experiment = "Bulgarelli", outPath = outPath, writeToFile = writeToFile)

