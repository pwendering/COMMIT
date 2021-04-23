# generate plots for the comparison of different gap filling solutions
library(scales)
library(pheatmap)
library(wesanderson)

##########################################################
# intersection of gap-filling solutions
##########################################################

# topDir <- "~/ComGapFill/figures/added_reactions/"
topDir = "C://Users/wende/MobaXterm/home/comgapfill/figures/added_reactions/"
writeToFile <- T
# par(family = "sans")
# habitat <- "Soil"
# experiments <- c("Schlaeppi", "Bulgarelli")
# spec <- c("", "_no_exc_rxns")
#  
# palettelength <- 100
# col1 <- "beige"
# col2 <- "firebrick4"
# pal <- colorRampPalette(c(col1, col2))(palettelength)
# mybreaks = c(seq(from = 0,to = 1, length.out =  palettelength))
# 
# for (E in experiments) {
#   for (S in spec){
#     data <- read.table(paste0(topDir, habitat, "_", E, "_dist_methods", S, ".txt"),
#                        header = T, na.strings = "NaN")
#     row.names(data) <- data[,1]
#     data = data[,-1]
#     data[is.na(data)] = 0
#     data <- 1 - data
#     if (writeToFile) {
#       png(paste0(topDir, "jaccard_index_solutions_",
#               habitat, '_', E, S, '.png'), height = 20, width = 20, units = "cm", res = 600)
#     }
#     labels <- c("KBase draft models cond.", "consensus - C cond.", "consensus cond.",
#                 "KBase draft models ind.", "consensus - C ind.", "consensus ind.")
#      c <- pheatmap(as.matrix(data),
#              cluster_cols = T,
#              cluster_rows = F,
#              labels_col = labels,
#              kmeans_k = 3,
#              main = "",
#              breaks = mybreaks,
#              fontsize = 18,
#              fontsize_col = 16,
#              family = "Arial",
#              border_color = NA,
#              color = pal,
#              angle_col = 315, 
#              gaps_row = c(1,2),
#              clustering_method = "average",
#              show_rownames = F,
#              width = 12,
#              height = 8
#              )
#     
#      if (writeToFile) dev.off()
#   }
# }

# ----------------- Alternative figure ----------------- #

library(igraph)

# store default graphical parameters
originalPar = par()
par(mar = c(3,5,3,0)+.1, family = "sans")

experiment = "Bulgarelli"
habitat = "Soil"

data <- read.table(paste0(topDir, habitat, "_",
                         experiment, "_dist_methods.txt"),
                   header = T, na.strings = "NaN", row.names = 1)
data[is.na(data)] = 0
data = 1 - data

# ----- HARDCODED -----
vertex_labels = rep(c("KBase", "-CarveMe", "consensus"),2)
row_labels = c("individual", "CompFill")
# ---------------------

# perform K-means clustering with 3 centers
clusters = kmeans(data,3)$cluster

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cluster_cols = cbp1[clusters]

# create a graph using the Jaccard similarity matrix as adjacency matrix
G = graph_from_adjacency_matrix(as.matrix(data), weighted = T, diag = F,
                                mode = "undirected")


# appearance of vertices
n = length(V(G))
V(G)$frame.color <- "grey40"
V(G)$color = cluster_cols
V(G)$size <- 30
V(G)$label <- vertex_labels
V(G)$label.dist <- rep(5,n)
V(G)$label.degree <- c(rep(.5*pi,n/2),rep(1.5*pi,n/2))
V(G)$label.color <- "black"
V(G)$label.font <- 4
V(G)$label.cex <- 1.2
V(G)$shape = "square"
V(G)$label.family <- "sans"



# appearance of edges
bins = cut(E(G)$weight,3, include.lowest = T)

lty_levels = c(3,5,1) #c(3,4,2,5,1)
edge_lty = as.numeric(paste(cut(E(G)$weight,3, include.lowest = T, labels = lty_levels)))

E(G)$lty <- edge_lty
E(G)$color <- "black"
E(G)$width <- 2
curvature = rep(0,length(E(G)))
curvature[c(2,14)] = c(.5,-.5)

if (writeToFile) {
  png(filename = paste0(topDir, "graph_jaccard_similarity_", habitat, "_",
                        experiment, ".png"),
      width = 12, height = 12, units = "cm", res = 600)
}

# plot graph
plot.igraph(G, layout = layout_on_grid(G, width = n/2), edge.curved = curvature)

# add labels for gap-filling type
text(y = c(1,-1),x = c(-1.8,-1.8), labels = row_labels, font = 2) 

# add legend
legend(y=0,x=-2.3,legend = levels(bins), cex = .8, lty = lty_levels,
       yjust = .5, bty = "n", lwd = 2)

if (writeToFile) dev.off()

# restore orignial graphical parameters
par(originalPar)