# Compare sets of reactions added during gap filling
library(scales)
library(pheatmap)
library(wesanderson)

##########################################################
# intersection of gap-filling solutions
##########################################################

topDir <- "figures/added_reactions/"
writeToFile <- T

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
row_labels = c("individual", "COMMIT")
# ---------------------

# perform K-means clustering with 3 centers
clusters = kmeans(data,3,nstart=100)$cluster

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

lty_levels = c(3,5,1)
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