# Comparison of structural distances to 16S rRNA sequence similarity
library(igraph, warn.conflicts = F, quietly = T)
library(plotrix, warn.conflicts = F, quietly = T)
library(ape, warn.conflicts = F, quietly = T)

# Function definitions
readDistMatrix <- function(file) {
  mt <- read.table(file, header = T)
  rownames(mt) = mt[,1];
  mt = mt[,-1]
  return(mt)
}

colorFunction <- function(x) {
  if (x < 0.25) {return(1)} 
  else if (0.25 <= x & x < 0.5) {return(6)}
  else if (0.5 <= x & x < 0.75) {return(2)}
  else return(8)
}

# Compare distance measures applied to consensus models to sequence similarity
topDir = "COMMIT/"
writeToFile = F
setwd(paste0(topDir, "figures/Comparison-of-methods/"))

# Input files
soil_file = "dist-phylo-Soil.txt"
root_file = "dist-phylo-Root.txt"
leaf_file = "dist-phylo-Leaf.txt"

# Read the values from respective files
soil_mt <- readDistMatrix(soil_file)
root_mt <- readDistMatrix(root_file)
leaf_mt <- readDistMatrix(leaf_file)

# Calculate the mean of all rows
data <- data.frame("Mantel_r" = (apply(cbind(soil_mt$Mantel_r, root_mt$Mantel_r, leaf_mt$Mantel_r), 1, mean)))

# Graph connecting each distance measure with the phylogeny
n <- nrow(data) + 1
adj_mat <- matrix(rep(0,n^2), ncol = n, nrow = n)
adj_mat[,1] = adj_mat[1,] = 1;
adj_mat[1,1] = 0

# Color gradient according to edge weight (or thickness)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
G <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")

# Edge attributes
E(G)$weight <- data$Mantel_r
E(G)$color <- cbp1[unlist(lapply(E(G)$weight, colorFunction))]
E(G)$width <- abs(5 * E(G)$weight/max(E(G)$weight)) + 0.5

# vertex attributes
vertex_names <- c("sequence similarity", "SVD distance", "reaction JD",
                  "metabolite JD", "dead-end metabolites",
                  "dead-end metabolite JD", "E.C. JD", "E.C. correlation",
                  "cofactor usage")
distances = c(8.5, 2.3, 2, 2.2, 2, 2.5, 4.5, 6, 1.5)
degrees = c(2*pi, -pi/3, -pi/2, pi/(3/4), pi/(5/2), pi/(3/4), pi/(9/10), pi/(9/10), pi/2)
order = c(1,2,3,4,6,7,8,5,9)

V(G)$size <- 8
V(G)$color <- "black"
V(G)$label <- vertex_names
V(G)$label.dist <- distances
V(G)$label.degree <- degrees
V(G)$label.color = "black"
V(G)$label.cex = 1.5

# Save figure as PNG
if (writeToFile) {
  png(paste("Comb-dist-seq-corr.png", sep = ""), height = 20, width = 20, units = "cm", res = 300)
}

# Graphical parameters
par(family = "Arial", mar = c(1, 4, 1, 10))

plot.igraph(G, layout=layout_in_circle(G, order = order), family = "Arial", margin = 0.2)

labels = c(expression(rho %in% "[0.75, 1]"),
           expression(rho %in%  "[0.5, 0.75)"),
           expression(rho %in% "[0.25, 0.5)"),
           expression(rho %in% "[0, 0.25)"))
legend(x = 1.2, y = 1.2, legend = labels,
       fill =  cbp1[c(8,2,6,1)],
       border = NA, y.intersp = 1, cex = 1.5, box.lwd = 0, bg = NA)

if (writeToFile) dev.off()