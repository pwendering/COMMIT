#!/usr/bin/Rscript
# Heatmap plots for distance matrices

library(pheatmap)
library(wesanderson)

args = commandArgs(trailingOnly = T)
inFile = args[1]
outFile = args[2]

# read file
mt <- read.table(inFile, header = T, sep = ",")
rownames(mt) = colnames(mt)
diag(mt) = NA

# Plotting parameters
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palettelength = 20
# mytitle = "Structural differences of the applied reconstruction methods"
# col1 = "#999999"#"gray92"
# col2 = "#E69F00"#"cornflowerblue"
# col3 = "#56B4E9"#"gold"
# col4 = "firebrick4"
# mycol = colorRampPalette(c(col1, col2, col3, col4))(palettelength)
pal <- wes_palette("Zissou1", palettelength, type = "continuous")
mybreaks = c(seq(from = 0,to = 1, length.out =  palettelength))

# plot pheatmap
pheatmap(mt,
         breaks = mybreaks,
         cluster_cols = F,
         cluster_rows = F,
         color = pal,
         fontsize = 18,
         fontsize_col = 16,
         fontsize_row = 16,
         width = 6, 
         height = 5,
         angle_col = 45,
         border_color = NA,
         na_col = "#999999",
         family = "Arial",
         filename = outFile
)
