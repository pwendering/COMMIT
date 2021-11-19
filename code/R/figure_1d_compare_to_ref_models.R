# Comparison of merged models to reference models

library(plotrix)
library(scales)

writeToFile = T
topDir = "figures/Comparison-to-reference-models/"

#####################################################################
#######  Data
#####################################################################
spec = ""
# scaling factor to plot precision
scale_factor_prec = 4;
# Bacillus megaterium
Bm <- read.table(paste0(topDir, spec, "B_megaterium.txt"),header = T,row.names = 1)
Bm_genus <- as.logical(Bm$genus)
Bm_species <- as.logical(Bm$species)
Bm$sensitivity <- Bm$sensitivity*Bm$relation
Bm$precision <- scale_factor_prec*Bm$precision*Bm$relation
Bm = Bm[,which(colnames(Bm)==c("precision", "sensitivity"))]

# Methylobacterium extorquens
Me <- read.table(paste0(topDir, spec, "M_extorquens.txt"), header = T,row.names = 1)
Me_genus <- as.logical(Me$genus)
Me_species <- as.logical(Me$species)
Me$sensitivity <- Me$sensitivity*Me$relation
Me$precision <- scale_factor_prec*Me$precision*Me$relation
Me = Me[,which(colnames(Me)==c("precision", "sensitivity"))]

#####################################################################
#######  Plotting
#####################################################################
if (writeToFile) {
  png(paste0(topDir, spec, "comparison_reference_models.png"),
      height = 800, width = 800)
}
par(xpd=T,mar = c(5.1,4.1,4.1,2.1)+c(2,2,4,6))
my_title = "Similarity of consensus models to reference models"
yaxis_label_1 = "Sensitivity x sequence similarity"
yaxis_label_2 = "Precision x sequence similarity"
cex = 2.3
col_genus = "black"
col_species = "red"
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###### plotrix violin plot

# precision data are multiplies by four to match the plotting range of sensitivity values
data = cbind(Bm,Me)
violin_plot(data,
            col = c("#999999", "#E69F00", "#999999", "#E69F00"),
            main = "",
            at = seq(1,4),
            axes = F,
            ylim = c(0,1.1))

# add points to the plot for species and genus values
# same genus
length(Bm$sensitivity[Bm_species])
length(Bm$precision[Bm_species])
length(Me$sensitivity[Me_species])
length(Me$precision[Me_species])
points(x = c(rep(1, length(Bm$sensitivity[Bm_genus])),
             rep(2, length(Bm$precision[Bm_genus])),
             rep(3, length(Me$sensitivity[Me_genus])),
             rep(4, length(Me$precision[Me_genus]))
),
y = c(Bm$sensitivity[Bm_genus],
      Bm$precision[Bm_genus],
      Me$sensitivity[Me_genus],
      Me$precision[Me_genus]
),
col = col_genus,
pch = 16
)

# same species
points(x = c(rep(1, length(Bm$sensitivity[Bm_species])),
             rep(2, length(Bm$precision[Bm_species])),
             rep(3, length(Me$sensitivity[Me_species])),
             rep(4, length(Me$precision[Me_species]))
),
y = c(Bm$sensitivity[Bm_species],
      Bm$precision[Bm_species],
      Me$sensitivity[Me_species],
      Me$precision[Me_species]
),
col = col_species,
pch = 16
)

# add left Y-axis
par(yaxt = "s")
axis(side=2,at = pretty(c(0,1.1)),
     col = "#666666", col.axis = "#666666", cex.axis = cex, cex.lab = cex,lwd=3)
title(ylab = yaxis_label_1, cex.lab = cex, col.lab = "#666666",
      family = "sans", font.lab = 2,line = 4)

# add right Y-axis
axis(side=4,at = pretty(c(0,1.1)),labels = pretty(c(0,1.1))/scale_factor_prec,padj = .5,
     col = "#E69F00", col.axis = "#E69F00", cex.axis = cex, cex.lab = cex,lwd=3)
mtext(yaxis_label_2, side = 4, col = "#E69F00", line = 5, cex = cex, font = 2)

# add X-axis labels
axis(side = 1, at = c(1.5,3.5), labels = c("B. megaterium", "M. extorquens"),
     col = "#999999", lwd = -1, font = 4, cex.axis = cex)

# add legend
legend("topright",c("same genus", "same species"),
       pch = 19, col = c(col_genus, col_species),
       bty = "n", cex = cex, pt.cex = 1.5)

if (writeToFile) dev.off()