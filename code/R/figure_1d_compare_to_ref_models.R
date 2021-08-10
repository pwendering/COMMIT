# Comparison of merged models to reference models

library(plotrix)
library(scales)

writeToFile = T
topDir = "figures/Comparison-to-reference-models/"

#####################################################################
#######  Data
#####################################################################
spec = ""
# Bacillus megaterium
Bm <- read.table(paste0(topDir, spec, "B_megaterium.txt"),header = T)
rownames(Bm) <- Bm[,1]
Bm <- Bm[,-1]
Bm_genus <- as.logical(Bm$genus)
Bm_species <- as.logical(Bm$species)
Bm$sensitivity <- Bm$sensitivity*Bm$relation
Bm$precision <- Bm$precision*Bm$relation
Bm = Bm[,which(colnames(Bm)==c("precision", "sensitivity"))]

# Methylobacterium extorquens
Me <- read.table(paste0(topDir, spec, "M_extorquens.txt"), header = T)
rownames(Me) <- Me[,1]
Me <- Me[,-1]
Me_genus <- as.logical(Me$genus)
Me_species <- as.logical(Me$species)
Me$sensitivity <- Me$sensitivity*Me$relation
Me$precision <- Me$precision*Me$relation
Me = Me[,which(colnames(Me)==c("precision", "sensitivity"))]

#####################################################################
#######  Plotting
#####################################################################
if (writeToFile) {
  png(paste0(topDir, spec, "comparison_reference_models.png"),
      height = 800, width = 800*1.1)
}


my_title = "Similarity of consensus models to reference models"
yaxis_label_1 = "Sensitivity x sequence similarity"
yaxis_label_2 = "Precision x sequence similarity"
cex = 2.3
my_xaxis_labels = c("B. megaterium", "B. megaterium", "M. extorquens", "M. extorquens")
col_genus = "black"
col_species = "red"
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###### plotrix violin plot
# precision data are multiplies by four to match the plotting range of sensitivity values
data = data.frame(Bm$sensitivity, Bm$precision*4, Me$sensitivity,  Me$precision*4, rep(-1, dim(Bm)[1]))
par(yaxt = "n", tck = 0, bty = "n", mar = c(2.5, 6, 4, 6) + 0.1,  mgp = c(3.5,1,0))
plot.new()
y_axis_limits = c(0,round(max(c(1, data$Bm.sensitivity, data$Me.sensitivity, data$Bm.precision, data$Me.precision)),1))#c(0, max(data))
lim_factor = y_axis_limits[2] / max(data$Bm.sensitivity, data$Me.sensitivity)
par(new = TRUE)

# Bacillus megaterium
violin_plot(X = as.data.frame(data[,c(1,5,3,5)]),
            at = c(1,2,3,4),
            col = c("#999999", "#999999", "#999999", "#999999"),
            main = "",
            x_axis_labels = rep("",4),
            axes = F,
            ylim = y_axis_limits
            )

# Methylobacterium extorquens
par(new = TRUE)
violin_plot(X = 
              data[,c(5,2,5,4)],
            at = c(1,2,3,4),
            main = "",
            col = c("#E69F00", "#E69F00", "#E69F00", "#E69F00"),
            ylim = y_axis_limits,
            x_axis_labels = rep("",4),
            axes = F
)

# add left Y-axis
par(yaxt = "s", tck = NA)
axis(side=2, at = pretty(c(0,0.9)),
     col = "#666666", col.axis = "#666666", cex.axis = cex, cex.lab = cex)
title(ylab = yaxis_label_1, cex.lab = cex, col.lab = "#666666",
      family = "sans", font.lab = 2)

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
            Bm$precision[Bm_genus]*4,
            Me$sensitivity[Me_genus],
            Me$precision[Me_genus]*4
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
            Bm$precision[Bm_species]*4,
            Me$sensitivity[Me_species],
            Me$precision[Me_species]*4
      ),
      col = col_species,
      pch = 16
)


# add transparent plot to set correct Y-limits for the precision
par(new = TRUE, yaxt = "n")
plot(c(data$Bm.precision...4/4, data$Me.precision...4/4), xlab = "", ylab = "", axes = F,
     ylim = c(0,max(data)/4),
     xlim = c(1,4),
     col = alpha("white", 0))

# add the box, secondary Y-axis and the X-axis
par(yaxt = "s", bty = "o")

axis(side=4, at = pretty(c(0, 0.3)),
     col = "#E69F00", col.axis = "#E69F00", cex.axis = cex, cex.lab = cex, mgp = c(3,1.5,0))
mtext(yaxis_label_2, side = 4, col = "#E69F00", line = 4, cex = cex, font = 2)

axis(side = 1, at = c(1.5,3.5), labels = c("B. megaterium", "M. extorquens"),
     col = "#999999", lwd = -1, font = 4, cex.axis = cex)

# add legend
legend("topright",c("same genus", "same species"),
       pch = 19, col = c(col_genus, col_species),
       bty = "n", cex = cex, pt.cex = 1.5)

if (writeToFile) dev.off()