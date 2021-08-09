# compare the numbers of added reactions between different gap-filling solutions

library(scales)
library(plotrix)

writeToFile = T

habitat <- "Soil"
experiment <- "Schlaeppi"

figureDir <- "figures/added_reactions/"

# specify file for figure
outFile <- paste(figureDir, habitat, "_", experiment, "_gf_sol_size.png", sep = "")

# read data from file
data <- read.table(paste(figureDir, habitat, "_",experiment, "_gf_sol_size.txt", sep = ""), header = T)
row.names(data) <- data[,1]
data <- data[,-1]

# iterative
idx_iter <- c(1,2,3)
# individual
idx_ind <- idx_iter + length(idx_iter)

# determine means and variances for each type of gap filling
av <- apply(X = data, MARGIN = 1, FUN = mean)
# error for error bars
error <- apply(X = data, MARGIN = 1, FUN = sd)

# test for difference in means
h <- rep(0, length(idx_iter))
alpha <- 0.05
for (i in 1:length(idx_iter)){
  d_1 <- as.numeric(data[idx_iter[i],])
  d_2 <- as.numeric(data[idx_ind[i],])
  n_1 <- shapiro.test(d_1)
  n_2 <- shapiro.test(d_2)
  if (n_1$p.value >= alpha & n_2$p.value >= alpha) {
    t <- t.test(d_1,d_2,alternative = "two.sided", paired = T, var.equal = F)
  } else {
    t <- wilcox.test(d_1,d_2,paired = T, alternative = "two.sided")
  }
  if (t$p.value < 0.001) {h[i] <- 3}
  else if (t$p.value < 0.01) {h[i] <- 2}
  else if (t$p.value < 0.05) {h[i] <- 1}
  print(t$p.value)
}

# limits
x_limits <- c(0,3)
# round y upper limit to next 10
y_max <- max(av+error) / 10
y_max <- ceiling(y_max) * 10 + 10 
y_limits <- c(0, y_max)

if (writeToFile) {
  png(filename = outFile, units = "cm", width = 12, height = 12, res = 600)
}

par(family = "sans", mgp = c(2.8,1,0), mar = c(4,4,4,2)+.1, xpd = T)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
violin_colors <- alpha(rep(cbp1[c(2,1)], 3), 0.8)
cex = 1.1
x_pos = c(1,2,3.5,4.5,6,7)
violin_plot(X = t(data)[,c(4,1,5,2,6,3)], at = x_pos, axes = F, main = "",
            col = violin_colors, ylim = y_limits)

# plot axes
x_pos <- c(1.5, 4, 6.5)
axis(side = 1, at = x_pos, labels = c("KBase", "-CarveMe", "consensus"),
     lty = 0, font = 4, tick = T, cex.axis = 1.2)
axis(side = 2, at = pretty(y_limits), cex.axis = cex, col = "gray30")
title(ylab = "Number of added reactions", cex.lab = cex)

# legend
legend(x = x_pos[2]*1.2, y = y_max*0.8, legend = c("individual", "COMMIT"),  bty = "n",
       cex = cex, col = violin_colors[1,3,5], fill = violin_colors[c(1,2)], text.col = "gray30")

box(lwd=2)

# draw significance stars
x_pos <- c(2, 4.5, 7)
for (i in 1:length(h)){
  if (h[i]==1) {
    points(x = x_pos[i], y = 3+max(data[idx_iter[i],]),
           pch = 8, cex = 0.4)
  } 
  else if (h[i]==2) {
    points(x = c(x_pos[i]-0.1,x_pos[i]+0.1), y = rep(3+max(data[idx_iter[i],]),2),
           pch = 8, cex = 0.4)
  }
  else if (h[i]==3) {
    points(x = c(x_pos[i]-0.1,x_pos[i],x_pos[i]+0.1), y = rep(3+max(data[idx_iter[i],]),3),
           pch = 8, cex = 0.4)
  }
    
}

if (writeToFile) dev.off()
