# compare the numbers of added reactions between different gap-filling solutions

library(scales)
library(plotrix)

par(family = "Arial", oma = c(0,0,0,0))

habitat <- "Soil"
experiment <- "Schlaeppi"

data <- read.table(paste("/stud/wendering/Masterthesis/FIGURES/added_reactions/", habitat, "_",
                         experiment, "_gf_sol_size.txt", sep = ""), header = T)
outFile <- paste("/stud/wendering/Masterthesis/FIGURES/added_reactions/", habitat, "_", experiment,
                 "_gf_sol_size.png", sep = "")
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
    # print(var.test(d_1, d_2)$p.value)
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

png(filename = outFile, units = "cm", width = 25, height = 20, res = 300)
# # create empty plot
# plot(1, xlim = x_limits, ylim = y_limits, axes = F, type = "n", xlab = "", ylab = "number of added reactions")
# 
# # point symbols
# ptypes <- c(17, 1)
# 
# # axes
# x_pos <- seq(from = 0.5, to = 2.5, by = 1)
# axis(side = 1, at = x_pos, labels = c("KBase draft", "no_CarveMe", "all"), lty = 0, font = 2, tick = T)
# axis(side = 2, at = pretty(y_limits))
# 
# # box
# box(which = "plot", lty = "solid")

# ~~~~~~~~~~~~~~ version 1 ~~~~~~~~~~~~~~ #
# # legend
# legend(x = 2.2, y = y_max*0.9, legend = c("individual", "iterative"), pch = ptypes,  bty = "n", cex = 1.3)
# 
# # error bars
# col_err = "gray40"
# x_err <- c(x_pos, x_pos+0.05)
# for (i in 1:length(error)) {
#   y_err_max <- av[i]+error[i]
#   y_err_min <- av[i]-error[i]
#   lines(x = c(x_err[i],x_err[i]), y = c(y_err_min, y_err_max), col = col_err)
# 
#   x_err_max <- x_err[i]+0.01
#   x_err_min <- x_err[i]-0.01
#   lines(x = c(x_err_min, x_err_max), y = c(y_err_min, y_err_min), col = col_err)
#   lines(x = c(x_err_min, x_err_max), y = c(y_err_max, y_err_max), col = col_err)
# }
# 
# # plot points
# # individual
# points(x = x_pos+0.05, y = av[idx_ind], pch = ptypes[1])
# # iterative
# points(x = x_pos, y = av[idx_iter], pch = ptypes[2])


# ~~~~~~~~~~~~~~ version 2 ~~~~~~~~~~~~~~ #
# # median
# med <- apply(X = data, MARGIN = 1, FUN = median)
# 
# # plotting parameters
# col_points <- c("lightblue", "firebrick")
# col_lines <- c("blue", "firebrick4")
# transp <- 0.5
# lwd_lines <- 4

# # legend
# legend(x = 2.2, y = y_max*0.9, legend = c("individual", "iterative"), pch = c(19,19),  bty = "n",
#        cex = 1.3, col = col_points)
# 
# for (i in 1:length(idx_ind)) {
#   points(x = rep(x_pos[i], ncol(data))+runif(ncol(data), min = 0, max = 0.05), y = data[idx_ind[i],],
#          pch = 19, col = col_points[1])
#   lines(x = c(x_pos[i]-0.05, x_pos[i]+0.05), y = rep(med[idx_ind[i]],2),
#         col = alpha(col_lines[1], transp), lwd = lwd_lines)
# }
# 
# for (i in 1:length(idx_iter)) {
#   points(x = rep(x_pos[i], ncol(data))+runif(ncol(data), min = 0, max = 0.05), y = data[idx_iter[i],],
#          pch = 19, col = col_points[2])
#   lines(x = c(x_pos[i]-0.05, x_pos[i]+0.05), y = rep(med[idx_iter[i]],2),
#         col = alpha(col_lines[2], transp), lwd = lwd_lines)
# }

# ~~~~~~~~~~~~~~ version 3 ~~~~~~~~~~~~~~ #
# boxplot(t(data)[,c(4,1,5,2,6,3)])

# ~~~~~~~~~~~~~~ version 4 ~~~~~~~~~~~~~~ #
par(family = "Arial", mgp = c(2.8,1,0))
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
violin_colors <- alpha(rep(cbp1[c(2,1)], 3), 0.8)
cex = 1.6
x_pos = c(1,2,3.5,4.5,6,7)
violin_plot(X = t(data)[,c(4,1,5,2,6,3)], at = x_pos, axes = F, main = "",
            col = violin_colors, ylim = y_limits)

# plot axes
x_pos <- c(1.5, 4, 6.5)
axis(side = 1, at = x_pos, labels = c("KBase draft", "consensus - CarveMe", "consensus"),
     lty = 0, font = 2, tick = T, cex.axis = cex)
axis(side = 2, at = pretty(y_limits), cex.axis = cex, col = "gray30")
title(ylab = "Number of added reactions", cex.lab = cex)
# legend
legend(x = x_pos[2]*1.3, y = y_max*0.7, legend = c("individual", "conditional"),  bty = "n",
       cex = 1.2*cex, col = violin_colors[1,3,5], fill = violin_colors[c(1,2)], text.col = "gray30")

# box
# box(which = "plot", lty = "solid")

# draw significance stars
x_pos <- c(2, 4.5, 7)
for (i in 1:length(h)){
  if (h[i]==1) {
    points(x = x_pos[i], y = 3+max(data[idx_iter[i],]), pch = 8, cex = 0.6)
  } 
  else if (h[i]==2) {
    points(x = c(x_pos[i]-0.05,x_pos[i]+0.05), y = rep(3+max(data[idx_iter[i],]),2), pch = 8, cex = 0.6)
  }
  else if (h[i]==3) {
    points(x = c(x_pos[i]-0.05,x_pos[i],x_pos[i]+0.05), y = rep(3+max(data[idx_iter[i],]),3), pch = 8, cex = 0.6)
  }
    
}

# experiment name
# text(x = x_pos[2], y = 0.9*y_limits[2], labels = experiment, cex = 1.5, font = 2)



dev.off()
