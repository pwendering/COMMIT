# plot characteristics of gap-filling solutions across all explored iterations
library(scales)

myDensityPlot <- function (X, col, nbins, xlim) {

  coord = binData(X, nbins = nbins, xlim = xlim)
  x = coord[[1]]
  y = coord[[2]]
  lines(x = x, y = y, col = alpha(col, 0.5), lwd = 6)
  y = c(0,y,0)
  x = c(x[1],x,x[length(x)])
  polygon(x,y, col = alpha(col,0.1), border = NA)#alpha(col, 0.4), lwd = 6)
}

binData <- function (X, nbins, xlim) {
  # range
  r = max(xlim) - min(xlim)
  # bin size
  s = r / nbins
  counts = rep(0, nbins)
  mids = rep(0, nbins)
  for (i in 1:nbins) {
    min = (i-1)*s+min(xlim)
    max = i*s+min(xlim)
    counts[i] = sum(as.numeric((X>=min)&(X<max)))
    mids[i] = mean(c(min, max))
  }
  res = list(mids, counts)
  return(res)
}

plot_gf_distribution <- function (wd, habitat, study, outFileBase, add_legend,
                                 y_axis, x_axis, writeToFile = F) {
  # number of exchanged metabolites
  exc <- read.table(paste(wd, study, "-exc.txt", sep = ""), header = T)
  exc = apply(exc, 1, sum)
  # number of added reactions
  gf <- read.table(paste(wd, study, "-gf.txt", sep = ""), header = T)
  gf = apply(gf, 1, sum)
  # sum of optimal biomass fluxes over all community members
  bio <- read.table(paste(wd, study, "-bio.txt", sep = ""), header = T)
  bio = apply(bio, 1, sum)
  # row index of the optimal ordering
  opt_idx <- read.table(paste(wd, study, "-opt.txt", sep = ""), header = T)
  opt_idx = opt_idx[1,1]
  
  # scale all values to the optimum determined above
  data = rbind(gf, bio, exc)
  data = apply(data,1,function(x, opt) {x=x-x[opt]; x=x/max(abs(x))}, opt=opt_idx)
  
  if (writeToFile) {
    png(paste(outFileBase, habitat, "-", study, "-", as.character(nrow(data)), ".png", sep = ""),
        units = "cm", height = 20, width = 20, res = 300)
  }
  
  # Plotting and graphical parameters
  cex_axis = 2
  col_axis = "gray30"
  par(family = "sans", mar = c(5, 7.8, 4, 2) + 0.1, mgp = c(3.5,1,0), xpd = T)
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  col <- cbp1[c(4,6,8)]
  
  x_limits <- c(-1, 1)
  n <- 10
    
  # y limits
  y_max <- c()
  for (i in 1:ncol(data)) {
    h = binData(data[,i], nbins = n, xlim = x_limits)
    y_max[i]=max(h[[2]])
  }
  
  y_max = max(y_max)/10
  y_max = ceiling(y_max)*10
  # y_limits <- c(0,y_max)
  y_limits <- c(0,70)
  # empty plot
  plot(1, xlim = x_limits, ylim = y_limits, type = "n", axes = F,
       xlab = "", ylab = "")

  # fill in density plots from histograms
  for (i in 1:ncol(data)) {
    myDensityPlot(data[,i], col = col[i], nbins = n, xlim = x_limits)
  }
  
  # axes and legend
  axis(side = 1, at = pretty(x_limits), col = col_axis, lwd = 1.5, cex.axis = cex_axis)
  if (x_axis) {
    title(xlab = "Scaled difference to optimal solution", cex.lab = cex_axis)
  }
  axis(side = 2, at = pretty(y_limits), col = col_axis, lwd = 1.5,
       cex.axis = cex_axis)
  if (y_axis){
    title(ylab = "Frequency", cex.lab = cex_axis)
  }
  if (add_legend){
    legend(x = -1, y = 95, legend = c("number of added reactions", "sum of biomass reaction fluxes", "number of exchanged metabolites"),
           col = alpha(col, 0.5), lwd = 6, bty = "n", cex = 1.5, bg = NA)
  }
  
  if (writeToFile) dev.off()
}

add_letter <- function (letter = "X", case = "upper") {
  if (is.numeric(letter)) {
    if (case == "upper") {
      letter = LETTERS[letter]
    }
    else
      letter = letters[letter]
  }
  
  tmp_usr = par("usr")
  par(usr = c(0,1,0,1), family = "sans")
  
  text(x = -0.1, y = 1.1, labels = letter, pos = 1, cex = 2.5, font = 2)
  
  par(usr = tmp_usr)
}
add_title <- function (title) {
  tmp_usr = par("usr")
  par(usr = c(0,1,0,1), family = "Arial")
  
  text(x = .5, y = 1.2, labels = title, cex = 2.5, font = 2)
  
  par(usr = tmp_usr)
}
add_annotation <- function (annotation) {
  tmp_usr = par("usr")
  par(usr = c(0,1,0,1), family = "sans")
  
  mtext(text = bquote(~underline(~bolditalic(.(annotation)))), side = 2,
        line = 6, cex = 1.5, adj = .5)
  par(usr = tmp_usr)
}
writeToFile = T
habitat = "Soil"
# topDir = "~/ComGapFill"
topDir = "C://Users/wende/MobaXterm/home/comgapfill/"
# store default graphical parameters
originalPar = par();



if (writeToFile) {
  png(filename = paste0(topDir, "/figures/gap-filling/Figure_3_gap_filling_iterations.png"), width = 20, height = 25, units = "cm",
      res = 600)
}
par(mfrow=c(3,2),xpd =T)

# Full consensus models
wd = paste(topDir, "/data/gap-filling/iterative/", habitat, "/all/", sep = "")
outFileBase <- paste(topDir, "/figures/gap-filling/all/", sep = "")
plot_gf_distribution(wd, habitat, "Schlaeppi", outFileBase, add_legend = F, y_axis = T, x_axis = F, writeToFile = F)
add_letter(1)
add_title("Schlaeppi")
add_annotation("consensus")
plot_gf_distribution(wd, habitat, "Bulgarelli", outFileBase, add_legend = F, y_axis = F, x_axis = F, writeToFile = F)
add_letter(2)
add_title("Bulgarelli")


# Consensus models without CarveMe models
wd = paste(topDir, "/data/gap-filling/iterative/", habitat, "/no_CarveMe/", sep = "")
outFileBase <- paste(topDir, "figures/gap-filling/no_CarveMe/", sep = "")
plot_gf_distribution(wd, habitat, "Schlaeppi", outFileBase, add_legend = T, y_axis = T, x_axis = F, writeToFile = F)
add_letter(3)
add_annotation("-CarveMe")
plot_gf_distribution(wd, habitat, "Bulgarelli", outFileBase, add_legend = F, y_axis = F, x_axis = F, writeToFile = F)
add_letter(4)

# Adapted KBase models
wd = paste(topDir, "/data/gap-filling/iterative/Soil/KBase/", sep = "")
outFileBase <- paste(topDir, "figures/gap-filling/iterative/KBase/kbase-", sep = "")
plot_gf_distribution(wd, habitat, "Schlaeppi", outFileBase, add_legend = F, y_axis = T, x_axis = F, writeToFile = F)
add_letter(5)
add_annotation("KBase")
plot_gf_distribution(wd, habitat, "Bulgarelli", outFileBase, add_legend = F, y_axis = F, x_axis = F, writeToFile = F)
add_letter(6)

# Add common x-axis
mtext("scaled difference to optimal solution", side = 1, outer = T, line = -1.5, cex = 1.2)

if (writeToFile) dev.off()

# reset parameters
par(originalPar)