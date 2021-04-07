# generate plots for the comparison of different gap filling solutions
library(scales)
library(pheatmap)
library(wesanderson)

par(family = "Arial")

##########################################################
# intersection of gap-filling solutions
##########################################################

topDir <- "~/ComGapFill"
writeToFile <- T
habitat <- "Soil"
experiments <- c("Schlaeppi", "Bulgarelli")
spec <- c("", "_no_exc_rxns")
 
palettelength <- 100
col1 <- "beige"
col2 <- "firebrick4"
pal <- colorRampPalette(c(col1, col2))(palettelength)
mybreaks = c(seq(from = 0,to = 1, length.out =  palettelength))

for (E in experiments) {
  for (S in spec){
    data <- read.table(paste(topDir, "/figures/added_reactions/", habitat, "_",
                             E, "_dist_methods", S, ".txt", sep = ""),
                       header = T, na.strings = "NaN")
    row.names(data) <- data[,1]
    data = data[,-1]
    data[is.na(data)] = 0
    data <- 1 - data
    if (writeToFile) {
      png(paste0(topDir, "/figures/added_reactions/jaccard_index_solutions_",
              habitat, '_', E, S, '.png'), height = 20, width = 20, units = "cm", res = 600)
    }
    labels <- c("KBase draft models cond.", "consensus - C cond.", "consensus cond.",
                "KBase draft models ind.", "consensus - C ind.", "consensus ind.")
     c <- pheatmap(as.matrix(data),
             cluster_cols = T,
             cluster_rows = F,
             labels_col = labels,
             kmeans_k = 3,
             main = "",
             breaks = mybreaks,
             fontsize = 18,
             fontsize_col = 16,
             family = "Arial",
             border_color = NA,
             color = pal,
             angle_col = 315, 
             gaps_row = c(1,2),
             clustering_method = "average",
             show_rownames = F,
             width = 12,
             height = 8
             )
    
     if (writeToFile) dev.off()
  print(paste(topDir, "/figures/added_reactions/", habitat, "_",
             E, "_dist_methods", S, ".txt", sep = ""))
  }
}

##########################################################
# compare exchanged metabolites and added reactions per model
##########################################################

# habitat = "Soil"
# modelType <- "KBase"
# 
# # ~~~~~~~~~~ reactions ~~~~~~~~~~ #
# 
# fileName <- paste("/stud/wendering/Masterthesis/FIGURES/added_reactions/", modelType, "_",
#                   habitat,"_comp_gf_exp.txt", sep = "")
# 
# 
# rxns <- read.table(fileName, header = T)
# row.names(rxns) = rxns[,1]
# rxns = rxns[,-1]
# 
# # find models that have a low intersection compared to the numbers of reactions specific to each experiment
# diff_rxns <- apply(X = rxns, MARGIN = 2, FUN = function(x) {x[3]-mean(x[c(1,2)])})
# JI_rxns <- apply(X = rxns, MARGIN = 2, FUN = function(x) {x[3] / sum(x)})
# 
# # ~~~~~~~~~~ imported metabolites ~~~~~~~~~~ #
# 
# fileName <- paste("/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/", modelType, "_",
#                   habitat,"_imported.txt", sep = "")
# 
# 
# imp <- read.table(fileName, header = T)
# row.names(imp) = imp[,1]
# imp = imp[,-1]
# 
# # find models that have a low intersection compared to the numbers of reactions specific to each experiment
# diff_imp <- apply(X = imp, MARGIN = 2, FUN = function(x) {x[3]-mean(x[c(1,2)])})
# JI_imp <- apply(X = imp, MARGIN = 2, FUN = function(x) {x[3] / sum(x)})
# # ~~~~~~~~~~ exported metabolites ~~~~~~~~~~ #
# 
# fileName <- paste("/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/", modelType, "_",
#                   habitat,"_exported.txt", sep = "")
# 
# 
# exp <- read.table(fileName, header = T)
# row.names(exp) = exp[,1]
# exp = exp[,-1]
# 
# # find models that have a low intersection compared to the numbers of reactions specific to each experiment
# diff_exp <- apply(X = exp, MARGIN = 2, FUN = function(x) {x[3]-mean(x[c(1,2)])})
# JI_exp <- apply(X = exp, MARGIN = 2, FUN = function(x) {x[3] / sum(x)})
# # ~~~~~~~~~~ plotting ~~~~~~~~~~ #
# # plot(x = c(rep(1, length(diff_rxns)), rep(2, length(diff_imp)), rep(3, length(diff_exp))),
# #      y = c(diff_rxns, diff_imp, diff_exp),
# #      axes = 1, bty = "n", xlab = "", ylab = "abs_diff")
# plot(x = c(rep(1, length(JI_rxns)), rep(2, length(JI_imp)), rep(3, length(JI_exp))),
#      y = c(JI_rxns, JI_imp, JI_exp),
#      axes = 1, bty = "n", xlab = "", ylab = "Jaccard Index")
# 
# plot(JI_rxns, JI_imp)
# fit <- line(JI_rxns, JI_imp)
# abline(fit)
# 
# Rsq = cor(JI_rxns, JI_imp)^2
# legend("topleft", legend = paste("R^2",deparse(round(Rsq, 2)), sep = " = "))
