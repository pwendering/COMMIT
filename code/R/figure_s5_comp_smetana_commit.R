# plot interaction overlap

library(pheatmap)
library(grid)

writeToFile = T

cbp1 <- c("#CC79A7", "#0072B2", "#009E73", "#D55E00",
          "#E69F00", "#56B4E9", "#F0E442", "#999999")

recov_mat = read.table("smetana-analysis/results/Soil_Schlaeppi_interaction_overlap.txt",
                       header = T, row.names = 1, sep = "\t")
recov_mat = as.matrix(recov_mat)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")

if (writeToFile) {
  png(filename = "figures/exchanged_metabolites/smetana/smetana_Soil_Schlaeppi_recovered.png",
      width = 20, height = 15, units = "cm", res = 600)
}

pheatmap(mat = recov_mat,
         cluster_rows = F,
         cluster_cols = F,
         color = cbp1[c(5,2,6,3)],
         legend = T,
         legend_breaks = c(0,1,2,3),
         legend_labels = c("none", "COMMIT", "SMETANA", "both"),
         fontsize = 16
)

setHook("grid.newpage", NULL, "replace")
grid.text("export", y=-0.05,x=.4,gp=gpar(fontsize=18,font=2))
grid.text("import", x=-0.07, rot=90, gp=gpar(fontsize=18,font=2))

if (writeToFile) dev.off()
