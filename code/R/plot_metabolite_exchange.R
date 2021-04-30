# plot import to and export from the common pool of exchanged metabolites per bacterial family

library(scales)

writeToFile = T
topDir = "figures/exchanged_metabolites/graph/"

modelType <- "all"
habitat = "Soil"
experiment = "Schlaeppi"

matrixFile <- paste0(topDir, habitat, "_", modelType, "_exchanged_metabolites_", experiment, ".txt")
outFile <- paste0(topDir, "exchange_", experiment, ".png")
# classFile <- paste0(topDir, modelType,"_brite_exchanged_", experiment, ".txt")
excFile <- paste0(topDir, habitat, "_", modelType, "_exchanged_metabolites_IDs_", experiment, ".txt")
dictFile <- paste0(topDir, habitat, "_", modelType, "_exchanged_metabolites_dict_", experiment, ".txt")
# better as separate file

ex_mt <- as.matrix(read.table(file = matrixFile, header = T, sep = "\t"))

brite_dict <- read.table(dictFile, header = F, sep = "\t", quote = "\"")
brite_order = order(brite_dict[,2])
brite_dict = brite_dict[brite_order,]
exc_ids = gsub("\\[.\\]","",brite_dict[,1])
exc_brite = brite_dict[,2]
exc_brite[which(exc_brite=="")] = "Other"
exc_names = brite_dict[,3]

exc_full <- read.table(excFile, header = T, row.names = 1, sep= "\t")
all_exported = unlist(strsplit(as.vector(exc_full$export_ID),","))
all_imported = unlist(strsplit(as.vector(exc_full$import_ID),","))
family_count = unlist(lapply(colnames(ex_mt),
                       function(x)length(which(exc_full$family==x))))

if (writeToFile) {
  png(filename = outFile, width = 20, height = 15,units = "cm",res = 600)
}

plot.new()

cbp1 <- c("#CC79A7", "#0072B2", "#009E73", "#D55E00",
          "#E69F00", "#56B4E9", "#F0E442", "#999999")
brite_colors = cbp1[unlist(lapply(exc_brite,
                                   function(x) which(unique(exc_brite)==x)))]
orignialPar = par()
par(usr = c(0,1,0,1), family = "sans", xpd = T)

y_min = -.1
y_max = 1

space = .2
width = (1 - 2*space)/3

cex_families = .8
cex_mets = .6
cex_legend = .7
x_pos = matrix(c(0, width,
                 width+space, 2*width+space,
                 2*width+2*space, 3*width+2*space
                 ),nrow = 3, byrow = T)

rect(x_pos[2,1],y_min-.05,x_pos[2,2],y_max+.05,col = "grey90",border = NA)

y_pos_families = seq(y_min,y_max,(y_max-y_min)/(ncol(ex_mt)-1))
y_pos_mets = seq(y_min,y_max,(y_max-y_min)/(length(exc_ids)-1))

text(x = rep(x_pos[1,1],ncol(ex_mt)), y = y_pos_families,
     labels = colnames(ex_mt), pos = 4, cex = cex_families, font= 2)
text(x = rep(.5,length(exc_ids)), y = y_pos_mets, labels = exc_names,adj = .5, cex = cex_mets,
     col = brite_colors, font = 2)
text(x = rep(x_pos[3,1],ncol(ex_mt)), y = y_pos_families,
     labels = colnames(ex_mt), pos = 4, cex = cex_families, font = 2)

t = 0

for (i in 1:ncol(ex_mt)) {
  
  current_familiy = colnames(ex_mt)[i]
  
  idx = which(exc_full$family==current_familiy)
  exported_ids = gsub("\\[.\\]","",unlist(strsplit(as.vector(exc_full$export_ID[idx]), ",")))
  exported_ids = exported_ids[exported_ids %in% exc_ids]
  exported_ids_uniq = unique(exported_ids)
  exported_ids_count = unlist(lapply(exported_ids_uniq,function (x) {length(which(exported_ids==x))}))
  exported_ids_count = exported_ids_count / family_count[i]
  
  if (any(exported_ids_count>t)) {
  
    exported_ids_uniq = exported_ids_uniq[which(exported_ids_count>t)]
    
    line_pos = unlist(lapply(exported_ids_uniq,function(x) {which(exc_ids==x)}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[1,2]+.01,
        y_pos_families[i],
        x_pos[2,1]-.01,
        y_pos_mets[line_pos[j]],
        lwd = exported_ids_count[j]/family_count[i])
    }
  }
  
  imported_ids = gsub("\\[.\\]","",unlist(strsplit(as.vector(exc_full$import_ID[idx]), ",")))
  imported_ids = imported_ids[imported_ids %in% exc_ids]
  imported_ids_uniq = unique(imported_ids)
  imported_ids_count = unlist(lapply(imported_ids_uniq,function (x) {length(which(imported_ids==x))}))
  
  imported_ids_count = imported_ids_count / family_count[i]
  
  if (any(imported_ids_count>t)) {
  
    imported_ids_uniq = imported_ids_uniq[which(imported_ids_count>t)]
    
    line_pos = unlist(lapply(imported_ids_uniq,function(x) {which(exc_ids==x)}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[2,2]+.01,
        y_pos_mets[line_pos[j]],
        x_pos[3,1]-.01,
        y_pos_families[i],
        lwd = imported_ids_count[j])
    }
  }

  
}

legend(x=.5,y=y_max+.15,xjust = .5, legend = unique(exc_brite), ncol = length(unique(exc_brite)),
       cex = cex_legend, fill = cbp1[1:length(unique(exc_brite))], bty = "n")

par(orignialPar)

if (writeToFile) dev.off()

export_abundances = sapply(exc_ids, function(x) {length(which(gsub("\\[.*$","",all_exported)==x))})
names(export_abundances) = brite_dict[exc_ids %in% gsub("\\[.*$","",brite_dict[,1]),3]
export_abundances = sort(export_abundances,decreasing = T)

import_abundances = sapply(exc_ids, function(x) {length(which(gsub("\\[.*$","",all_imported)==x))})
names(import_abundances) = brite_dict[exc_ids %in% gsub("\\[.*$","",brite_dict[,1]),3]
import_abundances = sort(import_abundances, decreasing = T)

