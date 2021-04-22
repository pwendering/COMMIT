library(scales)

writeToFile = F
topDir = "~/ComGapFill/"
# topDir = "C://Users/wende/MobaXterm/home/comgapfill/smetana-analysis/results/"
habitat = "Soil"
experiment = "Schlaeppi"

# read results of parsed smetana results
smetanaFile = paste0(topDir, "smetana-analysis/results/exchange_", habitat, "_", experiment, ".txt")
smetanaDat = read.table(smetanaFile,
           header = T, row.names = 1,  sep = "\t")
otus = rownames(smetanaDat)
# read BRITE dictionary file
brDictFile = paste0(topDir, "smetana-analysis/results/brite_dict_", habitat, "_", experiment, ".txt")
briteDict = read.table(brDictFile, header = T, sep = "\t")
# order by BRITE
brOrder = order(briteDict$BRITE)
briteDict = briteDict[brOrder,]

# color coding of metabolites according to BRITE class
cbp1 <- c("#CC79A7", "#0072B2", "#009E73", "#D55E00",
          "#E69F00", "#56B4E9", "#F0E442", "#999999")
col_brite = cbp1[unlist(lapply(briteDict$BRITE,
                                  function(x) which(unique(briteDict$BRITE)==x)))]

# read taxonomy file
# taxFile = paste0(topDir, "data/genomes/At-SPHERE-genera.txt")
taxFile = paste0(topDir, "data/genomes/At-SPHERE-families.txt")
taxTable = read.table(taxFile,header = T, sep = "\t")
taxClasses = taxTable[taxTable$isolate_ID %in% otus,2]


# get unique list of metabolites
metabolites = briteDict$NAME

if (writeToFile) {
  png(paste0(topDir, "figures/exchanged_metabolites/smetana/smetana_", habitat, "_", experiment, ".png"),
      width = 20, height = 12, units = "cm", res = 600)
}

originalPar = par()
par(usr = c(0,1,0,1), family = "sans", mar = c(0,0,0,0))

y_min = 0.12
y_max = 0.88

space = .18
width = (1 - 2*space)/3

cex_otus = .7
cex_mets = .8
cex_legend = .7
x_pos = matrix(c(0, width,
                 width+space, 2*width+space,
                 2*width+2*space, 3*width+2*space
),nrow = 3, byrow = T)

plot.new()

rect(x_pos[2,1],y_min-.05,x_pos[2,2],y_max+.05,col = "grey90",border = NA)

y_pos_families = seq(y_min,y_max,(y_max-y_min)/(nrow(smetanaDat)-1))
y_pos_mets = seq(y_min,y_max,(y_max-y_min)/(length(metabolites)-1))

text(x = rep(x_pos[1,1],nrow(smetanaDat)), y = y_pos_families,
     paste0(otus, " (", taxClasses, ")"), pos = 4, cex = cex_otus, font= 2)
# text(x = rep(x_pos[1,1],nrow(smetanaDat)), y = y_pos_families,
#      rownames(smetanaDat), pos = 4, cex = cex_otus, font= 2)
text(x = rep(.5,length(metabolites)), y = y_pos_mets, labels = metabolites,
     adj = .5, cex = cex_mets, font = 2, col = col_brite)
text(x = rep(x_pos[3,1],nrow(smetanaDat)), y = y_pos_families,
     otus, pos = 4, cex = cex_otus, font = 2)

t = 0

for (i in 1:nrow(smetanaDat)) {
  
  current_otu = otus[i]
  
  idx = which(otus==current_otu)
  exported_ids = unique(unlist(strsplit(paste(smetanaDat$export_NAME[idx], sep = ","),",")))
  

  if (length(exported_ids)) {
    line_pos = unlist(lapply(exported_ids,function(x) {which(metabolites==x)}))
    tmpBrCol = unlist(lapply(exported_ids,function(x) {col_brite[briteDict$NAME %in% x]}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[1,2]+.01,
        y_pos_families[i],
        x_pos[2,1]-.01,
        y_pos_mets[line_pos[j]],
        col = alpha(tmpBrCol[j],.5),
        lwd = 3)
    }
  }
  imported_ids = unique(unlist(strsplit(paste(smetanaDat$import_NAME[idx], sep = ","),",")))
  
  if (length(imported_ids)) {
    line_pos = unlist(lapply(imported_ids,function(x) {which(metabolites==x)}))
    tmpBrCol = unlist(lapply(imported_ids,function(x) {col_brite[briteDict$NAME %in% x]}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[2,2]+.01,
        y_pos_mets[line_pos[j]],
        x_pos[3,1]-.01,
        y_pos_families[i],
        col = alpha(tmpBrCol[j],.5),
        lwd = 3)
    }
  }
  
}

legend("top", legend = unique(briteDict$BRITE), ncol = length(unique(briteDict$BRITE)),
       cex = cex_legend, fill = cbp1[1:length(unique(briteDict$BRITE))], bty = "n")

par(originalPar)

if (writeToFile) dev.off()

