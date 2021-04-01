library(scales)

writeToFile = F
topDir = "~/ComGapFill/figures"

modelType <- "all"
habitat = "Soil"
experiment = "Schlaeppi"

matrixFile <- paste0(topDir, "/exchanged_metabolites/graph/", habitat, "_", modelType, "_exchanged_metabolites_", experiment, ".txt")
outFile <- paste0(topDir, "/exchanged_metabolites/graph/graph_", experiment, ".png")
classFile <- paste0(topDir, "/exchanged_metabolites/graph/", modelType,"_brite_exchanged_", experiment, ".txt")
excFile <- paste0(topDir, "/exchanged_metabolites/graph/", habitat, "_", modelType, "_exchanged_metabolites_IDs_", experiment, ".txt")


ex_mt <- as.matrix(read.table(file = matrixFile, header = T, sep = "\t"))

ex_classes <- read.table(file = classFile, header = T, row.names = 1, sep = "\t")
classes <- unique(as.vector(as.matrix(ex_classes)))
classes = setdiff(classes, "")

exc_full <- read.table(excFile, header = T, row.names = 1, sep= "\t")

plot.new()

tmp_usr = par("usr")
par(usr = c(0,1,0,1))

y_min = 0.2
y_max = 0.8

space = .2
width = (1 - 2*space)/3

x_pos = matrix(c(0, width,
                 width+space, 2*width+space,
                 2*width+2*space, 3*width+2*space
                 ),nrow = 3, byrow = T)

rect(x_pos[1,1],y_min-.05,x_pos[1,2],y_max+.05,col = "grey70")
rect(x_pos[2,1],y_min-.05,x_pos[2,2],y_max+.05,col = "grey70")
rect(x_pos[3,1],y_min-.05,x_pos[3,2],y_max+.05,col = "grey70")

y_pos_families = seq(y_min,y_max,(y_max-y_min)/ncol(ex_mt))
y_pos_classes = seq(y_min,y_max,(y_max-y_min)/length(classes))

text(x = rep(x_pos[1,1],ncol(ex_mt)), y = y_pos_families, labels = colnames(ex_mt), pos = 4, cex = .8)
text(x = rep(x_pos[2,1],length(classes)), y = y_pos_classes, labels = classes, pos = 4, cex = .8)
text(x = rep(x_pos[3,1],ncol(ex_mt)), y = y_pos_families, labels = colnames(ex_mt), pos = 4, cex = .8)

for (i in 1:ncol(ex_mt)) {
  current_familiy = colnames(ex_mt)[i]
  
  exported_uniq = unique(strsplit(exc))
  
}






par(usr = tmp_usr)


