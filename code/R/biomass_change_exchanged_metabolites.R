# relative biomass changes upon blocking uptake and export reactions
library(pheatmap)
library(wesanderson)

translateIDs2Name <- function(fname, topDir) {
  return(
    system(
    command = paste0("echo \"", topDir, "code/bash/mapIDs.sh ",
                    fname,
                    " 7 ", topDir, "data/tables/MNXref/", "MNXref-met-translation-table.csv\" | bash"),
    intern = T
    )
  )
}
cm2inches <- function(cm) return(cm*0.393701)
topDir <- "~/ComGapFill/"
habitat <- "Soil"
experiment <- "Schlaeppi"
modelType <- "all"
t_ratio <- 1-1E-3


# ~~~~~~~~~~~ General settings ~~~~~~~~~~~ #

palettelength = 100
col1 <- "darkorange4"
col2 <- "goldenrod"
col3 <- "gray50"
pal <- colorRampPalette(c(col1, col2, col3), bias = 0.1)(palettelength)
# pal <- wes_palette("Chevalier1", palettelength, type = "continuous")
mybreaks <- c(seq(from = 0,to = 1, length.out =  palettelength))

fileBaseName <- paste(topDir, "figures/exchanged_metabolites/reduction_biomass/", habitat, "_",
                     experiment, "_", modelType, "_biomass_", sep = "")
translationDir <- paste0(topDir, "data/tables/MNXref")
tmpFile <- paste(translationDir, "tmp", sep = "")

# ~~~~~~~~~~~ medium ~~~~~~~~~~~ #
medium <- translateIDs2Name(fname = paste0("<(cut -f 1 ", topDir, "data/media/minimal-medium.csv)"), topDir = topDir)

# ~~~~~~~~~~~ Import ~~~~~~~~~~~ #

import <- read.table(paste(fileBaseName, "import.txt", sep = ""), header = T, sep = "\t", row.names = 1)
ids_import <- colnames(import)
write.table(x = ids_import, file = tmpFile, row.names = F,
            col.names = F, quote = F)
cnames <- translateIDs2Name(tmpFile, topDir = topDir)

colnames(import) <- cnames
# apply a threshold of 10E-6 to the ratio
import <- apply(X = import, MARGIN = 2, FUN = function(x) {x[x>=t_ratio]=1; return(x)})
# select only columns with an observed difference
ids_import = ids_import[which(colSums(import)<nrow(import))]
import = import[,which(colSums(import)<nrow(import))]

# ~~~~~~~~~~~ Export ~~~~~~~~~~~ #
export <- read.table(paste(fileBaseName, "export.txt", sep = ""), header = T, sep = "\t", row.names = 1)
ids_export <- colnames(export)
write.table(x = ids_export, file = tmpFile, row.names = F,
            col.names = F, quote = F)
cnames <- translateIDs2Name(tmpFile, topDir = topDir)
colnames(export) <- cnames

# apply a threshold of 10E-6 to the ratio
export <- apply(X = export, MARGIN = 2, FUN = function(x) {x[x>=t_ratio]=1; return(x)})
# select only columns with an observed difference
ids_export = ids_export[which(colSums(export)<nrow(export))]
export = export[,which(colSums(export)<nrow(export))]

# ~~~~~~~~~~~ Annotation ~~~~~~~~~~~ #

# Rows
families <- read.csv(paste0(topDir, "data/genomes/At-SPHERE-families.txt"), header = T,
                       sep = "\t")
classes <- read.csv(paste0(topDir, "data/genomes/At-SPHERE-classes.txt"), header = T,
                       sep = "\t")

families = families[families[,1] %in% row.names(import),]
classes = classes[classes[,1] %in% row.names(import),]

unique_families = unique(families[,2])
unique_classes = unique(classes[,2])

p_length_families = length(unique_families)
p_length_classes = length(unique_classes)

p_families = wes_palette("Moonrise2", p_length_families, type = "continuous")
p_classes = wes_palette("Cavalcanti1", p_length_classes, type = "continuous")

ann_colors <- list(Class = c(Bacilli = "navyblue",
                             Actinobacteria = "darkolivegreen",
                             Gammaproteobacteria = "firebrick4"),
                   Family = c(Paenibacillaceae = p_families[1],
                             Bacillaceae = p_families[2],
                             Microbacteriaceae = p_families[3],
                             Mycobacteriaceae = p_families[4],
                             Micrococcaceae = p_families[5],
                             Intrasporangiaceae = p_families[6],
                             Xanthomonadaceae = p_families[7],
                             Nocardioidaceae = p_families[8]))

ann_colors$Family <- p_families
names(ann_colors$Family) <- unique_families

ann_colors$Class <- p_classes
names(ann_colors$Class) <- unique_classes

annotation_r <- data.frame(Class = factor(classes[,2]),
                           #Family = factor(families[,2]),
                           row.names = row.names(import))

# Columns

# mapping of MNXref identifiers to KEGG br08001
brite_import <- read.csv(paste(fileBaseName, "import_brite.txt", sep = ""),
                           header = T, sep = "\t", na.strings = "NA",
                         colClasses = c("character", "character"))
# remove rows that are not in imported
brite_import <- brite_import[which(brite_import[,1] %in% ids_import), ]
# empty/ non-classified compounds are grouped as "Other"
brite_import[which(brite_import[,2]==""),2] <- "Other"
annotation_c_import <- data.frame(BRITE = factor(brite_import[,2]),
                                  row.names = colnames(import))

brite_export <- read.csv(paste(fileBaseName, "export_brite.txt", sep = ""),
                         header = T, sep = "\t", na.strings = "NA",
                         colClasses = c("character", "character"))
# remove rows that are not in exported
brite_export <- brite_export[which(brite_export[,1] %in% ids_export),]

# empty/ non-classified compounds are grouped as "Other"
brite_export[which(brite_export[,2]==""),2] <- "Other"
annotation_c_export <- data.frame(BRITE = factor(brite_export[,2]),
                                  row.names = colnames(export))

# re-arrange the dataframes clustered by annotation

# Import
brite_classes <- sort(unique(annotation_c_import$BRITE))
brite_classes <- setdiff(brite_classes, "Other")
brite_classes <- c(brite_classes, "Other")
import_re <- import[,1]
cnames <- c()
for (i in 1:length(brite_classes)) {
  
  # find matching indices of current class in brite classification
  idx <- annotation_c_import$BRITE %in% brite_classes[i]
  # append these columns to the new data frame 
  import_re <- cbind(import_re, import[, idx])
  # same for column names
  cnames <- c(cnames, colnames(import)[idx])
}
import_re <- import_re[,-1]
colnames(import_re) <- cnames

# Export
brite_classes <- sort(unique(annotation_c_export$BRITE))
brite_classes <- setdiff(brite_classes, "Other")
brite_classes <- c(brite_classes, "Other")
export_re = export[,1]
cnames <- c()
for (i in 1:length(brite_classes)) {
  idx <- annotation_c_export$BRITE %in% brite_classes[i]
  export_re <- cbind(export_re, export[, idx])
  cnames <- c(cnames, colnames(export)[idx])
  
}
export_re <- export_re[,-1]
colnames(export_re) <- cnames

# ~~~~~~~~~~~ Plotting ~~~~~~~~~~~ #

# trim the vairable Names
colnames(import) <- strtrim(colnames(import), 30)
colnames(export) <- strtrim(colnames(export), 30)

# metabolite class colors
p_length_metclass <- max(c(length(unique(annotation_c_import[,1])),
                           length(unique(annotation_c_export[,1]))))

palette_col <- c("deeppink4", "darkgreen", "indianred", "darkseagreen4", "darkslateblue", "cornflowerblue", "firebrick")
ann_colors$BRITE <- c(palette_col[1:length(brite_classes)-1], "white")
names(ann_colors$BRITE) <- c(brite_classes)

# Import
font_size <- 8
c_import<- pheatmap(t(import_re),
                    cluster_cols = T,
                    cluster_rows = F,
                    clustering_method = "average",
                    clustering_distance_cols = "euclidean",
                    breaks = mybreaks,
                    color = pal,
                    annotation_row = annotation_c_import,
                    fontsize_row = font_size,
                    fontsize_col = font_size,
                    annotation_col = annotation_r,
                    annotation_colors = ann_colors,
                    border_color = "gray40",
                    filename = paste(fileBaseName, "import.png", sep = ""),
                    cellwidth = 10,
                    cellheight = 8,
                    height = cm2inches(10),
                    width = cm2inches(18),
                    fontsize = font_size,
                    angle_col = 315
)


# Export
c_export <- pheatmap(t(export_re),
                     cluster_cols = T,
                     cluster_rows = F,
                     clustering_method = "average",
                     clustering_distance_cols = "euclidean",
                     #clustering_distance_rows = "euclidean",
                     breaks = mybreaks,
                     color = pal,
                     annotation_row = annotation_c_export,
                     fontsize_row = font_size,
                     fontsize_col = font_size,
                     annotation_col = annotation_r,
                     annotation_colors = ann_colors,
                     border_color = "gray40",
                     filename = paste(fileBaseName, "export.png", sep = ""),
                     cellwidth = 20,
                     cellheight = 12,
                     fontsize = font_size,
                     angle_col = 315
)
