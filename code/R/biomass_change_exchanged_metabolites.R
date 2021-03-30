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
cnames <- colnames(import)
write.table(x = ids_import, file = tmpFile, row.names = F,
            col.names = F, quote = F)
cnames <- translateIDs2Name(tmpFile, topDir = topDir)
#cnames <- c("Fe(2+)", "beta-D-galactose", "PPi", "L-histidine", "keto-D-fructose",
#            "propanoyl P", "glutathione disulfide", "acetoacetate", "3'-AMP",
#            "alpha,alpha-trehalose", "H(+)", "3,4-dihydroxybenzoate", "3'-GMP", "3'-UMP",
#            "2-hydroxy-3-oxobutyl P", "indole-3-acetaldehyde", "shikimate","hydrogencarbonate",
#            "maltopentaose", "maltotriose", "allantoin", "D-mannitol", "folate", "sulfur", "Pi")
# 
#  cnames <- c("sulfite", "Fe(2+)", "PPi", "L-asparagine", "propanoyl P", "glutathione disulfite",
#              "3'-AMP", "L-valine", "3'-CMP", "indole-3-acetaldehyde", "sulfate", "keto-L-sorbose",
#              "D-xylulose", "iminosuccinate", "shikimate", "hydrogencarbonate", "maltohexaose",
#              "maltopentaose", "maltotriose", "allantoin", "D-glucopyranose 1-P(2-)", "L-ornithine",
#              "sulfur", "Pi")
colnames(import) <- cnames
# apply a threshold of 10E-6 to the ratio
import <- apply(X = import, MARGIN = 2, FUN = function(x) {x[x>=t_ratio]=1; return(x)})
# select only columns with an observed difference
ids_import = ids_import[which(colSums(import)<nrow(import))]
import = import[,which(colSums(import)<nrow(import))]
# colnames(import) <- cnames


# ~~~~~~~~~~~ Export ~~~~~~~~~~~ #
export <- read.table(paste(fileBaseName, "export.txt", sep = ""), header = T, sep = "\t", row.names = 1)
ids_export <- colnames(export)
write.table(x = ids_export, file = tmpFile, row.names = F,
            col.names = F, quote = F)
cnames <- translateIDs2Name(tmpFile, topDir = topDir)

# Schlaeppi
#cnames <- c("CO",
#    "L-proline",
#    "diphosphate",
#    "N-formyl-L-Glu",
#    "thiamine tri-P",
#    "keto-D-Fructose",
#    "CO2",
#    "melibiose",
#    "propanoyl P",
#    "acetoacetate",
#    "uracil",
#    "thiamine(1+) chloride",
#    "propionate (n-C3:0)",
#    "methyl-D-erythritol P",
#    "sucrose",
#    "L-proline amide",
#    "adenine",
#    "phosphoethanolamine",
#    "mannosylfructose 6-P",
#    "alpha,alpha-trehalose",
#    "D-ribose",
#    "(S)-4,5-dihydroxypentane-2,3-dione",
#    "2-deoxy-D-ribose",
#    "heme b",
#    "thiamine diphosphate",
#    "succinate",
#    "acetate",
#    "acetyl phosphate",
#    "oxalate",
#    "glycine",
#    "triphosphate",
#    "alpha-D-galactose 1-P",
#    "D-glucono-1,5-lactone",
#    "propanoate",
#    "L-aspartate 4-semialdehyde",
#    "D-pantetheine 4'-phosphate",
#    "alpha-D-galactose",
#    "formate",
#    "autoinducer-2",
#    "L-aspartate",
#    "5-deoxy-D-ribose",
#    "butanoate",
#    "dihydroxyacetone",
#    "malonate",
#    "(S)-2,3,4,5-tetrahydrodipicolinate",
#    "alpha-maltose 1-P",
#   "sulfate",
#    "L-erythrulose",
#    "L-threitol",
#    "D-xylulose",
#    "iminosuccinate",
#    "hydrogencarbonate",
#    "D-mannitol",
#    "O-phosphorylhomoserine",
#    "beta-D-glucose 1-P",
#    "amylotriose",
#    "carbamate",
#    "cis-aconitate",
#    "O-succinyl-L-homoserine",
#    "gly-asn",
#    "3-hydroxypropanoate",
#    "hydroxymethylpyrimidine",
#    "alpha-D-ribose 1,5-bisphosphate",
#    "H2S",
#    "D-glucopyranose 1-P(2-)",
#    "isocitrate",
#    "beta-maltose",
#    "(S)-methylmalonate semialdehyde",
#    "fumarate",
#    "Thiamine thiazole",
#    "Pi")

colnames(export) <- cnames

# apply a threshold of 10E-6 to the ratio
export <- apply(X = export, MARGIN = 2, FUN = function(x) {x[x>=t_ratio]=1; return(x)})
# select only columns with an observed difference
ids_export = ids_export[which(colSums(export)<nrow(export))]
export = export[,which(colSums(export)<nrow(export))]
# colnames(export) <- cnames
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
brite_export <- brite_export[which(brite_export[,1] %in% ids_export), ]
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
# palette_col <- wes_palette("GrandBudapest1", p_length_metclass, type = "continuous")
# palette_col <- colorRampPalette(c("cornsilk2", "dodgerblue4"))(p_length_metclass)
palette_col <- c("deeppink4", "darkgreen", "indianred", "darkseagreen4", "darkslateblue", "cornflowerblue", "firebrick")
ann_colors$BRITE <- c(palette_col[1:length(brite_classes)-1], "white")
names(ann_colors$BRITE) <- c(brite_classes)
# ann_colors$BRITE <- c("Carbohydrates" = palette_col[1],
#                                "Minerals" = palette_col[2],
#                                "Nucleic acids" = palette_col[3],
#                                "Organic acids" = palette_col[4],
#                                "Peptides" = palette_col [5],
#                                "Vitamins and Cofactors" = palette_col[6],
#                                "Lipids" = palette_col[7],
#                                "Other" = "white")

# Import
font_size <- 14
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
                     # height = 14,
                     # width = 14,
                     cellwidth = 20,
                    cellheight = 12,
                     fontsize = font_size,
                     angle_col = 315
)
# c_import <- pheatmap(import_re,
#                      cluster_cols = F,
#                      cluster_rows = T,
#                      clustering_method = "average",
#                      breaks = mybreaks,
#                      color = pal,
#                      annotation_row = annotation_r,
#                      fontsize_row = 12,
#                      fontsize_col = 10,
#                      annotation_col = annotation_c_import,
#                      annotation_colors = ann_colors,
#                      border_color = NA,
#                      filename = paste(fileBaseName, "import.png", sep = ""),
#                      height = 6,
#                      width = 12,
#                      angle_col = 315
#                      )
# cutree(c_import$tree_row, k = 8)

# Export
# c_export <- pheatmap(export_re,
#                      cluster_cols = F,
#                      cluster_rows = T,
#                      clustering_method = "average",
#                      breaks = mybreaks,
#                      color = pal,
#                      annotation_row = annotation_r,
#                      fontsize_row = 12,
#                      fontsize_col = 6,
#                      annotation_col = annotation_c_export,
#                      annotation_colors = ann_colors,
#                      border_color = NA,
#                      filename = paste(fileBaseName, "export.png", sep = ""),
#                      height = 5,
#                      width = 10,
#                      fontsize = 8,
#                      angle_col = 315
#                      )
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
                     # height = 14,
                     # width = 14,
                     cellwidth = 20,
                     cellheight = 12,
                     fontsize = font_size,
                     angle_col = 315
)
