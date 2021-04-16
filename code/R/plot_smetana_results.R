smetanaDat = read.table("C://Users/wende/MobaXterm/home/comgapfill/smetana-analysis/results/exchange_Soil_Schlaeppi.txt",
           header = T, row.names = 1,  sep = "\t")
metabolites = unique(unlist(strsplit(c(smetanaDat$import_NAME,
                                smetanaDat$export_NAME), ",")))
tmp_usr = par("usr")
par(usr = c(0,1,0,1), family = "sans")

y_min = 0.12
y_max = 0.88

space = .2
width = (1 - 2*space)/3

cex_otus = .5
cex_mets = .4
cex_legend = .5
x_pos = matrix(c(0, width,
                 width+space, 2*width+space,
                 2*width+2*space, 3*width+2*space
),nrow = 3, byrow = T)

plot.new()

rect(x_pos[2,1],y_min-.05,x_pos[2,2],y_max+.05,col = "grey90",border = NA)

y_pos_families = seq(y_min,y_max,(y_max-y_min)/(nrow(smetanaDat)-1))
y_pos_mets = seq(y_min,y_max,(y_max-y_min)/(length(metabolites)-1))

text(x = rep(x_pos[1,1],nrow(smetanaDat)), y = y_pos_families,
     rownames(smetanaDat), pos = 4, cex = cex_otus, font= 2)
text(x = rep(.5,length(metabolites)), y = y_pos_mets, labels = metabolites,
     adj = .5, cex = cex_mets, font = 2)
text(x = rep(x_pos[3,1],nrow(smetanaDat)), y = y_pos_families,
     rownames(smetanaDat), pos = 4, cex = cex_otus, font = 2)

t = 0

for (i in 1:nrow(smetanaDat)) {
  
  current_familiy = rownames(smetanaDat)[i]
  
  idx = which(rownames(smetanaDat)==current_familiy)
  exported_ids = unique(unlist(strsplit(smetanaDat$export_NAME[idx],",")))
  

  if (length(exported_ids)) {
    line_pos = unlist(lapply(exported_ids,function(x) {which(metabolites==x)}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[1,2]+.01,
        y_pos_families[i],
        x_pos[2,1]-.01,
        y_pos_mets[line_pos[j]])
    }
  }
  imported_ids = unique(unlist(strsplit(smetanaDat$import_NAME[idx],",")))
  
  if (length(imported_ids)) {
    line_pos = unlist(lapply(imported_ids,function(x) {which(metabolites==x)}))
    
    for (j in 1:length(line_pos)) {
      segments(
        x_pos[2,2]+.01,
        y_pos_mets[line_pos[j]],
        x_pos[3,1]-.01,
        y_pos_families[i])
    }
  }
  
  
}


par(usr = tmp_usr)
