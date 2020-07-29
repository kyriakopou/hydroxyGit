maha4.clustering =  read.table(file = "Desktop/DeepResults/WT/clustering/kerror/maha4.txt",
                               header = T, sep = "\t")


est.eff.levels = read.table(file = "Desktop/DeepResults/WT/T.compact.filt.txt",
                                         header = T, sep = "\t")

est.eff.levels = read.table(file = "Desktop/DeepResults/WT/T_compact_BI_WΤ.csv",
                            header = T, sep = ",")


eff.d0 = data.frame(est.eff.levels[, c(17:19)], day = "d0")
eff.d3 = data.frame(est.eff.levels[, c(20:22)], day = "d3")
eff.d6 = data.frame(est.eff.levels[, c(23:25)], day = "d6")

std.d0 = data.frame(est.eff.levels[, c(27:29)], day = "d0" )
std.d3 = data.frame(est.eff.levels[, c(30:32)], day = "d3")
std.d6 = data.frame(est.eff.levels[, c(33:35)], day = "d6")


colnames(eff.d0)[1:3] = c("maint", "deNovo", "hydroxy")
colnames(eff.d3)[1:3] = c("maint", "deNovo", "hydroxy")
colnames(eff.d6)[1:3] = c("maint", "deNovo", "hydroxy")
colnames(std.d0)[1:3] = c("maint", "deNovo", "hydroxy")
colnames(std.d3)[1:3] = c("maint", "deNovo", "hydroxy")
colnames(std.d6)[1:3] = c("maint", "deNovo", "hydroxy")


p.eff = list()
p.std = list()

library(ggplot2)
library(reshape)
library(gridExtra)
for (i in unique(maha4.clustering$Cluster)){

  rows = (maha4.clustering$Cluster == i)
  cluster.eff = rbind(eff.d0[rows, ], eff.d3[rows, ], eff.d6[rows, ])
  cluster.eff = melt(cluster.eff, id = c("day"))
  
  cluster.std = rbind(std.d0[rows, ], std.d3[rows, ], std.d6[rows, ])
  cluster.std = melt(cluster.std, id = c("day"))
  
  p.eff[[i]] = ggplot(cluster.eff, aes(x=day, y = value, fill = variable)) + 
                geom_boxplot() + scale_fill_manual(values = c("firebrick3", "blue3", "goldenrod1"))
  
  p.std[[i]] = ggplot(cluster.std, aes(x=day, y = value, fill = variable)) + 
    geom_boxplot() + scale_fill_manual(values = c("firebrick3", "blue3", "goldenrod1"))
  
  
}


do.call(grid.arrange, p.eff)
do.call(grid.arrange, p.std)



ossta.gene.info = read.table("~/Desktop/OsttaV2_all.gtf", sep = "\t", fill = T)
