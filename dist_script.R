# A. JAFFE CAROLINENSIS
# MARCH '15
# Distance Helper Script

setwd("") #current directory
setwd("./Data")

# load in data
load("data_dists.Rdata")
final2=read.csv("localities_final.csv", header=T)
final2 = final2[,2:5]
final2$type = as.character(final2$type)

thresh = c(75, 100, 125,150, 175,200, 225,250, 275,300)
# see coverage by distance threshold
for (dist in thresh){
  nos = c()
  for (i in 1:nrow(dists)){
    nos = cbind(nos,length(dists[i,][dists[i,] < dist]))}
  print(paste(dist,"km threshold: ", length(nos[nos>0]), "/14 populations covered", sep=" "))
}
  
# find closest genetic population for each morph populatioj
for (i in 1:nrow(dists)){
  loc = colnames(dists)[which(dists[i,] == min(dists[i,]))]
  print(paste(rownames(dists)[i], loc, final2[final2$City == loc,]$type, sep=":"))
  }
# sort distances to find closest populations
# closest population does not always have data for all genes
sort(dists[1,])

# using this we can then manually collect genetic data from genBank
