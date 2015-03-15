# A. JAFFE CAROLINENSIS
# MARCH '15
# Main Analysis Script

setwd("") #current directory
setwd("./Data")

# load packages
library(ggplot2)
library(ggsubplot)
library(plyr)
library(mapproj)
library(maps)
library(RColorBrewer)
library(reshape)
library(stringr)
library(ICC)
library(MODISTools)
library(MCMCglmm)
library(gridExtra)
library(RColorBrewer)
library(plotrix)

#colors for multiclass plots
# from http://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25=c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black",
  "gold1","skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F",
  "gray70","khaki2","maroon","orchid1","deeppink1","blue1",
  "steelblue4","darkturquoise","green1","yellow4","yellow3","darkorange4", "brown")

# read in raw data
morph <- read.csv("Acar_morph_data_M_final.csv", header = TRUE)
# remove shoulder width and front parietal measurement
morph <- cbind(morph[,1:13], morph[,15:21], morph[,23:length(morph)])
# change SVL to mm to match other data
morph$Snout.Vent.Length = morph$Snout.Vent.Length*10
# remove miami individs
morph = morph[morph$Site != "Miami",]

# extract mean and sd for each trait
means = c()
sds = c()
for (i in 4:30){
  means = cbind(means, mean(morph[,i]))
  sds = cbind(sds, sd(morph[,i]))}
as.character(sds)
as.character(means)

# log transform data
logmorph = log(morph[,4:30])
logmorph = cbind(morph[,1:3], logmorph)
names(logmorph)[1] = "Observation.Number"

## SECTION 3.1 - GEOGRAPHIC VARIATION IN MORPHOLOGY

# Figure 4ab - snout vent length and exemplar traits
line = coef(lm(Femur ~ Snout.Vent.Length, data = logmorph))
f = ggplot(data=logmorph,aes(x = Snout.Vent.Length, y=Femur, color=Site)) + 
  geom_point(size=4) + geom_abline(intercept = line[1], slope=line[2], color='black') + 
  scale_color_manual(values=c25) + theme(axis.text = element_text(size=13), 
  axis.title = element_text(size=14), legend.text=element_text(size=12), 
  legend.title = element_text(size=12)) + xlab("Snout Vent Length (ln mm)") + 
  ylab("Femur Length (ln mm)")
line2 = coef(lm(Braincase.Width ~ Snout.Vent.Length, data = logmorph))
b = ggplot(data=logmorph,aes(x = Snout.Vent.Length, y=Braincase.Width, color=Site)) + 
  geom_point(size=4) + geom_abline(intercept = line2[1], slope=line2[2], color='black') + 
  scale_color_manual(values=c25) + theme(axis.text = element_text(size=13), 
  axis.title = element_text(size=14), legend.text=element_text(size=12), 
  legend.title = element_text(size=12)) + xlab("Snout Vent Length (ln mm)") + 
  ylab("Braincase Width (ln mm)")
grid.arrange(f,b, ncol=1)

# calculate residuals
resids = logmorph
for (i in 4:29){
  resids[i] = resid(lm(resids[,i] ~ resids$Snout.Vent.Length))
}
#write.csv(resids, "Acar_log_resids_XXX.csv")

#Create full morphological distance matrix from raw data
male_morph=resids
names(male_morph)[2] = "population"
# add in mass data
mass = read.csv("mass_data.csv")
# trim mass data
mass = mass[1:201,]
# change ID to match morph
mass$Observation.Number = as.factor(paste("A",mass$ID,sep=''))
# subselect data
mass = mass[, 5:6]
male_morph = join(male_morph, mass)

# FIGURE 5 - plot residuals by trait
res = cbind(male_morph$population, male_morph[,4:29])
resid.melt= melt(res, id=c("male_morph$population"))
names(resid.melt)[1] = "Population"
ggplot(data=resid.melt, aes(x=Population, y=value)) + geom_boxplot() + geom_point(size=1) + 
  theme(legend.position = "none") + facet_wrap(~ variable, scale="free_y", ncol=7) + 
  theme(axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust=0.4), 
  strip.text = element_text(size=11), axis.title = element_text(size=13)) + 
  ylab("Residual Trait Value")

## run principal componenent analysis on full data
male_morph_pc <- prcomp(male_morph[,4:29], scale=T)
summary(male_morph_pc)
# examine loadings
full.pc = male_morph_pc$rotation
sort(full.pc[,1])
as.character(round(full.pc[,1], 3))
# aggregate by population
pop_pca_means<-aggregate(male_morph_pc$x ~ male_morph$population, data=male_morph_pc$x, mean)
# use all pc's in generation of distance matrix
male_full_morph_matrix <- dist(pop_pca_means[,2:length(pop_pca_means)], method = "euclidean", diag = TRUE, upper = TRUE)

## variation regressions
#generate pc matrix
male.pc = male_morph_pc$x
male.pc = data.frame(male.pc[,1:6])
male.pc = cbind(male.pc, male_morph$population)
names(male.pc)[length(male.pc)] = "Population"
#anova
summary(aov(male.pc$PC1 ~ male.pc$Population))
summary(aov(male.pc$PC2 ~ male.pc$Population))
summary(aov(male.pc$PC3 ~ male.pc$Population))

# hierarchical clustering
labs = pop_pca_means[,1]
dataset = male_full_morph_matrix
morph.clust = hclust(dataset, method = "average")
## FIGURE 7
plot(morph.clust,labels=as.character(labs))

# allen's rule PC analysis
# pelvis width not incuded
male_limb_pc <- prcomp(cbind(male_morph[,4:8], male_morph[,10:13]), scale=T)
limb.pc = male_limb_pc$rotation
sort(limb.pc[,1])
# average by population
pop_pca_limb_means<-aggregate(male_limb_pc$x ~ male_morph$population, data=male_limb_pc$x, mean)

# calculate average SVL
male_svl_avg = aggregate(male_morph$Snout.Vent.Length ~ male_morph$population, data = male_morph, mean)
names(male_svl_avg) = c("City", "SVL")

# bergmann's rule PC analysis
# remove na's
male_morph_nan = na.omit(male_morph)
male_size_pc <- prcomp(cbind(male_morph_nan[,30],male_morph_nan[,31]), scale=T)
size.pc = male_size_pc$rotation
sort(size.pc[,1])
# average by population
pop_pca_size_means<-aggregate(male_size_pc$x ~ male_morph_nan$population, data=male_size_pc$x, mean)

## FIGURE 6
## generate box plots for PC by population
pc.full = data.frame(cbind(male_morph_pc$x[,1:2],
  male_limb_pc$x[,1], male_size_pc$x[,1]))
names(pc.full)[1:length(pc.full)] = c("Full PC1", "Full PC2", "Limb PC1", "Size PC1")
pc.full = cbind(pc.full, male_morph$population)
names(pc.full)[length(pc.full)] = "Population"
pc.full.melt = melt(pc.full, id=c("Population"))
box.pc = ggplot(data=pc.full.melt, aes(x=Population, y=value, fill=Population)) + geom_boxplot(alpha=0.65) + 
  geom_point() + facet_wrap(~ variable, scale="free_y")  + theme(legend.position = "none") + 
  scale_fill_manual(values=c25)
box.pc + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.4), 
               strip.text=element_text(size=15), axis.title = element_text(size=15),
               axis.text = element_text(size=14)) + ylab("Principal Component Score")

# Figure 8ab - plot PC space
# put data in correct format
pcs = cbind(as.character(male_morph$population),male_morph_pc$x[,1:3])
pcs.df = data.frame(pcs)
names(pcs.df) = c("site", "pc1", "pc2", "pc3")
pcs.df$pc1 = as.numeric(as.character(pcs.df$pc1))
pcs.df$pc2 = as.numeric(as.character(pcs.df$pc2))
pcs.df$pc3 = as.numeric(as.character(pcs.df$pc3))
# aggregate to find mean and standard error
ag = aggregate(. ~ site,data = pcs.df,FUN=function(x) c(mn =mean(x), n=std.error(x) ) )
# now plot
pc12 = ggplot(data=ag, aes(x=pc1[,1], y=pc2[,1], colour=site)) +
  scale_color_manual(values=c25) + 
  geom_errorbar(aes(ymax = pc2[,1] + pc2[,2], ymin=pc2[,1] - pc2[,2]), color="black", alpha=.5, width=.15) +
  geom_errorbarh(aes(xmax = pc1[,1] + pc1[,2], xmin=pc1[,1] - pc1[,2]), color="black", alpha=.5, height=0.15) +
  geom_point(alpha=0.5, size=20) + geom_text(aes(label = ag$site), color="black", size=5) + 
  xlab("Full Morph PC1 Score") + ylab("Full Morph PC2 Score") +
  theme(legend.position = "none", axis.text = element_text(size=15), axis.title =element_text(size=15))
pc13 = ggplot(data=ag, aes(x=pc1[,1], y=pc3[,1], colour=site)) +
  scale_color_manual(values=c25) + 
  geom_errorbar(aes(ymax = pc3[,1] + pc3[,2], ymin=pc3[,1] - pc3[,2]), color="black", alpha=.5, width=.15) +
  geom_errorbarh(aes(xmax = pc1[,1] + pc1[,2], xmin=pc1[,1] - pc1[,2]), color="black", alpha=.5, height=0.15) +
  geom_point(alpha=0.5, size=20) + geom_text(aes(label = ag$site), color="black", size=5) + 
  xlab("Full Morph PC1 Score") + ylab("Full Morph PC3 Score") +
  theme(legend.position = "none", axis.text = element_text(size=15), axis.title = element_text(size=15))
grid.arrange(pc12, pc13)

## SECTION 3.2 - ENVIRONMENTAL CORRELATION

#create temperature distance matrix from raw data
male_temp <- read.csv("climate_data.csv")
male_temp = male_temp[order(male_temp$City),]
male_temp$Temp.Seasonality = male_temp$Temp.Seasonality/100
# remove Miami anoles
male_temp  = male_temp[male_temp$City != "Miami",]

# run temperature PC analysis
male_temp_pc <- prcomp(male_temp[,5:15], scale=T)
temp_pca_means<-aggregate(male_temp_pc$x ~ male_temp$City, data=male_temp_pc$x, mean)

#run temperature PC analysis
male_precip_pc <- prcomp(male_temp[,16:23], scale=T)
precip_pca_means<-aggregate(male_precip_pc$x ~ male_temp$City, data=male_precip_pc$x, mean)

# retrieve vegetation data
modis.subset = male_temp[,2:4]
# re-pull NDVI data from MODIS, or use provided vegetation_data.csv below
#modis.subset$start.date <- rep(2003, nrow(modis.subset))
#modis.subset$end.date <- rep(2013, nrow(modis.subset))
#GetDates(Product = "MOD13Q1", Lat = modis.subset$Latitude[1], Long = modis.subset$Longitude[1])
names(modis.subset)[2:3] = c("long","lat")
#MODISSubsets(LoadDat = modis.subset, Product = "MOD13Q1", Bands = c("250m_16_days_NDVI"),Size = c(0,0), StartDate=TRUE)

# vegetation data, parse to combine with geo data  
veggie = read.csv("vegetation_data.csv")
names(veggie)[1:2] = c("Latitude","Longitude")
veggie$Longitude = as.numeric(str_sub(veggie$Longitude, start=1, end = 8))
veggie$Latitude = as.numeric(str_sub(veggie$Latitude, start=1, end = 7))
modis.subset$Latitude = as.numeric(str_sub(modis.subset$lat, start=1, end = 7))
modis.subset$Longitude = as.numeric(str_sub(modis.subset$long, start=1, end = 8))
veggie.combine = join(veggie, modis.subset)
#subselect populations to match other data
pops = male_temp$City
veggie.combine = veggie.combine[veggie.combine$City %in% pops,]
veg_subset = data.frame(cbind(veggie.combine$NDVI, as.character(veggie.combine$City)))
names(veg_subset) = c("NDVI", "City")
veg_subset$NDVI = as.numeric(as.character(veg_subset$NDVI)
                             
## create a final dataframe for regressions
male.master = cbind(male_temp[,2:4], male_svl_avg[,2], pop_pca_means[,2:4],
   pop_pca_size_means[,2], pop_pca_limb_means[,2:3], temp_pca_means[,2:3],
   precip_pca_means[,2:3])
names(male.master) = c("City", "Longitude","Latitude","SVL",
  rep(paste("Full", rep(1:3), sep='')),"Size1",rep(paste("Limb", rep(1:2), sep='')), 
  rep(paste("Temp", rep(1:2), sep='')), rep(paste("Precip", rep(1:2), sep='')))
                             
# add annual temp
ann = cbind(as.character(male_temp$City), male_temp$Annual.Mean.Temp)
ann = data.frame(ann)
names(ann) = c("City", "AnnTemp")
ann$AnnTemp = as.numeric(as.character(ann$AnnTemp))
male.master = join(male.master, ann)
# add unPCAd veg
male.master = join(male.master, veg_subset)
# now scale non pca variables
male.master$AnnTemp = scale(male.master$AnnTemp, center=T, scale=T)
male.master = join(male.master, veg_subset)
male.master$NDVI = as.numeric(as.character(male.master$NDVI))
male.master$NDVI = scale(male.master$NDVI, center=T, scale=T)

# read in mtdna and ndna genetic data and format
mtdna = read.table("ND2_distances.txt", header=F)
ndna = read.table("concat_distances.txt", header=F)
rownames(mtdna) = mtdna$V1
rownames(ndna) = ndna$V1
mtdna = mtdna[,2:length(mtdna)]
ndna = ndna[,2:length(ndna)]
colnames(mtdna) = rownames(mtdna)
colnames(ndna) = rownames(ndna)
mtdna.ord = mtdna[order(rownames(mtdna))]
ndna.ord = ndna[order(rownames(ndna))]
mtdna.t = data.frame(t(mtdna.ord))
ndna.t = data.frame(t(ndna.ord))
mtdna.final = mtdna.t[order(rownames(mtdna))]
ndna.final = ndna.t[order(rownames(ndna))]

# decompose genetic matrix
Msvd = svd(mtdna.final)
#if running ndna, subset male.master nuclear data
#pops = c("Austin", "Brownsville", "Sinton")
#male.master = male.master[!male.master$City %in% pops,]
#Msvd = svd(ndna.final)
# format and append decomposed matrix
Msvd = Msvd$v%*%(t(Msvd$u)*sqrt(Msvd$d))
male.master$Msvd = Msvd

# run regressions
# use prior if necessary
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
# change variables here as necessary
M1 = MCMCglmm(Full2~Temp1, random = ~idv(Msvd), data = male.master, verbose=F)
#M1 = MCMCglmm(Limb1~AnnTemp, random = ~idv(Msvd), data = male.master, verbose=F, prior=prior)

# FIGURE 9 - plot the regression from above
coefs = summary(M1)$solutions
reg = ggplot(data=male.master, aes(x=Temp1, Full2)) + geom_point(size=3) + 
  geom_abline(intercept=coefs[1], slope=coefs[2]) + 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=15)) +
  ylab("Full Morph PC2") + xlab("Temperature PC1")
# run above two sections repeatedly for variables of choice
#grid.arrange(reg, ...)


## INTRODUCTION, DISCUSSION, SUPPLEMENTARY PLOTS

# set base map for geographic plots
all_states <- map_data("state")
p = ggplot()
p <- p + geom_polygon(data=all_states, aes(x=long, y=lat, group = group),colour="grey60", fill="grey80")
p = p + theme(panel.background = element_rect(fill="white"), panel.grid.major = 
  element_line(color="black", size=.1)) + 
  coord_map(xlim=c(-102,-75), ylim=c(25,37)) + xlab("") + ylab("")

#FIGURE 2 - sample distribution/density
# get population counts
n = data.frame(summary(male_morph$population))
n$City = rownames(n)
names(n) = c("n", "City")
male.master.n = join(male.master, n)
sample_density = p + geom_point(data=male.master.n, aes(x=Longitude, y=Latitude), fill = "lightblue", color="black", size=13, pch=21)
sample_density = sample_density + geom_text(data=male.master.n, aes(x = Longitude, y=(Latitude-.6), label = male.master.n$City), size=4) + theme(legend.position = "none")
sample_density + geom_text(data=male.master.n,aes(x = Longitude, y=(Latitude), label = male.master.n$n))

## FIGURE 11 - PC values with geography
# subselect variables we want to plot
trim = data.frame(cbind(as.character(male.master[,1]), male.master[,6], male.master[,9], male.master[,11]))
names(trim) = c("City", "Full2", "Limb1", "Temp1")
# get geo information
male.geo = male.master[,1:3]
# change data format
trim.melt = melt(trim, id="City")
trim.geo = join(trim.melt, male.geo)
trim.geo$value = as.numeric(as.character(trim.geo$value))
# plot
var = ggplot(trim.geo)+ geom_polygon(data=all_states, 
  aes(x=long, y=lat, group = group),colour="grey60", fill="grey80") + 
  theme(panel.background = element_rect(fill="white"), 
  panel.grid.major = element_line(color="black", size=.1)) + 
  coord_map(xlim=c(-102,-75), ylim=c(25,37)) + xlab("") + ylab("") + 
  geom_subplot(height=rel(1.5), ref=ref_box(fill="white", color="black", alpha=0.5), 
  aes(Longitude, Latitude, group=City, subplot = geom_bar(aes(x = variable, y = value, 
  fill=variable), color = "black", stat="identity"))) + 
  geom_text(aes(x = Longitude, y=(Latitude+.83), label = trim.geo$City), size=3)

# FIGURE 12 = bergmann's rule plot
# regenerate state locality for some plots
short = data.frame(Locality = unique(morph$Locality), state = c("MS", "TN", "NC", "GA", "FL", "LA", "TX", "OK"))
sites = data.frame(unique(morph[c("Site", "Locality")]))
final = join(sites, short)
names(final)[1] = "City"
# join dataframes
berg = join(male.master, final)
# divide populations into floridian and non-floridian
berg$inFlorida = (berg$state == "FL")
lat = ggplot(data=berg, aes(Latitude, Size1)) +geom_point(size=6, color="black") + 
  geom_point(size=5, aes(color=berg$inFlorida)) + ylab("Size PC1")
temp = ggplot(data=berg, aes(Temp1, Size1)) +geom_point(size=6, color="black") + 
  geom_point(size=5, aes(color=berg$inFlorida)) + theme(legend.position="none") +
  xlab("Temperature PC1") + ylab("Size PC1")
grid.arrange(temp, lat, ncol=2)

# FIGURE 13 - allen's rule and interaction plots
# calculate coefficients for these two regressions with MCMCglmm
coefs_full = summary(MCMCglmm(Full2~Latitude, random = ~idv(Msvd), 
  data = male.master, verbose=F))$solutions
coefs_limb = summary(MCMCglmm(Limb1~Latitude, random = ~idv(Msvd), 
  data = male.master, verbose=F))$solutions
# now plot
full = ggplot(data=male.master, aes(x=Latitude, y=Full2)) + 
  geom_point(data = male.master, aes(color=Temp1), size=15, alpha=0.6) + 
  scale_color_gradient(low = "blue", high = "red") + 
  geom_text(aes(x=Latitude, y=Full2,label=City), size=4) + theme_minimal() + 
  geom_abline(intercept = coefs_full[1], slope=coefs_full[2], color='black') +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=17), 
  legend.title=element_text(size=15), legend.text=element_text(size=15)) +
  ylab("Full Morphology PC 2")
limb = ggplot(data=male.master, aes(x=Latitude, y=Limb1)) + 
  geom_point(data = male.master, aes(color=Temp1), size=15, alpha=0.6) + 
  scale_color_gradient(low = "blue", high = "red") + 
  geom_text(aes(x=Latitude, y=Limb1,label=City), size=4) + theme_minimal() + 
  geom_abline(intercept = coefs_limb[1], slope=coefs_limb[2], color='black') +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=17), 
  legend.title=element_text(size=15), legend.text=element_text(size=15))
grid.arrange(full, limb, ncol=1)

# SUPP FIGURE 1 - genetic data localities
final2=read.csv("localities_final.csv", header=T)
final2 = final2[,2:5]
final2$type = as.character(final2$type)
# insert full paper names
final2$type[final2$type == "SCS 2012"] <- "Campbell-Staton et al. 2012"
final2$type[final2$type == "Tollis 2014"] <- "Tollis & Boissinot 2014"
final2$type[final2$type == "Tollis 2012"] <- "Tollis et al. 2012"
names(final2)[4] = "Data_Source"
# plot on base map
p + geom_point(data=final2, aes(x=Longitude, y=Latitude, fill=Data_Source), color="black", pch=21, size=3) +
  theme(legend.title=element_text(size=15), legend.text = element_text(size=14), axis.text = element_text(size=12))                
