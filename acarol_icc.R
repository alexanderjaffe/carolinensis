# A. JAFFE CAROLINENSIS
# MARCH '15
# ICC Calculation Script

# set working directory
# setwd("./Data")

library(stringr)
library(ICC)
library(reshape)

# import data
data = read.csv("acarol_raw_morph.csv")
# get into proper format
data[data == 0] = NA
data[data == "#REF!"] = NA
data = na.omit(data)
data = data[2:nrow(data),]
colnames(data) = as.matrix(data[1,])
data = data[2:(nrow(data)-1),]
names(data)[1] = 'Observation.Number'
# create new names column
data$name.strip = str_sub(data$Observation.Number, start=0, end=6)
# load in full dataset
morph <- read.csv("Acar_morph_data_M_final.csv", header = TRUE)
# check overlap with actual dataset
keep = intersect(unique(male_morph$Observation.Number), unique(data$name.strip))
# limit to these individuals
data = data[data$name.strip %in% keep,]
# cast as numeric for each column
for (i in 2:29){
  data[,i] = as.numeric(as.character(data[,i]))}

# run ICC tests on each trait
options(warn=-1)
ICCbare(name.strip, Hindlimb.Phalanx, data=data)$ICC
ICCbare(name.strip, Hindlimb.Metatarsal, data=data)$ICC
ICCbare(name.strip, Tibia, data=data)$ICC
ICCbare(name.strip, Fibula, data=data)$ICC
ICCbare(name.strip, Femur, data=data)$ICC
ICCbare(name.strip, Pelvis.Width, data=data)$ICC
ICCbare(name.strip, Metacarpal, data=data)$ICC
ICCbare(name.strip, Ulna, data=data)$ICC
ICCbare(name.strip, Radius, data=data)$ICC
ICCbare(name.strip, Humerus, data=data)$ICC
ICCbare(name.strip, Opening.Inlever, data=data)$ICC
ICCbare(name.strip, Closing.Inlever, data=data)$ICC
ICCbare(name.strip, Whole.Head, data=data)$ICC
ICCbare(name.strip, Outlever, data=data)$ICC
ICCbare(name.strip, Snout.Plus.Eye, data=data)$ICC
ICCbare(name.strip, Snout, data=data)$ICC
ICCbare(name.strip, Real.Eye.Length, data=data)$ICC
ICCbare(name.strip, Braincase.Width, data=data)$ICC
ICCbare(name.strip, Head.Width.Retroarticulars, data=data)$ICC
ICCbare(name.strip, Head.Width.Jugals, data=data)$ICC
ICCbare(name.strip, Head.Width.Quadrates, data=data)$ICC
ICCbare(name.strip, Snout.Width.front.of.eye, data=data)$ICC
ICCbare(name.strip, Real.Lower.Jaw.Length, data=data)$ICC
ICCbare(name.strip, Real.Quadrate.to.Symphasis, data=data)$ICC
ICCbare(name.strip, Real.Jugal.to.Symphysis, data=data)$ICC
ICCbare(name.strip, Real.Orbit.to.Symphysis, data=data)$ICC
options(warn=0)

## ICC for SVL
svl = read.csv("svl.csv")
svl = svl[,1:4]
svl = svl[svl$Specimen.No. %in% keep,]
svl.melt = melt(svl, id.vars = "Specimen.No.")
ICCbare(Specimen.No., value, data=svl.melt)$ICC
