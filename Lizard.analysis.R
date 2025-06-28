#load packages needed
library(PCDimension)
library(geiger)
library(mvMORPH)
library(geomorph)
library(car)
library(ips)
library(scatterplot3d)
library(phytools) 

# set the directory to read in the *.nts files
setwd("/Users/tree/Desktop/Lizards/Pod 1")

filelist <- list.files(pattern=".nts")
dat <- readmulti.nts(filelist)

################################################
### Run a Generalized Procrustes Analysis
###   This creates a 3D description of shape from the raw coordinate data
################################################

Y.gpa <- gpagen(dat, Proj=TRUE, ProcD=TRUE, curves=NULL, surface=NULL)
summary(Y.gpa)
plot(Y.gpa)
y <- two.d.array(Y.gpa$coords)

attributes(Y.gpa$coords)$dimnames[[3]] <- 
  c("Sceloporus_acanthinus", "Sceloporus_adleri", "Sceloporus_aeneus", "Sceloporus_angustus", "Sceloporus_arenicolus", "Sceloporus_bicanthalis",  "Sceloporus_carinatus", "Sceloporus_chrysostictus", "Sceloporus_clarkii", "Sceloporus_consobrinus", "Sceloporus_couchii", "Sceloporus_cowlesi", "Sceloporus_cozumelae", "Sceloporus_cyanogenys", "Sceloporus_dugesii", "Sceloporus_edwardtaylori", "Sceloporus_formosus", "Sceloporus_graciosus", "Sceloporus_grammicus", "Sceloporus_grandaevus", "Sceloporus_hartwegi", "Sceloporus_horridus", "Sceloporus_internasalis", "Sceloporus_jalapae", "Sceloporus_jarrovii", "Sceloporus_magister", "Sceloporus_megalepidurus", "Sceloporus_merriami", "Sceloporus_mucronatus", "Sceloporus_occidentalis", "Sceloporus_olivaceus", "Sceloporus_orcutti", "Sceloporus_poinsettii", "Sceloporus_siniferus", "Sceloporus_spinosus", "Sceloporus_tristichus", "Sceloporus_variabilis", "Sceloporus_virgatus", "Sceloporus_woodi")
attributes(Y.gpa$Csize)$names <- attributes(Y.gpa$coords)$dimnames[[3]]


# Many objects are generated
attributes(Y.gpa) 
# 3D array of Procrustes coordinates
Y.gpa$coords 
# Vector of centroid sizes (these are a measure of the overall size of the structure
Y.gpa$Csize 
# The number of GPA iterations until convergence was found
Y.gpa$iter
# Variance-covariance matrix among landmark coordinates
Y.gpa$points.VCV
# Variances of landmark points
Y.gpa$points.var
# The consensus (mean) configuration; akin to using mshape
Y.gpa$consensus
# Data frame with an n x (pk) matrix of Procrustes residuals and centroid size
Y.gpa$data
# Final convergence criterion value
Y.gpa$Q

##creating a geomorph dataframe for downstream analyses
gfd <- geomorph.data.frame(Y.gpa)

## Checking outliers for possible landmarking problems
outliers <- plotOutliers(Y.gpa$coords) 
M <- mshape(Y.gpa$coords)   
## Specimens falling above the upper quartile are plotted in red and their 
## address returned, for inspection by plotRefToTarget

##Outlier Code Should not be needed anymore, but is here just in case##

# Example (for the first outlier)
#plotRefToTarget(M,Y.gpa$coords[,,outliers[1]], method="TPS", label = T)
#plotRefToTarget(M,Y.gpa$coords[,,outliers[1]], method="vector", label = T)

# Example (for the second outlier)
#plotRefToTarget(M,Y.gpa$coords[,,outliers[2]], method="vector", label = T)

# Example (for the second outlier)
#plotRefToTarget(M,Y.gpa$coords[,,outliers[3]], method="vector", label = T)

##This code allows to you see the deformation of all the species
#plotRefToTarget(M, Y.gpa$coords[,,4], method="vector", label = F, verbose=T)


##########################################################################
######## identify comparison groups ########## Update for Size!!! ########
##########################################################################

#Note that SVL size needs to be updated from Round 3 to be exact, and categorical size should also be updated

group <- read.csv("scelinfoV2.csv")
### species names here need to match those above exactly
group$species

head(group)

View(group)

### need to get data only for species for which we have CT scans 
filelist2 <- gsub(".nts", "", attributes(Y.gpa$Csize)$names)
vars <- group$species %in% filelist2
group2 <- group$species[vars]
group2 <- as.data.frame(group2)
names(group2) <- "species" 
group2 <- merge(group, group2)

arb <- group2$arb <- as.factor(group2$arb)
species <- group2$species <- as.character(group2$species)
size <- group2$SVL

## associating the categorical data with species it belongs to
names(arb) <- group2$species
arb
names(size) <- group2$species
size

gfd$arb <- arb
gfd$svl <- size

gfd$species <- attributes(gfd$Csize)$name

plot(gfd$Csize ~ size)


####################################################
# 2D distances
####################################################

#A = 37,41 (snout length)
#B = 41,45 (height at coronoid)
#C = 45,46 (snout width)
#D = 41,53 (closing-in lever)
#E = 51,53 (opening-in lever)
#F = 53,54 (jaw width)
#G = A + D (lower jaw length)

lmks <- matrix(c(37,41, 41,45, 45,46, 41,53, 51,53, 53,54), ncol=2, byrow=TRUE, 
               dimnames = list(c("A", "B", "C", "D", "E", "F"), c("start", "end")))
TwoD <- interlmkdist(dat, lmks)
TwoD <- data.frame(TwoD)

rownames(TwoD)<-c("Sceloporus_acanthinus", "Sceloporus_adleri", "Sceloporus_aeneus", "Sceloporus_angustus", "Sceloporus_arenicolus", "Sceloporus_bicanthalis",  "Sceloporus_carinatus", "Sceloporus_chrysostictus", "Sceloporus_clarkii", "Sceloporus_consobrinus", "Sceloporus_couchii", "Sceloporus_cowlesi", "Sceloporus_cozumelae", "Sceloporus_cyanogenys", "Sceloporus_dugesii", "Sceloporus_edwardtaylori", "Sceloporus_formosus", "Sceloporus_graciosus", "Sceloporus_grammicus", "Sceloporus_grandaevus", "Sceloporus_hartwegi", "Sceloporus_horridus", "Sceloporus_internasalis", "Sceloporus_jalapae", "Sceloporus_jarrovii", "Sceloporus_magister", "Sceloporus_megalepidurus", "Sceloporus_merriami", "Sceloporus_mucronatus", "Sceloporus_occidentalis", "Sceloporus_olivaceus", "Sceloporus_orcutti", "Sceloporus_poinsettii", "Sceloporus_siniferus", "Sceloporus_spinosus", "Sceloporus_tristichus", "Sceloporus_variabilis", "Sceloporus_virgatus", "Sceloporus_woodi")

#Lengths
gfd$A <- TwoD$A
gfd$B <- TwoD$B
gfd$C <- TwoD$C
gfd$D <- TwoD$D
gfd$E <- TwoD$E
gfd$F <- TwoD$F
gfd$G <- (gfd$A+gfd$D)

species

names(gfd$A) <- group2$species
names(gfd$B) <- group2$species
names(gfd$C) <- group2$species
names(gfd$D) <- group2$species
names(gfd$E) <- group2$species
names(gfd$F) <- group2$species
names(gfd$G) <- group2$species

gfd$A
gfd$B
gfd$C
gfd$D
gfd$E
gfd$F
gfd$G

##########################################
### 2D Adjustments for PC Graphs Below ###
##########################################

AS <- (gfd$A/size)
BS <- (gfd$B/size)
CS <- (gfd$C/size)
DS <- (gfd$D/size)
ES <- (gfd$E/size)
FS <- (gfd$F/size)
GS <- (gfd$G/size)

AS
BS
CS
DS
ES
FS
GS

#quartile analysis

plot(AS ~ arb)

#t test

ASt <- t.test(AS ~ arb)
ASt

##########################################
########## IGNORE ##########
##########################################

plot(AS ~ arb)

# Create a new dataset including the measurements with the species name and arb variable
nosize <- data.frame(
  species = gfd$species,
  arb = gfd$arb,  # Adding the arb variable here
  A = gfd$A,
  B = gfd$B,
  C = gfd$C,
  D = gfd$D,
  E = gfd$E,
  F = gfd$F,
  G = gfd$G
)

# Write the new dataset to a CSV file
write.csv(nosize, row.names = FALSE)

View(nosize)

wsize <- data.frame(
  species = gfd$species,
  arb = gfd$arb,  # Adding the arb variable here
  A = AS,
  B = BS,
  C = CS,
  D = DS,
  E = ES,
  F = FS,
  G = GS
)

# Write the new dataset to a CSV file
write.csv(wsize, row.names = FALSE)

View(wsize)

number_of_rows <- nrow(nosize)

# Print the number of rows
print(number_of_rows)


# Example t-test
result <- t.test(AS ~ arb)
p_value <- result$p.value

  
  
  
  

