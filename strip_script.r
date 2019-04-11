#----------------------------------------------------#
# # load required packages
# install.packages("backports")
# install.packages("rioja")
# install.packages("dplyr")
# install.packages("vegan")
# install.packages("tidyverse")
# install_github("richardjtelford/ggpalaeo")
# install.packages("devtools")
# install.packages("ggpalaeo")

library(rioja)
library(dplyr)
library(vegan)
library(tidyverse)
library(devtools)
library(backports)
library(ggpalaeo)
# citation("rioja")
# citation("tidyverse")
# citation("vegan")
# citation("dplyr")
# citation("devtools")
# citation("backports")
# citation("RStudio")
# citation()
# RStudio.Version()
#setRepositories()
#ap <- available.packages()

setwd("C:/Users/krist/Dropbox/Master/kristina_strip/data")
pol <- read.csv("stripage.csv", header = FALSE, row.names = 1, sep = ",", dec = ".")
#-------------------------------------------------#
# START! 

pol <- as.data.frame(t(pol)) #transpose matrix
head(pol)

pol <- pol[nrow(pol):1, ]
pol[is.na(pol)] <- 0 # turn no data values into 0s
depth <- as.numeric(pol[,"Depth"])
row.names(pol) <- depth
head(pol)

#----------------------------------------------------#
# Age modelling
setwd("C:/Users/krist/Dropbox/Master/kristina_strip/clam")
source("clam.R")
clam("Strip", type=4, sep="," , depthseq = seq(0, 310, 0.5))
agesTab <- read.table("Cores/Strip/Strip_smooth_spline_ages.txt", header=TRUE, sep="\t")

# match up depths with ages
wantAge <- agesTab$depth %in% depth
age <- agesTab$best[wantAge]
accrate <- agesTab$accrate[wantAge]
plot(age, depth)
plot(age, accrate)
setwd("C:/Users/krist/Dropbox/Master/kristina_strip/data") # set the working directory back to the main directory
pol$Age <- age # overwrites the age file you had on the original dataset with the new ages so you don't have any confusion

# Now you have the age model calculations built into the script !

#----------------------------------------------------#

# Influx data. Calculated for later.
strip <- read.csv("Strip_influx_counts.csv", sep = ",", dec = ".")
head(strip)

#data.frame av kun pollentellingene
pollen_counts <- strip %>% select(Alnus:Vaccinium.sp.)
strip$accrate <- accrate # overwriting acc.rate in this object so we can use those calculated from the age model
strip$Age <- age # overwriting Age in this object so we can use those calculated from the age model
head(strip)

# Funksjon som gir kalkuleringene av pollenkonsentrasjonene og influx verdier basert på sedimentasjonsraten
#fra dybde-aldersmodellen, f.eks CLAM. Viktig: Er det deposition_time må du endre koden under slik at du deler
#pollenkonsentrasjonene på deposition_time (dvs influx_rate = pollen_concentration/deposition_time); lurt å skifte ut navnet
#på variabelen sedimentation_rate med deposition_time.  ### FRA CLAM ER ACCRATE LIK DEPOSITION TIME!
conc_influx <- function(df, Marker, Marker_added, Sed_vol, accrate){
  pollen_concentration <- (df/Marker) * (Marker_added/Sed_vol)   
  pollen_concentration <- round(pollen_concentration, digits = 1)
  influx_rate <- round(pollen_concentration/accrate, digits = 1)
  results <- list(pollen_concentration = pollen_concentration, influx_rate = influx_rate)}

# returnerer en liste som inneholder både datasettet med pollenkonsentrasjoner og influx-verdier-
influx_strip <- conc_influx(df = pollen_counts, Marker = strip$Marker, Marker_added = strip$Marker_added, Sed_vol = strip$Sed_vol,
                            accrate = strip$accrate)

#data
influx_strip$pollen_concentration
influx_strip$influx_rate
influx_strip$sedimentation_rate ##Stod originalt influx_heim$deposition_time, men null der også, regner med dette er siden det var de verdiene jeg allerede hadde (accrate)

#write csv file e.g.
#write.csv(influx_strip$pollen_concentration, "pollen_concentration_strip.csv", row.names = depth)
#write.csv(influx_strip$influx_rate, "influx_rate_strip.csv", row.names = depth)

#-------------------------------------------------------#

#LAGER GRUPPER:

TREES <- c("Betula", "Alnus", "Carpinus", "Corylus", "Fagus", "Picea", "Pinus", "Populus", "Quercus", "Sorbus", "Tilia", "Ulmus")

HERBS <- c("Apiaceae", "Artemisia", "Asteraceae", "Beckwithia.glacialis", "Chenopodiaceae", "Cirsium.sp.",
           "Compositae.cich.",
           "Cyperaceae", "Drosera", "Filipendula", "Geranium", "Humulus", "Linnaea.borealis", "Lotus.sp.", "Melampyrum",
           "Onagraceae", "Oxyria.sp.", "Poaceae", "Potentilla.sp.", "Ranunculus.acris", "Ranunculus.flammula",
           "Ranunculus.sp.", "Rosaceae.undiff.",
           "Rubiaceae", "Rubus.chamaemorus", "Rumex.acetosa", "Rumex.sp.", "Saussurea.alpina", "Sedum", "Solidago.sp.",
           "Stachys.sp.", "Thalictrum")

SHRUBS <- c("Hippophae", "Juniperus", "Salix")

DSHRUBS <- c("Calluna", "Empetrum", "Ericaceae.undiff.", "Vaccinium.sp.", "Betula.nana")

UNKNOWN <- c("Unknown")

WATER <- c("Menyanthes", "Potamogeton.eu.", "Sparganium", "Hippuris")

FAS <- c("Botryococcus", "Dryopteris.sp.", "Equisetum", "Gymnocarpium.dryopteris", "Lycopodium.annotinum", 
         "Lycopodium.clavatum", "Pteridium", "Selaginella", "Pediastrum", "Arnium.sp.", "Scenedesmus")

MOSS <- c("Sphagnum")

CHC <- c("Charcoal.small", "Charcoal.big")

MARKER <- c("Marker")

STOMATA <- c("Stomata.Pinus", "Stomata.unknown")

# Regner ut prosent for alle gruppene jeg vil ha videre:
# TREES, SHRUBS, HERBS, UNKNOWN:
sumT <- rowSums(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN)])
#polPercent <- (pol/ sumT)*100

# PLUSS FAS:
sumTFAS <- rowSums(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS)])
#polpercent1 <- (pol/ sumTFAS)*100

# PLUSS WATER:

sumTFASW <- rowSums(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER)])
#polpercent2 <- (pol/ sumTFASW)*100

# PLUSS MOSER:

sumTFASWM <- rowSums(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER, MOSS)])
#polpercent3 <- (pol/ sumTFASWM)*100

# PLUSS KULL:

sumTFASWMC <- rowSums(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER, MOSS, CHC)])

Newpol <- cbind(pol[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN)]/sumT,
                pol[,c(FAS)]/sumTFAS,
                pol[,c(WATER)]/sumTFASW,
                pol[,c(CHC)]/sumTFASWMC,
                pol[,c(MOSS), drop = FALSE]/sumTFASWM) * 100

wanted <- c("Alnus", "Betula", "Corylus", "Picea", "Pinus", "Ulmus", "Juniperus", "Salix", "Betula.nana",
            "Vaccinium.sp.", "Cyperaceae", "Poaceae", "Dryopteris.sp.", "Equisetum",
            "Gymnocarpium.dryopteris", "Botryococcus", "Pediastrum", "Sphagnum")


Newpola <- Newpol[, wanted]
head(Newpola)
colnames(Newpola)
Newpola
#write.csv(Newpol, "pollen_prosen_strip.csv")
LOI <- pol$LOI
setdiff(names(pol), c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, WATER, FAS, MOSS, CHC, MARKER, STOMATA, LOI))
setdiff(wanted, names(Newpola)) #Sjekke om noen av navnene/taxaene ikke er med i noen gruppe

#----------------------------------------------------#
#Preparing the influx data:

polinf <- read.csv("influx_rate_strip.csv", header = TRUE, sep = ",", dec = ".")
head(polinf)
polinf[is.na(polinf)] <- 0 # turn no data values into 0s
#plot(age, polinf$Alnus)

wanted <- c("Alnus", "Betula", "Corylus", "Picea", "Pinus", "Ulmus", "Juniperus", "Salix", "Betula.nana",
            "Vaccinium.sp.", "Cyperaceae", "Poaceae", "Dryopteris.sp.", "Equisetum",
            "Gymnocarpium.dryopteris", "Sphagnum", "Total")
wanted

Newpolinf <- polinf[, wanted]/1000
Newpolinf <- rename(Newpolinf, `Gymnocarpium` = `Gymnocarpium.dryopteris`, `Betula nana` = `Betula.nana`, 
                    `Vaccinium` = `Vaccinium.sp.`, `Dryopteris` = `Dryopteris.sp.`)
setdiff(wanted, names(Newpolinf))

#----------------------------------------------------#
#Rarify richness
polR.df<- pol[,c(TREES, HERBS, SHRUBS, DSHRUBS, UNKNOWN)]
countTotals <- rowSums(polR.df) # You counted 2000 grains in one sample! Woah! Is that correct?
countTotals
rarefy.n <- 432

#rarefy.n <- sample(countTotals) # maybe you should standardise this between the two sites (e.g. rarify to the lowest count sum across both sites?)

strip_richness <- rarefy(polR.df, rarefy.n) 
strip_richness
plot(age, strip_richness, type ="l")
#write.csv(strip_richness, "Richness_Strip.csv")

#------------------------------------------#
# Preparing the percentage pollen data for the plotting and the stats
Newpol2 <- rename(Newpola, `Gymnocarpium` = `Gymnocarpium.dryopteris`, `Betula nana` = `Betula.nana`, 
                   `Vaccinium` = `Vaccinium.sp.`, `Dryopteris` = `Dryopteris.sp.`)

Newpol2
# Percentages 
pollenstrip <- pol[,c(TREES, SHRUBS, DSHRUBS, HERBS)]
pollenstrip2 <- rename(pollenstrip, `Betula nana` = `Betula.nana`, `Rubus` = `Rubus.chamaemorus`)
pollenstrip_perc <- pollenstrip2/rowSums(pollenstrip2) *100
pollenstrip_perc <- pollenstrip_perc[, colSums(pollenstrip_perc>0)>2]

# Statistical analysis
#1: Cluster og soner:
ma.dist <- vegdist(sqrt(pollenstrip_perc), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")
plot(ma.chclust)
b2 <- bstick(ma.chclust, 10)   #MEN KJØR FOR ALLE TERR. TAXA

cbind(cutree(ma.chclust, k = 4), age)

# 2. ordination 

# DCA
pollenstrip.dca <- decorana(sqrt(pollenstrip_perc))
scor_s1 <- scores(pollenstrip.dca, choices = 1:2, display = "sites")
summary(pollenstrip.dca)
pollenstrip.dca
#plot(pollenstrip.dca)

# PCA 
pollenstrip.pca <- rda(sqrt(pollenstrip_perc))
scor_s2 <- scores(pollenstrip.pca, choices = 1:2, display = "sites")

#RICHNESS
scor_s3 <- scores(strip_richness, choices = 1:2, display = "sites")

# Sorting of data ready for plotting
var_s1 <- as.matrix(cbind(LOI, scor_s2, scor_s3))
colnames(var_s1) <- c("LOI", "PCA1", "PCA2", "Species richness")
var_s2 <- as.data.frame(cbind(Newpol[,CHC], pol[,STOMATA]))  
var_s2 <- var_s2 %>%
  mutate(Charcoal = Charcoal.big + Charcoal.small) %>% 
  select(-Stomata.unknown, -Charcoal.big, -Charcoal.small)

#----------------------------------------------------#
####### PLOTTING THE PERCENTAGE DATA
cc <- c(rep("chartreuse1", 6), rep("blue", 2), rep("cyan", 2), rep("mediumorchid1", 2), rep("brown1", 5), rep
        ("gold", 1)) #Velger farger
par(omi = c(0.3, 0.5, 0.6, 0.6)) #Fikse marginen til venstre. http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/#comment-173
par(mgp = c(0, 2, 1))

#Plott med cluster
pol.plot <- strat.plot(Newpol2,
                       yvar = age, y.rev=TRUE, scale.percent=TRUE, plot.poly=TRUE, col.poly=cc, 
                       plot.bar=TRUE, col.bar="grey69", col.poly.line="grey40", exag=TRUE, col.exag="auto", 
                       exag.alpha = 0.4, exag.mult=10, x.pc.inc = 5, xSpace=0.01, x.pc.lab=TRUE, 
                       x.pc.omit0=FALSE, srt.xlabel=60, ylabel=NA, cex.ylabel=1.05,
                       cex.yaxis=1, cex.xlabel=1, xRight = 0.65, xLeft = 0.06, yTop=0.805, yBottom=0.01, 
                       y.tks= c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500,
                                              6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000),
                       title="", graph.widths = 1)

# ADDING STOMATA/CHARCOAL TO THE PLOT
char.plots <- strat.plot(rename(var_s2, `Stomata Pinus` = `Stomata.Pinus`), xRight = 0.85, xLeft = 0.65, yTop=0.805, 
                         yBottom=0.01, 
                         yvar = age, y.rev=TRUE, plot.poly = FALSE,  plot.line = FALSE,  
                         y.axis = FALSE, plot.bar = TRUE, lwd.bar = 1.8, col.bar = "firebrick3", 
                         plot.symb = FALSE, srt.xlabel = 60, cex.xlabel = 1, 
                         xSpace = 0.01, add = TRUE, graph.widths = 0.4) 

# ADDING LOI/PCA/RICHNESS TO THE PLOT
pc.plots <- strat.plot(var_s1, xRight = 1, xLeft = 0.832, yTop=0.805, 
                         yBottom=0.01,
                         yvar = age, y.rev= TRUE, plot.poly = FALSE, plot.line = TRUE, col.line = "forestgreen", lwd.line = 1.8,
                         y.axis = FALSE,
                         plot.bar = FALSE, plot.symb = FALSE, srt.xlabel = 60, 
                         cex.xlabel = 1, xSpace = 0.01, add = TRUE, graph.widths = 0.4)

secondary.scale(yvar = rev(c(0,age)), yvar2 = rev(c(0,depth)), n = 6, y.rev = TRUE, xLeft = 0.02, yBottom = 0.01,
                yTop=0.812, ylabel2 = NA, cex.ylabel2 = 1)

#Add cluster zones
addClustZone(pol.plot, ma.chclust, nZone=4, lwd=1.5, lty=2, col="grey25")
addClustZone(pc.plots, ma.chclust, nZone=4, lwd=1.5, lty=2, col="grey25")
addClustZone(char.plots, ma.chclust, nZone=4, lwd=1.5, lty=2, col="grey25")

text(x=0.61, y=1080, labels="Age (cal. yrs BP)", srt=60, xpd=NA)
text(x=0.58, y=1320, labels="Depth (cm)", srt=60, xpd=NA)
text(x=0.62, y=1220, labels="Zone names", srt=60, xpd=NA)
text(x=0.91, y=25, labels="Stripåtmyrin 990 m a.s.l.", cex=1.5, xpd=NA)

# text(x=0.62, y=9, labels = "ST-1", xpd=NA, cex=0.9)
# text(x=0.62, y=8700, labels = "ST-2", xpd=NA, cex=0.9)
# text(x=0.62, y=6712, labels = "ST-\n3.1", xpd=NA, cex=0.9)
# text(x=0.62, y=8000, labels = "ST-\n3.2", xpd=NA, cex=0.9)
# text(x=0.62, y=4950, labels = "ST-\n4.1", xpd=NA, cex=0.9)
# text(x=0.62, y=2500, labels = "ST-\n4.2", xpd=NA, cex=0.9)
# 
# text(x=0.617, y=9453, labels = "Limus, \nsilt, \ndark \n(2-)", xpd=NA, cex=0.8)
# text(x=0.617, y=6000, labels = "Limus, \ndetritus, \npeat, \ndark \n(2)", xpd=NA, cex=0.8)
# text(x=0.617, y=4000, labels = "Limus, \ndetritus, \npeat, \ndark \n(2)", xpd=NA, cex=0.8)

#----------------------------------------------------#
# TOTAL DIAGRAM PROSENT:
#GROUPS:

TREES <- c("Betula", "Alnus", "Carpinus", "Corylus", "Fagus", "Picea", "Pinus", "Populus", "Quercus", "Sorbus", "Tilia", "Ulmus")

HERBS <- c("Apiaceae", "Artemisia", "Asteraceae", "Beckwithia.glacialis", "Chenopodiaceae", "Cirsium.sp.",
           "Compositae.cich.",
           "Cyperaceae", "Drosera", "Filipendula", "Geranium", "Humulus", "Linnaea.borealis", "Lotus.sp.", "Melampyrum",
           "Onagraceae", "Oxyria.sp.", "Poaceae", "Potentilla.sp.", "Ranunculus.acris", "Ranunculus.flammula",
           "Ranunculus.sp.", "Rosaceae.undiff.",
           "Rubiaceae", "Rubus.chamaemorus", "Rumex.acetosa", "Rumex.sp.", "Saussurea.alpina", "Sedum", "Solidago.sp.",
           "Stachys.sp.", "Thalictrum")

SHRUBS <- c("Hippophae", "Juniperus", "Salix")

DSHRUBS <- c("Calluna", "Empetrum", "Ericaceae.undiff.", "Vaccinium.sp.", "Betula.nana")

# POLLEN PERCENTAGES LIFEFORM GROUPS

trees_p <- rowSums(pollen_counts[,c(TREES)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

shrubs_p <- rowSums(pollen_counts[,c(SHRUBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

dshrubs_p <- rowSums(pollen_counts[,c(DSHRUBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

herbs_p <- rowSums(pollen_counts[,c(HERBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

# THREES;SHRUBS; HERBS; DSHRUBS
plot(strip$Age, trees_p,  type = "n", ylim = c(0, 100), xlim = c(0, 9700), ylab = "Groups pollen %", xlab = "Age (cal. yrs BP)", axes = F, cex.lab = 1.4)
par(mgp = c(0, 1, 0))
axis(1, at = seq(0, 9700, 500))
axis(2, at = seq(0, 100, 20), las = 1, tcl = -0.2)
xx <- c(strip$Age, rev(strip$Age))
yy1 <- c(rep(0, 37), rev(trees_p))
polygon(xx, yy1, col = "chartreuse1")
#Hvis man vil fjerne sorte linjer -> border = NA inni hvert polygon

yy2 <- c(trees_p, rev(trees_p) + rev(shrubs_p))
polygon(xx, yy2, col = "blue")

yy3 <- c(trees_p + shrubs_p, rev(trees_p)+ rev(shrubs_p)+ rev(dshrubs_p))
polygon(xx, yy3, col = "cyan")

yy4 <- c(trees_p + shrubs_p + dshrubs_p, rev(trees_p) + rev(shrubs_p) + rev(dshrubs_p) + rev(herbs_p))
polygon(xx, yy4, col = "mediumorchid1")

legend(x = 7620, y = 20, c("TREES","SHRUBS", "DSHRUBS", "HERBS"), fill= c("chartreuse1", "blue", "cyan", "mediumorchid1"), cex = 0.7)
box()
mtext("Stripåtmyrin, 990 m a.s.l.", side = 3, cex = 1.7)

#abline(v = ) ---> "v" for vertikale linjer, sett på en for hver sone. 
dev.off()
# # TOTAL POLLEN DIAGRAM; TRANSPOSED
par(mgp = c(3, 0.5, 0))
par(omi = c(.1, .2, 0, 0))
xx <- c(strip$Age, rev(strip$Age))
yy1 <- c(rep(0, 37), rev(trees_p))

plot(trees_p, strip$Age, type ="n", ylim = c(9700, 0), xlim = c(0, 100), xlab = "Groups %", axes = FALSE, ylab = "Age (cal. yrs BP)", cex.lab = 1.5)
axis(1, at = seq(0, 100, 10))
axis(2, at = seq(0, 10000, 500), las = 1, tcl = -0.2)
#axis(3, labels = FALSE)
polygon(yy1, xx, col = "chartreuse1")

yy2 <- c(trees_p, rev(trees_p) + rev(shrubs_p))
polygon(yy2, xx, col = "blue")

yy3 <- c(trees_p + shrubs_p, rev(trees_p)+ rev(shrubs_p)+ rev(dshrubs_p))
polygon(yy3,xx, col = "cyan")

yy4 <- c(trees_p + shrubs_p + dshrubs_p, rev(trees_p) + rev(shrubs_p) + rev(dshrubs_p) + rev(herbs_p))
polygon(yy4, xx, col = "mediumorchid1")
legend(x = 0.7, y = 180, c("TREES","SHRUBS", "DSHRUBS", "HERBS"), fill= c("chartreuse1", "blue", "cyan", "mediumorchid1"), cex = 0.8)
box()
mtext("Stripåtmyrin, 990 m a.s.l.", side = 3, cex = 1.7)
abline(h = c(5300, 8750, 9300), lty=4, lwd=1, col="black")
abline(h = c(3100, 7000), lty=3, lwd=1, col="grey47")
text(x=-2, y=9600, labels = "ST-1", xpd=NA, cex=0.9)
text(x=-2, y=9000, labels = "ST-2", xpd=NA, cex=0.9)
text(x=-2, y=7900, labels = "ST-\n3.1", xpd=NA, cex=0.9)
text(x=-2, y=6200, labels = "ST-\n3.2", xpd=NA, cex=0.9)
text(x=-2, y=4300, labels = "ST-\n4.1", xpd=NA, cex=0.9)
text(x=-2, y=1500, labels = "ST-\n4.2", xpd=NA, cex=0.9)

#----------------------------------------------------------#
# PCA plot
dev.off()
par(mfrow=c(1,1), mar= c(5,4,4,2))

loadPC7 = order(abs(summary(pollenstrip.pca)$species[,1]), decreasing=T)[1:9]
loadPC8 = order(abs(summary(pollenstrip.pca)$species[,2]), decreasing=T)[1:9]

plot(summary(pollenstrip.pca)$sites[,1], summary(pollenstrip.pca)$sites[,2], cex=1, xlab="PC1", ylab="PC2", type="n",
     ylim=c(-1,2), xlim=c(-1.7,1.5), main= "PCA Stripåtmyrin (990 m a.s.l.)")

abline(h=0, v=0, lty=2)

specPC5 = summary(pollenstrip.pca)$sites[,1]
specPC6 = summary(pollenstrip.pca)$sites[,2]
spec1 = c(loadPC7, loadPC8)

text(x=c(specPC5), y=c(specPC6), labels=names(specPC5), cex=1)

specPC7 = summary(pollenstrip.pca)$species[,1]
specPC8 = summary(pollenstrip.pca)$species[,2]

arrows(x0=rep(0,24), y0=rep(0,24), x1=c(specPC7[spec1]), y1=c(specPC8[spec1]), col="dark red", length=0.07)
text(x=c(specPC7[spec1]), y=c(specPC8[spec1]), labels=names(specPC7[spec1]), cex=1, col="dark red", font=3)

# alternative
#plot(pollenstrip.pca) # see here https://www.fromthebottomoftheheap.net/2013/01/13/decluttering-ordination-plots-in-vegan-part-2-orditorp/

#----------------------------------------------------#
#Plotting the influx diagram
cc <- c(rep("chartreuse1", 6), rep("blue", 2), rep("cyan", 2), rep("mediumorchid1", 2), rep("brown1", 3),
        rep("gold", 1), rep("grey40", 1)) #Velger farger
par(omi = c(0.3, 0.3, 0.6, 0.6)) #Fikse marginen til venstre. http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/#comment-173
par(mgp = c(0, 2, 1))

polinf.plot <- strat.plot(Newpolinf, 
                          yvar = age, y.rev=TRUE, scale.percent=FALSE, scale.minmax = TRUE, plot.poly=TRUE, col.poly=cc, 
                          plot.bar= TRUE, col.bar="grey50", col.poly.line="grey69", exag=TRUE, col.exag="auto",
                          exag.alpha = 0.4, exag.mult=10, x.pc.inc = 5000, x.pc.lab=TRUE, x.pc.omit0=FALSE,
                          srt.xlabel=60, ylabel="", cex.ylabel=1, cex.yaxis=1, cex.axis = 0.6,
                          cex.xlabel=1.2, title="", xLeft = 0.09, yBottom = 0.05, 
                          y.tks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 
                                    4000, 4500, 5000, 5500, 6000, 6500, 7000, 
                                    7500, 8000, 8500, 9000, 9500, 10000), graph.widths=1)

secondary.scale(yvar = rev(c(0,age)), yvar2 = rev(c(0,depth)), n = 6, y.rev = TRUE,  cex.ylabel2 = 1, 
                xLeft = 0.05, yBottom = 0.05, yTop=0.81)  

text(x=0.63, y=1080, labels="Age (cal. yrs BP)", srt=60, xpd=NA)
text(x=0.60, y=1320, labels="Depth (cm)", srt=60, xpd=NA)
text(x=0.65, y=1200, labels="Zone names", srt=60, xpd=NA)
text(x=0.90, y=25, labels="Stripåtmyrin 990 m a.s.l.", cex=1.5, xpd=NA)

addClustZone(polinf.plot, ma.chclust, nZone=4, lwd=1.2, lty=2, col="grey25")

# text(x=0.6432, y=9020, labels = "ST-1", xpd=NA, cex=0.9)
# text(x=0.6432, y=8700, labels = "ST-2", xpd=NA, cex=0.9)
# text(x=0.643, y=6700, labels = "ST-3.1", xpd=NA, cex=0.9)
# text(x=0.643, y=8000, labels = "ST-3.2", xpd=NA, cex=0.9)
# text(x=0.643, y=4950, labels = "ST-4.1", xpd=NA, cex=0.9)
# text(x=0.643, y=2500, labels = "ST-4.2", xpd=NA, cex=0.9)
# 
# text(x=0.6425, y=9300, labels = "Limus, \nsilt, dark \n(2-)", xpd=NA, cex=0.8)
# text(x=0.6425, y=6000, labels = "Limus, \ndetritus, \npeat, \ndark \n(2)", xpd=NA, cex=0.8)
# text(x=0.6425, y=4000, labels = "Limus, \ndetritus, \npeat, \ndark \n(2)", xpd=NA, cex=0.8)

#abline(h=c(2750, 8150, 9075), col="grey50")
Newpolinf
#---------------------------------------------------------------#
#TOTAL DIAGRAM INFLUX:
library(tidyverse)
#GROUPS:

TREES <- c("Betula", "Alnus", "Carpinus", "Corylus", "Fagus", "Picea", "Pinus", "Populus", "Quercus", "Sorbus", "Tilia", "Ulmus")

HERBS <- c("Apiaceae", "Artemisia", "Asteraceae", "Beckwithia.glacialis", "Chenopodiaceae", "Cirsium.sp.", "Compositae.cich.", "Cyperaceae", "Drosera", "Filipendula", "Geranium", "Humulus", "Linnaea.borealis", "Lotus.sp.", "Melampyrum", "Onagraceae", "Oxyria.sp.", "Poaceae", "Potentilla.sp.", "Ranunculus.acris", "Ranunculus.flammula", "Ranunculus.sp.", "Rosaceae.undiff.", "Rubiaceae", "Rubus.chamaemorus", "Rumex.acetosa", "Rumex.sp.", "Saussurea.alpina", "Sedum", "Solidago.sp.", "Stachys.sp.", "Thalictrum")

SHRUBS <- c("Hippophae", "Juniperus", "Salix")

DSHRUBS <- c("Calluna", "Empetrum", "Ericaceae.undiff.", "Vaccinium.sp.", "Betula.nana")

# # TOTAL POLLEN DIAGRAM; TRANSPOSED
par(mgp = c(3, 0.5, 0))
par(omi = c(.1, .2, 0, 0))
xx <- c(strip$Age, rev(strip$Age))
yy1 <- c(rep(0, 37), rev(TREES))

plot(TREES, strip$Age, type ="n", ylim = c(20, 0), xlim = c(0, 100), xlab = "Groups influx", axes = FALSE, ylab = "Age (cal. yrs BP)", cex.lab = 1.5)
axis(1)
axis(2, at = seq(0, 9700, 1000), las = 1, tcl = -0.2)
#axis(3, labels = FALSE)
polygon(yy1, xx, col = "chartreuse1")

yy2 <- c(TREES, rev(TREES) + rev(SHRUBS))
polygon(yy2, xx, col = "blue")

yy3 <- c(TREES + SHRUBS, rev(TREES)+ rev(SHRUBS)+ rev(DSHRUBS))
polygon(yy3,xx, col = "cyan")

yy4 <- c(trees_p + shrubs_p + dshrubs_p, rev(trees_p) + rev(shrubs_p) + rev(dshrubs_p) + rev(herbs_p))
polygon(yy4, xx, col = "mediumorchid1")
legend(x = 0.4, y = 180, c("TREES","SHRUBS", "DSHRUBS", "HERBS"), fill= c("chartreuse1", "blue", "cyan", "mediumorchid1"), cex = 0.8)
box()
mtext("Stripåtmyrin, 990 m a.s.l.", side = 3, cex = 1.7)
abline(h = c(2714, 8143, 9071), lty=4, lwd=1, col="black")

#mtext("Age cal. yrs BP", side = 2,line = 3)
#box(bty="c", lwd = 0.5)
#mtext("TREES", side = 3, line = 1.5)
# 
# ?mtext
# #need to add names on top - not done.
# 
# 
# 
# 

















