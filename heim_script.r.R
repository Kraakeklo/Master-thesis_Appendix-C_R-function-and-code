#----------------------------------------------------#
# load required packages
library(rioja)
library(dplyr)
library(vegan)
library(tidyverse)
#devtools::install_github("richardjtelford/ggpalaeo")
library(devtools)
library(backports)
library(ggpalaeo)

setwd("C:/Users/krist/Dropbox/Master/kristina_heim/data") # see edits to make this work
# START! 
polh <- read.csv("heim_u191.csv", header = FALSE, row.names = 1, sep=",", dec=".")
polh <- as.data.frame(t(polh)) #transpose matrix
head(polh)

# EDIT: Far better to get samples in order here, otherwise it causes confusion when doing the plotting later.
polh <- polh[nrow(polh):1, ]
polh[is.na(polh)] <- 0 #turn no data values into 0s
depth <- as.numeric(polh[,"Depth"]) #extract depth variable
row.names(polh) <- depth
head(polh)

#----------------------------------------------------#
# Age modelling
setwd("C:/Users/krist/Dropbox/Master/kristina_heim/clam")
source("clam.R")
clam("Heim", type=4, sep="," , depthseq = seq(0, 202, 0.5))
agesTab <- read.table("Cores/Heim/Heim_smooth_spline_ages.txt", header=TRUE, sep="\t")

# match up depths with ages
wantAge <- agesTab$depth %in% depth
age <- agesTab$best[wantAge]
accrate <- agesTab$accrate[wantAge]
plot(age, depth)
accrate
setwd("C:/Users/krist/Dropbox/Master/kristina_heim/data") # set the working directory back to the main directory
polh$Age <- age # overwrites the age file you had on the original dataset with the new ages so you don't have any confusion

# Now you have the age model calculations built into the script !

#----------------------------------------------------#

# Influx data. Calculated for later.

heim <- read.csv("Heim_influx_counts.csv", sep = ",", dec = ".")
head(heim)

#data.frame av kun pollentellingene
pollen_counts <- heim %>% select(Achillea.t.:Vaccinium.sp.)
heim$accrate <- accrate # overwriting acc.rate in this object so we can use those calculated from the age model
heim$Age <- age # overwriting Age in this object so we can use those calculated from the age model
head(heim)


# Funksjon som gir kalkuleringene av pollenkonsentrasjonene og influx verdier basert pÃ‚ sedimentasjonsraten
#fra dybde-aldersmodellen, f.eks CLAM. Viktig: Er det deposition_time mÃ‚ du endre koden under slik at du deler
#pollenkonsentrasjonene pÃ‚ deposition_time (dvs influx_rate = pollen_concentration/deposition_time); lurt Ã‚ skifte ut navnet
#pÃ‚ variabelen sedimentation_rate med deposition_time.  ### FRA CLAM ER ACCRATE LIK DEPOSITION TIME!
conc_influx <- function(df, Marker, Marker_added, Sed_vol, accrate){
   pollen_concentration <- (df/Marker) * (Marker_added/Sed_vol)   
   pollen_concentration <- round(pollen_concentration, digits = 1)
   influx_rate <- round(pollen_concentration/accrate, digits = 1)
   results <- list(pollen_concentration = pollen_concentration, influx_rate = influx_rate)}

# returnerer en liste som inneholder bÃ‚de datasettet med pollenkonsentrasjoner og influx-verdier-
influx_heim <- conc_influx(df = pollen_counts, Marker = heim$Marker, Marker_added = heim$Marker_added, Sed_vol = heim$Sed_vol,
                      accrate = heim$accrate)
  
#data
influx_heim$pollen_concentration
influx_heim$influx_rate
influx_heim$sedimentation_rate ##Stod originalt influx_heim$deposition_time, men null der ogsÃ‚, regner med dette er siden det var de verdiene jeg allerede hadde (accrate)

#write.csv(influx_heim$pollen_concentration, "pollen_concentration_heim_feb_u191.csv", row.names= depth)
#write.csv(influx_heim$influx_rate, "influx_rate_heim_feb_u1912.csv", row.names = depth)


#----------------------------------------------------#
#LAGER GRUPPER:
TREES <- c("Betula", "Alnus", "Carpinus", "Corylus", "Fagus", "Picea", "Pinus", "Populus", "Quercus", "Sorbus", "Tilia", 
           "Ulmus")

HERBS <- c("Achillea.t.", "Alchemilla", "Apiaceae", "Arenaria", "Artemisia", "Asteraceae", "Beckwithia.glacialis",
           "Caryophyllaceae", "Chenopodiaceae", "Cirsium.sp.", "Compositae.cich.", "Cyperaceae", "Drosera", "Epilobium", 
           "Filipendula", "Geranium", "Humulus", "Jasione.montana", "Linnaea.borealis", "Lotus.sp.", "Melampyrum",
           "Onagraceae", "Oxyria.sp.", "Plantago.media", "Poaceae", "Potentilla.sp.", "Ranunculus.acris",
           "Ranunculus.flammula", "Ranunculus.sp.", "Rosaceae.undiff.", "Rubiaceae", "Rubus.chamaemorus",
           "Rumex.acetosa", "Rumex.sp.", "Saussurea.alpina", "Sedum",
           "Solidago.sp.", "Stachys.sp.", "Thalictrum", "Trientalis")

SHRUBS <- c("Hippophae", "Juniperus", "Salix")

DSHRUBS <- c("Calluna", "Empetrum", "Ericaceae.undiff.", "Vaccinium.sp.", "Betula.nana")

UNKNOWN <- c("Indet:Corroded", "Unknown")

WATER <- c("Menyanthes", "Potamogeton.eu.", "Sparganium", "Hippuris")

FAS <- c("Botryococcus", "Dryopteris.sp.", "Equisetum", "Gymnocarpium.dryopteris", "Huperzia.selago", "Lycopodium.annotinum",
         "Lycopodium.clavatum", "Pteridium", "Selaginella", "Pediastrum")

MOSS <- c("Sphagnum")

CHC <- c("Charcoal.small", "Charcoal.big")

MARKER <- c("Marker")

STOMATA <- c("Stomata.Pinus", "Stomata.Unknown")

# Regner ut pollenprosent for alle gruppene jeg vil ha videre:
# TREES, SHRUBS, HERBS, UNKNOWN:
sumT <- rowSums(polh[,c(TREES, HERBS, SHRUBS, DSHRUBS, UNKNOWN)])

# PLUSS FAS:
sumTFAS <- rowSums(polh[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS)])

# PLUSS WATER:
sumTFASW <- rowSums(polh[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER)])

# PLUSS MOSER:
sumTFASWM <- rowSums(polh[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER, MOSS)])

# PLUSS KULL:
sumTFASWMC <- rowSums(polh[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, FAS, WATER, MOSS, CHC)])

Newpolh <- cbind(polh[,c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN)]/sumT,
                 polh[,c(FAS)]/sumTFAS,
                 polh[,c(WATER)]/sumTFASW,
                 polh[,c(CHC)]/sumTFASWMC,
                 polh[,c(MOSS), drop = FALSE]/sumTFASWM) * 100

wanted <- c("Alnus", "Betula", "Corylus", "Picea", "Pinus", "Ulmus", "Salix", "Juniperus", "Betula.nana", 
            "Vaccinium.sp.", "Cyperaceae", "Poaceae", "Dryopteris.sp.", "Equisetum", 
             "Gymnocarpium.dryopteris", "Sphagnum")

Newpolha <- Newpolh[, wanted]
head(Newpolha)
colnames(Newpolha)
Newpolha
#write.csv(Newpolh, "pollen_prosen_heim.csv")

LOI <- polh$LOI
setdiff(names(polh), c(TREES, SHRUBS, DSHRUBS, HERBS, UNKNOWN, WATER, FAS, MOSS, CHC, MARKER, STOMATA, LOI))
setdiff(wanted, names(Newpolha)) #Sjekke om noen av navnene/taxaene ikke er med i noen gruppe

#----------------------------------------------------#
# Preparing the influx data

polinf1 <- read.csv("influx_rate_heim_feb_u191.csv", header = TRUE, sep = ",", dec = ".")
head(polinf1)
polinf1[is.na(polinf1)] <- 0 # turn no data values into 0s
plot(age, polinf1$Alnus)
wanted <- c("Alnus", "Betula", "Corylus", "Picea", "Pinus", "Ulmus", "Salix", "Juniperus", "Betula.nana", 
            "Vaccinium.sp.", "Cyperaceae", "Poaceae", "Dryopteris.sp.", "Equisetum", 
            "Gymnocarpium.dryopteris", "Sphagnum", "Total")

Newpolinf1 <- polinf1[, wanted]/1000
Newpolinf1 <- rename(Newpolinf1, `Gymnocarpium` = `Gymnocarpium.dryopteris`,
                            `Betula nana` = `Betula.nana`, `Vaccinium` = `Vaccinium.sp.`, `Dryopteris` = `Dryopteris.sp.`)
setdiff(wanted, names(Newpolinf1))
#----------------------------------------------------#
#----------------------------------------------------#
#Rarify richness
polhR.df<- polh[,c(TREES, HERBS, SHRUBS, DSHRUBS, UNKNOWN)]
countTotals <- rowSums(polhR.df) # You counted 2000 grains in one sample! Woah! Is that correct?
countTotals
rarefy.n <- 432

#rarefy.n <- sample(countTotals) # maybe you should standardise this between the two sites (e.g. rarify to the lowest count sum across both sites?)

heim_richness <- rarefy(polhR.df, rarefy.n) 
heim_richness
plot(age, heim_richness, type ="l")
#write.csv(heim_richness, "Richness_Heim.csv")

#------------------------------------------#
# Preparing the percentage pollen data for the plotting and the stats
Newpolh1 <- rename(Newpolha, `Gymnocarpium` = `Gymnocarpium.dryopteris`, `Betula nana` = `Betula.nana`, 
                   `Vaccinium` = `Vaccinium.sp.`, `Dryopteris` = `Dryopteris.sp.`)

# Percentages
pollenh <- polh[,c(TREES, SHRUBS, DSHRUBS, HERBS)]
pollenh2 <- rename(pollenh, `Vaccinium` = `Vaccinium.sp.`, `Rumex` = `Rumex.sp.`, `Rubus` = `Rubus.chamaemorus`, `Betula nana` = `Betula.nana`)
pollenh_perc <- pollenh2/rowSums(pollenh2) *100
pollenh_perc <- pollenh_perc[, colSums(pollenh_perc>0)>2]

# Statistical analysis
#1: Cluster og soner:
ma.dist <- vegdist(sqrt(pollenh_perc), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")
plot(ma.chclust)
b1 <- bstick(ma.chclust, 10)

cbind(cutree(ma.chclust, k = 4), age)

# 2. ordination # QUESTION- why not do ordination on Newpolh? Because they just the "wanted" - alchemilla, r. chamaemorus etc might give information about climate gradients etc. --> Spør Aage og Anne. (Altså, Alistair lurer på hvorfor jeg ikke bare kjører ordinasjon på "wanted".)

# DCA
pollenh.dca <- decorana(sqrt(pollenh_perc)) # Note- changed percentage data to 
scor_h1 <- scores(pollenh.dca, choices = 1:2, display = "sites")
summary(pollenh.dca)
pollenh.dca
plot(pollenh.dca)

# PCA 
pollenh.pca <- rda(sqrt(pollenh_perc))
scor_h2 <- scores(pollenh.pca, choices = 1:2, display = "sites")

#RICHNESS
scor_h3 <- scores(heim_richness, choices = 1:2, display = "sites")

# Sorting of data ready for plotting
var_h1 <- as.matrix(cbind(LOI, scor_h2, scor_h3))
colnames(var_h1) <- c("LOI", "PCA1", "PCA2", "Species richness")
var_h2 <- as.data.frame(cbind(Newpolh[,CHC], polh[,STOMATA]))         
var_h2 <- var_h2 %>%
  mutate(Charcoal = Charcoal.big + Charcoal.small) %>% 
  select(-Stomata.Unknown, -Charcoal.big, -Charcoal.small)

#----------------------------------------------------#
####### PLOTTING THE PERCENTAGE DATA
cc <- c(rep("chartreuse1", 6), rep("blue", 2), rep("cyan", 1), rep("mediumorchid1", 3), rep("brown1", 3),
        rep("gold", 1)) #Velger farger
par(omi = c(0.3, 0.5, 0.6, 0.6)) #Fikse marginen til venstre. http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/#comment-173
par(mgp = c(0, 2, 1))

#med cluster
polh.plot <- strat.plot(Newpolh1,
                        yvar = age, y.rev = TRUE, scale.percent=TRUE, plot.poly=TRUE, col.poly=cc,
                        plot.bar=TRUE, col.bar ="grey69", col.poly.line="grey40", exag=TRUE, col.exag="auto",
                        exag.alpha = 0.4, exag.mult=10, x.pc.inc = 5, xSpace=0.01, x.pc.lab=TRUE,
                        x.pc.omit0=FALSE, srt.xlabel=60, ylabel=NA, cex.ylabel=1.05,
                        cex.yaxis=1, cex.xlabel=1, xRight = 0.65, xLeft = 0.06, yTop=0.805, yBottom=0.01,
                        y.tks= c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500
                                 , 7000, 7500, 8000, 8500, 9000, 9500, 10000), 
                        title="", graph.widths = 1)

# ADDING STOMATA/CHARCOAL TO THE PLOT
char.plot <- strat.plot(rename(var_h2, `Stomata Pinus` = `Stomata.Pinus`), xRight = 0.85, xLeft = 0.65, yTop=0.805, 
                         yBottom=0.01, 
                         yvar = age, y.rev=TRUE, plot.poly = FALSE,  plot.line = FALSE,  
                         y.axis = FALSE, plot.bar = TRUE, lwd.bar = 1.8, col.bar = "firebrick3", 
                         plot.symb = FALSE, srt.xlabel = 60, cex.xlabel = 1, 
                         xSpace = 0.01, add = TRUE, graph.widths = 0.4) 

# ADDING LOI/PCA/RICHNESS TO THE PLOT
pc.plot <- strat.plot(var_h1, xRight = 1, xLeft = 0.832, yTop=0.805, 
                       yBottom=0.01,
                       yvar = age, y.rev= TRUE, plot.poly = FALSE, plot.line = TRUE, col.line = "forestgreen", lwd.line = 1.8,
                       y.axis = FALSE,
                       plot.bar = FALSE, plot.symb = FALSE, srt.xlabel = 60, 
                       cex.xlabel = 1, xSpace = 0.01, add = TRUE, graph.widths = 0.4)

secondary.scale(yvar = rev(age), yvar2 = rev(depth), n = 6, y.rev = TRUE, xLeft = 0.02, yBottom = 0.01, ylabel2 = NA, cex.ylabel2 = 1)

# mtext("Depth (m)", side=3, line=-8, adj = 0.2)
# text(x=0.650, y=-100, labels = "Depth (m)", srt=60)

# Add cluster zones
addClustZone(polh.plot, ma.chclust, nZone=4, lwd=1.2, lty=2, col="grey25")
addClustZone(pc.plot, ma.chclust, nZone=4, lwd=1.2, lty=2, col="grey25")
addClustZone(char.plot, ma.chclust, nZone=4, lwd=1.2, lty=2, col="grey25")

text(x=0.61, y=1080, labels="Age (cal. yrs BP)", srt=60, xpd=NA)
text(x=0.58, y=1320, labels="Depth (cm)", srt=60, xpd=NA)
text(x=0.62, y=1220, labels="Zone names", srt=60, xpd=NA)
text(x=0.91, y=25, labels="Heimfjellsmyren 1097 m a.s.l.", cex=1.5, xpd=NA)


# text(x=0.6432, y=8800, labels = "HM-1", xpd=NA, cex=0.9)
# text(x=0.643, y=8000, labels = "HM-2.1", xpd=NA, cex=0.9)
# text(x=0.643, y=6000, labels = "HM-2.2", xpd=NA, cex=0.9)
# text(x=0.643, y=4000, labels = "HM-2.3", xpd=NA, cex=0.9)
# text(x=0.6432, y=3200, labels = "HM-3", xpd=NA, cex=0.9)
# text(x=0.6432, y=2500, labels = "HM-4", xpd=NA, cex=0.9)
# 
# text(x=0.6425, y=9120, labels = "Silt, \ndark \n(3+)", xpd=NA, cex=0.8)
# text(x=0.6425, y=6700, labels = "Detritus, \n dark \n(2)", xpd=NA, cex=0.8)
# text(x=0.6425, y=4700, labels = "Limus\n & \nsilt, \ndark \n(3)", xpd=NA, cex=0.8)
# text(x=0.6425, y=2000, labels = "Peat", xpd=NA, cex=0.8)
#----------------------------------------------------#
# TOTAL DIAGRAM:

TREES <- c("Betula", "Alnus", "Carpinus", "Corylus", "Fagus", "Picea", "Pinus", "Populus", "Quercus", "Sorbus", "Tilia", 
           "Ulmus")
HERBS <- c("Achillea.t.", "Alchemilla", "Apiaceae", "Arenaria", "Artemisia", "Asteraceae", "Beckwithia.glacialis",
           "Caryophyllaceae", "Chenopodiaceae", "Cirsium.sp.", "Compositae.cich.", "Cyperaceae", "Drosera", "Epilobium", 
           "Filipendula", "Geranium", "Humulus", "Jasione.montana", "Linnaea.borealis", "Lotus.sp.", "Melampyrum",
           "Onagraceae", "Oxyria.sp.", "Plantago.media", "Poaceae", "Potentilla.sp.", "Ranunculus.acris",
           "Ranunculus.flammula", "Ranunculus.sp.", "Rosaceae.undiff.", "Rubiaceae", "Rubus.chamaemorus",
           "Rumex.acetosa", "Rumex.sp.", "Saussurea.alpina", "Sedum",
           "Solidago.sp.", "Stachys.sp.", "Thalictrum", "Trientalis")
SHRUBS <- c("Hippophae", "Juniperus", "Salix")
DSHRUBS <- c("Calluna", "Empetrum", "Ericaceae.undiff.", "Vaccinium.sp.", "Betula.nana")

# POLLEN PERCENTAGES LIFEFORM GROUPS

trees_p <- rowSums(pollen_counts[,c(TREES)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

shrubs_p <- rowSums(pollen_counts[,c(SHRUBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

dshrubs_p <- rowSums(pollen_counts[,c(DSHRUBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

herbs_p <- rowSums(pollen_counts[,c(HERBS)])/rowSums(pollen_counts[, c(TREES, SHRUBS, DSHRUBS, HERBS) ]) *100

# TOTAL DIAGRAM THREES;SHRUBS; HERBS; DSHRUBS
plot(heim$Age, trees_p,  type = "n", ylim = c(0, 100), xlim = c(-67, 9700), ylab = "Groups pollen %", xlab = "Age (cal. yrs BP)", axes = F, cex.lab = 1.4)
par(mgp = c(0, 1, 0))
par(omi = c(0, 0.2, 0, 0))
axis(1, at = seq(0, 10000, 500))
axis(2, at = seq(0, 100, 10), las = 1, tcl = -0.2)
xx <- c(heim$Age, rev(heim$Age))
yy1 <- c(rep(0, 29), rev(trees_p))
polygon(xx, yy1, col = "chartreuse1")

#Hvis man vil fjerne sorte linjer -> border = NA inni hvert polygon

yy2 <- c(trees_p, rev(trees_p) + rev(shrubs_p))
polygon(xx, yy2, col = "blue")

yy3 <- c(trees_p + shrubs_p, rev(trees_p)+ rev(shrubs_p)+ rev(dshrubs_p))
polygon(xx, yy3, col = "cyan")

yy4 <- c(trees_p + shrubs_p + dshrubs_p, rev(trees_p) + rev(shrubs_p) + rev(dshrubs_p) + rev(herbs_p))
polygon(xx, yy4, col = "mediumorchid1")

legend(x = 150, y = 20, c("TREES","SHRUBS", "DSHRUBS", "HERBS"), fill= c("chartreuse1", "blue", "cyan", "mediumorchid1"), cex = 0.8)
box()
mtext("Heimfjellsmyren, 1097 m a.s.l.", side = 3, cex = 1.7)

# TOTAL POLLEN DIAGRAM; TRANSPOSED
par(mgp = c(3, 0.5, 0))
par(omi = c(.2, .2, 0, 0))
xx <- c(heim$Age, rev(heim$Age))
yy1 <- c(rep(0, 29), rev(trees_p))

plot(trees_p, heim$Age, type ="n", ylim = c(9637, -67), xlim = c(0, 100), xlab = "Groups pollen %", axes = FALSE, ylab = "Age (cal. yrs BP)", cex.lab = 1.5)
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
legend(x = 0.7, y = 1, c("TREES","SHRUBS", "DSHRUBS", "HERBS"), fill= c("chartreuse1", "blue", "cyan", "mediumorchid1"), cex = 0.8)
box()
mtext("Heimfjellsmyren, 1097 m a.s.l.", side = 3, cex = 1.7)
abline(h = c(1400, 1800, 8900), lty=4, lwd=1.8, col="black")
abline(h = c(4000, 6500), lty=3, lwd=1.8, col="grey47")
text(x=-2, y=9200, labels = "HM-1", xpd=NA, cex=0.9)
text(x=-2, y=7750, labels = "HM-\n2.1", xpd=NA, cex=0.9)
text(x=-2, y=5250, labels = "HM-\n2.2", xpd=NA, cex=0.9)
text(x=-2, y=2750, labels = "HM-\n2.3", xpd=NA, cex=0.9)
text(x=-2, y=1600, labels = "HM-3", xpd=NA, cex=0.9)
text(x=-2, y=500, labels = "HM-4", xpd=NA, cex=0.9)
#----------------------------------------------------#

# PCA plot
dev.off()
par(mfrow=c(1,1), mar= c(5,4,4,2))

loadPC1 = order(abs(summary(pollenh.pca)$species[,1]), decreasing=T)[1:9]
loadPC2 = order(abs(summary(pollenh.pca)$species[,2]), decreasing=T)[1:9]

plot(summary(pollenh.pca)$sites[,1], summary(pollenh.pca)$sites[,2], cex=1, xlab="PC1", 
     ylab="PC2", type="n", main= "PCA Heimfjellsmyren (1097 m a.s.l.)", xlim = c(-2.1,2.1), ylim=c(-2,1.7))

abline(h=0, v=0, lty=2)

#points(summary(pollenh.pca)$sites[,1], summary(pollenh.pca)$sites[,2], cex=1, bg="dark grey", pch=21)

specPC1 = summary(pollenh.pca)$sites[,1]
specPC2 = summary(pollenh.pca)$sites[,2]
spec = c(loadPC1, loadPC2)

text(x=c(specPC1), y=c(specPC2), labels=names(specPC1), cex=1)

specPC3 = summary(pollenh.pca)$species[,1]
specPC4 = summary(pollenh.pca)$species[,2]

arrows(x0=rep(0,24), y0=rep(0,24), x1=c(specPC3[spec]), y1=c(specPC4[spec]), col="dark red", length=0.07)
text(x=c(specPC3[spec]), y=c(specPC4[spec]), labels=names(specPC3[spec]), cex=1, col="dark red", font=3)

# alternative
#plot(pollenh.pca) # see here https://www.fromthebottomoftheheap.net/2013/01/13/decluttering-ordination-plots-in-vegan-part-2-orditorp/

#----------------------------------------------------#
#Plotting the influx diagram
cc <- c(rep("chartreuse1", 6), rep("blue", 2), rep("cyan", 2), rep("mediumorchid1", 2), rep("brown1", 3),
        rep("gold", 1), rep("grey40", 1)) #Velger farger
par(omi = c(0.3, 0.3, 0.6, 0.6)) #Fikse marginen til venstre. http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/#comment-173
par(mgp = c(0, 2, 1))

polinf1.plot <- strat.plot(Newpolinf1,
                           yvar = age, y.rev=TRUE, scale.percent=FALSE, scale.minmax = TRUE, plot.poly=TRUE, 
                           col.poly=cc, plot.bar= TRUE, col.bar="grey50", col.poly.line="grey69", exag=TRUE, 
                           col.exag="auto", exag.alpha = 0.4, exag.mult=10, x.pc.inc = 5000, x.pc.lab=TRUE, 
                           x.pc.omit0=FALSE, srt.xlabel=60, ylabel="",
                           cex.ylabel=1, cex.yaxis=1, cex.axis = 0.6, cex.xlabel=1.2,
                           title="", xLeft = 0.09, 
                           yBottom = 0.05, y.tks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000,
                                                     4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500,
                                                     10000), graph.widths=1)

poinf1.plot <- secondary.scale(yvar = rev(age), yvar2 = rev(depth), n = 6, y.rev = TRUE,  
                               cex.ylabel2 = 1, xLeft = 0.05, yBottom = 0.05, yTop = 0.795)  

text(x=0.63, y=1080, labels="Age (cal. yrs BP)", srt=60, xpd=NA)
text(x=0.60, y=1320, labels="Depth (cm)", srt=60, xpd=NA)
text(x=0.645, y=1200, labels="Zone names", srt=60, xpd=NA)
text(x=0.90, y=25, labels="Heimfjellsmyren 1097 m a.s.l.", cex=1.5, xpd=NA)

addClustZone(polinf1.plot, ma.chclust, nZone=4, lwd=1.2, lty=2, col="grey25")

# text(x=0.6432, y=8800, labels = "HM-1", xpd=NA, cex=0.9)
# text(x=0.643, y=8000, labels = "HM-2.1", xpd=NA, cex=0.9)
# text(x=0.643, y=6000, labels = "HM-2.2", xpd=NA, cex=0.9)
# text(x=0.643, y=4000, labels = "HM-2.3", xpd=NA, cex=0.9)
# text(x=0.6432, y=3200, labels = "HM-3", xpd=NA, cex=0.9)
# text(x=0.6432, y=2500, labels = "HM-4", xpd=NA, cex=0.9)
# 
# text(x=0.6425, y=9120, labels = "Silt, \ndark \n(3+)", xpd=NA, cex=0.8)
# text(x=0.6425, y=6700, labels = "Detritus, \n dark \n(2)", xpd=NA, cex=0.8)
# text(x=0.6425, y=4700, labels = "Limus\n & \nsilt, \ndark \n(3)", xpd=NA, cex=0.8)
# text(x=0.6425, y=2000, labels = "Peat", xpd=NA, cex=0.8)
#----------------------------------------------------#
#Rarify richness
polhR.df<- polh[,c(TREES, HERBS, SHRUBS, DSHRUBS, UNKNOWN)]
countTotals <- rowSums(polhR.df) # You counted 2000 grains in one sample! Woah! Is that correct?
countTotals
rarefy.n <- min(countTotals) # maybe you should standardise this between the two sites (e.g. rarify to the lowest count sum across both sites?)

heim_richness <- rarefy(polhR.df, rarefy.n) 
heim_richness
plot(age, heim_richness, type ="l")
#write.csv(heim_richness, "Richness_Heim.csv")

