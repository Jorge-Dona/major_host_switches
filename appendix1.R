
########### Host-specificity, infrequent major host-switching and the diversification of highly host-specific symbionts: the case of feather mites ##############################
# Jorge Doña 12/09/2017


require(ape)
require(Hmisc)
require(rworldmap)
require(RColorBrewer)
require(vegan)
require(betapart)
require(dplyr)
require(caper)
require(diversitree)
require(picante)
require(stringi)

########## Map of the dataset, Figure 1 ###########################################################################

# Database
data<- read.csv("global_fmbird.csv", sep="\t")

# Prunning to only high quality data.
data<- droplevels(subset(data, data_quality == 2))

pivot<- as.data.frame(table(data$country_collected))

#Plotting into a map

sPDF <- joinCountryData2Map(pivot,
                            joinCode = "NAME",
                            nameJoinColumn = "Var1",
                            verbose = TRUE)

pdf("Figure1.pdf", width= 18.5,height= 12)

mapParams<-mapCountryData(sPDF, numCats =13,nameColumnToPlot="Freq",catMethod="logFixedWidth", colourPalette="topo",addLegend=FALSE, mapTitle="")

#Modifying the legend

do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=1
                          , legendIntervals = "page"
                          ,digits =0.8,labelFontSize=2,legendMar=10.1
))

dev.off()

png("Figure1.png", width= 18,height= 12, units= "in",res=300)

mapParams<-mapCountryData(sPDF, numCats =13,nameColumnToPlot="Freq",catMethod="logFixedWidth", colourPalette="topo",addLegend=FALSE, mapTitle="")

#Modifying the legend

do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=1
                          , legendIntervals = "page"
                          ,digits =0.8,labelFontSize=2,legendMar=10.1
))

dev.off()

#### Host range ####

#data
data<- read.csv("global_fmbird.csv", sep="\t")

# Quality trimming
dataHN<- droplevels(subset(data, data_quality == 2))


# Removing author column
dataHN<- dataHN[,-6]


# Adding a bird genus column

inter<-with(dataHN, table (dataHN$mite_species,dataHN$host_species))
Bnames=strsplit(as.character(dataHN$host_species), split=" ")
hostgenus<-NULL
#bird genus names
for (i in 1:length((Bnames)))
{
  hostgenus[i]<-sapply(Bnames[i], function(x) paste((x[1])))
}

#Adding the new column to matrix
dataHN$host_genus <- hostgenus



##Building matrix to stablish the plot order (i.e. mite_species vs bird_species, mite_species vs host_species, etc.)
mforloop_birds= matrix(c(5,5,5,5,4,4,4,4,3,3,3,3,2,2,2,2,8,15,7,6,8,15,7,6,8,15,7,6,8,15,7,6),nrow=16,ncol=2) #The order of interaction to be plotted for dataHN dataset. For example: The first one is 5(column five from dataHN, "mite_species" and 8 column "host_species").
mforloop_mites= matrix(c(8,8,8,8,15,15,15,15,7,7,7,7,6,6,6,6,5,4,3,2,5,4,3,2,5,4,3,2,5,4,3,2),nrow=16,ncol=2) 


##This function creates a matrix with zeros and ones and then sum rows and columns to obtain the degree (i.e. the numer of link for each node).
##mitecolumn and hostcolumn comes from mforloop_birds (birds), mforloop_mites (mites) 


#This is to know the number of mite species having X hosts. It is necesary to run manually plotting_degrees  with these parameters from below:
# then it is necessary to manually run scaling barplots.

# mitecolumn = 5
# hostcolumn = 8

# number of mites on different families

#mitecolumn = 5
#hostcolumn = 7

plotting_degrees_mite <- function(mitecolumn, hostcolumn) {
  {
    mitesp_birdsp <-with(dataHN, table (dataHN[,mitecolumn], dataHN[,hostcolumn]))#Create a interaction matrix from the two columns selected, the number will represent the number of times that that interaction is repeated
    mitesp_birdsp <-ifelse(mitesp_birdsp>0,1,0) #Transform all values above 1 in 1, that is, a presence-absence matrix not weighted.
    matrix= as.data.frame(mitesp_birdsp)
    mites_degree <- rowSums(matrix) # Bird number of species where a mite is present
    mites_degree= sort(mites_degree,decreasing=TRUE)
    birds_degree <- colSums(matrix) # Mite number of species of a single bird species
    birds_degree <- sort(birds_degree,decreasing=TRUE)#Sorting data for plot
  }
  result2<<-mites_degree
}


### This function adjust data to power of two intervals. See Jovani 2006 for the rationale behind.

# To estimate the percentage of mite species with X hosts.

#object= result2
# 
# These are the percentages extracted from each letter  for example for mite species in 3 hosts (see below)
# 494/1876*100= 26%
# 334/1876*100= 18%
# 160/1876*100 = 9%
# 
# 57+26= 83%

scaling_barplots <- function(object) {
  {
    c= as.data.frame(object)
    a=nrow(subset(c, object == 1)) #Subset of all the values equal to one (from plotting_degrees_mite output)
    # a=nrow(subset(c, object == 3)) # For example, to estimate the number of mite specie with 3 hosts
    # a=nrow(subset(c, object == 2))
    #a=nrow(subset(c, object >1)) # To know the mites in more than one family
    
    b<- nrow(subset(c, object >= 2 & object <= 3)) #Subset of all the values among 2 and 3 (from plotting_degrees_mite output), the following, the same.
    c1=nrow(subset(c, object >= 4 & object <= 7))
    d1=nrow(subset(c, object >= 8 & object <= 15))
    e=nrow(subset(c, object >= 16 & object <= 31))
    f=nrow(subset(c, object >= 32 & object <= 63))
    g=nrow(subset(c, object >= 64 & object <= 127))
    h=nrow(subset(c, object >= 128 & object <= 255))
    i=nrow(subset(c, object >= 256 & object <= 511))
    j=nrow(subset(c, object >= 512 & object <= 1023))
    k=nrow(subset(c, object >= 1024 & object <= 2047))
    z=rbind(a,b,c1,d1,e,f,g,h,i,j,k) # A vector with all the values for the power of two interval
    maxz=max(z)*1.1 # To include space in relation with the maximum in the y-axis
    minz=0-(maxz-max(z)) # To include the same space but in the lower part
    z=as.vector(z)
    z[z == 0] <- NA  #To not plot the 0 value
  }
  m=barplot(z,xlim=c(-0.25,13.75), ylim=c(minz,maxz), axisnames=FALSE, col="black",lwd=2) #Final barplot, it requires the posterior foor loop for axes and graphical anotations
}




##Graphical parameters 

coordenates=1:11 #Setting coordenates to x axes lab; i.e. number of categories to legend
mids <- barplot(coordenates, xlab="",plot=FALSE)#Setting coordenates to x axes lab

#### Host range (Figure 2) #####

pdf("Figure2.pdf", width= 18.5,height= 12)

par(mfrow=c(4,4)) #vertical way; 4 columns by 4 rows
par(xaxs = "i",yaxs="i",par(pin=c(1.4, 1.4), mai=c(1.02,0.42,0.32,0.42)))#Modifying plot space parameters


for (i in (1:4))#The number of plots
{
  a=mforloop_birds[i,1] #The position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  b=mforloop_birds[i,2] #The second position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  plotting_degrees_mite(a,b) #The input for plotting_degrees_mite
  scaling_barplots(result2) #The input (coming from plotting_degrees_mite) for scaling_barplots
  #if(i>= 13) {#To only plot the x-axes labels in the las row of plots
  
  axis(1, at=mids,labels=c("1","2-3","4-7", "8-15","16-31","32-63","64-127","128-255", "256-511","512-1023","1024-2047"),cex.axis=1,las=2)
  box(,lwd=2)
  #}
  # else{
  #   axis(1, at=mids,cex.axis=1,labels=FALSE)
  #   box(,lwd=2)
}


dev.off()

# Figure S1

#### Host range (Figure S1)

pdf("FigureS1.pdf", width= 18.5,height= 12)

par(mfrow=c(4,4)) #vertical way; 4 columns by 4 rows
par(xaxs = "i",yaxs="i",par(pin=c(1.4, 1.4), mai=c(1.02,0.42,0.32,0.42)))#Modifying plot space parameters

# Proctophyllodidae

data<- read.csv("global_fmbird.csv", sep="\t")

# Data
dataHN<- droplevels(subset(data, data_quality == 2))

# Removing author column
dataHN<- dataHN[,-6]

# Adding bird genus column

inter<-with(dataHN, table (dataHN$mite_species,dataHN$host_species))
Bnames=strsplit(as.character(dataHN$host_species), split=" ")
hostgenus<-NULL
#bird genus names
for (i in 1:length((Bnames)))
{
  hostgenus[i]<-sapply(Bnames[i], function(x) paste((x[1])))
}

#Adding the new column to matrix
dataHN$host_genus <- hostgenus


dataHN <-droplevels(filter(dataHN, mite_family == "Proctophyllodidae")) 

for (i in (1:4))#The number of plots
{
  a=mforloop_birds[i,1] #The position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  b=mforloop_birds[i,2] #The second position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  plotting_degrees_mite(a,b) #The input for plotting_degrees_mite
  scaling_barplots(result2) #The input (coming from plotting_degrees_mite) for scaling_barplots
  #if(i>= 13) {#To only plot the x-axes labels in the las row of plots
  
  axis(1, at=mids,labels=c("1","2-3","4-7", "8-15","16-31","32-63","64-127","128-255", "256-511","512-1023","1024-2047"),cex.axis=1,las=2)
  box(,lwd=2)
  #}
  # else{
  #   axis(1, at=mids,cex.axis=1,labels=FALSE)
  #   box(,lwd=2)
}


# Trouessartiidae
# Data
data<- read.csv("global_fmbird.csv", sep="\t")

#
dataHN<- droplevels(subset(data, data_quality == 2))

# Removing author column
dataHN<- dataHN[,-6]


# Adding bird genus column

inter<-with(dataHN, table (dataHN$mite_species,dataHN$host_species))
Bnames=strsplit(as.character(dataHN$host_species), split=" ")
hostgenus<-NULL
#bird genus names
for (i in 1:length((Bnames)))
{
  hostgenus[i]<-sapply(Bnames[i], function(x) paste((x[1])))
}

# Adding the new column to matrix
dataHN$host_genus <- hostgenus

dataHN <-droplevels(filter(dataHN, mite_family == "Trouessartiidae")) 
for (i in (1:4))#The number of plots
{
  a=mforloop_birds[i,1] #The position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  b=mforloop_birds[i,2] #The second position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  plotting_degrees_mite(a,b) #The input for plotting_degrees_mite
  scaling_barplots(result2) #The input (coming from plotting_degrees_mite) for scaling_barplots
  #if(i>= 13) {#To only plot the x-axes labels in the las row of plots
  
  axis(1, at=mids,labels=c("1","2-3","4-7", "8-15","16-31","32-63","64-127","128-255", "256-511","512-1023","1024-2047"),cex.axis=1,las=2)
  box(,lwd=2)
  #}
  # else{
  #   axis(1, at=mids,cex.axis=1,labels=FALSE)
  #   box(,lwd=2)
}

# Pteronyssidae

# Data
data<- read.csv("global_fmbird.csv", sep="\t")

# Quality trimming
dataHN<- droplevels(subset(data, data_quality == 2))

# Removing author column
dataHN<- dataHN[,-6]


# Adding bird genus column

inter<-with(dataHN, table (dataHN$mite_species,dataHN$host_species))
Bnames=strsplit(as.character(dataHN$host_species), split=" ")
hostgenus<-NULL
#bird genus names
for (i in 1:length((Bnames)))
{
  hostgenus[i]<-sapply(Bnames[i], function(x) paste((x[1])))
}

# Adding the new column to matrix
dataHN$host_genus <- hostgenus

dataHN <-droplevels(filter(dataHN, mite_family == "Pteronyssidae")) 
for (i in (1:4))#The number of plots
{
  a=mforloop_birds[i,1] #The position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  b=mforloop_birds[i,2] #The second position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  plotting_degrees_mite(a,b) #The input for plotting_degrees_mite
  scaling_barplots(result2) #The input (coming from plotting_degrees_mite) for scaling_barplots
  #if(i>= 13) {#To only plot the x-axes labels in the las row of plots
  
  axis(1, at=mids,labels=c("1","2-3","4-7", "8-15","16-31","32-63","64-127","128-255", "256-511","512-1023","1024-2047"),cex.axis=1,las=2)
  box(,lwd=2)
  #}
  # else{
  #   axis(1, at=mids,cex.axis=1,labels=FALSE)
  #   box(,lwd=2)
}


# Alloptidae

# Data
data<- read.csv("global_fmbird.csv", sep="\t")

# Quality trimming
dataHN<- droplevels(subset(data, data_quality == 2))

# Removing author column
dataHN<- dataHN[,-6]


# Adding bird genus column

inter<-with(dataHN, table (dataHN$mite_species,dataHN$host_species))
Bnames=strsplit(as.character(dataHN$host_species), split=" ")
hostgenus<-NULL
#bird genus names
for (i in 1:length((Bnames)))
{
  hostgenus[i]<-sapply(Bnames[i], function(x) paste((x[1])))
}

#Adding the new column to matrix
dataHN$host_genus <- hostgenus


dataHN <-droplevels(filter(dataHN, mite_family == "Alloptidae")) # Prueba solo proctophyllodidae
for (i in (1:4))#The number of plots
{
  a=mforloop_birds[i,1] #The position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  b=mforloop_birds[i,2] #The second position to be selected in the mforloop_birds matrix and the input for plotting_degrees_mite
  plotting_degrees_mite(a,b) #The input for plotting_degrees_mite
  scaling_barplots(result2) #The input (coming from plotting_degrees_mite) for scaling_barplots
  #if(i>= 13) {#To only plot the x-axes labels in the las row of plots
  
  axis(1, at=mids,labels=c("1","2-3","4-7", "8-15","16-31","32-63","64-127","128-255", "256-511","512-1023","1024-2047"),cex.axis=1,las=2)
  box(,lwd=2)
  #}
  # else{
  #   axis(1, at=mids,cex.axis=1,labels=FALSE)
  #   box(,lwd=2)
}
dev.off()

##### Estimating phylgenetic host-specificity ####

# Database
data<- read.csv("global_fmbird.csv", sep="\t")

# Prunning to only high quality data.
data<- droplevels(subset(data, data_quality == 2))

# Enabling the correspondence of bird taxonomy in the Doña et al. (2016) database to that of Jetz et al. (2014). Requires Table S1 ()
data2<- data
data1<- read.csv("Table_S1.csv", h=TRUE)
old=as.vector(data1$IOC_5.4)
n= as.vector(data1$Birdtree)
a= data2$host_species
for (i in 1:length(old))
{
  a=gsub(old[i],n[i], a)
}
data2$host_species = a
data<- as.data.frame( data2)

# Adding underscore in bird species names. It is require to work with the bird phylogenetic distances matrix.
Bnames= strsplit(as.vector(data$host_species), split=" ")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
birdnames2<- sapply(Bnames, function(x) paste(x[2]))
data[["host_species"]] <- as.factor(paste(birdnames, birdnames2, sep="_"))


# To create a matrix with a single entry per bird-mite association
data[["concatenate"]] <- as.factor(paste(data$mite_species, data$host_species, sep="_"))
data<- droplevels(data[!duplicated(data$concatenate), ])


# Bird tree 
fmtree <- read.nexus("fm_bird_tree.tre")


#13) Preparing the tree
fmtree$node.label <- NULL # This is required to use comparative.data object to check the congruence among phylogeny and data
fmtree$edge.length[which(fmtree$edge.length <= 0)] <- 0.0000001 # It is a problem with zeroes with comparative.data. I find this solution in a Liam Revell post.


# Presence-absence matrices

mitesp_birdsp2 <-with(data, table (droplevels(data[,9]), droplevels(data[,5])))
mitesp_birdsp2 <-as.data.frame(ifelse(mitesp_birdsp2>0,1,0))
species <- as.vector(row.names(mitesp_birdsp2))
mitesp_birdsp2["species"]<- species


# Now selecting species with more than one host
a<- droplevels(mitesp_birdsp2[,1:1876])
z<-  colSums(a)
at<- z>1.5
a<-rbind(a,at)
def<- a[,a[2145,] == 1] 
def<- def[-2145,]
species <- as.vector(row.names(def))
def["species"]<- species


# 
fmcompsp <- comparative.data(fmtree, def, species) #Species
removez<- fmcompsp$dropped$tips # Object, after create the compdata object and to delete the taxa that are not in phylogeny or in data. In this case, these are the phylogeny tips not present in data.

# Re-loading the bird tree 
fmtree <- read.nexus("fm_bird_tree.tre")
fmtree1<- drop.tip(fmtree, removez) # Now, removing the tips not present in data
rmdata<-fmcompsp$dropped$unmatched.rows # Object,the same but with the data rows not present in phylogeny.
mitesp_birdsp<-def[ !(rownames(def) %in% rmdata), ] #

#
mitesp_birdsp= t(mitesp_birdsp)
SES=ses.pd(mitesp_birdsp, fmtree1, null.model="taxa.labels", runs=1000, iterations=1000) 


#NET PD

NET<- SES$pd.obs-SES$pd.rand.mean
SES_NET<- cbind(SES, net)

write.csv(SES_NET, "Table_S2.csv")


#Plot ses.pd versus nº of species

data= read.csv("Table_S2.csv")
data=data[-809,]

model=lm(data$ntaxa ~ data$pd.obs.z)
summary(model)
plot(data$ntaxa,data$pd.obs.z, xlab="Number of taxa", ylab="SES.PD")



###################### Geographic specificity ##########################################

### Geographic specificity (Continent approach) ####

#Database
data_last<- read.csv("global_fmbird.csv", sep="\t")

# Prunning associations without continent information and Antarctica (because it has only 1 association).
data_last<- subset(data_last, continent != "")
data_last<- subset(data_last, continent != "Antarctica")

# Prunning the database to retain only high quality data 
data_last<- subset(data_last, data_quality == 2)

# Creating supercontinents (i.e. EurAfrica = Europe + Africa, AsiaOce = Asia + Oceania, America = North America + South America)
a= data_last$continent
a = gsub("Africa","EurAfrica",a)
a = gsub("Europe","EurAfrica",a)
a = gsub("South America","America",a)
a = gsub("North America","America",a)
a = gsub("Asia","AsiaOce",a)
a = gsub("Oceania","AsiaOce",a)
data_last$supracontinent =as.factor(a)


# Continents:
# Creating a column "concatenate" with a concatenation of species and continent
data_last[["concatenate"]] <- as.factor(paste(data_last$host_species, data_last$continent, sep="_"))

# "matrix" is a matrix where rows are the unique combinations of bird-continent and columns unique species of mites
# Inside the matrix it is the number of registers of each bird-mite combination
matrix<- with(data_last, table (droplevels(data_last[,17]),droplevels (data_last[,5])))

# data is like matrix but with presence-absence data to use Jaccard
data <-as.data.frame(ifelse(matrix>0,1,0))

# continent is a matrix with as many rows as bird-continent combinations (as matrix)
# continent only has two columns, one for bird species and the other one for continent 
a<-row.names(droplevels(data))
Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continent<- sapply(Bnames, function(x) paste(x[2]))
continent<- as.data.frame(cbind(birdnames, continent))


#  Bird species in more than one continent
# Creating a matrix with only those species in more than one continent 
matrixduplicados <- cbind (data,continent)
# 1774 is the bird species
matrixduplicados <- subset(matrixduplicados,duplicated(matrixduplicados[,1774]) | duplicated(matrixduplicados[,1774], fromLast=TRUE))
matrixduplicados <- matrixduplicados[,1:1773]
# continent is a matrix with as many rows as bird-continent combinations (as matrix)
# continent only has two columns, one for bird species and the other one for continent 

a<-row.names(droplevels(matrixduplicados))
Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continent <- sapply(Bnames, function(x) paste(x[2]))
continent <- as.data.frame(cbind(birdnames, continent))



#### PERMANOVA:

birddel<-adonis(matrixduplicados ~ continent$birdnames + continent$continent, permutations=100, method="jaccard" )
matrizdistancia= vegdist(matrixduplicados, method="jaccard")

# betadisper
beta <- betadisper(matrizdistancia, continent$birdnames)
permutest(beta)

beta <- betadisper(matrizdistancia, continent$continent)
permutest(beta)



#Partitioning Jaccard dissimilarity into nestedness and replacement matrices following Baselga 2010, 2012.

cont.dist<-beta.pair(matrixduplicados, index.family="jac")

#beta.jtu
betaturn<- as.matrix(cont.dist$beta.jtu)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=100 )

#beta.jne

betaturn<- as.matrix(cont.dist$beta.jne)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=1000)



### Geographic specificity (Continent approach, high quality data) ######

# To retain only birds with at least three mites in at least two continents

matrixduplicados<- cbind (matrixduplicados, continent)
#number of mites per row (bird-continent combination)
matrixduplicados <- cbind (matrixduplicados, rowSums(matrixduplicados[,1:1773]))
# Removing rows with less than three mites 
matrixduplicados<- subset(matrixduplicados, matrixduplicados[,1776] > 1 )
# Removing birds with a single row (i.e., only one continent)
matrixduplicados <- subset(matrixduplicados, duplicated(matrixduplicados[,1774]) | duplicated(matrixduplicados[,1774], fromLast=TRUE))
# Removing columns with a single species of birds
matrixduplicados <- matrixduplicados[,which(!apply(matrixduplicados,2,FUN = function(x){all(x == 0)}))]

a<-row.names(droplevels(matrixduplicados))
Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
supracontinent<- sapply(Bnames, function(x) paste(x[2]))
supracontinent<- as.data.frame(cbind(birdnames, supracontinent))

birddel<-adonis(matrixduplicados[1:273] ~ supracontinent$birdnames + supracontinent$supracontinent, permutations=100,method="jaccard")

#Partitioning Jaccard dissimilarity into nestedness and replacement matrices following Baselga 2010, 2012.

cont.dist<-beta.pair(matrixduplicados[,1:273], index.family="jac")

#beta.jtu
betaturn<- as.matrix(cont.dist$beta.jtu)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=100 )

#beta.jne

betaturn<- as.matrix(cont.dist$beta.jne)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=1000)


########## Geographic specificity (Supercontinent approach) ###############

#Database
data_last<- read.csv("global_fmbird.csv", sep="\t")

# Prunning associations without continent information and Antarctica (because it has only 1 association).
data_last<- subset(data_last, continent != "")
data_last<- subset(data_last, continent != "Antarctica")

# Prunning the database to retain only high quality data 
data_last<- subset(data_last, data_quality == 2)

# Creating supercontinents (i.e. EurAfrica = Europe + Africa, AsiaOce = Asia + Oceania, America = North America + South America)
a= data_last$continent
a = gsub("Africa","EurAfrica",a)
a = gsub("Europe","EurAfrica",a)
a = gsub("South America","America",a)
a = gsub("North America","America",a)
a = gsub("Asia","AsiaOce",a)
a = gsub("Oceania","AsiaOce",a)
data_last$supracontinent =as.factor(a)

#Super continents:
# Creating a column "concatenate" with a concatenation of species and continent
data_last[["concatenate"]] <- as.factor(paste(data_last$host_species, data_last$supracontinent, sep="_"))

# "matrix" is a matrix where rows are the unique combinations of bird-supercontinent and columns unique species of mites
# Inside the matrix it is the number of registers of each bird-mite combination
matrix<- with(data_last, table (droplevels(data_last[,17]),droplevels (data_last[,5])))

# data is like matrix but with presence-absence data to use Jaccard
data <-as.data.frame(ifelse(matrix>0,1,0))

# birds_supracontinents is a matrix with as many rows as bird-supercontinent combinations (as matrix)
# birds_supracontinents only has two columns, one for bird species and the other one for continent 

a <-row.names(droplevels(data))
Bnames <- strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
supracontinents <- sapply(Bnames, function(x) paste(x[2]))
birds_supracontinents <- as.data.frame(cbind(birdnames,supracontinents))


#  Bird species in more than one supercontinent
# Creating a matrix with only those species in more than one supercontinent 
matrixduplicados <- cbind (data,birds_supracontinents)
# 1774 is the bird species
matrixduplicados <- subset(matrixduplicados,duplicated(matrixduplicados[,1774]) | duplicated(matrixduplicados[,1774], fromLast=TRUE))
matrixduplicados <- matrixduplicados[,1:1773]
# birds_supracontinents is a matrix with as many rows as bird-supercontinent combinations (as matrix)
# birds_supracontinents only has two columns, one for bird species and the other one for supercontinent 
a<-row.names(droplevels(matrixduplicados))

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
supracontinents<- sapply(Bnames, function(x) paste(x[2]))
birds_supracontinents<- as.data.frame(cbind(birdnames, supracontinents))


# PERMANOVA 

birddel<-adonis(matrixduplicados ~ birds_supracontinents$birdnames + birds_supracontinents$supracontinents, permutations=100,method="jaccard")

matrizdistancia= vegdist(matrixduplicados, method="jaccard")

# betadisper
beta <- betadisper(matrizdistancia, birds_supracontinents$birdnames)
permutest(beta)

beta <- betadisper(matrizdistancia, birds_supracontinents$supracontinents)
permutest(beta)


#Partitioning Jaccard dissimilarity into nestedness and replacement matrices following Baselga 2010, 2012.

cont.dist<-beta.pair(matrixduplicados, index.family="jac")


#beta.jtu
betaturn<- as.matrix(cont.dist$beta.jtu)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=100 )

#beta.jne

betaturn<- as.matrix(cont.dist$beta.jne)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=1000)


#### Geographic specificity (Super continent approach, Figure S3) ####

#Using all species
#Building the dendrogram
dist.mat<-vegdist(matrixduplicados,method="jaccard")
x=hclust(dist.mat)
phy<- as.phylo(x)

# Color vectors
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color2=  sample(color, 227)
color3=  sample(color, 4)

#Plotting
pdf("FigureS3.pdf", width= 18,height= 12)

plot(phy, type = "fan",tip.color= color2[birds_supracontinents$birdnames], label.offset=0.02,cex = 0.5 )
tiplabels(bg=color3[as.numeric(birds_supracontinents$supracontinent)],pch=21,  cex = 0.8)


dev.off()


png("FigureS3.png", width= 18,height= 12, units= "in",res=300)


dev.off()



### Geographic specificity (Supercontinent approach, high quality data) #####

# To retain only birds with at least three mites in at least two supercontinents

matrixduplicados<- cbind (matrixduplicados, birds_supracontinents)
#number of mites per row (bird-supercontinent combination)
matrixduplicados <- cbind (matrixduplicados, rowSums(matrixduplicados[,1:1773]))
# Removing rows with less than three mites 
matrixduplicados<- subset(matrixduplicados, matrixduplicados[,1776] > 1 )
# Removing birds with a single row (i.e., only one continent)
matrixduplicados <- subset(matrixduplicados, duplicated(matrixduplicados[,1774]) | duplicated(matrixduplicados[,1774], fromLast=TRUE))
# Removing columns with a single species of birds
matrixduplicados <- matrixduplicados[,which(!apply(matrixduplicados,2,FUN = function(x){all(x == 0)}))]

matrixduplicados<- matrixduplicados[, 1: 181]

a<-row.names(droplevels(matrixduplicados))
Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
supracontinent<- sapply(Bnames, function(x) paste(x[2]))
supracontinent<- as.data.frame(cbind(birdnames, supracontinent))


birddel<-adonis(matrixduplicados ~ supracontinent$birdnames + supracontinent$supracontinent, permutations=100,method="jaccard")

matrizdistancia= vegdist(matrixduplicados, method="jaccard")


#Partitioning Jaccard dissimilarity into nestedness and replacement matrices following Baselga 2010, 2012.

cont.dist<-beta.pair(matrixduplicados, index.family="jac")


#beta.jtu
betaturn<- as.matrix(cont.dist$beta.jtu)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=100 )

#beta.jne

betaturn<- as.matrix(cont.dist$beta.jne)
a<-row.names(betaturn)

Bnames= strsplit(as.vector(a), split="_")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
continents<- sapply(Bnames, function(x) paste(x[2]))
birds_continents<- as.data.frame(cbind(birdnames, continents))

birddel<-adonis(betaturn ~ birds_continents$birdnames + birds_continents$continents, permutations=1000)



### Geographic specificity (Supercontinent approach, high quality data, Figure 3) #####

#Subset of better sampled species
#Building the dendrogram
dist.mat<-vegdist(matrixduplicados,method="jaccard")
x=hclust(dist.mat)
phy<- as.phylo(x)

# Color vectors
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color2=  sample(color, 227)
color3=  sample(color, 4)

#Plotting
pdf("Figure3.pdf", width= 9,height= 18)

plot(phy, type = "fan",tip.color= color2[birds_supracontinents$birdnames], label.offset=0.02,cex = 0.5 )
tiplabels(bg=color3[as.numeric(birds_supracontinents$supracontinent)],pch=21,  cex = 0.8)

dev.off()


png("Figure3.png", width= 9,height= 18, units= "in",res=300)


dev.off()



############# Host switching ecological patterns, Figure 3 ####################

# Database
data<- read.csv("global_fmbird.csv", sep="\t")

# Prunning to only high quality data.
data<- droplevels(subset(data, data_quality == 2))

# Enabling the correspondence of bird taxonomy in the Doña et al. (2016) database to that of Jetz et al. (2014). Requires Table S1 ()
data2<- data
data1<- read.csv("Table_S1.csv", h=TRUE)
old=as.vector(data1$IOC_5.4)
n= as.vector(data1$Birdtree)
a= data2$host_species
for (i in 1:length(old))
{
  a=gsub(old[i],n[i], a)
}
data2$host_species = a
data<- as.data.frame( data2)

# Adding underscore in bird species names. It is require to work with the bird phylogenetic distances matrix.
Bnames= strsplit(as.vector(data$host_species), split=" ")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
birdnames2<- sapply(Bnames, function(x) paste(x[2]))
data[["host_species"]] <- as.factor(paste(birdnames, birdnames2, sep="_"))


# To create a matrix with a single entry per bird-mite association
data[["concatenate"]] <- as.factor(paste(data$mite_species, data$host_species, sep="_"))
data<- droplevels(data[!duplicated(data$concatenate), ])


# Bird tree 

fmtree <- read.nexus("fm_bird_tree.tre")


# Creating bird phylogenetic distance matrix

dist<-cophenetic.phylo(fmtree) 
dist<- as.data.frame(dist)


# Creating a matrix with the the number of columns required by the interval size (currently by 20, so ten rows is enough).

matrizdef<- matrix(0 , nrow = 10, ncol = 2)
matrizdef[,1] = seq(from = 10, to = 190, by = 20)

npairwise<-numeric()

for (i in levels(data$mite_species)){
  
  focalmite<- droplevels(subset(data, data$mite_species == i )) # For each iteration i equals mite species. 
  host<- as.vector(levels(focalmite$host_species)) # A vector with all the host of that mite species
  number_interactions<- (length(host))# Vector to calculate the number of pairwise comparisons
  final<-(number_interactions)*(number_interactions-1)/2 # Formule to extract the number of comparisons given a number of hosts
  npairwise<-c(npairwise,final)# Saving the number of pairwise comparisons. It is necessary to subset below.
  
  distance<- as.data.frame(dist[rownames(dist) %in% host , colnames(dist)%in% host]) # It creates a cophenetic phylo matrix between all the host of that mite species. It is done subsetting the bigger matrix (dist)
  if (nrow(distance) == 0) # This loop is because there is one species with a host not included in the matrix (Guira guira)
    distance = -1
  
  #Intervals  
  for (i in 1:10) { # 10 because is the maximum numbe of rows expected to matrizdef
    if (any(distance < (i*20) & distance >= (i-1)*20 & distance > 0)) { # Here it is specified the length of the interval. 20 here. The first one is from 0-20 then 21-40, etc.
      matrizdef[i,2] <- matrizdef[i,2] + 1 # Each time there is a distance fullfilling the requirements it adds one to the second column of matrizdef.
    }
  }
}
number_interactions<- sum(subset(npairwise, npairwise != 0)) # This is the number of comparisons. The zero ones are those mite species having only a single hosts (and are not included in this calculation).

#Total number of associations, 1305

#9) To estimate the likelihood
matrizdef[,2] <- matrizdef[,2] / sum(matrizdef[,2]) # Divided by the total number of associations.


# 10) Plotting
par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(5, 1, 0))

#Function to estimate CI from a proportion. Z=1.96 equales to alpha=0.05.
simpasym <- function(n, p, z=1.96, cc=TRUE){
  out <- list()
  if(cc){
    out$lb <- p - z*sqrt((p*(1-p))/n) - 0.5/n
    out$ub <- p + z*sqrt((p*(1-p))/n) + 0.5/n
  } else {
    out$lb <- p - z*sqrt((p*(1-p))/n)
    out$ub <- p + z*sqrt((p*(1-p))/n)
  }
  out
}

# Calculating CI
CI<-simpasym(n=1305, matrizdef[,2], z=1.96, cc=TRUE)

pdf("Figure4.pdf", width= 9,height= 6)

#Plotting CI

errbar(matrizdef[,1],matrizdef[,2],yplus= CI$ub, yminus= CI$lb, pch=16,  xaxt='n',xlab="", ylab="",las=2 ) # Este es el plot de la prob versus las distancias filogenéticas de los hospedadores.

box(lwd=2)
# axis(1, at=c(0,10,30,50,70,90,110,130,150,170,190), labels=c("","0-20","21-40","41-60","61-80", "81-100","101-120","121-140","141-160","161-180","181-200")) # Eje antiguo

axis(1, at=c(0,8,10,30,47,50,70,90,92,110,130,150,168,170,190,200), las=2,labels=c("0","","10","30","","50","70", "90","","110","130","150","","170","190","200"))#, cex.axis=0.65) # Así pinta también, las comparaciones para P.anthi

title(xlab="Phylogenetic distance between bird species sharing mite species", ylab= "Probability")


dev.off()


png("Figure4.png", width= 9,height= 6, units= "in",res=300)

#Plotting CI

errbar(matrizdef[,1],matrizdef[,2],yplus= CI$ub, yminus= CI$lb, pch=16,  xaxt='n',xlab="", ylab="",las=2 ) # Este es el plot de la prob versus las distancias filogenéticas de los hospedadores.

box(lwd=2)
# axis(1, at=c(0,10,30,50,70,90,110,130,150,170,190), labels=c("","0-20","21-40","41-60","61-80", "81-100","101-120","121-140","141-160","161-180","181-200")) # Eje antiguo

axis(1, at=c(0,8,10,30,47,50,70,90,92,110,130,150,168,170,190,200), las=2,labels=c("0","","10","30","","50","70", "90","","110","130","150","","170","190","200"))#, cex.axis=0.65) # Así pinta también, las comparaciones para P.anthi

title(xlab="Phylogenetic distance between bird species sharing mite species", ylab= "Probability")


dev.off()

#### Host-switching evolutionary consequences (Percentage of species consequence of major host-switches) #####

# Overall measure

m <-read.csv("global_fmbird.csv", sep="\t", header=TRUE)
m <- droplevels(filter(m, data_quality == 2, mite_superfamily == "Analgoidea"))
m <- droplevels(filter(m, mite_family == "Proctophyllodidae" | mite_family == "Trouessartiidae"| mite_family == "Pteronyssidae"))

# Counting the number of species and genera from these families out of passerines 
nspeciesa<- length(levels(m$mite_species)) # Total number of species
ngeneraa<-length(levels(m$mite_genus)) # Total number of genera
z <- paste(m$mite_species, m$host_species)
m$z <- z 
p <- stri_duplicated(m$z) 
data_f <- cbind(z,p)
m$rep <- data_f[,2]
m <- m[m$rep==FALSE,]
mpass <- droplevels(filter(m, host_order=="Passeriformes"))
genpass <- as.data.frame(with(mpass, table (mite_genus)))
genpass <- droplevels(filter(genpass,Freq>0))
m <- droplevels(anti_join(m, genpass, by = "mite_genus"))
p <- stri_duplicated(m$mite_species) 
data_f <- cbind(m,p)
m$rep <- data_f[,18]
m <- m[m$rep==FALSE,]
nspecies_outa<- length(levels(m$mite_species)) # Number of species from these families only present out of passerines
percentagesp<-(nspecies_outa/nspeciesa)*100
ngenus_outa<-length(levels(m$mite_genus)) # Number of genera of these families only present out of passerines
percentagenera<- (ngenus_outa/ngeneraa)*100

#Results
percentagesp
percentagenera


#### Proctophyllodidae

m <-read.csv("global_fmbird.csv", sep="\t", header=TRUE)
m <- droplevels(filter(m, data_quality == 2, mite_superfamily == "Analgoidea"))


m <- droplevels(filter(m, mite_family == "Proctophyllodidae"))

# Counting the number of species and genera from Proctophyllodidae out of passerines 
nspeciesa<- length(levels(m$mite_species)) # Total number of species
ngeneraa<-length(levels(m$mite_genus)) # Total number of genera
z <- paste(m$mite_species, m$host_species)
m$z <- z 
p <- stri_duplicated(m$z) 
data_f <- cbind(z,p)
m$rep <- data_f[,2]
m <- m[m$rep==FALSE,]
mpass <- droplevels(filter(m, host_order=="Passeriformes"))
genpass <- as.data.frame(with(mpass, table (mite_genus)))
genpass <- droplevels(filter(genpass,Freq>0))
m <- droplevels(anti_join(m, genpass, by = "mite_genus"))
p <- stri_duplicated(m$mite_species) 
data_f <- cbind(m,p)
m$rep <- data_f[,18]
m <- m[m$rep==FALSE,]


nspecies_outa<- length(levels(m$mite_species)) # Number of species from Proctophyllodidae only present out of passerines
percentagesp<-(nspecies_outa/nspeciesa)*100


ngenus_outa<-length(levels(m$mite_genus)) # Number of genera of Proctophyllodidae only present out of passerines
percentagenera<- (ngenus_outa/ngeneraa)*100

# Results
percentagesp
percentagenera


#### Trouessartiidae

m <-read.csv("global_fmbird.csv", sep="\t", header=TRUE)
m <- droplevels(filter(m, data_quality == 2, mite_superfamily == "Analgoidea"))
m <- droplevels(filter(m, mite_family == "Trouessartiidae"))


nspeciesa<- length(levels(m$mite_species)) # Total number of species
ngeneraa<-length(levels(m$mite_genus)) # Total number of genera

# Counting the number of species and genera from Trouessartiidae out of passerines 
z <- paste(m$mite_species, m$host_species)
m$z <- z 
p <- stri_duplicated(m$z) 
data_f <- cbind(z,p)
m$rep <- data_f[,2]
m <- m[m$rep==FALSE,]
mpass <- droplevels(filter(m, host_order=="Passeriformes"))
genpass <- as.data.frame(with(mpass, table (mite_genus)))
genpass <- droplevels(filter(genpass,Freq>0))
m <- droplevels(anti_join(m, genpass, by = "mite_genus"))
p <- stri_duplicated(m$mite_species) 
data_f <- cbind(m,p)
m$rep <- data_f[,18]
m <- m[m$rep==FALSE,]


nspecies_outa<- length(levels(m$mite_species)) # Number of species from Trouessartiidae only present out of passerines
percentagesp<-(nspecies_outa/nspeciesa)*100


ngenus_outa<-length(levels(m$mite_genus)) # Number of genera of Trouessartiidae
percentagenera<- (ngenus_outa/ngeneraa)*100

# Results 

percentagesp
percentagenera

#### Pteronyssidae

# Paseriformes origin scenario


m <-read.csv("global_fmbird.csv", sep="\t", header=TRUE)
m <- droplevels(filter(m, data_quality == 2, mite_superfamily == "Analgoidea"))
m <- droplevels(filter(m, mite_family == "Pteronyssidae"))


nspeciesa<- length(levels(m$mite_species)) # Total number of species
ngeneraa<-length(levels(m$mite_genus)) # Total number of genera

# Counting the number of species and genera from Pteronyssidae out of passerines 
z <- paste(m$mite_species, m$host_species)
m$z <- z 
p <- stri_duplicated(m$z) 
data_f <- cbind(z,p)
m$rep <- data_f[,2]
m <- m[m$rep==FALSE,]
mpass <- droplevels(filter(m, host_order=="Passeriformes"))
genpass <- as.data.frame(with(mpass, table (mite_genus)))
genpass <- droplevels(filter(genpass,Freq>0))
m <- droplevels(anti_join(m, genpass, by = "mite_genus"))
p <- stri_duplicated(m$mite_species) 
data_f <- cbind(m,p)
m$rep <- data_f[,18]
m <- m[m$rep==FALSE,]


nspecies_outa<- length(levels(m$mite_species)) # Number of species from Pteronyssidae only present out of passerines
percentagesp<-(nspecies_outa/nspeciesa)*100


ngenus_outa<-length(levels(m$mite_genus)) # Number of genera from Pteronyssidae only present out of passerines
percentagenera<- (ngenus_outa/ngeneraa)*100

# Results

percentagesp
percentagenera


# Piciformes origin scenario


m <-read.csv("global_fmbird.csv", sep="\t", header=TRUE)
m <- droplevels(filter(m, data_quality == 2, mite_superfamily == "Analgoidea"))
m <- droplevels(filter(m, mite_family == "Pteronyssidae"))


nspeciesa<- length(levels(m$mite_species)) # Total number of species
ngeneraa<-length(levels(m$mite_genus)) # Total number of genera


# Counting the number of species and genera from Pteronyssidae out of passerines 
z <- paste(m$mite_species, m$host_species)
m$z <- z 
p <- stri_duplicated(m$z) 
data_f <- cbind(z,p)
m$rep <- data_f[,2]
m <- m[m$rep==FALSE,]
mpass <- droplevels(filter(m, host_order=="Piciformes"))
genpass <- as.data.frame(with(mpass, table (mite_genus)))
genpass <- droplevels(filter(genpass,Freq>0))
m <- droplevels(anti_join(m, genpass, by = "mite_genus"))
p <- stri_duplicated(m$mite_species) 
data_f <- cbind(m,p)
m$rep <- data_f[,18]
m <- m[m$rep==FALSE,]


nspecies_outa<- length(levels(m$mite_species)) # Number of species from Pteronyssidae only present out of Piciformes
percentagesp<-(nspecies_outa/nspeciesa)*100


ngenus_outa<-length(levels(m$mite_genus)) # Number of genera from Pteronyssidae only present out of passerines
percentagenera<- (ngenus_outa/ngeneraa)*100


# Results

percentagesp
percentagenera




#### Host-switching evolutionary consequences (Figure 5) #####


# Database
data<- read.csv("global_fmbird.csv", sep="\t")

# Retaining only high-quality data
data<- droplevels(subset(data, data_quality == 2))

# Retaining only species from Proctophyllodidae, Trouessartiidae and Pteronyssidae
data <-droplevels(filter(data, mite_family == "Proctophyllodidae" | mite_family == "Trouessartiidae"| mite_family == "Pteronyssidae"))

# Enabling the correspondence of bird taxonomy in the Doña et al. (2016) database to that of Jetz et al. (2014). Requires Table S1 ()
data2<- data
data1<- read.csv("bird_species_names.csv", h=TRUE)
old=as.vector(data1$IOC_5.4)
n= as.vector(data1$Birdtree)
a= data2$host_species
for (i in 1:length(old))
{
  a=gsub(old[i],n[i], a)
}
data2$host_species = a
data<- as.data.frame( data2)

# Adding underscore in bird species names. It is require to work with the bird phylogenetic distances matrix.
Bnames= strsplit(as.vector(data$host_species), split=" ")
birdnames <- sapply(Bnames, function(x) paste(x[1]))
birdnames2<- sapply(Bnames, function(x) paste(x[2]))
data[["host_species"]] <- as.factor(paste(birdnames, birdnames2, sep="_"))

# To create a matrix with a single entry per bird-mite association
data[["concatenate"]] <- as.factor(paste(data$mite_species, data$host_species, sep="_"))
data<- droplevels(data[!duplicated(data$concatenate), ])

# Bird phylogeny 
fmtree <- read.nexus("fm_bird_tree.tre")

#
mitesp_birdsp2 <-with(data, table (droplevels(data[,9]), droplevels(data[,4])))
mitesp_birdsp2 <-as.data.frame(ifelse(mitesp_birdsp2>0,1,0))
species <- as.vector(row.names(mitesp_birdsp2))
mitesp_birdsp2["species"]<- species


# Creating the bird order (for plotting purposes)
birdorder <-with(data, table (droplevels(data[,9]), droplevels(data[,7])))
birdorder <-as.data.frame(ifelse(birdorder>0,1,0))
species <- as.vector(row.names(birdorder))
birdorder["species"]<- species
index<- as.vector (c(1:14))
birdorder2<- as.matrix(t(birdorder[1:1110,1:14]) * index)
birdorder2<- t(birdorder2)
test<- rowSums(birdorder2) 
test<- as.vector (test)
birdorder$order<- test
mitesp_birdsp2[,89]<- birdorder$order

# Preparing the tree
fmtree$node.label <- NULL # This is required to use comparative.data object to check the congruence among phylogeny and data
fmtree$edge.length[which(fmtree$edge.length <= 0)] <- 0.0000001 # It is a problem with zeroes with comparative.data. I found this solution in a Liam Revell post.

# Checking the comparative and tree data
fmcompf <- comparative.data(fmtree, mitesp_birdsp2, species) #Family; This is to check the congruence among phylogeny and database, i.e. which species are not present in both sites. 
fmtree <- read.nexus("fm_bird_tree.tre") # It is necessary to re-load the tree to avoid a problem when plotting, showing the root.

# Trimming the comparative and tree data
removez<- fmcompf$dropped$tips # Object, after create the compdata object and to delete the taxa that are not in phylogeny or in data. In this case, these are the phylogeny tips not present in data.

fmtree1<- drop.tip(fmtree, removez) # Now, removing the tips not present in data

rmdata<-fmcompf$dropped$unmatched.rows # Object,the same but with the data rows not present in phylogeny.

mitesp_birdsp<-mitesp_birdsp2[ !(rownames(mitesp_birdsp2) %in% rmdata), ] #Now removing data from dataset

mitesp_birdsp= mitesp_birdsp[,-88]



# Creating the color vector for plotting
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color2=  sample(color, 100)

#### Preparing data for trait.plot (creating a trait value for inhabiting only Passerines, only non-passerines or both)

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Proctophyllodidae"))
generosconmassaltos= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order != "Passeriformes")) 
levels(generosconmassaltos$host_order)
vec=sort(table(generosconmassaltos$mite_genus))
vec= dimnames(vec)
vec= as.vector(vec)
vec=(unlist(vec))

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Proctophyllodidae"))
generosconmassaltos_1= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order == "Passeriformes")) 
vec_1=sort(table(generosconmassaltos_1$mite_genus))
vec_1= dimnames(vec_1)
vec_1= as.vector(vec_1)
vec_1=(unlist(vec_1))

proctophyllodidaegenus= levels(proctophyllodidaegenus$mite_genus)
generosenlosdos=vec[(vec %in% vec_1)]
generersolopaser= proctophyllodidaegenus[!(proctophyllodidaegenus %in% vec)]
generossolonopaser= vec[!(vec %in% generosenlosdos)]
generossolonopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generossolonopaser] )
generossolonopaserorden=sort(table(generossolonopaserorden))
vecb= dimnames(generossolonopaserorden)
vecb= as.vector(vecb)
vecb=unlist(vecb)
generosenlosdossorden= droplevels(data$mite_genus[data$mite_genus %in% generosenlosdos] )
generosenlosdossorden=sort(table(generosenlosdossorden))
vecc= dimnames(generosenlosdossorden)
vecc= as.vector(vecc)
vecc=unlist(vecc)
generossolopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generersolopaser] )
generossolopaserorden=sort(table(generossolopaserorden))
vecd= dimnames(generossolopaserorden)
vecd= as.vector(vecd)
vecd=unlist(vecd)


# Vector to plot in black non-passerines orders 

test <- c("black","black","black","black","black","black","black","black","black","gray",
          "black","black","black","black")
test2<- "white"


# Trait plot matrix

# Proctophyllodidae

a=mitesp_birdsp2[,c(proctophyllodidaegenus, "V89")]
test3 <- droplevels(a[, c(vecb)])
test6<-a[,"V89"]
test7 <-droplevels(a[, c(vecc)])
test8 <-droplevels(a[, c(vecd)])
test5<- cbind(test8,test7,test3,test6)
ordeninfoparacolorear=colnames(test5)

#Genera only in passerines
vecd
#Genera only in non-passerines
vecb
#Generos en los dos
vecc

par(mfrow=c(3,1))
pdf("Figure5.pdf", width= 18,height= 12)

trait.plot(fmtree1, test5, cex.lab = .001, legend = FALSE, type="f", margin= 2,w=1/40, root.edge=FALSE, use.edge.length = TRUE, cex.legend= 0.15,
           
           cols =  list(Afroproterothrix= c("white", "dodgerblue2"), Berladectes= c("white", "dodgerblue2")
                        ,Bradyphyllodes= c("white", "dodgerblue2"),Favettea= c("white", "dodgerblue2"),Hemipterodectes= c("white", "dodgerblue2"), Lamelladectes= c("white", "dodgerblue2")
                        ,Neodectes= c("white", "dodgerblue2"),Philepittalges= c("white", "dodgerblue2"), Ptaismatophyllodes= c("white", "dodgerblue2"),Tyranniphyllodes= c("white", "dodgerblue2")
                        , Vireodectes= c("white", "dodgerblue2"), Atrichophyllodes= c("white", "dodgerblue2"), Cotingodectes= c("white", "dodgerblue2"),Rupicolacarus= c("white", "dodgerblue2")
                        ,Tanyphyllodes= c("white", "dodgerblue2"), Hemitriccodectes= c("white", "dodgerblue2"), Mimicalges= c("white", "dodgerblue2"),Exochojoubertia= c("white", "dodgerblue2")
                        , Metapterodectes= c("white", "dodgerblue2"),Nanopterodectes= c("white", "dodgerblue2"), Anisophyllodes= c("white", "dodgerblue2"),Pterodectes= c("white", "dodgerblue2")
                        ,Diproctophyllodes= c("white", "dodgerblue2"),Pedanodectes= c("white", "dodgerblue2"), Joubertophyllodes= c("white", "dodgerblue2"),Nycteridocaulus= c("white", "dodgerblue2")
                        ,Anisodiscus= c("white", "dodgerblue2"),Hadrophyllodes= c("white", "dodgerblue2"), Tyrannidectes= c("white", "dodgerblue2"), Alaudicola= c("white", "dodgerblue2"),Alaudicola= c("white", "dodgerblue2")
                        , Monojoubertia= c("white", "dodgerblue2"),Dolichodectes= c("white", "dodgerblue2"), Platyacarus= c("white","dodgerblue2"),Amerodectes= c("white", "dodgerblue2")
                        ,
                        Proterothrix= c("white", "darkorange2"), Montesauria= c("white", "darkorange2"), Proctophyllodes= c("white", "darkorange2")
                        ,
                        Anorthalloptes= c("white", "darkolivegreen3"), Picipterodectes= c("white", "darkolivegreen3"), Pteralloptes= c("white", "darkolivegreen3"),Schizodectes= c("white", "darkolivegreen3")
                        ,Sclerodectes= c("white", "darkolivegreen3"),Syntomodectes= c("white", "darkolivegreen3"), Ptyctophyllodes= c("white", "darkolivegreen3"),Rhamphocaulus= c("white", "darkolivegreen3")
                        ,Xynonodectes= c("white", "darkolivegreen3"),Allodectes= c("white", "darkolivegreen3"), Toxerodectes= c("white", "darkolivegreen3"),Trochilodectes= c("white","darkolivegreen3")
                        , 
                        test6= c(test)))


#Trouessartidae

####Preparing data for trait.plot (creating a trait value for inhabiting only Passerines, only non-passerines or both)

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Trouessartiidae"))
generosconmassaltos= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order != "Passeriformes")) 
levels(generosconmassaltos$host_order)
vec=sort(table(generosconmassaltos$mite_genus))
vec= dimnames(vec)
vec= as.vector(vec)
vec=(unlist(vec))

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Trouessartiidae"))
generosconmassaltos_1= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order == "Passeriformes")) 
vec_1=sort(table(generosconmassaltos_1$mite_genus))
vec_1= dimnames(vec_1)
vec_1= as.vector(vec_1)
vec_1=(unlist(vec_1))
proctophyllodidaegenus= levels(proctophyllodidaegenus$mite_genus)
generosenlosdos=vec[(vec %in% vec_1)]
generersolopaser= proctophyllodidaegenus[!(proctophyllodidaegenus %in% vec)]
generossolonopaser= vec[!(vec %in% generosenlosdos)]
generossolonopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generossolonopaser] )
generossolonopaserorden=sort(table(generossolonopaserorden))
vecb= dimnames(generossolonopaserorden)
vecb= as.vector(vecb)
vecb=unlist(vecb)
generosenlosdossorden= droplevels(data$mite_genus[data$mite_genus %in% generosenlosdos] )
generosenlosdossorden=sort(table(generosenlosdossorden))
vecc= dimnames(generosenlosdossorden)
vecc= as.vector(vecc)
vecc=unlist(vecc)
generossolopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generersolopaser] )
generossolopaserorden=sort(table(generossolopaserorden))
vecd= dimnames(generossolopaserorden)
vecd= as.vector(vecd)
vecd=unlist(vecd)

# Vector to plot in black non-passerines orders 

test <- c("black","black","black","black","black","black","black","black","black","gray",
          "black","black","black","black")
test2<- "white"


# Trait plot matrix

a=mitesp_birdsp2[,c(proctophyllodidaegenus, "V89")]
test3 <- droplevels(a[, c(vecb)])
test6<-a[,"V89"]
test7 <-droplevels(a[, c(vecc)])
test8 <-droplevels(a[, c(vecd)])
test5<- cbind(test8,test7,test3,test6)

ordeninfoparacolorear=colnames(test5)

#Genera only in passerines
vecd
#Genera only in non-passerines
vecb
#Genera in both
vecc


trait.plot(fmtree1, test5, cex.lab = .001, legend = FALSE, type="f", margin= 2,w=1/40, root.edge=FALSE, use.edge.length = TRUE, cex.legend= 0.15,
           , cols =  list(Hemicalcealges= c("white", "dodgerblue2"), Neocalcealges= c("white", "dodgerblue2"),  Bicentralges= c("white", "dodgerblue2")
                          ,
                          Arthrogynalges = c("white", "darkorange2"),Calcealges = c("white", "darkorange2"),  Trouessartia= c("white", "darkorange2")
                          , 
                          Anisanchus = c("white", "darkolivegreen3"), Ptyctalloptes= c("white", "darkolivegreen3")
                          ,Steatacarus = c("white", "darkolivegreen3"), Analloptes= c("white", "darkolivegreen3"), Allanalges= c("white", "darkolivegreen3"),Proterocaulus= c("white", "darkolivegreen3")
                          ,Pseudalges= c("white", "darkolivegreen3"),Uniscutalges= c("white", "darkolivegreen3"), test6= c(test),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2), test6= c(test2),test6= c(test2)))



########### Pteronyssidae

#### Preparing data for trait.plot (creating a trait value for inhabiting only Passerines, only non-passerines or both)

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Pteronyssidae"))
generosconmassaltos= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order != "Passeriformes")) 
levels(generosconmassaltos$host_order)
vec=sort(table(generosconmassaltos$mite_genus))
vec= dimnames(vec)
vec= as.vector(vec)
vec=(unlist(vec))

proctophyllodidaegenus= droplevels(subset(data, data$mite_family == "Pteronyssidae"))
generosconmassaltos_1= droplevels(subset(proctophyllodidaegenus, proctophyllodidaegenus$host_order == "Passeriformes")) 
vec_1=sort(table(generosconmassaltos_1$mite_genus))
vec_1= dimnames(vec_1)
vec_1= as.vector(vec_1)
vec_1=(unlist(vec_1))
proctophyllodidaegenus= levels(proctophyllodidaegenus$mite_genus)
generosenlosdos=vec[(vec %in% vec_1)]
generersolopaser= proctophyllodidaegenus[!(proctophyllodidaegenus %in% vec)]
generossolonopaser= vec[!(vec %in% generosenlosdos)]
generossolonopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generossolonopaser] )
generossolonopaserorden=sort(table(generossolonopaserorden))
vecb= dimnames(generossolonopaserorden)
vecb= as.vector(vecb)
vecb=unlist(vecb)
generosenlosdossorden= droplevels(data$mite_genus[data$mite_genus %in% generosenlosdos] )
generosenlosdossorden=sort(table(generosenlosdossorden))
vecc= dimnames(generosenlosdossorden)
vecc= as.vector(vecc)
vecc=unlist(vecc)
generossolopaserorden= droplevels(data$mite_genus[data$mite_genus %in% generersolopaser] )
generossolopaserorden=sort(table(generossolopaserorden))
vecd= dimnames(generossolopaserorden)
vecd= as.vector(vecd)
vecd=unlist(vecd)


# Vector to plot in black non-passerines orders 

test <- c("black","black","black","black","black","black","black","black","black","gray",
          "black","black","black","black")
test2<- "white"


# Trait plot matrix

a=mitesp_birdsp2[,c(proctophyllodidaegenus, "V89")]
test3 <- droplevels(a[, c(vecb)])
test6<-a[,"V89"]
test7 <-droplevels(a[, c(vecc)])
test8 <-droplevels(a[, c(vecd)])
test5<- cbind(test8,test7,test3,test6)

ordeninfoparacolorear=colnames(test5)

#Genera only in passerines
vecd
#Genera only in non-passerines
vecb
#Genera in both
vecc


trait.plot(fmtree1, test5, cex.lab = .001, legend = FALSE, type="f", margin= 2,w=1/40, root.edge=FALSE, use.edge.length = TRUE, cex.legend= 0.15,
           , cols =  list(Mouchetia= c("white", "dodgerblue2"), Micropteroherpus= c("white", "dodgerblue2"),  Dicrurobius= c("white", "dodgerblue2"),  Vanginyssus = c("white", "dodgerblue2"), 
                          Timalinyssus = c("white", "dodgerblue2"),  Sturnotrogus= c("white", "dodgerblue2"), Metapteronyssus = c("white", "dodgerblue2"), Pteroherpus= c("white", "dodgerblue2")
                          ,Pteronyssoides = c("white", "dodgerblue2"), Scutulanyssus= c("white", "dodgerblue2")
                          , 
                          Parapteronyssus= c("white", "darkolivegreen3"),Stenopteronyssus= c("white", "darkolivegreen3")
                          ,Megalaimobius= c("white", "darkolivegreen3"),Zygepigynia= c("white", "darkolivegreen3"),    Pegopteronyssus= c("white", "darkolivegreen3"), Cleyastobius= c("white","darkolivegreen3")
                          ,Pterotrogus= c("white", "darkolivegreen3"), Monapsidus= c("white", "darkolivegreen3"),Pteronyssus= c("white","darkolivegreen3"), Ramphastobius= c("white", "darkolivegreen3")
                          ,Neopteronyssus= c("white", "darkolivegreen3"),Anephippius= c("white", "darkolivegreen3"),Conomerus= c("white", "darkolivegreen3"), 
                          test6= c(test),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),
                          test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2),test6= c(test2), test6= c(test2)))

dev.off()



