library(ggplot2)

#Install sensiphy from github
library(devtools)
devtools::install_github("paternogbc/sensiPhy")



#1. tree_phylm: tree method for phylogenetic linear regression


# Load data
alien <- data

# This analysis needs a multiphylo file:
class(alien.phy$phy) #Checking it's multiphylo
alien.phy$phy

#Load Tree phylum (tree_phylum doc)
#Run PGLS accounting for phylogenetic uncertain
tree <- tree_phylm(log(gestaLen) ~ log(adultMass), phy = data$phy, 
                   data = alien.phy$data, n.tree = 30)
# To check summary results:
summary(tree)
# Visual diagnostics
sensi_plot(tree)

View(alien)


#2. tree_physig(): tree method for phylogenetic signal

# Load data:
data(alien)
# Logtransform data
alien.data$logMass <- log(alien.data$adultMass) 
# Run sensitivity analysis:
tree <- tree_physig(trait.col = "logMass", data = alien.data, phy = alien.phy)
summary(tree)
sensi_plot(tree)
sensi_plot(tree, graphs = 1)
sensi_plot(tree, graphs = 2)


#3. tree_continuous() and tree_discrete: tree method for trait evolution of continous and discrete characters

# Load data:
primates

# Model trait evolution accounting for phylogenetic uncertainty
adultMass<-primates$data$adultMass
names(adultMass)<-rownames(primates$data)
tree_cont<-tree_continuous(data = adultMass,phy = primates$phy,
                           model = "OU",n.tree=30,n.cores = 2,track = TRUE)
# Print summary statistics for the transitions rates, aic-values and (if applicable) 
# optimisation parameter
summary(tree_cont)
## Visual diagnostics
sensi_plot(tree_cont,graphs="sigsq")
sensi_plot(tree_cont,graphs="optpar")


#Discrete characters 
# Load data:
primates
# Create a binary trait factor 
adultMass_binary<-ifelse(primates$data$adultMass > 7350, "big", "small")
adultMass_binary<-as.factor(as.factor(adultMass_binary))
names(adultMass_binary)<-rownames(primates$data)
# Model trait evolution accounting for phylogenetic uncertainty
tree_binary<-tree_discrete(data = adultMass_binary,phy = primates$phy,
                           model = "ARD",transform = "none",n.tree = 30,n.cores = 2,track = TRUE)
# Print summary statistics for the transitions rates, aic-values and (if applicable) 
# optimisation parameter
summary(tree_binary)
## Visual diagnostics
sensi_plot(tree_binary,graphs="q12")
sensi_plot(tree_binary,graphs="q21")

sensi_plot(tree_binary,graphs="q11")




#Continuous characters 
# Load data:
primates
# Model trait evolution accounting for phylogenetic uncertainty
adultMass<-primates$data$adultMass
names(adultMass)<-rownames(primates$data)
tree_cont<-tree_continuous(data = adultMass,phy = primates$phy,
                           model = "OU",n.tree=30,n.cores = 2,track = TRUE)
# Print summary statistics for the transitions rates, aic-values and (if applicable) 
# optimisation parameter
summary(tree_cont)
## Visual diagnostics
sensi_plot(tree_cont,graphs="sigsq")
sensi_plot(tree_cont,graphs="optpar")
#Plot
sensi_plot(tree_cont) 


# Use a different evolutionary model or transformation (e.g. Pagel's lambda)
tree_binary_lambda<-tree_discrete(data = adultMass_binary,phy = primates$phy,
                                  model = "SYM",transform = "lambda",n.tree = 30,n.cores = 2,track = TRUE)
summary(tree_binary_lambda)

sensi_plot(tree_binary_lambda) 


#Paper example 
data <- primates.data
data.frame <- primates.data
fit2 <- tree_clade_phylm(log(sexMaturity) ~ log(adultMass))

is.data.frame(primates.data)




