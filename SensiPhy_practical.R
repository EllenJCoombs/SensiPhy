library(ggplot2)

#Install sensiphy from github
library(devtools)
devtools::install_github("paternogbc/sensiPhy")

# Load data
alien <- data
# This analysis needs a multiphylo file:
class(alien$phy) #Checking it's multiphylo
alien$phy

#Load Tree phylum (tree_phylum doc)
#
#Run PGLS accounting for phylogenetic uncertain
tree <- tree_phylm(log(gestaLen) ~ log(adultMass), phy = data$phy, 
                   data = alien$data, n.tree = 30)
# To check summary results:
summary(tree)
# Visual diagnostics
sensi_plot(tree)

View(alien)


