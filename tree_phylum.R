

tree_phylm <- function(formula,data,phy,
                       n.tree=2,model="lambda",track=TRUE,...){
  #Error check
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if ( (model == "trend") && (ape::is.ultrametric(phy)))
    stop("Trend is unidentifiable for ultrametric trees., see ?phylolm for details")
  else
    
    #Matching tree and phylogeny using utils.R
    datphy<-match_dataphy(formula,data,phy, ...)
  full.data<-datphy[[1]]
  phy<-datphy[[2]]
  
  # If the class of tree is multiphylo pick n=n.tree random trees
  trees<-sample(length(phy),n.tree,replace=F)
  
  #Create the results data.frame
  sensi.estimates<-data.frame("n.tree"=numeric(),"intercept"=numeric(),"se.intercept"=numeric(),
                              "pval.intercept"=numeric(),"estimate"=numeric(),"se.estimate"=numeric(),
                              "pval.estimate"=numeric(),"aic"=numeric(),"optpar"=numeric())
  
  #Model calculation
  counter=1
  errors <- NULL
  c.data<-list()
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  for (j in trees){
    
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #phylolm model
    mod = try(phylolm::phylolm(formula, data=full.data,model=model,phy=phy[[j]]),FALSE)
    
    
    if(isTRUE(class(mod)=="try-error")) {
      error <- j
      names(error) <- rownames(c.data$full.data)[j]
      errors <- c(errors,error)
      next }
    
    
    else{
      intercept            <- phylolm::summary.phylolm(mod)$coefficients[[1,1]]
      se.intercept         <- phylolm::summary.phylolm(mod)$coefficients[[1,2]]
      estimate                <- phylolm::summary.phylolm(mod)$coefficients[[2,1]]
      se.estimate             <- phylolm::summary.phylolm(mod)$coefficients[[2,2]]
      pval.intercept       <- phylolm::summary.phylolm(mod)$coefficients[[1,4]]
      pval.estimate           <- phylolm::summary.phylolm(mod)$coefficients[[2,4]]
      aic.mod              <- mod$aic
      n                    <- mod$n
      d                    <- mod$d
      
      if (model == "BM"){
        optpar <- NA
      }
      if (model != "BM"){
        optpar               <- mod$optpar
      }
      
      if(track==TRUE) utils::setTxtProgressBar(pb, counter)
      
      #write in a table
      estim.simu <- data.frame(j, intercept, se.intercept, pval.intercept,
                               estimate, se.estimate, pval.estimate, aic.mod, optpar,
                               stringsAsFactors = F)
      sensi.estimates[counter, ]  <- estim.simu
      counter=counter+1
      
    }
  }
  if(track==TRUE) on.exit(close(pb))
  #calculate mean and sd for each parameter
  #variation due to tree choice
  statresults<-data.frame(min=apply(sensi.estimates,2,min),
                          max=apply(sensi.estimates,2,max),
                          mean=apply(sensi.estimates,2,mean),
                          sd_tree=apply(sensi.estimates,2,stats::sd))
  
  
  statresults$CI_low  <- statresults$mean - qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
  statresults$CI_high <- statresults$mean + qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
  
  res <- list(   call = match.call(),
                 formula=formula,
                 data=full.data,
                 sensi.estimates=sensi.estimates,N.obs=n,
                 stats = round(statresults[c(1:6),c(3,5,6)],digits=3),
                 all.stats = statresults)
  class(res) <- "sensiTree"
  return(res)
}

