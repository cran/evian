#
#  Part of the evian R package, http://strug.ccb.sickkids.ca/evian
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details. The license describes
#  your right to use this software.
#
#  evian is Copyrighted and the ownership of intellectual rights belongs
#  to Dr. Lisa J. Strug at the Hospital for Sick Children.
#
#  Copyright © 2010 by Lisa J Strug. Hospital for Sick Children, Toronto.
#
glr <-
  function (data, bim, xcols=NULL, ycol=NULL, covariateCol=NULL, formula=NULL, c, robust=FALSE, model="additive",m=200, bse=5, lolim=NULL, hilim=NULL, multiThread=1,family='gaussian',plinkCC=F)
  {
    tryCatch(
      {
        if (multiThread>1){
          cl=parallel::makeCluster(as.numeric(multiThread))
          registerDoParallel(cl)
        }
      }, error = function(e) {
        multiThread=1
        warning('Not supported multi-threads, only 1 thread being used')
      }
    )

    if (model != "additive" & model != "dominant" & model != "recessive" & model != "overdominance" & model != "2df") { stop ("Model needs to be specified as one of the followings: additive, dominant, recessive, or overdominance") }
    #print ("Computing likelihood intervals...")

    if ((is.null(xcols)|is.null(ycol))&is.null(formula)){stop ("No X/Y specified, please use either xcols, ycols argument or use formula argument") }
    if (!(is.null(xcols)&is.null(ycol)|is.null(formula))){stop ("Please use either xcols and ycols argument or use formula argument, don't provide both") }

    if (!is.null(formula)){ # we just need to extract SNP columns and fit to profileLikehood function
      # to fit in the profileLikelihood, we need to have a list of x input as 'profile.theta' option, and the rest as the formula inputted
      formula_type=as.formula(formula)
      xandCov=strsplit(as.character(formula_type[3]),split =' \\+ ')[[1]] # Formula format has separator as '[space]+[space]'
      x_list=xandCov[grepl('^rs',xandCov)]; covList=paste(xandCov[!grepl('^rs',xandCov)],collapse = '+')
      # we don't want to place a different representation for the formula input
      # However, the profileLikelihood only does scalar calculation
      # so if multiple SNPs were input from the formula option, we can only work on one
      # we shall assume to be the first one, and output a warning message
      if (length(x_list)>1){
        # indicating more than 1 rs SNP found
        # provide a warning message and only using the first SNP
        warnings('More than 1 rsSNP detected in the formula specification. Evian only handles one parameter of interests at a time.
                 Only the first rsSNP input will be considered as the parameter of interest and the rest will be
                 considered as nuisance parameters. If you are looking for running associations of multiple SNPs
                 individually with the same phenotype, please use the xcol, ycol options to specify.')
        x_list=x_list[1]
      }

      if (covList==''){ covList='1'} # so this means no covariate in the original model
      #this is the formula that will be sent to profileLikelihood function
      formula_tofit=as.formula(paste0(as.character(formula_type[2]),'~',covList))
    } else { # so Xcols and ycol are inputted, maybe covariateCol are also inputted, we can generate new formula option to fit profileLikelihood function
      y_name=names(data)[as.numeric(ycol)]
      x_list_raw=names(data)[as.numeric(xcols)]
      x_list=mapply(FUN=function(t) strsplit(x_list_raw[t],split='_')[[1]][1],1:length(x_list_raw))
      if (is.null(covariateCol)){ covList='1' } else { covList=paste(names(data)[as.numeric(covariateCol)],collapse ='+')}
      formula_tofit=as.formula(paste0(y_name,'~',covList))
    }

    # decide which model to use
    if (multiThread != 1){
      result=foreach(snpIndex=1:length(x_list),.combine = rbind) %dopar% {
        calculateGLR(snp=x_list[snpIndex],formula_tofit=formula_tofit,model=model,data=data,bim=bim,lolim=lolim,hilim=hilim,m=m,bse=bse,c=c,family=family,plinkCC=plinkCC)
      }
      parallel::stopCluster(cl)
    } else {
      result=NULL
      for (snpIndex in 1:length(x_list)){
        result=rbind(result,calculateGLR(snp=x_list[snpIndex],formula_tofit=formula_tofit,model=model,data=data,bim=bim,lolim=lolim,hilim=hilim,m=m,bse=bse,c=c,family=family,plinkCC=plinkCC))
      }
    }
    return(result)
  }
