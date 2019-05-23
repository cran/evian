calculateGLR=function(snp,formula_tofit,model,data,bim,lolim,hilim,m,bse,family,c,plinkCC){
  #same as evian linear/logit in cleaning data, picking models, and calculate likelihood
  # GLR does not have robust correction
  data_subset=subsetData(snp,formula_tofit,data,plinkCC) # subset so that only contains 1 snp
  af=sum(as.numeric(data_subset$X))/(2*nrow(data_subset)) #it is wrt to ref at this stage
  data_nomiss=adjustModel(data_subset,model) #adjust genotypes based on model
  names(bim)=c('chr','snp','physicalDist','pos','ref','alt') #ref is the allele in effect
  ref=bim$ref[bim$snp==snp]; alt=bim$alt[bim$snp==snp]; pos=bim$pos[bim$snp==snp]

  if (model=='overdominance'){ # parameter of interest is the deviation from additive
    formula_tofit=as.formula(paste(as.character(formula_tofit)[2],'~X1+',as.character(formula_tofit)[3]))
  }
  robustFactor=1
  if (is.null(lolim)| is.null(hilim)){ # if bse is provided, it will ignore the input lower and upper bound.
    #running this feature will significant increase calculation time
    bounds=getGridBound(formula=formula_tofit,data=data_nomiss,bse=bse,k=1,m=m,family=family,robust=robustFactor)
    lolim=bounds[1]; hilim=bounds[2]
  }
  # the therotical derivation of GLR require a symmetrical region, i.e. H0 -c<=theta<=c
  c=abs(c)
  # we have to fix the theta estimation region (lolim and hilim) such that it contains (-c,c) and symmetrical
  maxBound=max(abs(lolim),abs(hilim),c)
  if (maxBound==c){
    stop(paste0('paramater c input (|c|=',c,') exceeds the current grid search range (',lolim,', ',hilim,') for theta estimation. Please provide a suitable search range or c value.'))
  }
  lolim=-maxBound; hilim=maxBound

  rst=profilelike.glm(formula=formula_tofit,profile.theta ='X',data=data_nomiss,lo.theta=lolim,hi.theta = hilim,length=m,family=family)

  alternative=rst$profile.lik.norm[rst$theta < -c | rst$theta > c]
  null=rst$profile.lik.norm[rst$theta >= -c & rst$theta <= c]
  glr=max(alternative, na.rm=T)/max(null, na.rm=T)
  # summarize all results to output
  summaryStats=data.frame(GLR=glr,boundary=c,AF=af,SNP=snp,bp=pos,effect=ref,ref=alt,stringsAsFactors = F)
  return(summaryStats)
}
