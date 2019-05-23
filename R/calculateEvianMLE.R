calculateEvianMLE=function(snp,formula_tofit,model,data,bim,lolim,hilim,m,bse,k,robust,family,plinkCC){
  data_subset=subsetData(snp=snp,formula_tofit=formula_tofit,data=data,plinkCC=plinkCC) # subset so that only contains 1 snp
  af=sum(as.numeric(data_subset$X))/(2*nrow(data_subset)) #it is wrt to ref at this stage
  data_nomiss=adjustModel(data_subset,model) #adjust genotypes based on model
  names(bim)=c('chr','snp','physicalDist','pos','ref','alt') #ref is the allele in effect
  ref=bim$ref[bim$snp==snp]; alt=bim$alt[bim$snp==snp]; pos=bim$pos[bim$snp==snp]

  if (model=='overdominance'){ # parameter of interest is the deviation from additive
    formula_tofit=as.formula(paste(as.character(formula_tofit)[2],'~X1+',as.character(formula_tofit)[3]))
  }
  # Robust factor correction
  if (robust){
    robustFactor=as.numeric(robust_forCluster(data=data_nomiss,formula=formula_tofit,family=family))
  } else {
    robustFactor=1
  }
  if (is.null(lolim)| is.null(hilim)){ # we used the MLE to find the lo/high limit
    # Directly from the original evian package, linear_plotSNP function
    bounds=getGridBound(formula=formula_tofit,data=data_nomiss,bse=bse,k=k,m=m,family=family,robust=robustFactor)
    lolim=bounds[1]; hilim=bounds[2]
  }


  rst=profilelike.glm(formula=formula_tofit,profile.theta = 'X',data=data_nomiss,lo.theta=lolim,hi.theta = hilim,length=m,family=family)
  rst_new=list(theta=rst$theta,profile.lik.norm=(rst$profile.lik.norm)^robustFactor)

  #calculate the intervals here:
  k_df=data.frame(temp=1); k_ord=k[order(k)]
  for (kcut in k_ord){
    cutoff=1/as.numeric(kcut)
    k_df=cbind(k_df,data.frame(min(rst_new$theta[rst_new$profile.lik.norm >=cutoff]),
                               max(rst_new$theta[rst_new$profile.lik.norm >=cutoff])))
  }
  k_df$temp=NULL;
  colnames(k_df)=c(rbind(paste0('lo_',1:length(k_ord)),paste0('hi_',1:length(k_ord))))
  ##might need a flip feature here

  # summarize tghe result
  rst_new$k_cutoff=k_ord
  rst_new$SummaryStats=cbind(data.frame(mle=as.numeric(rst_new$theta[rst_new$profile.lik.norm==1][1]),
                                        maxlr=1/rst_new$profile.lik.norm[which.min(abs(rst_new$theta))],
                                        AF=af,SNP=snp,bp=pos,effect=ref,ref=alt,robustFactor=robustFactor,stringsAsFactors = F),k_df)


  return(rst_new)
}
