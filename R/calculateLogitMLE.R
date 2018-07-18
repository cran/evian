calculateLogitMLE=function(snp,formula_tofit,model,data,bim,lolim,hilim,m,bse,k,robust,family){
  data_subset=subsetData(snp,formula_tofit,data) # subset so that only contains 1 snp
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
  
  if (!is.null(bse)| missing(bse)){ # if bse is provided, it will ignore the input lower and upper bound.
    #running this feature will significant increase calculation time
    bounds=getGridBound(formula=formula_tofit,data=data_nomiss,bse=bse,k=k,m=m,family=family,robust=robustFactor)
    lolim=bounds[1]; hilim=bounds[2]
  }

  rst=profilelike.glm(formula=formula_tofit,profile.theta ='X',data=data_nomiss,lo.theta=lolim,hi.theta = hilim,length=m,family=family)

  # if there are multiple theta values that correspondingly to the highest PL, use the first one (we dont distinguish them anyways)
  betaMLE=as.numeric(rst$theta[rst$profile.lik.norm==1][1])
  #FLip the OR depending on the sign of beta
  if (betaMLE<0){
    mleOR=1/exp(betaMLE); flip=T; theta_new=1/exp(rst$theta);effectAlelle=alt; protectiveAllele=ref
    af=1-af #We want af to represent the effective allele
  } else {
    mleOR=exp(betaMLE); flip=F; theta_new=exp(rst$theta);effectAlelle=ref; protectiveAllele=alt
  }
  ## Here we need to worry about dominant/recessive/2df model, since they cannot be simply flipped
  # We need to recalculate it if we need to flip the genotype
  if ((model %in% c("dominant","recessive",'2df'))&flip){
    #inside this if statement, we already flipped the alleles, so we can keep it unchanged
    # Since the allele is swapped, the allele frequency should be changed too.
    #But we need to recalculate mle since data is completely different
    data_subset$X=2-data_subset$X; af=1-af
    data_nomiss=adjustModel(data_subset,model)
    rst=profilelike.glm(formula=formula_tofit,profile.theta ='X',data=data_nomiss,lo.theta=lolim,hi.theta = hilim,length=m,family=family)
    mleOR=exp(as.numeric(rst$theta[rst$profile.lik.norm==1][1]))
    theta_new=exp(rst$theta)
  }

  ### From here, we start to prepare for output
  # For logistic regression, we are more interested in OR rather than original beta. So we output ORs
  rst_new=list(theta=theta_new,profile.lik.norm=(rst$profile.lik.norm)^robustFactor)
  #calculate summary stats at each k cut-offs
  k_df=data.frame(temp=1); k_ord=k[order(k)]
  for (kcut in k_ord){
    cutoff=1/as.numeric(kcut)
    #Bound for beta(with robust factor)
    loBound=min(rst$theta[rst_new$profile.lik.norm >=cutoff]);
    hiBound=max(rst$theta[rst_new$profile.lik.norm >=cutoff])
    # rst_new and rst has the same order, can be used interchangably when searching index
    # the bound defined as (normlized PL)^robustFactor >= 1/k
    #Output bound for OR (depending on the flip)
    if (flip){
      k_df=cbind(k_df,data.frame(exp(-hiBound), exp(-loBound)))
    } else {
      k_df=cbind(k_df,data.frame(exp(loBound), exp(hiBound)))
    }
  }
  k_df$temp=NULL;
  colnames(k_df)=c(rbind(paste0('lo_',1:length(k_ord)),paste0('hi_',1:length(k_ord))))

  # summarize all results to output

  rst_new$k_cutoff=k_ord
  rst_new$SummaryStats=cbind(data.frame(mle=mleOR,
                                        maxlr=1/rst_new$profile.lik.norm[which.min(abs(rst$theta))],
                                        AF=af,SNP=snp,bp=pos,effect=effectAlelle, ref=protectiveAllele,robustFactor=robustFactor,stringsAsFactors = F),k_df)

  return(rst_new)
}
