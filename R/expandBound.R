expandBound=function(data,bse,parameters,formula,m,k,family){
  # This is essentially a while loop; This is using Lehang's method to provide a proper bound
  bhat=parameters[1]; bhat_se=parameters[2]; ratio_c=parameters[3]
  lolim=bhat-bse*bhat_se/min(ratio_c,1); hilim=bhat+bse*bhat_se/min(ratio_c,1)
  rst=profilelike.glm(formula=formula,profile.theta = 'X',data=data,lo.theta=lolim,hi.theta = hilim,length=m,round=round,family=family)
  # we need to make sure there is value that is less than the min 1/k value, i.e. the distribution had more values in all bounds
  kcut=1/max(k)
  if (sum(rst$profile.lik.norm<kcut)>0){
    return(c(lolim,hilim))
  } else {
    expandBound(data=data,bse=bse+1,parameters,formula,m,k,family)
  }
}
