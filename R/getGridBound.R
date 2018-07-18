getGridBound=function(formula,data,bse,k,m,family,robust){
  formula_withX=paste0(as.character(formula)[2],'~X+',as.character(formula)[3])
  reg=glm(formula_withX,data=data,family=family)
  bhat=summary(reg)$coefficient['X',1] ; bhat.se=summary(reg)$coefficient['X',2]
  mleParameters=c(bhat,bhat.se,robust)
  finalRange=expandBound(data,bse,parameters=mleParameters,formula=formula,m=m,k=k,family=family)
  return(finalRange)
}
