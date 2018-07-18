robust_forCluster=function(formula,data,family){
  # Robust function was newly developed by Zeynep Baskurt - Mar 23, 2018.
  # The robust factor will adjust the likelihood llik1[] for indivs who are related (ie. in the same FamID, FID).

  if (is.null(data$FID)){
    warning('no FID detected in the data, automatically assume each data is independent')
    data$FID=seq(1,nrow(data))
  }
  covList=as.character(formula)[3]
  data_sorted=data[order(data$FID),]
  reg=glm(paste0(formula[2],'~X+',covList),data=data_sorted,family=family)
  famID=data_sorted$FID
  if (family=='gaussian'){
    robustness_factor=(vcov(reg)/vcovCL(reg,famID,type="HC1"))[2,2]
  } else{
    robustness_factor=(vcov(reg)/vcovCL(reg,famID,type="HC0"))[2,2]
  }

  return(robustness_factor)
}


