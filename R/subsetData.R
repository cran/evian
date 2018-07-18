subsetData=function(snp,formula_tofit,data){
  # We subset the dataframe, and change the genotype accordingly based on choice of the model
  data_Xcol=data[,names(data)==snp]
  covList=as.character(formula_tofit)[3]
  data_FID=data[,names(data)=='FID']
  data_cov=data[,names(data) %in% strsplit(covList,split =' \\+ ')[[1]]] # get the covariance columns
  data_y=data[,names(data)==as.character(formula_tofit)[2]]
  data_subList=cbind(data_Xcol,data_y,data_FID=data[,names(data)=='FID'],data_cov,stringsAsFactors=F)
  names(data_subList)[1:3]=c('X',as.character(formula_tofit)[2],'FID')
  data_nomiss=na.omit(data_subList)
  return(data_nomiss)
}
