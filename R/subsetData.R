subsetData=function(snp,formula_tofit,data,plinkCC){
  # We subset the dataframe, and change the genotype accordingly based on choice of the model
  # the raw format for each SNP is marked as rsxxx_[ATCG] or rsxxx_[1234]
  # we substring just the first part, left out the ATCG (this info will be recovered in the bim file anyways)
  snpName=mapply(FUN=function(t) strsplit(names(data)[t],split='_')[[1]][1],1:ncol(data))
  data_Xcol=data[,which(snpName==snp)]
  covList=as.character(formula_tofit)[3]
  data_FID=data[,names(data)=='FID']
  data_cov=data[,names(data) %in% strsplit(covList,split =' \\+ ')[[1]]] # get the covariance columns
  data_y=data[,names(data)==as.character(formula_tofit)[2]]
  if (plinkCC){
    #the case/control output from plink is 1/2 instead of 0/1, change it
    data_y=as.numeric(data_y)-1
  }
  data_subList=cbind(data_Xcol,data_y,data_FID=data[,names(data)=='FID'],data_cov,stringsAsFactors=F)
  names(data_subList)[1:3]=c('X',as.character(formula_tofit)[2],'FID')
  data_nomiss=na.omit(data_subList)
  return(data_nomiss)
}
