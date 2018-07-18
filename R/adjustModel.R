adjustModel=function(data_nomiss,model){
  # change genotype coding based on different choices of models
  #here model can only be chosen from c("additive","dominant","recessive","overdominance")
  if (model == 'dominant'){
    data_nomiss$X[data_nomiss$X==2]=1
  } else if (model == 'recessive'){
    data_nomiss$X[data_nomiss$X==1]=0
    data_nomiss$X[data_nomiss$X==2]=1
  } else if (model == "overdominance"){
    data_nomiss$X1=data_nomiss$X;
    data_nomiss$X=0;
    data_nomiss$X[data_nomiss$X1==1]=1
 # } else if (model == "2df"){
#    print('Not yet supported...')
  } # leave additive as it is
  return(data_nomiss)
}
