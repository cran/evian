densityPlot = function(dList,snpName,kcut=NULL,pl='linear',xlim=NULL,color=c('red','orange','green','blue'),round=2,legend_cex=1){
  snpRowIndex=which(do.call(rbind,dList[,'SummaryStats'])$SNP==snpName)
  dataToPlot=dList[snpRowIndex,]
  if (pl=='logit'){ xlab='Odds Ratio';nullLine=1} else {xlab='Beta'; nullLine=0}
  if (is.null(xlim)){ # We just plot a bit more than the largest LI
    loBound=min(nullLine-0.1,dataToPlot$SummaryStats[,ncol(dataToPlot$SummaryStats)-1]) #we want to show the null
    hiBound=max(nullLine+0.1,dataToPlot$SummaryStats[,ncol(dataToPlot$SummaryStats)])
    xlim=c(loBound,hiBound)
  }
  plot(dataToPlot$theta,dataToPlot$profile.lik.norm,xlim=xlim,type='l',xlab=xlab,ylab='standardized profile likelihood')
  abline(v=nullLine,col='grey',lty=2)
  #plot the LIs
  colorVector=NULL; cutoffVector=NULL
  for (i in 1:length(dataToPlot$k_cutoff)) {
    segments(x0=dataToPlot$SummaryStats[,paste0('lo_',i)],y0=1/dataToPlot$k_cutoff[i]
             ,x1=dataToPlot$SummaryStats[,paste0('hi_',i)],y1=1/dataToPlot$k_cutoff[i],col=color[i])
    colorVector=c(colorVector,color[i])
    cutoffVector=c(cutoffVector,paste0('1/',dataToPlot$k_cutoff[i],' LI (',round(dataToPlot$SummaryStats[,paste0('lo_',i)],round)
                                       ,', ',round(dataToPlot$SummaryStats[,paste0('hi_',i)],round),')'))
  }
  mle=round(dataToPlot$SummaryStats$mle,round)
  if (mle > nullLine){posLegend='topleft';pos=loBound+0.1} else {posLegend='topright';pos=hiBound-0.1}
  # write out the legend and some summary stats
  legend(posLegend,legend=cutoffVector,col = colorVector,lty=1,cex=legend_cex)
  text(pos,0.8,paste0('MLE at ',mle))
  text(pos,0.75,paste0('Max LR is ',round(dataToPlot$SummaryStats$maxlr,round)))
  #output robustFactor if it is not 1
  if (dataToPlot$SummaryStats$robustFactor!=1){
    text(pos,0.7,bquote(paste(hat(a)/hat(b),' = ',.(round(dataToPlot$SummaryStats$robustFactor,round)))))
  }
}
