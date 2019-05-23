glr_plot=function(glr,snpName=glr$SNP[which.max(glr$GLR)],kcut=c(1/8,8,100,1000),col=c('red','red','green','blue'),legend_cex=1,...){
  plot(as.numeric(glr$bp),log10(as.numeric(glr$GLR)),xlab='Location (bp)',ylab=expression('log'[10]*'(GLR)'),...)
   #mark different k threshold
  for (k_inx in 1:length(kcut)){
    k=as.numeric(kcut[k_inx]);col_k=col[k_inx]
    abline(h=log10(k),col=col_k)
  }
  #mark top SNP, add legend
  if (!is.null(snpName)){
    for (snp in snpName){
      text(glr$bp[glr$SNP==snp],log10(as.numeric(glr$GLR[glr$SNP==snp]))+0.1,snp)
    }
  }
  #legend
  legend_txt=c(paste0('1/',kcut,' LI'))
  legend("topleft",legend_txt,col=col,lty=c(1.5),bty="n",cex=legend_cex)
}
