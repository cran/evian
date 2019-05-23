multiLine_plot <-
  function (bpstart=0, bpend=1000000000, dList,  title=NULL, showmaxlr=3, kcut=NULL, pl='linear'
            ,ylim=c(-0.5,10),color=c('violet','green','red','blue'),markSNP=NULL,round=2,legend_cex=1)
  {
    #dList need to be a list with at least one column named 'SummaryStats'
    # In this column should also be a list that contains at least the following columns: (and in such order)
    #mle, maxlr, AF,SNP, bp, effect, ref,lower and upper limits
    #markSNP need to be a set of SNP ID that matched the snp column in dList
    # We only need the summaryStats column to do the mulit-line plot
    # we allow input of kcut, as a single numeric value or a vector as threshold. If it is a vector, the smallest value will be used as cutoffs, but the other values will be color-coded
    # if kcut is not provided, it will use the smallest cut-off from the output of the previous logit/lm function
    summaryStats=do.call(rbind,dList[,'SummaryStats'])
    if (is.null(kcut)){
      kcut=dList[,'k_cutoff'][[1]] # The cut-off are the same for all SNPs
    }
    kcutCol=which(dList[,'k_cutoff'][[1]]==min(kcut)) # So this finds the position for the designed k-threshold

    # If a vector is inputted, then the smallest value shall be the threshold for k
    # ---defaults:  maxlr=3, k=32, model="all"
    #    if (showmaxlr < 1) { showmaxlr = 3 }

    dframe_inRange=summaryStats[summaryStats$bp>=bpstart&summaryStats$bp<=bpend,]
    # we plot (-0.5,10) for Beta, (0,10) for OR
    if (pl=='logit'){
	    ylim=c(0,10); ylab='Odds Ratio';threshold=1
    } else {
    	ylab='Beta';threshold=0
    }
    plot(0,0,main=title,xlab='Location (bp)',ylab=ylab, xlim=range(dframe_inRange$bp),ylim=ylim)
    abline(h=threshold)
    skipCol=which(names(dframe_inRange)=='lo_1')-1 #skip the other summary Stats before the LIs
    # Plot from the highest k-value to the smallest
    for (i in 1:length(kcut)){ # we plot it grey if the chosen cut-off
      signColor=color[length(kcut)-i+1]
      boundCol=(length(kcut)-i+1)*2+skipCol+c(-1,0)
      lineColor=c('grey',signColor)[as.numeric(dframe_inRange[,(2*kcutCol-1+skipCol)]>threshold | dframe_inRange[,(2*kcutCol+skipCol)]<threshold)+1]
      segments(dframe_inRange$bp,dframe_inRange[,boundCol[1]],dframe_inRange$bp,dframe_inRange[,boundCol[2]],col=lineColor)
    }

    signAtKcutoff = dframe_inRange[dframe_inRange[,(2*kcutCol-1+skipCol)]>threshold | dframe_inRange[,(2*kcutCol+skipCol)]<threshold,]
    signAtKcutoff=signAtKcutoff[order(signAtKcutoff$maxlr,decreasing = T),]
    points(signAtKcutoff$bp,signAtKcutoff$mle,pch='-') #label the mle

    #if there's no significant SNP on all boundaries
    # plot all in grey and don't mark them
    num_of_hits = nrow(signAtKcutoff)
      if (num_of_hits>0 & is.null(markSNP)){
        text(signAtKcutoff$bp,signAtKcutoff[,ncol(signAtKcutoff)],signAtKcutoff$SNP,srt=90,pos=4,offset=0.3,cex=0.75)
      } else if (!is.null(markSNP)) {
        snpIndex=which(dframe_inRange$SNP %in% markSNP)
        text(dframe_inRange$bp[snpIndex],dframe_inRange[snpIndex,ncol(dframe_inRange)],dframe_inRange$SNP[snpIndex],srt=90,pos=4,offset=0.3,cex=0.75)
      }


    legend_txt=c(paste0('> 1/',kcut[1],' LI'),paste0('1/',kcut,' LI'))

    if (num_of_hits > 0 & showmaxlr > 0) {
      if (num_of_hits < showmaxlr) { showmaxlr = num_of_hits }
      legend.txt2=''
      for (i in 1:showmaxlr) {
        hit = paste ("max LR =", round(signAtKcutoff$maxlr[i],round), "at", signAtKcutoff$snp[i], "(", signAtKcutoff$bp[i] ,")")
        legend.txt2 = c(legend.txt2, hit)
      }
      legend("topright",8,legend.txt2,col=par("col"),lty=0,bty="n",cex=legend_cex)
    }
    legend("topleft",legend_txt,col=c("grey",color),lty=c(1.5),bty="n",cex=legend_cex)

  }
