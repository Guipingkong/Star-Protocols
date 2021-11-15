


####============================================================================================================
##  signed v unsigned GTEX7 beta10 ----------------------
#╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╗╔╦╗╔═╗╔═╦╗╔═╦═╦╦╦╦╗╔═╗╔╗═╦╗╔═╦╗╔╦╦╗
options(stringsAsFactors=F);library(colorout);rm(list=ls());ls()#╚═╣║ ╚╣║¯\_(•_•)_/¯║╚╣╔╣╔╣║║║║╚╣
options(menu.graphics=F)  #═╦╩╬╝╚╬╬╚╗═╚╩╗═╩╠╬╝╔╚╗╬╚║╣══╣╦╬╬╗╠╔╗╔╣╣╝╝╣╠╔╠╚╔╔═╦╩╬╣╦╣╔╚╬╦╣╩╬╚╩╗╣╝╚╠╣
library('adds') #╔╦╣═╩╚╣║╔╔╣╦═║║╔╗║╔╚╔╣╩╚╚╦╣║╩╔╦║║ ╚╩╣╚╚╣║╣╚╩╔╦╩╚╦╚╩╣╬╝╚╗╔╝╬╚╝ ╔╣═╦╦╦╩╠╔╠╗╔╝╚═╗╩║
# devtools::install_github("ks4471/addR") #╗╣╠═╩╠╣╠╬═╚╬╗╩╩═║╚╝╝╣╠╗╗╠╔║╩╬╠╝╣╬╔╬╬╚╦╝╔╗╩╠╚╝╠═╝╝╦╔═╚╠╝╣║
#╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩═╩╝╚═╝╩═╩╝╚═╩══╩═╩═╩═╩╝╩═╩╝╚═╩═╩╩═╝


projdir='~/Dropbox/PROJ/icl/giovanni.llaria/proteo/'

Load(paste0(projdir,'/dtb/ref/deconvolution.results.MEDIAN.Rdata'))


####===================================================================================================
##   log2 plots   -------------------------------------------


# icon=condi[2]
for(icon in condi){
   lprogr(icon,condi,T)
  holder=mranks[,grepl(icon,colnames(mranks))]
  dummy=geno[,grepl('fdr|log2[.]ratio',colnames(geno))&grepl(icon,colnames(geno)),drop=F]
  cutoff_log2r=min(abs(dummy[dummy[grepl('fdr',colnames(dummy))]<0.05,grepl('log2[.]ratio',colnames(dummy)) ]),na.rm=T)
  holder=rmerge(holder, dummy )

  holder=holder[apply(holder,1,function(x){sum(is.na(x))<3}),]
    str(holder)

  minmaxx=c(floor(min(holder[,!grepl('fdr|log2[.]ratio',colnames(holder))],na.rm=T)),ceiling(max(holder[,!grepl('fdr|log2[.]ratio',colnames(holder))],na.rm=T)))
  minmaxy=c(0,ceiling(-log10(min(
    holder[,grepl('fdr',colnames(holder))]
    ,na.rm=T))))

  min(abs(holder[,4]),na.rm=T)   ## mean log2 ratio is 4th column "by construction" (sort of)

    for(irep in 1:3){
      pdat=holder[,c(irep,ncol(holder))]
      pdat=pdat[complete.cases(pdat),]

		if(irep==1){
	        # plot(psig[,1], -log10(psig[,2]), pch = 16, col=rgb(0.4,0,0,alpha=0.8), frame.plot = F, xlim=minmaxx, ylim=minmaxy,xlab='log2 Ratio Heavy/Light',ylab='-log10 FDR',main=paste0(icon))
	        # points(potr[,1], -log10(potr[,2]), pch = 16, col=rgb(0,0,0,alpha=0.5))

		      psig=pdat[pdat[,2]<=0.05 & abs(pdat[,1])>=cutoff_log2r,]
		      colnames(psig)=c('a','b')
		        summary(psig)
		      potr=pdat[!rownames(pdat)%in%rownames(psig),]
		      colnames(potr)=c('a','b')
		        summary(potr)

		}
		if(irep!=1){
	        # points(psig[,1], -log10(psig[,2]), pch = 16, col=rgb(0.4,0,0,alpha=0.8))
	        # points(potr[,1], -log10(potr[,2]), pch = 16, col=rgb(0,0,0,alpha=0.5))
	        ander=pdat[pdat[,2]<=0.05 & abs(pdat[,1])>=cutoff_log2r,]
	        	colnames(ander)=c('a','b')
	        psig=rbind(psig,ander)

	        ander=pdat[!rownames(pdat)%in%rownames(psig),]
	        	colnames(ander)=c('a','b')
	        potr=rbind(potr,ander)

		}


    }

	if(icon=='SNA.vs.Sham'){
	pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/pretty/log10_pval.MEDIAN.',icon,'.logYscale.pdf'),height=8,width=8)
	    plot(psig[,1], -log10(psig[,2]),pch=16,col=rgb(0.4,0,0,alpha=0.8),xlab='log2 Ratio',ylab='-log10 FDR',main=icon,log='y',ylim=c(2e-1,90),xlim=c(-6.5,6.5))#,ylim=c(0.1,ceiling(max(-log10(psig[,2])))))
	    points(potr[,1], -log10(potr[,2]),pch=16,col=rgb(0,0,0,alpha=0.5))
	    abline(v=0.58,col='dodgerblue',lty=2,lwd=2)
	    abline(v=-0.58,col='dodgerblue',lty=2,lwd=2)
	dev.off()
	}

	if(icon=='DCA.vs.Lam'){
	pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/pretty/log10_pval.MEDIAN.',icon,'.logYscale.pdf'),height=8,width=8)
	    plot(psig[,1], -log10(psig[,2]),pch=16,col=rgb(0.4,0,0,alpha=0.8),xlab='log2 Ratio',ylab='-log10 FDR',main=icon,log='y',ylim=c(2e-1,90),xlim=c(-5,5))#,ylim=c(0.1,ceiling(max(-log10(psig[,2])))))
	    points(potr[,1], -log10(potr[,2]),pch=16,col=rgb(0,0,0,alpha=0.5))
	    abline(v=0.58,col='dodgerblue',lty=2,lwd=2)
	    abline(v=-0.58,col='dodgerblue',lty=2,lwd=2)
	dev.off()
	}

	if(icon=='Central.vs.Peripheral'){
	pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/pretty/log10_pval.MEDIAN.',icon,'.logYscale.pdf'),height=8,width=8)
	    plot(psig[,1], -log10(psig[,2]),pch=16,col=rgb(0.4,0,0,alpha=0.8),xlab='log2 Ratio',ylab='-log10 FDR',main=icon,log='y',ylim=c(2e-1,90),xlim=c(-10.5,10.5))#,ylim=c(0.1,ceiling(max(-log10(psig[,2])))))
	    points(potr[,1], -log10(potr[,2]),pch=16,col=rgb(0,0,0,alpha=0.5))
	    abline(v=0.58,col='dodgerblue',lty=2,lwd=2)
	    abline(v=-0.58,col='dodgerblue',lty=2,lwd=2)
	dev.off()
	}
}

	 #  if(icon=='SNA.vs.Sham'){
	 #  	ygap=c(22,39)


		# Library('plotrix')
		# pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/log10_pval.',icon,'.breaks.pdf'),height=12,width=7)
		# 	    gap.plot(psig[,1], -log10(psig[,2]), gap=ygap,pch = 16, col=rgb(0.4,0,0,alpha=0.8),breakcol=rgb(1,1,1),ylim=c(0,ceiling(max(-log10(psig[,2])))),xlab='log2 Ratio',ylab='-log10 FDR',main=icon)
		# 	    gap.plot(potr[,1], -log10(potr[,2]), gap=ygap,pch = 16, col=rgb(0,0,0,alpha=0.5),add=T,breakcol=rgb(1,1,1))

		# 		axis.break(2, ygap[1], breakcol="black", style="slash")
		# 		axis.break(4, ygap[1], breakcol="black", style="slash")

		# 		# axis(2, at=ygap[1]-6)
		# 		# axis(2, at=ygap[2])


		# dev.off()

	 #  }

	 #  if(icon=='DCA.vs.Lam'){
	 #  	ygap=c(51,78)

		# Library('plotrix')
		# pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/log10_pval.',icon,'.breaks.pdf'),height=12,width=7)

		# 	 #    gap.plot(psig[,1], -log10(psig[,2]), gap=ygap,pch = 16, col=rgb(0.4,0,0,alpha=0.8),breakcol=rgb(1,1,1),ylim=c(1e-5,ceiling(max(-log10(psig[,2])))),xlab='log2 Ratio',ylab='-log10 FDR',main=icon,log='y')
		# 	 #    gap.plot(potr[,1], -log10(potr[,2]), gap=ygap,pch = 16, col=rgb(0,0,0,alpha=0.5),add=T,breakcol=rgb(1,1,1))

		# 		# axis.break(2, ygap[1], breakcol="black", style="slash")
		# 		# axis.break(4, ygap[1], breakcol="black", style="slash")

		# 	    plot(psig[,1], -log10(psig[,2]),pch=16,col=rgb(0.4,0,0,alpha=0.8),xlab='log2 Ratio',ylab='-log10 FDR',main=icon,log='y',ylim=c(0.0001,90))#,ylim=c(0.1,ceiling(max(-log10(psig[,2])))))
		# 	    points(potr[,1], -log10(potr[,2]),pch=16,col=rgb(0,0,0,alpha=0.5))


		# dev.off()
	 #  }

	 #  if(icon=='Central.vs.Peripheral'){
	 #  	ygap=c(5.8,9.5)

		# pdf(paste0('~/Dropbox/PROJ/icl/giovanni.llaria/proteo/out/img/log10_pval.',icon,'.breaks.pdf'),height=12,width=7)
		# 	    gap.plot(psig[,1], -log10(psig[,2]), gap=ygap,pch = 16, col=rgb(0.4,0,0,alpha=0.8),breakcol=rgb(1,1,1),ylim=c(0,ceiling(max(-log10(psig[,2])))),xlab='log2 Ratio',ylab='-log10 FDR',main=icon)
		# 	    gap.plot(potr[,1], -log10(potr[,2]), gap=ygap,pch = 16, col=rgb(0,0,0,alpha=0.5),add=T,breakcol=rgb(1,1,1))

		# 		axis.break(2, ygap[1], breakcol="black", style="slash")
		# 		axis.break(4, ygap[1], breakcol="black", style="slash")

		# dev.off()
	 #  }









Library('plotrix')
x <- c(1:5, 6.9, 7)
y <- 2^x
from <- 33
to <- 110































