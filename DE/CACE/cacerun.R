#
# 
# Author: nmatasci
# v 2.1.0
###############################################################################
# v 2.1.0
# Added '-a' option to process all columns
#
#v2.0.0 Revision
# TODO Change output format
# Nexus support input
# TODO Change output to png

#v1.4.0 Refactored fel loading table function
#
#v1.3.4
#Fixed (yet again) Phylip/Fel import support
#
#v 1.2
#Added dataout.csv file
#Fixed Phylip/Fel import support
#
#v 1.1
#Added Phylip/Fel import support
#
#Runs the R function ace (package ape) to estimate continuous ancestral characters.
#Returns 5 files for each character analyzed.
#character_plot.pdf: original phylogenetic tree with colored nodes and node names superimposed
#character_nodes.txt: estimated character values for the internal nodes. Names (numbers) are the same as in treeplot.pdf
#character_tree.txt: phylogenetic tree in newick format with named internal nodes
#character_out.txt: Analysis results: likelihood and sigma2 estimate
#dataout.csv: Observed and estimated values


#LIBRARIES
library(ape)

infile="infile"
intree="intree"

traits=1

#ARGUMENTS
args<-commandArgs(trailingOnly=TRUE)

	
if(length(args)>=2){
	infile<-args[1]
	intree<-args[2]
	args<-args[-1:-2]
	if(length(args)>0){
		if(args=='-a'){
			traits='-a'
		} else{
			traits<-make.names(args)
		}
	}else{
		traits<-1
	}
}


#VARIABLES
#tr: Phylogenetic tree
#ct: Characters (traits) table

#FUNCTIONS

#saveplot()
#writes the tree to the file treeplot.pdf

saveplot<-function(results,tt, cv,trait){
	pdf(file=paste(trait,"_plot.pdf",sep=""))
			
	pc<-function (x){
		return(rgb(colorRamp(c("blue","red"),space="Lab")(x),maxColorValue=255))
	}	
	minmax<-range(c(cv,results$ace), na.rm=TRUE, finite=TRUE)
	lm<-pretty(minmax,5)
	tip.color<-(cv-minmax[1])/(minmax[2]-minmax[1])
	node.color<-(results$ace-minmax[1])/(minmax[2]-minmax[1])
	
	tt$tip.label<-gsub("."," ",tt$tip.label,fixed=TRUE)
	plot(tt,type="p",use.edge.length=FALSE,label.offset=.2)
	title(main=trait)
	tiplabels(NULL,pch=16,col=pc(tip.color),cex=1.5)
	nodelabels(NULL,pch=16,col=pc(node.color),cex=1.5,frame="none")
	nodelabels(text=tt$node.label,frame="none",bg=NULL,adj=c(1.5,-1),pch="",cex=0.5)
	
	#plots the color key using 20 levels and ~5 tick marks
	par(new=TRUE)
	par(fig=c(0,0.2,0.5,1))
	image(matrix(1:20/20,1,20),axes=FALSE,col=pc(1:20/20))
	axis(2,at=(0:(length(lm)-1)/(length(lm)-1)),labels=lm, par(cex=0.7,las=1),lwd=0.7,tcl=-0.3)
	dev.off()
}

#Print the results
print.res<-function(x,i) {
	cat(
		paste(c("\tAncestral Character Reconstruction: ",i,"\n\nCall: ",x$call,"\n\n\tLog-likelihood: ",round(x$loglik,6),"\n\n\tSigma2: ",(tt<-round(x$sigma2,6))[1]," ",tt[2],"\n\nSee \"",i,"_nodes.txt\" for ancestral character estimates.\n"),collapse=""),
    	file=paste(i,"_out.txt",sep="")
	)
	
	
}
make.table<-function(results,tt,cv,trait){
	outtable<-as.data.frame(results$ace,row.names=row.names(results))
	colnames(outtable)='val'
	intable<-as.data.frame(cv,row.names=names(cv))
	colnames(intable)='val'
	outtable<-rbind(intable,outtable)
	return(outtable)	
}


fel<-function(infile=infile){
	ctt<-readLines(con=infile)
	nrowsncols<-as.numeric(strsplit(ctt[1],split="\\s+", perl=TRUE)[[1]][-1])
	ct<-t(as.data.frame(lapply(ctt[-1],function(x) splitrecord(x, ncols=nrowsncols[2]))))
	row.names(ct)<-gsub(" +$","",ct[,1])
	ct<-ct[,-1]
	colnames(ct)<-paste('V',1:nrowsncols[2],sep="")
	return(ct)
}

splitrecord<-function(row,ncols){
	tname<-substr(row,1,10)
	field<-strsplit(substring(row,11),split="  ")[[1]]
    md<-ncols-length(field)
 	field<-c(field,rep("",md))
	field[field==""]<-NA
	return(c(tname,field))
}





run<-function(i) {
	if(sum(is.na(as.numeric(ct[,i]) ) ) > sum( is.na(ct[,i]) ) ) {
		wm<-paste("Character ",i," is not a continuous trait.",sep="")
		warning(wm," Skipping.")
		cat(wm,"Skipped.\n",file=paste(i,"_warn.txt",sep=""))
		return(FALSE)
	}
	cv<-as.double(ct[,i])
	names(cv)<-row.names(ct)
	tt<-tr
	if(any(is.na(cv))) {
		tt<-drop.tip(tr,names(cv[is.na(cv)]))
		cv<-cv[!is.na(cv)]
	}

	cv<-cv[match(tt$tip.label,names(cv))]
	results<-ace(cv,tt,type="c",CI=TRUE)	
	
	

	saveplot(results,tt,cv,i)
	
	write.tree(tt,file=paste(i,"_tree.txt",sep="")) 
	out<-data.frame(cbind(results$ace, results$CI95))
	colnames(out)<-c("Estimated value","lower 95CI","upper 95CI")
	row.names(out)<-tt$node.label
	write.csv(out,file=paste(i,"_nodes.txt",sep=""))
	write.csv(make.table(results,tt,cv,i),file="dataout.csv")
	print.res(results,i)
	return(results)
	
}




#INPUT
#FILE FORMATS
#Tree: newick
#Data: csv with header
#or Felsenstein


#TREE
t<-readLines(intree,n=1,warn=FALSE)

if( length(grep('#nexus',t,ignore.case = TRUE) ) ){
	tr<-read.nexus(file=intree)
} else {
	tr<-read.tree(file=intree)
}
tr$tip.label<-make.names(tr$tip.label,allow_=FALSE,unique=TRUE)

#If internal node names are not present, numbers are assigned
if(is.null(tr$node.label) || any(is.na(tr$node.label) || any(tr$node.label==""))) 	tr$node.label<-c((length(tr$tip.label)+1):((length(tr$tip.label)+1+tr$Nnode-1)))

#DATA
#
#DATA
h<-scan(file=infile,nlines=1,what='raw',sep=',',quiet=TRUE)
if(length(h)==1){	
	ct<-fel(infile)
} else {
	ct<-read.csv(file=infile,row.names=1,head=TRUE,as.is=FALSE)
}



row.names(ct)<-make.names(row.names(ct),allow_=FALSE,unique=TRUE)

if(length(traits)==1){
	
 if(is.numeric(traits)){
	traits<-colnames(ct)[traits]
}
else if(traits=='-a'){
	traits=colnames(ct)
}
}

#FUNCTION CALL	
graphics.off()
r<-unlist(apply(as.array(traits),1,function(x) run(x)))
