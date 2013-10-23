#
# 
# Author: nmatasci
# v 1.7.0
###############################################################################
#v.1.7.0
# Refactored fel loading table function
#
#v1.6.4
#Fixed Phylip/Fel import support
#
#v1.5
#Added dataout.csv file
#Added Phylip/Fel import support
#v1.4
#Runs the R function ace (package ape) to estimate discrete ancestral characters.
#Returns 4 files for each character analyzed.
#character_plot.pdf: original phylogenetic tree with likelihood pies and node names superimposed
#character_nodes.txt: likelihood values for the internal nodes. Names (numbers) are the same as in treeplot.pdf
#character_tree.txt: phylogenetic tree in newick format with named internal nodes
#character_out.txt: likelihood value and rate estimate
#dataout.csv: Observed and estimated values


#LIBRARIES
library(ape)


##VARIABLES
#tr: Phylogenetic tree
#ct: Characters (traits) table
IP=0.1
infile="infile"
intree="intree"
traits=1

#ARGUMENTS
args<-commandArgs(trailingOnly=TRUE)


if(length(args)==0){

} else{
	infile<-args[1]
	intree<-args[2]
	args<-args[-1:-2]
	iin<-which(args=='-i')
	if(length(iin)>0) {
		IP=args[iin+1] #removes the IP value
		traits<-args[-iin:-(iin+1)]
	}
	
	if(length(traits)>0){
		traits<-make.names(traits)
	}else{
		traits<-1
	}
}


##FUNCTIONS

#traitfreq(x)
#x: named vector of different character states
#returns: data frame with characters as colums and frequency of 1
#Used to plot pies at the leaves of the tree
traitfreq<-function(x){
	df<-matrix(0, length(x), nlevels(x))
	df[cbind(1:length(x),x)]<-1
	dimnames(df)<-list(names(x),sort(levels(x)))
	
	return(df)
}


#saveplot()
#writes the tree to the file treeplot.pdf
saveplot<-function(results,tt, cv,trait){
	pdf(file=paste(trait,"_plot.pdf",sep=""))
	pc<-rainbow(nlevels(cv))
	tt$tip.label<-gsub("."," ",tt$tip.label,fixed=TRUE)
	plot(tt,type="p",use.edge.length=FALSE,show.tip.label = TRUE,label.offset=1)
	title(main=trait)
	tiplabels(pie=as.matrix(traitfreq(cv)),piecol=pc,cex=.5)
	
	
	nodelabels(pie=results$lik.anc,piecol=pc,cex=.5,frame="none")
	nodelabels(text=tt$node.label,frame="none",cex=.5,bg=NULL,adj=c(1.5,-1),pch="")
	legend("topleft",sort(levels(cv)),fill=pc)
	dev.off()
}

make.table<-function(results,tt,cv){
	outtable<-as.data.frame(results$lik,row.names=as.character(tt$node.label))
	outtable<-rbind(as.data.frame(traitfreq(cv)),outtable)
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


#Print the results
print.res<-function(x,i) {
	ratemat <- x$index.matrix
	names<-dimnames(x$lik.anc)[[2]]
	ratemat<-cbind(names,ratemat)
	ratemat<-rbind(c(" ",names),ratemat)
	ratemat[is.na(ratemat)]<-"."
	ratemat<-rbind(matrix(paste("",ratemat,""),dim(ratemat),byrow=TRUE),"\n")
	estim <- data.frame(1:length(x$rates), round(x$rates, 6), round(x$se, 6))
	names(estim) <- c("rate index", "estimate", "std-err")
	
	cat(
			paste(c("\tAncestral Character Reconstruction: ",i,"\n\nCall: ",x$call,"\n\n\tLog-likelihood: ",round(x$loglik,6),"\n\nRate index matrix:\n",ratemat,"\nParameter estimates:\n\n",paste(names(estim),collapse="\t"),"\n",paste(estim, collapse="\t"),"\n\nSee \"",i,"_nodes.txt\" for ancestral character states likelihoods.\n"),collapse="")
			,file=paste(i,"_out.txt",sep="")
	)
	
	
}

run<-function(i) {
	
	cv<-factor(ct[,i])
	if(nlevels(cv)<=1) {
		wm<-paste("Character ",i," is uniform across all taxa.",sep="")
		warning(wm," Skipping.")
		cat(wm,"Skipped.\n",file=paste(i,"_warn.txt",sep=""))
		return(FALSE)
	}
	names(cv)<-row.names(ct)
	tt<-tr
	
	if(any(is.na(cv))) {
		tt<-drop.tip(tr,names(cv[is.na(cv)]))
		cv<-cv[!is.na(cv)]
	}
	cv<-cv[match(tt$tip.label,names(cv))]
	
	results<-ace(cv,tt,type="d",CI=TRUE,ip=IP)
	
	graphics.off()
	saveplot(results,tt,cv,i)
	write.tree(tt,file=paste(i,"_tree.txt",sep="")) 
	write.csv(file=paste(i,"_nodes.txt",sep=""),as.data.frame(results$lik.anc, row.names=tt$node.label))
	write.csv(make.table(results,tt,cv),file="dataout.csv")
	print.res(results,i)
	return(results)
	
}


#INPUT
#FILE FORMATS
#Tree: newick
#Data: csv with header


#TREE
t<-readLines(intree,n=1,warn=FALSE)
if( length(grep('#nexus',t,ignore.case = TRUE))  ){
	tr<-read.nexus(file=intree)
} else {
	tr<-read.tree(file=intree)
}
tr$tip.label<-make.names(tr$tip.label,allow_=FALSE,unique=TRUE)

#If internal node names are not present, numbers are assigned
if(is.null(tr$node.label) || any(is.na(tr$node.label) || any(tr$node.label==""))) tr$node.label<-c((length(tr$tip.label)+1):((length(tr$tip.label)+1+tr$Nnode-1)))


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

if(length(traits)==1 && is.numeric(traits)){
	traits<-colnames(ct)[traits]
}

#FUNCTION CALL
r<-unlist(apply(as.array(traits),1,function(x) run(x)))