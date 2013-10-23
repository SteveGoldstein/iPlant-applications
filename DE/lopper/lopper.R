#!/usr/bin/env Rscript
# lopper.R: A tool to subset phylogenetic trees, trait matrices and taxa lists.
# Author: Naim Matasci (nmatasci@iplantcollaborative.org), Barb Banbury (darwininthesun@gmail.com)
#
# The contents of this file are subject to the terms listed in the LICENSE file you received with this code.
# Copyright (c) 2011, The Arizona Board of Regents on behalf of 
# The University of Arizona
#
###############################################################################

VERSION=1.0
tnrs.service<-'http://tnrs.iplantc.org/tnrsm-svc'
#Load libs. Used to prevent library loading if program cannot run.
load.libs<-function() {
	library(geiger)
	library(rjson)	
    library(RCurl) #Required to call TNRS	
}#End loab.libs


#Loads the input files
load.file<-function(file.path,file.format){	
#TO DO:  should add some automagic to guess the file type using count.fields
	file<-switch(file.format,
			csv = read.csv(file = file.path,row.names=1),
			csvx = read.csv(file=file.path,row.names=1, header=FALSE),
			nwk = read.tree(file=file.path),
			nex = read.nexus(file=file.path),
			lst = {tn_temp<-readLines(file.path);return(data.frame(tn_temp,row.names=make.names(tn_temp)))},
			)
	return(file)
} #End load.file


#Identifies whether an object is a tree
get.type<-function(file){
	if(class(file) == 'phylo') return('t') else return('d')
}#End get.type


#Extracts lists of taxa from a data.frame or a phylo object
extract.names<-function(file) {
	taxa.names<-switch(class(file),
			
			'phylo' = file$tip.label,
			'data.frame' = row.names(file),
	)

	return(as.character(taxa.names))
	
}#End extract.names


#Transforms _ . into blankspaces
fix.names<-function(taxa.names){
	taxa.names<-sapply(list(taxa.names),function(x) gsub("\\.|_",' ',x,fixed=FALSE))
	return(as.vector(taxa.names))
}#End fix.taxa

#TO DO: Identify abbreviations (optional)
#Some awesome code goes here
#that looks for clusters of taxa with the same initial

#Send to TNRS
send.TNRS<-function(taxa.names,authority="tropicos"){
	uri<-tnrs.service
	maxn=50 #maximum number of names to send at one time
	taxa.matrix<-matrix(c(taxa.names,rep(NA,maxn-length(taxa.names)%%maxn)),nrow=maxn) #the names are stored in a matrix, the last column is padded with NAs
	names<-unlist(apply(taxa.matrix,2,function(x) submit.TNRS(x,authority=authority,uri=uri))) #One column at the time, names are sent to TNRS
	if(any(is.na(names))){
		return(taxa.names)
	} else{ 
		return(names)
	}
	
}#End send.TNRS

submit.TNRS<-function(taxa.names,authority,uri){

	taxa<-paste(taxa.names[!is.na(taxa.names)],collapse=',')
	taxa<-gsub('&',' ',taxa) 
	url<-paste(uri,'/matchNames?retrieve=best&names=',curlEscape(taxa),sep='')
#	url<-paste(uri,'/matchNames?retrieve=best&names=',curlPercentEncode(taxa),sep='')
#	print(url)
	results<-fromJSON(getURL(url))
#	return(taxa.names)
	names<-sapply(results[[1]], function(x) c(x$nameSubmitted,x$acceptedName))
	names<-as.data.frame(t(names),stringsAsFactors=FALSE)
	names[names[,2]=="",2]<-names[names[,2]=="",1] #If TNRS did not return any accepted name (no match, or name is already accepted), the submitted nam eis returned
	return(names[,2])
	
}

#Substitute taxa name in the original data with the corrected names
subst.names<-function(file,taxa.names){
	if(class(file) == 'phylo'){
		file$tip.label<-taxa.names
	} else  {
		row.names(file)<-taxa.names
	} 
	return(file)
}#End subst names
	



#Matches the taxa names in file1 and file2
#Value: A list composed of $file1, $file2, $file1.not.file2, $file2.not.file1 and $matched
match.names<-function(file1,file2){
	type<-paste(get.type(file1),get.type(file2),sep="")
	results<-switch(type,
				tt = compare.trees(file1,file2),
				dt = swap.labels(compare.treedata(file2,file1)),
				td = compare.treedata(file1,file2),
				dd = compare.data(file1,file2)		
			)
	return(results)
}#End match.names

#Compares two trees
compare.trees<-function(phy1,phy2){
	
	a<-TreeNameCheck(phy1, phy2)
	
	if (a[1] == "OK") {
		a<-list()
		a$phy1.not.phy2<-NULL
		a$phy2.not.phy1<-NULL
	}
	r<-list(phy1, phy2, a$phy1.not.phy2,a$phy2.not.phy1)

	names(r)<-c("file1","file2","file1.not.file2","file2.not.file1")

	if (length(a$phy1.not.phy2) > 0) {
		r$file1<-drop.tip(phy1, a$phy1.not.phy2)
	} 
	if (length(a$phy2.not.phy1) > 0) {
		r$file2<-drop.tip(phy2, a$phy2.not.phy1)
	}

	r$matched<-r$file1$tip.label

	return(r)
}#End compare.trees

#Compares data
compare.data<-function(data1,data2){
	a<-DataNameCheck(data1, data2)
	if (a[1] == "OK") {
		a<-list()
		a$data1.not.data2<-NULL
		a$data2.not.data1<-NULL
	} 

	r<-list(data1, data2, a$data1.not.data2,a$data2.not.data1)
	names(r)<-c("file1","file2","file1.not.file2","file2.not.file1")
	
	if (length(a$data1.not.data2) > 0) {
		dn<-dimnames(data1)
		to.remove<-which(rownames(data1) %in% a$data1.not.data2)
		r$file1<-as.data.frame(data1[-to.remove,],row.names=dn[[1]][-to.remove])
		colnames(r$file1)<-dn[[2]]
		

	}	
	if (length(a$data2.not.data1) > 0) {
		dn<-dimnames(data2)
		to.remove<-which(rownames(data2) %in% a$data2.not.data1)
		r$file2<-as.data.frame(data2[-to.remove,],row.names=dn[[1]][-to.remove])
		colnames(r$file2)<-dn[[2]]
	}	
	
	r$matched<-row.names(r$file1)
	
	return(r)
}#End compare.data

#compare tree to data
compare.treedata<-function(tree,dt){

	a<-name.check(tree,dt)

	if (a[1] == "OK") {
		a<-list()
		a$Tree.not.data<-NULL
		a$Data.not.tree<-NULL
		
	}

	r<-list(tree, dt, a$Tree.not.data,a$Data.not.tree)
	names(r)<-c("file1","file2","file1.not.file2","file2.not.file1")
	
	if (length(a$Tree.not.data) > 0) {
		r$file1<-drop.tip(tree, a$Tree.not.data)
	}
	if (length(a$Data.not.tree) > 0) {
		dn<-dimnames(dt)
		to.remove<-which(rownames(dt) %in% a$Data.not.tree)
		r$file2<-as.data.frame(dt[-to.remove,],row.names=dn[[1]][-to.remove])
		colnames(r$file2)<-dn[[2]]
	}

	r$matched<-row.names(r$file2)
	return(r)
}#End compare.treedata


#Swaps the labels file1 and file2 if file1 is data nd file2 is a tree
#is needed so that tree-data and data-tree can use the same function
swap.labels<-function(x){
	r<-x
	r$file1<-x$file2
	r$file2<-x$file1
	r$file1.not.file2<-x$file2.not.file1
	r$file2.not.file1<-x$file1.not.file2
	return(r)
}#End swap.labels


#checks name concordance among data files
DataNameCheck<-function(data1, data2) {
	t1<-rownames(data1)[which(is.na(match(rownames(data1), rownames(data2))))]
	t2<-rownames(data2)[which(is.na(match(rownames(data2), rownames(data1))))]

	t<-list(t1, t2)
	names(t)<-c("data1.not.data2", "data2.not.data1")
	if (length(t1) == 0 && length(t2) == 0) {
		return("OK")
	}
	else {
		return(t)
	}
}#End DataNameCheck


#checks name concordance among tree files
TreeNameCheck<-function(phy1, phy2) {
	
	t1<-phy1$tip.label[which(is.na(match(phy1$tip.label, phy2$tip.label)))]
	t2<-phy2$tip.label[which(is.na(match(phy2$tip.label, phy1$tip.label)))]
	t<-list(t1, t2)
	names(t)<-c("phy1.not.phy2", "phy2.not.phy1")
	if (length(t1) ==0 && length(t2) == 0) {
		return("OK")
	}
	else {
		return(t)
	}	
}#End TreeNameCheck


#Renames taxa in a data structure
rename.taxa<-function(file,useTNRS){
	file.names<-fix.names(extract.names(file))
	if(useTNRS) file.names<-send.TNRS(file.names)
	return(subst.names(file,file.names))
}#End rename.taxa


#Write files to disk
write.file<-function(results,format,label){
	switch(format,
			txt = write(results,file=label),
			csv = write.csv(results,file=label),
			csvx = write.csv(results,file=label,col.names=FALSE),
			nwk = write.tree(add.dist(results),file=label),
			nex = write.nexus(results,file=label),
			lst = write(row.names(results),file=label),
			
			)
	
}#End write.file

add.dist<-function(tr){
	if(is.null(tr$edge.length)){
		tr<-compute.brlen(tr,1)
	}
	return(tr)
}

#Formats results for JSON output
format.json<-function(object){

	if(get.type(object) == 't'){
		return(write.tree(object,file=""))
	}
	else{
		if(class(object)=='character'){
			return(object)
		} else if(colnames(object)[1]=='tn_temp' && all(fix.names(object[,1]) == row.names(object) ) ) {
			return(row.names(object))		
		}else{
			object<-cbind(Species=row.names(object),object)
			return(as.data.frame( t(object) ) ) 
		}
	}
	
}#End format.json

#Initialization function
init<-function(...){
	#Default values. Custom output filename is disabled	
	auth="tropicos"
	mismatches=0
	output=NULL
	use.json=FALSE
	formats=c()
	output=c() #Disabled
	use.tnrs=0
	file1.useTNRS=FALSE
	file2.useTNRS=FALSE	
	file1.saveMM=FALSE
	file2.saveMM=FALSE
	
	
	sep = alist(
			'--help'=usage(),
			'-h'=usage(), #prints help
			'--version'=version(),
			'-v'=version(), #prints version
			'--formats'=formats<-get.args(cline,'--formats'),
			'-f'=formats<-get.args(cline,'-f'), #formats of file 1 and file 2
			'--mismatches'=mismatches<-cline[which(cline == '--mismatches')+1],
			'-m'=mismatches<-cline[which(cline == '-m')+1], #which mismatch file to save: 0=none, 1=of file 1, 2=of file 2, 3=both files
			'--auth'=auth<-cline[which(cline == '--auth')+1],
			'-a'=auth<-cline[which(cline == '-a')+1],
			'--output'=output<-get.args(cline,'--output'),
			'-o'= get.args(cline,'-o'),
			'--tnrs'= use.tnrs<-cline[which(cline == '-tnrs')+1],
			'-t'= use.tnrs<-cline[which(cline == '-t')+1],
			'--json'= use.json<-TRUE,
			'-j'= use.json<-TRUE #output is JSON string sent to stdout				
	)
	
	
	usage <- function(){ 
		text=paste("\n Usage: Rscript lopper.R -f FORMATS [OPTIONS] FILE1 FILE2\n",
				"Reduces a pair of files which contain taxa names (tree, trait matrix, list) to the common set of taxa.\n\n",
				"-h --help\tPrints this help.\n",
				"-v --version\tPrints the program's version.\n",
				"-f FORMAT --format FORMAT\tFormat of FILE1 and FILE2. If only one format is provided it is used for both files.\n",
				"\tCurrently supported formats are:\n",
				"\t\tcsv: table of traits in comma separated value format\n",
				"\t\tcsvx: table of traits in comma separated value without header format\n",
				"\t\tnwk: phylogenetic tree in newick format\n", # TODO: should think about multiple trees support
				"\t\tnex: phylogenetic tree in nexus format\n",
				"\t\tlst: list of taxa names in plain text.\n",
				"-m NUM --mismatches NUM: lists of mismatched taxa to save:\n\t0 = none (default)\n\t1 = FILE1\n\t2 = FILE2\n\t3 = both FILES.\n",
				"-a AUTHORITY --auth AUTHORITY: Authority for resolution of taxa name. Ignored if TNRS is not used\n",
				"-t NUM --tnrs NUM\tUse Taxonomic Resolution Service for:\n\t0 = none (default)\n\t1 = FILE1\n\t2 = FILE2\n\t3 = both FILES.\n",
				"\tCurrently supported authorities\n",
				"\t\t\"tropicos\" (default): Tropicos, botanical information system at the Missouri Botanical Garden\n",
#				"-o --output:\tName of the output files (default: FILE1_FILE2 and FILE2_FILE1). If only one name is given the suffixes \"1\" and \"2\" will be appended to the name.\n",
				"-j --json\tThe output is sent as a JSON string to STDOUT and no files are generated (options -m and -o are ignored).\n\n") # TODO: JSON export
		cat(text)
		q(status=0)
	}
	
	version<-function(){
		text=paste(
				" lopper 0.9\n",
				"Copyright (c) 2011, The Arizona Board of Regents on behalf of The University of Arizona\n",
				"This software is released under a BSD 3-clause license. For details see http://iplantcollaborative.org/opensource\n")
		cat(text)
		q(status=0)
	}
	
	#Helper function to retrun 2 duplicate arguments if only 1 is provided
	get.args<-function(cline,opt){
		kp<-which(cline == opt)
		if(is.na(cline[kp+2]) || cline[kp+2] %in% sep){
			return(c(cline[kp+1], cline[kp+1]))
		} else {
			return(c(cline[kp+1], cline[kp+2]))
		}
	}
	
	
	
	#Reads in the CL input
	cline<-commandArgs(TRUE)
	
	#Returns the version
	if (any(c('-v','--version') %in% cline)){
		version()	
	}
	
	#Returns the usage instruction in case it has been requested, no arguments werre given or no file format has been specified 
	if(is.null(cline) || length(cline)==0 ||  any(c('-h','--help') %in% cline) || !any(c('-f','-format') %in% cline) ){			
		usage()
	}
	
	#Gets the file names (last 2 CLI input)
	file1<-tail(cline,2)[1]
	file2<-tail(cline,1)
	if(file1 %in% sep || file2 %in% sep){
		cat("Please specify valid file names.\n")
		usage()
	}
	cline<-head(cline,-2)	
	
	#Evaluates the options
	for(i in cline[cline %in% names(sep)]){
		eval(sep[[i]])
	}
	
	#All is fine: the program can run.
	#Loads the libraries
	load.libs()
	
	if(use.tnrs == 1 || use.tnrs == 3){
		file1.useTNRS<-TRUE
	} 
	if(use.tnrs == 2 || use.tnrs == 3){
		file2.useTNRS<-TRUE
	}
	
	if(mismatches == 1 || mismatches == 3){
		file1.saveMM<-TRUE
	} 
	if(mismatches == 2 || mismatches == 3){
		file2.saveMM<-TRUE
	}
	
	
	
	return(
			list(
					file1=file1 ,
					file1.format= formats[1],
					file2 = file2,
					file2.format = formats[2],
					file1.useTNRS = file1.useTNRS,
					file2.useTNRS=file2.useTNRS,
					TNRS.auth =auth,
					file1.saveMM=file1.saveMM,
					file2.saveMM=file2.saveMM,
					file1.newname =NULL, #
					file2.newname=NULL,
					use.json=use.json
			)
	)
	
	
	
}#End init

main<-function(file1,file1.format,file2,file2.format,
		file1.useTNRS=TRUE,file2.useTNRS=TRUE,TNRS.auth="tropicos",
		file1.saveMM=FALSE,file2.saveMM=FALSE,
		file1.newname=NULL,file2.newname=NULL,
		use.json=FALSE) {
	
	# TODO: Add support for output filename
	file1.original<-load.file(file1,file1.format)
	file2.original<-load.file(file2,file2.format)
	file1.name<-strsplit(tail(strsplit(file1,'\\/')[[1]],1),'\\.')[[1]]
	file2.name<-strsplit(tail(strsplit(file2,'\\/')[[1]],1),'\\.')[[1]]
	
	file1<-rename.taxa(file1.original,file1.useTNRS)
	file2<-rename.taxa(file2.original,file2.useTNRS)
	
	results<-match.names(file1,file2)

	
	if(use.json){
		r<-lapply(results, function(x) format.json(x))
		return(cat(toJSON(r)))
	} else {
	
#		if(is.na(file1.newname) || is.na(file2.newname)){
#			file1.newname<-paste(head(file1.name,-1),"_",head(file2.name,-1),".",tail(file1.name,1),sep="")
#			file2.newname<-paste(head(file2.name,-1),"_",head(file1.name,-1),".",tail(file2.name,1),sep="")
#		} else {
#			type<-paste(get.type(file1),get.type(file2),sep="")
#			switch(type,
#					tt = {
#						file1.newname<-"tree1.txt"
#						file2.newname<-"tree2.txt"
#					},
#					dt = {
#						file1.newname<-"data.txt"
#						file2.newname<-"tree.txt"
#					},
#					td = {
#						file1.newname<-"tree.txt"
#						file2.newname<-"data.txt"
#					},
#					dd = {
#						file1.newname<-"data1.txt"
#						file2.newname<-"data2.txt"
#					}	
#			)
#			
#		}
		file1.newname<-"newfile1.txt"
		file2.newname<-"newfile2.txt"
		write.file(results$file1,file1.format,file1.newname)
		write.file(results$file2,file2.format,file2.newname)
		write.file(results$matched,"txt","matched.txt")
		if(file1.saveMM){
			write.file(results$file1.not.file2,"txt",paste(head(file1.name,-1),"_mismatches.txt",sep=""))
		}
		if(file2.saveMM){
			write.file(results$file2.not.file1,"txt",paste(head(file2.name,-1),"_mismatches.txt",sep=""))
		}
		return(results)
	}
	
}


do.call(main,init())
q(status=0)

