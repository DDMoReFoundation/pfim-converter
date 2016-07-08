
rm(list=ls(all=TRUE))

under.unix<-!(version$os=='Microsoft Windows' ||
			version$os=='Win32' || version$os=='mingw32')

library(nlme)
library(deSolve)
library(numDeriv)
options(expressions=10000,max.deparse.lines=10000,max.deparse.length=10000, deparse.max.length=10000)
if(!under.unix) windows.options(record=TRUE) 
subjects<-NULL

directory<-"DIRECTORY_REPLACE"
directory.program<-"DIRECTORY_PROGRAM_REPLACE"

#############################

PFIM<-function(model.file="stdin.r") {
	cat("Directory for PFIM4.0 programs:",directory.program,"\n")
	cat("Directory for project:",directory,"\n")
	source(file.path(directory,model.file))
	source(file.path(directory,file.model))
	# download the file "modelcreated.r" linked to the R function create_formED()
	
	if(length(grep("CreateModel_PKPDdesign.r",readLines(file.path(directory,file.model),warn = F)))!=0) {
		source(file.path(directory,"model_created.r"))}
	
	#test to accept stdin file of version PFIM 3.0
	
	#names.datax and names.datay objects
	ngraph<-2
	vec<-c("x","y")
	err1<-tryCatch(names.data_test<-lapply(1:ngraph,function(ngraph,x,vec) get(x), x=paste("names.data",vec[ngraph],sep=""), vec=vec), error=function(e) 4)
	
	if(length(which(is.null(unlist(err1))))>=1 || length(which(is.null(unlist(err1)==4)))>=1)
	{
		names.datax<-rep("Time",nr)
		names.datay<-rep("Concentration",nr)
	}
	
	
	#option object
	err2<-tryCatch(get("option"), error=function(e) 4)
	if(err2==4) option<-1
	
	#covariate.model
	err3<-tryCatch(get("covariate.model"), error=function(e) 4)
	if(err3==4) covariate.model<-F
	
	#IOV
	err4<-tryCatch(get("n_occ"), error=function(e) -4)
	if(err4==-4) n_occ<-0
	
	#covariate_occ.model
	err5<-tryCatch(get("covariate_occ.model"), error=function(e) 4)
	if(err5==4) covariate_occ.model<-F
	
	#choice of population, idividual or bayesian Fisher information matrix
	err7<-tryCatch(get("FIM"), error=function(e) 4)
	if(err7==4) FIM<-"P"
	
	if (run=="EVAL"){
		source(file.path(directory,file.model))
		if (modelform=="DE")
		{
			if(option==2)
			{ if (FIM!="P") {stop("Individual and bayesian matrices are only computed using option 1")}
				else{
					if (covariate.model==T || n_occ>1) stop("You can not use option 2 with covariate or inter-occasion variability currently")
					else source(file.path(directory.program,"EQPfim4.0op2.r"))
					out()
				}
			}
			else
			{
				source(file.path(directory.program,"EQPfim4.0op1.r"))
				out()
				
			}
			
		}
		else #modelform =AF
		{
			if(option==2)
			{ if (FIM!="P") {stop("Individual and bayesian matrices are only computed using option 1")}
				else
				{if (covariate.model==T || n_occ>1) {stop("You can not use option 2 with covariate or inter-occasion variability currently")}
					else {source(file.path(directory.program,"Pfim4.0op2.r"))}
					out()}
			}
			else
			{
				source(file.path(directory.program,"Pfim4.0op1.r"))
				out()
				
			}
		}
	}
	
	else {
		
		#FEDOROV-WYNN
		if (algo.option=="FW") {
			if (FIM!="P" && option==2) {stop("Individual and bayesian matrices are only computed using option 1")}
			source(file.path(directory.program,"algofedorov4.0.r"))
			
			if (previous.FIM!="")
			{cat("\n",'Warning: Previous information was taken into account only for the best one group protocol and not for the Fedorov-Wynn algorithm',"\n \n")}
			out.fedorov(modelfile=model.file,directory=directory,
					directory.program=directory.program)
		}
		
		else #SIMPLEX
		{
			
			
			if (modelform=="DE")
			{
				if (option==1)
					source(file.path(directory.program,"EQPfimOPT4.0op1.r"))
				else {
					if (FIM!="P") {stop("Individual and bayesian matrices are only computed using option 1")}
					else{
						if (covariate.model==T || n_occ>1) {stop("You can not use option 2 for Simplex algorithm with covariate or inter-occasion variability currently")}
						else { source(file.path(directory.program,"EQPfimOPT4.0op2.r")) }
					}
				}
			}
			else
			{
				if (option==1)
					source(file.path(directory.program,"PfimOPT4.0op1.r"))
				else {
					if (FIM!="P") {stop("Individual and bayesian matrices are only computed using option 1")}
					else{
						if (covariate.model==T || n_occ>1) {stop("You can not use option 2 for Simplex algorithm with covariate or inter-occasion variability currently")}
						else {source(file.path(directory.program,"PfimOPT4.0op2.r")) }
					}
				}
			}
			
			if (algo.option=="SIMP")
			{
				source(file.path(directory,file.model))
				source(file.path(directory.program,"algosimplex4.0.r"))
				x<-out.simplex()
				print(x)
				#remove(subjects,inherits = FALSE)
			}
			else
				cat("Error in the specification of the optimisation algorithm.","\n")
			
		}
	}
}
