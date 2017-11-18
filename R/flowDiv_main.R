#' @title Cytometric Diversity Indices from Gated Data
#'
#' @description Uses multidimensional contingency tables to calculate ecological diversity indices from gated populations.
#' @author Bruno M.S. Wanderley, Maria Victoria Quiroga, Andre M. Amado, Fernando Unrein
#'
#' @usage flowDiv(myworkspaces, gate.name=NULL, ialpha="invsimpson", ibeta="bray", do.plot=FALSE,
#'
#' static=FALSE, dilutions=NULL,  pmax=0.05, env=NULL, use.beads=FALSE, beads=NULL,
#'
#' transform.hell=FALSE, dimension=20, psize=3, autotrans=FALSE)
#'
#'
#' @export flowDiv
#' @return A list containing alpha index, beta matrix and Pielou's indices for each cytogram.
#'
#' @import flowWorkspace
#' @import flowCore
#' @import vegan
#' @import ggplot2
#' @importFrom grDevices adjustcolor rainbow
#' @importFrom graphics points text
#' @importFrom methods new
#' @importFrom stats IQR kmeans na.omit
#' @importFrom utils capture.output select.list
#' @importFrom gridExtra arrangeGrob
#' @param myworkspaces  A list containing the paths to FlowJo workspaces or GatingSet which are meant to be analyzed. More than one workspace can be analyzed at the same time. Workspaces should contain .fcs files (versions 2.0 or 3.0) with its original names.
#' @param gate.name Name of the gate to be analyzed. Must be a single-valued string. The gate should be named exactly the same in all samples from the workspaces.
#' @param ialpha Method used to calculate alpha diversity index of cytograms (i.e. the degree of similarity between two cytograms). Should be one of "shannon", "simpson" or "invsimpson", as for vegan::diversity function. Default is "invsimpson".
#' @param ibeta Method used to calculate beta diversity index of cytograms. Should be one of  the sixteeen avaiable for vegan::vegdist function. Default is “bray”.
#' @param do.plot Logical. It plots scatterplots of gated populations and displays grid lines corresponding to the limits of the bins.
#' @param static Logical. If TRUE, one must set minimum and maximum range for each channel individually.
#' @param dilutions A vector containing dilution factors for each sample in the workspaces. Its lenght is the same as the total number of samples meant to be analyzed and must follow the exact same order in which samples are presented to the function (the order of samples corresponds both to the order of workspace described in "myworkspaces" parameter and their order in each of these FlowJo workspaces). By default, it is assumed that all samples have no dilutions whatsoever (i.e. all dilutions factors equal 1).
#' @param pmax Maximum estimated p value for displayed variables on NMDS plot, as required by vegan::plot.envifit
#' @param env Enviromental matrix associated to the set of cytograms under analysis.
#' @param use.beads Logical. If “TRUE” , it centers all cytograms under a common point, based on the arithmetic mean of a standard region for all cytograms (usually beads regions), before proceeding to analysis. It is recommended to proceed this way only if samples were analyzed with different settings.
#' @param beads Name of the gate describing the standard region to be used (usually bead's regions). Necessary only if use.beads is set to “TRUE”. The beads gate should be named exactly the same in all samples from the workspaces.
#' @param transform.hell Logical. If TRUE, Hellinger transformation is applyed to data before analysis.
#' @param dimension Dimensions of plot devices, in centimeters.
#' @param psize Defines the plotting size of bins.
#' @param autotrans Logical. Argument to be passed to vegan::metaMDS function.

#' @references Li, W.K.W. (1997). Cytometric diversity in marine ultraphytoplankton. Limnology and Oceanography 42, 874–880.

#' @examples \dontrun{
#'
#'
#'# Analyzing a .xml FlowJo workspace.
#'
#' indexes <- flowDiv(myworkspaces = “my_workspace.xml”,gate.name = “my_gate_name”)
#' #assume that the FlowJo workspace is below the current directory
#'
#'
#'}


flowDiv<- function(myworkspaces, gate.name=NULL, ialpha="invsimpson", ibeta="bray", do.plot=FALSE, static=FALSE, dilutions=NULL,  pmax=0.05, env=NULL, use.beads=FALSE, beads=NULL, transform.hell=FALSE, dimension=20, psize=3, autotrans=FALSE){

  classes=sapply(myworkspaces, function(x) return(class(x)), simplify = T)
  gsets=which(classes=="GatingSet")
  if(length(gsets)>0){
    wksp1=myworkspaces[-gsets]
    gtsets=myworkspaces[gsets]
    wksp2<-opc(wksp1)
    wksp<-append(wksp2, gtsets)

  } else (wksp=opc(myworkspaces))
  nnn=nn(wksp, nod=gate.name,use.beads=use.beads, nod2=beads)

  ops<-unique(unlist(lapply(unlist(nnn$nodesample), colnames)))
  selection<-select.list(ops, c(1:length(ops)), multiple = T, title="Please select channels for use")

  fixed=NULL
  if (static==TRUE){
    fixeds=lapply(selection, function(x){
      message(paste("Please enter minimum value for", x))
      min=scan(nmax=1, quiet=T)

      message(paste("Please enter maximum value for", x))
      max=scan(nmax=1, quiet=T)

      return(c(minimum=min, maximum=max))
    } )

    fixed=do.call(rbind, fixeds)
    rownames(fixed)<-selection
  }


  binss=lapply(nnn$nodesample, function(x) lapply(x, function(y)lapply(selection, function(x)Freedman.Diaconis(exprs(y)[,x]))))
  suggested.bins=round(mean(unlist(binss)))

  message(paste("Suggested number of bins:", suggested.bins, "\n"), "How many bins do you want to use?")
  nbins=scan(nmax=1, quiet=T)

  myseqs <- seq_fun(nnn$nodesample, selection, nbins, fixed)

  hits<-lapply(nnn$nodesample, function(x) lapply(x, function(y)sapply(selection, function(x){
    cuts<-cut(exprs(y)[,x], myseqs$cuts[,x],include.lowest = T, labels=F)
    return(cuts)
  }, simplify=F)))
  hits=unlist(hits, recursive=FALSE)
  bin.levels=lapply(hits, function(x)lapply(x, function(y)factor(as.numeric(y), levels=1:nbins)))
  tables=lapply(bin.levels, function(x)table(x))

  if (!is.null(dilutions))tables=mapply(function(x,y)x*y, tables, dilutions, SIMPLIFY = F)

  matrices<-do.call(mapply, c(cbind, tables))

  colnames(matrices)<-as.character(c(1:nbins^length(selection)))
  rownames(matrices)<-unlist(nnn$names, recursive = F)


  if (do.plot==TRUE){


    expressions=lapply(unlist(nnn$nodesample), function(x)exprs(x)[,selection])
    tubed=lapply(expressions, function(x)data.frame(x))

    breaks=apply(myseqs$cuts, 2, function(x){
      bind=cbind(x[1:length(x)-1], x[2:length(x)]) # Bindinfg for processing on next statge
      brk=apply(bind, 1, mean)
      return(brk)
    }
    )

    set.seed(1)
    nmds=metaMDS(matrices,autotransform = autotrans, try=100)
    specs=na.omit(nmds$species)


    ccKM=cascadeKM(specs, 2 ,10)
    ngroups=names(which(ccKM$results["SSE",]==min(ccKM$results["SSE",])))
    ngroups=gsub(" groups","", ngroups)
    message(paste("Suggested number of clusters:", ngroups, "\n"), "Please enter the number of clusters to use:")
    nclusters=scan(nmax=1, quiet=T)

    km=kmeans(specs, nclusters)

    cbins=as.numeric(names(km$cluster))
    bb=sapply(cbins, function(x)binpos(tables[[1]], x), simplify = F)
    bb=do.call(rbind, bb)
    clusters=lapply(selection, function(x)breaks[bb[,x],x])
    clusters=t(do.call(mapply, c(rbind, clusters)))
    colnames(clusters)<-selection
    clusters<-data.frame(clusters)


    invisible(capture.output(mapply(graphs, tubed=tubed, names=nnn$names, MoreArgs = list(NBINS=nbins,clusterscols=km$cluster, selection=selection, myseqs=myseqs, breaks=breaks, clusters=clusters, psize=psize, dimension=dimension))))

    plot(nmds, type="n")
    points(specs, pch=19, col=adjustcolor(km$cluster, alpha.f = 0.3))
    text(nmds$points, labels = c(1:nrow(matrices)))
    if (!is.null(env)){
      fit <- envfit(nmds, env, perm = 999)
      plot(fit, p.max = pmax, col = "red")
    }





  }

  alpha<-apply(matrices, 1, function(x)diversity(x, index=ialpha))

  pielou<-apply(matrices, 1, function(x)diversity(x)/log(vegan::specnumber(x)))

  if(transform.hell){
    beta<-vegdist(decostand(matrices, "hell"), method=ibeta)
  }

  else {
    beta<-as.matrix(vegdist(matrices, method=ibeta))
  }

  return(list(Alpha=alpha, Pielou=pielou, Beta=beta, Bins=nbins, Channels=selection, Matrices=matrices))
}
