
################################################
# OPEN, PARSE AND CLOSE (OPC) WORKSPACES
# This function is used for pasrsing all the workspaces used in analysis.

opc <- function(y) {

  lapply(y, function(x){
    wksp<-openWorkspace (x)
    return(parseWorkspace(wksp, name=1))
    closeWorkspace(wksp)
  }
  )
}

############ meantans for transformation of flowframes

meantrans <- function(transformId, a=0){
  t = new("transform", .Data = function(x) x+a)
  t@transformationId = "transformId"
  t
}
############
############ logtans for transformation of flowframes
logtrans <- function(transformId){
  t = new("transform", .Data = function(x) log10(x))
  t@transformationId = "transformId"
  t
}


################################################
# NAMES AND NODES (NN)
# Creates a new object to store the names & nodes used in flowDiv

nn<- function(x, nod, nod2, use.beads)
{

  names=unlist(lapply(x, function(y)sampleNames(y)))
  nodesample=lapply(x, function(y) getData(y, nod))

  if (use.beads){
    beads=lapply(x, function(y) getData(y, nod2))
    summa=lapply(nodesample, function(x)summary(x))
    means=rapply(summa, function(x)x["Mean",], how="list")
    all.means=data.frame(lapply(means, function(x) do.call(mapply, c(rbind, x))))
    grand.mean=apply(all.means, 2, mean)
    diff.mean=sweep(all.means, 2, grand.mean)
    my=diff.mean
  }

  else {
    cols=unique(unlist(lapply(nodesample, colnames)))
    my=matrix(ncol=length(cols), nrow=length(unlist(names)))
    my=ifelse(is.na(my),0, 0)
  }
  # Mean trans
  selec=gsub("\\.","-",cols)
  transVec=apply(my, 1, function(x)lapply(x, function(w)meantrans(a=w)))
  meansCorrection=lapply(transVec, function(x)transformList(selec, x))
  transdata=lapply(nodesample, function(k)mapply(function(x,y)transform(x, y), k, meansCorrection))
  nodesample=transdata
  # Log trans
  transVeclog=apply(my, 1, function(x)lapply(x, function(w)logtrans(w)))
  logTransformation=lapply(transVeclog, function(x)transformList(selec, x))
  transdata2=lapply(nodesample, function(k)mapply(function(x,y)transform(x, y), k, logTransformation))
  nodesample2=transdata2
  ####
  return(list(names=names, nodesample=nodesample2))

}


################################################
# SEQ FUNCTION
# Gets thte sequences of each selected channel
seq_fun<-function(object1, selection, NBINS, fixed){

  # Maximum
  d=lapply(object1, function(x) lapply(x, function(y)lapply(selection, function(x)max(exprs(y)[,x]))))
  e=lapply(d, function(x) do.call(mapply, c(rbind, x)))
  e= do.call(rbind, e)
  maximum=apply(e, 2, max)
  # Minimum
  a=lapply(object1, function(x) lapply(x, function(y)lapply(selection, function(x)min(exprs(y)[,x]))))
  b=lapply(a, function(x) do.call(mapply, c(rbind, x)))
  b= do.call(rbind, b)
  minimum=apply(b, 2, min)
  # Tabled max and min
  max.min.table=cbind(minimum,maximum)
  rownames(max.min.table)<-selection
  if (!is.null(fixed)) max.min.table <- fixed
  # Cuts
  cuts<- apply(max.min.table, 1, function(x) seq(from = x[1], to = x[2], length = NBINS+1))

  return(list(cuts=cuts, max.min.table=max.min.table))
}
################################################
# BINPOS FUNCTION
# Returns coordinates of a specific bin

binpos<-function(table, binposition){

  return(which(table==table[binposition], arr.ind = T)[which(which(table==table[binposition], arr.ind = F)==binposition),])

}
################################################
# GRAPHS FUNCTION
# Plots
graphs <- function(tubed, names, NBINS, selection, myseqs, breaks, clusterscols, clusters, binsize, psize, dimension){


  newframe=expand.grid(selection, selection)
  newframe2=as.matrix(newframe)
  newframe2= gsub("-",".", newframe2)

  breaks=data.frame(breaks)
  myseqs=lapply(myseqs, function(x)data.frame(x))
  rownames(myseqs$max.min.table)<-gsub( "-", ".", rownames(myseqs$max.min.table))

  gplots=mapply(function(x,y){
    clusters2=cbind(clusters, clust=clusterscols)
    clusters2=clusters2[!duplicated(clusters2[,c(x,y)]),]

    if(x==y){
      g1=ggplot(tubed, aes_string(x=x)) + geom_density(colour="blue", fill="blue", alpha=0.2)+
        theme_bw() +
        theme(panel.border = element_blank(),
              axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+

        scale_x_continuous(breaks=breaks[,x], limits= c(myseqs$max.min.table[x, "minimum"], myseqs$max.min.table[x, "maximum"]))

    }


    else{

      g2=ggplot(tubed)+
        aes_string(y=y, x=x)+
        geom_hex(bins = 100, na.rm = T) +
        scale_fill_gradientn("", colours = rev(rainbow(3, end = 4/6)))+


        scale_x_continuous(breaks=breaks[,x],  limits= c(myseqs$max.min.table[x, "minimum"], myseqs$max.min.table[x, "maximum"])) +
        scale_y_continuous(breaks=breaks[,y],  limits= c(myseqs$max.min.table[y, "minimum"], myseqs$max.min.table[y, "maximum"])) +
        theme_bw() +
        theme(panel.border = element_blank(),
              axis.line = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        theme(legend.position="none")+
        geom_point(data=clusters2, aes_string(y=y, x=x), colour=clusters2$clust, alpha=.5, shape=15, size=psize)

    }

  }, newframe2[,1], newframe2[,2], SIMPLIFY = F)

  grid=arrangeGrob(gplots, nrow=length(selection), ncol=length(selection))
  ggsave(filename=names, plot=grid, device="jpeg", units="cm", width=dimension, height=dimension)

}


################################################
# Best number of bins (Freedman.Diaconis)
Freedman.Diaconis<-function(x){
  h = 2*IQR(x)*(length(x)^(-(1/3)))
  bins=(max(x)-min(x))/h
  return(bins)
}
