


### Comparing clusters TODO
compare.dendrograms<-function(dend1,dend2){
    dend1<-as.dendrogram(dend1)
    dend2<-as.dendrogram(dend2)
    cf1<-get.cluster.factor(dend1,k=4)
    cf2.1<-get.cluster.factor(dend2,k=4)[names(cf1)]
}


clust.sim.matrix<-function(clustmat){
    same.clust.info<-lapply(1:NCOL(clustmat),function(i){
     ot<-outer(clustmat[,i],clustmat[,i],"==")
     storage.mode(ot)<-"integer"
     ot
    })    
    same.clust.info<-Reduce(function(a,b)a+b,same.clust.info)/(NCOL(clustmat)-1)
    same.clust.info
}

heatmap.clust.summary<-function(iteration.states,k,n=8,same.clust.info=NULL,...){
  pname<-plotName(paste("patient_cluster_analysis.png",sep=""))    
  
  if (is.null(same.clust.info)){  
    clustmat<-sapply(head(tail(iteration.states,n),n-1),function(it)cutree(as.hclust(it$hmplot.last$colDendrogram),k=k))
    same.clust.info<-clust.sim.matrix(clustmat)
  }

  itstate<-last(iteration.states)                            
  cnames<-itstate$sampleList
                            
  colnames(same.clust.info)<-cnames                            
                            
  sample.labels<-itstate$sample.labels

  heatmap.clust(pname,same.clust.info,
                   cnames,
                   pca.type="none",
                   plot.only.raw=TRUE,
                   dend.both=TRUE,
                   #nv=length,
                   #nv=function(d){d2<-d^2
                   #               cs<-cumsum(d2/sum(d2))
                   #               head(which(cs>=0.999),1)
                   #               },
                   main.title="Similarity between clusters",
                   sample.labels=c("violet","lightblue","green","purple")[sample.labels],
                   #sample.labels.marks=orig.sample.labels=="Micropap",
                   gene.colors=c("violet","lightblue","green","purple")[sample.labels],
                   centralize=FALSE, 
                   centralize.for.plot=FALSE,
                   normalize=FALSE,
                   normalize.for.plot=FALSE,                   
                   do.plot=TRUE,
                   use.L1.for.clustering=FALSE,
                   absolute.gene.clustering=FALSE,
                   rect.column=NULL,
                   rect.row=NULL,
                   sample.method="ward",
                   genes.method="ward",
                   plot.extreme.quantile=0.0,
                   col.panel=colorpanel(50, "transparent", "red", "blue"),
                   row.cluster.order.by.vars=FALSE,
                   density.info="histogram",
                   #hm.width=15*NCOL(same.clust.info),
                   #term.delta=10^(-5),
                   #hm.width=1200,
                   #hm.height=1000,
                   #trace=FALSE,
                   #max.iter=99
                   ...
             )
     invisible(same.clust.info)
    
}

