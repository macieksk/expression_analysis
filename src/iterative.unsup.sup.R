library(statmod)
library(limma)
library(cluster)


perform.iteration.fit<-function(selDataMatrix, labels, design, 
        iteration.state=list(),
        recreate.genes.data=FALSE,
        dont.fit=FALSE
){
 within(iteration.state,{
    sample.labels<-labels
    design <- design
    #dataMatrix <- exprs(all.lumi.vst.norm)
    #To speed up the processing and reduce false positives, remove the unexpressed and un-annotated genes
    #presentCount <- detectionCall(all.lumi.vst.norm)
    #selDataMatrix <- dataMatrix[presentCount > 0,-to.remove]
    if (recreate.genes.data||!exists("probeList",inherits=FALSE)){
       probeList <- rownames(selDataMatrix)
       sampleList <- colnames(selDataMatrix)
       geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
       geneName <- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
    }
    if (!dont.fit){
      print(c(fit.time=system.time(
          fit <- lmFit(selDataMatrix, design,method="robust",maxit=80)
      )))
      #fit <- eBayes(fit,robust=TRUE)	
      fit$genes <- data.frame(ID=probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)
    }
  })
}

last<-function(l)l[[length(l)]]

ith.from.end<-function(l,i)tail(l,i)[[1]]

col.clustering.from.hmplot<-function(hmplot,k)cutree(as.hclust(hmplot$colDendrogram),k=k)
row.clustering.from.hmplot<-function(hmplot,k)cutree(as.hclust(hmplot$rowDendrogram),k=k)

last.clustering<-function(iterstates,k,i=1)col.clustering.from.hmplot(ith.from.end(iterstates,i)$hmplot.last,k)

hclust.model.labels<-function(dendrogram,k=4){
    hc<-as.hclust(dendrogram)
    ct<-cutree(hc,k=k)
    #or<-order(as.numeric(unlist(dendrapply(dendrogram,function(x)attr(x,"label")))))
    ct
}

hclust.model.matrix.from.dendrogram<-function(dendrogram,k=4){
    stopifnot(k>1)
    hc<-as.hclust(dendrogram)
    ct<-cutree(hc,k=1:k)
    hclust.model.matrix.from.dendmatrix(ct)
}

hclust.model.matrix.from.dendmatrix<-function(ct){
    k<-NCOL(ct)
    #Minimal new cluster computation
    min.cl<-1
    for (i in 2:k){
       ft<-factor(interaction(ct[,i-1],ct[,i]))    
       sl<-sub("^([^.]+).*","\\1",levels(ft))
       duplev<-which(duplicated(sl))
       duplev<-sort(unique(c(duplev,which(duplicated(sl,fromLast=TRUE)))))
       dupvec<-ft%in%levels(ft)[duplev]
       tb<-table(factor(ft[dupvec]))   
       mc<-sub("^.*[.](.*)$","\\1",names(tb)[which.min(tb)])
       min.cl<-c(min.cl,as.integer(mc))    
    }
    
    #return(ct)
    #print(min.cl)
    res<-ct
    for(i in 1:k){
                  res[ct[,i]==min.cl[i],i]<-1
                  res[ct[,i]!=min.cl[i],i]<-0
    }
    colnames(res)<-paste("cluster.",colnames(res),sep="")    
    #or<-order(as.numeric(unlist(dendrapply(dendrogram,function(x)attr(x,"label")))))
    res
}



hclust.model.matrix.from.factor<-function(labels){
    f<-factor(labels)   
    dm<-sapply(1:length(levels(f)),function(l)as.integer(f)>=l)
    colnames(dm)<-1:length(levels(f))
    hclust.model.matrix.from.dendmatrix(dm)    
}


iterative.clustering<-function(...){
    iteration.states<-list()
    iteration.states[[1]]<-perform.first.iteration.clustering(...)
    perform.iterative.clustering(iteration.states,startit=2,...)    
}

iterative.clustering.from.ttab<-function(ttab,selDataMatrix,sample.labels=NULL,...){    
    #stopifnot(!is.null(sample.labels))
    iteration.state<-perform.iteration.fit(selDataMatrix,sample.labels,design=NULL,
                     recreate.genes.data=TRUE,dont.fit=TRUE)
    
    iteration.state$ttab<-ttab
    iteration.states<-list()
    iteration.states[[1]]<-iteration.state
    perform.iterative.clustering(iteration.states,selDataMatrix,startit=2,...)    
}


perform.first.iteration.clustering<-function(selDataMatrix,
    sample.labels=NULL,   
    design=hclust.model.matrix.from.factor(sample.labels),
    N=1500,
    gene.order=c("abscoef","fstat","pamlike","coef.alltop","coef.fstat"),
    p.abscoef=2,
    coef.fstat.pvweight=300,
    pval.threshold=0.001,
    lfc.threshold=log2(2),
    obligatory.genes=NULL,
    ...)
{    
    stopifnot(!is.null(sample.labels))
    iteration.state<-perform.iteration.fit(selDataMatrix,sample.labels,design)
    iteration.state<-within(iteration.state,{
        k<-length(levels(factor(sample.labels)))
        ttab<-select.top.table.from.fit(fit,N,k,pval.threshold,lfc.threshold,gene.order,
                    p.abscoef, coef.fstat.pvweight,
                    obligatory.genes=obligatory.genes)
        print(c(selected.genes=NROW(ttab)))
        ttab<-cbind(ttab,entrezID=unlist(lookUp(ttab$ID,'lumiHumanAll.db','ENTREZID')))
        print("log(Adj.PVal) summary:")
        print(summary(log10(ttab$adj.P.Val)))
     })
    iteration.state
}



perform.iterative.clustering<-function(iteration.states,selDataMatrix,
    N=1500, k=6, 
    gene.order=c("abscoef","fstat","pamlike","coef.alltop","coef.fstat"),
    p.abscoef=2,
    coef.fstat.pvweight=300,
    plotname.core=match.arg(gene.order),
    maxit=14,
    startit=2,
    avail.colors=c("violet","lightblue","green","purple","red","yellow"),
    sample.labels.marks=NULL,
    pval.threshold=0.001,
    lfc.threshold=log2(2),
    explained.variance.for.clustering.thr=1,
    obligatory.genes=NULL,
    ...)
{

    stopifnot(startit>=2)
    gene.order<-match.arg(gene.order)
    iteration.states<-iteration.states[1:(startit-1)]
    
    pname.core<-plotNameOut(paste("heatmap_top",N,"_k",k,"_",plotname.core,sep=""))    
    
    if(N=="ALL"||is.null(N)) N<-Inf
    
    heatmap.clust.expression<-expression({              
         itnum.str<-sub(" ",0,format(itnum-1,flag="0",width=2))
         pname<-paste(pname.core,"_iterations_",itnum.str,".png",sep="")
         res<-with(ttab[1:min(N,NROW(ttab)),],{
             #coefs<-unclass(fit)["coefficients"]$coefficients[ID,]
             genedata<-selDataMatrix[ID,]#-coefs[,1]
             print(paste("Plotting",length(ID),"genes X",NCOL(genedata),"samples"))
             
             heatmap.clust(pname,genedata,geneSymbol,dend.both=TRUE,
                           #nv=length,
                           nv=function(d){d2<-d^2
                                          cs<-cumsum(d2/sum(d2))
                                          head(which(cs>=explained.variance.for.clustering.thr),1)
                                          },
                           main.title="log(FoldChange) relative to median gene expression",
                           sample.labels= if (!is.null(sample.labels)) avail.colors[sample.labels]
                                          else NULL,
                           sample.labels.marks=sample.labels.marks,
                           gene.colors=get.color.for.gene(geneSymbol),
                           normalize.for.plot=FALSE,
                           normalize=FALSE,
                           do.plot=TRUE,
                           rect.column=k, rect.row=2*k+2)
        })
    })
    
   for (itnum in startit:maxit){
    ttab<-iteration.states[[itnum-1]]$ttab
    sample.labels<-iteration.states[[itnum-1]]$sample.labels
    
    eval(heatmap.clust.expression)
    
    if (!("hmplot"%in%names(iteration.states[[itnum-1]]))){
          iteration.states[[itnum-1]]$hmplot<-list()
    }
    iteration.states[[itnum-1]]$hmplot<-res$hmplot
    iteration.states[[itnum-1]]$hmplot.last<-res$hmplot      
        
    design.matrix<-hclust.model.matrix.from.dendrogram(iteration.states[[itnum-1]]$hmplot$colDendrogram,k=k)
    sample.labels<-hclust.model.labels(iteration.states[[itnum-1]]$hmplot$colDendrogram,k=k)

    if (itnum>startit)
        if (all(sample.labels==iteration.states[[itnum-1]]$sample.labels)){
          print(paste("All sample.labels equal. Iteration STOP at itnum:",itnum))    
          return(iteration.states)
        }
        
    iteration.states[[itnum]]<-perform.iteration.fit(selDataMatrix,sample.labels,design.matrix,
                                                     iteration.state=iteration.states[[itnum-1]])

    iteration.states[[itnum]]<-within(iteration.states[[itnum]],{

        ttab<-select.top.table.from.fit(fit,N,k,pval.threshold,lfc.threshold,gene.order,
                            p.abscoef, coef.fstat.pvweight,
                            obligatory.genes=obligatory.genes)

        print(NROW(ttab))
        #ttab<-cbind(ttab,entrezID=unlist(lookUp(ttab$ID,'lumiHumanAll.db','ENTREZID')))
        #print("AbsLogFoldChange summary:")
        #print(summary(abs(ttab$logFC)))
        #print("LogFoldChange summary:")
        #print(summary(ttab$logFC))
        print("log(Adj.PVal) summary:")
        print(summary(log10(ttab$adj.P.Val)))
     })
        
    }
    itnum<-itnum+1
    
    ttab<-iteration.states[[itnum-1]]$ttab
    sample.labels<-iteration.states[[itnum-1]]$sample.labels  
    
    eval(heatmap.clust.expression)
    
    if (!("hmplot"%in%names(iteration.states[[itnum-1]]))){
           iteration.states[[itnum-1]]$hmplot<-list()
    }
    iteration.states[[itnum-1]]$hmplot<-res$hmplot
    iteration.states[[itnum-1]]$hmplot.last<-res$hmplot    
        
    iteration.states
}
    

select.top.table.from.fit<-function(fit,N,k,
                            pval.threshold=0.001,lfc.threshold=log2(2),
                            gene.order=c("abscoef","fstat","pamlike","coef.alltop","coef.fstat"),
                            p.abscoef=2,
                            coef.fstat.pvweight=300,
                            adjust="fdr",
                            robust=TRUE,
                            obligatory.genes=NULL
                            ){
        stopifnot(k>=2)
        gene.order<-match.arg(gene.order)  
        if (is.infinite(N)) stopifnot(gene.order=="abscoef")                     
        
        cfnum=2:k        
        ttab<-modifTopTableF(eBayes(fit[, cfnum],robust=robust), 
                       adjust.method=adjust, 
                       p.value=pval.threshold,
                       lfc=lfc.threshold,
                       number=NROW(fit$coefficients),
                       sort.by="F", #c("t","F")[(k>2)+1]
                       obligatory.genes=obligatory.genes
                       )
                       
        #ttab<-ttab[ttab$adj.P.Val<pval.threshold,]
        #print(NROW(ttab))
        
        take.top.best.equally.for.each.column<-function(ttab,
                modt.pvals.sortednames, N,
                tocheck.seq=seq(N/k,N,by=30)){
            first.best<-function(nn)unique(as.vector(t(sapply(modt.pvals.sortednames,function(sn)head(sn,nn)))))
            
            checked.sizes<-sapply(tocheck.seq,function(nn)length(first.best(nn)))
            ind<-which.max(checked.sizes>=N)
            if (checked.sizes[ind]<N) 
                ind<-length(checked.sizes)
            sel.set<-first.best(tocheck.seq[ind])
            sel.id<-which(ttab$ID%in%sel.set)
            if(length(sel.id)>0)
                ttab<-rbind(ttab[sel.id,],ttab[-sel.id,])  
            ttab
        }
        
        #if (k<=2){
        #       coefs<-ttab[,"logFC",drop=FALSE]
        #} else {
        coefs<-ttab[,paste("cluster",cfnum,sep="."),drop=FALSE]
        #}
        if (gene.order=="abscoef"){
            or<-order(rowSums(abs(coefs)^p.abscoef),decreasing=TRUE)
            ttab<-ttab[or,]    
        } else 
        if (gene.order=="pamlike"){
            #modt.pvals<-fit$p.value[ttab$ID,2:k,drop=FALSE]
            modt.pvals<-sapply(cfnum,function(cf){
                eBayes(fit[ttab$ID, cf],robust=robust)$F.p.value
                })
            rownames(modt.pvals)<-rownames(ttab)
            modt.pvals.sortednames<-lapply(1:NCOL(modt.pvals),function(i)rownames(modt.pvals)[order(modt.pvals[,i])])
            ttab<-take.top.best.equally.for.each.column(ttab,modt.pvals.sortednames,N)
          
        } else 
        if (gene.order=="coef.fstat"){
            #modt.pvals<-fit$p.value[ttab$ID,2:k,drop=FALSE]
            coefs<-abs(coefs)
            modt.pvals<-sapply(cfnum,function(cf){
                coefs[,cf-1]-eBayes(fit[ttab$ID, cf],robust=robust)$F.p.value*coef.fstat.pvweight
                })
            rownames(modt.pvals)<-rownames(ttab)
            modt.pvals.sortednames<-lapply(1:NCOL(modt.pvals),function(i)rownames(modt.pvals)[order(modt.pvals[,i])])
            ttab<-take.top.best.equally.for.each.column(ttab,modt.pvals.sortednames,N)
          
        } else 
        if (gene.order=="coef.alltop"){
            coefs<-abs(coefs)
            coefs.sortednames<-lapply(1:NCOL(coefs),function(i)rownames(coefs)[order(coefs[,i],decreasing=TRUE)])

            ttab<-take.top.best.equally.for.each.column(ttab,coefs.sortednames,N)
        }
        if (!is.null(obligatory.genes)){       
            obl.id<-which(ttab$ID%in%obligatory.genes)
            if(length(obl.id)>0)
                ttab<-rbind(ttab[obl.id,],ttab[-obl.id,])
        }
        ttab  
}



### Genes Table generation

get.cluster.factor<-function(rowDendrogram,k=4){
# here returned order is the same as in original data, not as in dendrogram
    ctr<-cutree(as.hclust(rowDendrogram),k=2:k)    
    ret<-if(k>2) apply(ctr,1,function(row)paste(row,collapse=".")) else ctr
    #stats<-addmargins(table(ret))
    #names(stats)<-1:length(stats)
    #print(t(stats))
    ret<-factor(ret)
    print(summary(ret))
    ret
}


give.intensity.design.matrix<-function(fit){
#From hierarchical desing matrix to a matrix with indicator columns
 tree.nodes<-unique(apply(fit$design,1,function(r)paste(r,collapse="")))
 clust.design<-sapply(1:NCOL(fit$design),
     function(i)sort(tree.nodes[substr(tree.nodes,i,i)=="1"])[1]
 )
 df<-as.data.frame(do.call(cbind,strsplit(clust.design,"")),stringsAsFactors=FALSE)
 colnames(df)<-colnames(fit$design)
 sapply(df,function(x)as.logical(as.integer(x)))
}



genes.cluster.table<-function(iter,N=NULL,k=3,k.genes=2*k+2,ttab=NULL,hmplot=NULL,
        foldchange.formulas=NULL,
        refseqid=TRUE){
 fit<-iter$fit
 if (is.null(ttab))
    ttab<-iter$ttab
 if (is.null(hmplot))
     hmplot<-iter$hmplot.last
 if (is.null(N)) N<-NROW(ttab)
 print(N)
 #fname.l[[N]]<-fname
 with(ttab[1:N,],{
  coefs<-unclass(fit)$coefficients[ID,]
  indes<-give.intensity.design.matrix(fit)
  avgs<-sapply(1:NCOL(coefs),function(i)rowSums(coefs[,indes[,i],drop=FALSE]))
  #avgs[,2]<-avgs[,1]+avgs[,2]
  avg.intensity<-2^avgs  
  colnames(avg.intensity)<-paste("Avg.intensity.",colnames(coefs),sep="")
  #fc<-ifelse(logFC>0,2^logFC,-2^(-logFC))
  rowclf<-get.cluster.factor(hmplot$rowDendrogram,k=k.genes)
  #colclf<-get.cluster.factor(hmplot$colDendrogram,k=k)

  RefSeqID<-rep(NA,N)
  if (refseqid)
     RefSeqID<-nuID2RefSeqID(ID, lib.mapping='lumiHumanIDMapping')

  df<-data.frame(Official.Gene.Symbol=geneSymbol,RefSeqID=RefSeqID,
                Marker.Set=get.marker.type.for.gene(geneSymbol),
               Gene.Cluster=rowclf,
               #Fold.Change=fc,
               avg.intensity,
               FDR.q.value=adj.P.Val)
 
  fcstats<-NULL
  if (!is.null(foldchange.formulas)){
    fcs<-lapply(foldchange.formulas,function(fcform){
         fc<-rowMeans(avg.intensity[,fcform[[1]],drop=FALSE])/rowMeans(avg.intensity[,fcform[[2]],drop=FALSE])
         ifelse(fc<1,-1/fc,fc)
    })
         
    df<-cbind(df,as.data.frame(fcs))
         
    #Add indicator columns
    quantiles<-c(0.1,0.5,0.9)
    fcstats<-lapply(fcs,function(fcss){
       tapply(fcss,rowclf,function(fc){
         #fivenum(fc)
         quantile(fc,p=quantiles)
       })
    })

    #Make it data.frame
    fcstats<-lapply(names(fcstats),function(fcsname){
       fcs<-fcstats[[fcsname]]
       fcs<-as.data.frame(do.call(rbind,fcs))
       colnames(fcs)<-paste("q",quantiles,sep="")#c("Min","1st","Median","3rd","Max")       
       #print(fcs)  
       cbind(cluster=rownames(fcs),fcs,FCset=fcsname)
    })
    names(fcstats)<-names(fcs)

    #Add indicator columns to df
    indic.df<-do.call(cbind,lapply(names(fcstats),function(fcsname){
        indf<-as.data.frame(matrix(0,ncol=3,nrow=length(rowclf)))
        colnames(indf)<-paste(fcsname,c("UP","DOWN","UPorDOWN"),sep=".")
        fcs<-fcstats[[fcsname]]
        indf[rowclf%in%fcs$cluster[fcs$q0.1<(-2)],2]<-1
        indf[rowclf%in%fcs$cluster[fcs$q0.9>( 2)],1]<-1
        indf[,3]<-indf[,1]+indf[,2]
        indf
    }))
    df<-cbind(df,indic.df)
        
    fcstats<-do.call(rbind,fcstats)    
    rownames(fcstats)<-NULL
    

  }
  
  df<-df[rev(hmplot$rowInd),]
  rownames(df)<-NULL
  df<-data.frame(RowNum=1:N,df)
  #print(all(order(df$Gene.Cluster)==hmplot[[300]]$rowInd))
  list(df=df,row.clusters=rowclf,col.clusters=colclf,fcstats=fcstats)
})
}
      
genes.cluster.table.save<-function(fname.core,iter,N=NULL,k=3,k.genes=2*k+2,...){
    fname.core<-paste(fname.core,"_top",N,"_k",k,sep="")
    df<-genes.cluster.table(iter,N=N,k=k,k.genes=k.genes,...)
    write.xlsx2(df$df,file=paste(fname.core,".xls",sep=""),quote=FALSE,row.names=FALSE,showNA=FALSE)   
    if (!is.null(df$fcstats))
       write.xlsx2(df$fcstats,file=paste(fname.core,"_fcstats.xls",sep=""),row.names=FALSE,showNA=FALSE)        
    invisible(df)
}

#Patient table
### Patients table generation
patient.cluster.table<-function(iter=NULL,hmplot=iter$hmplot,N,k,k.genes=2*k+2){
    N<-300
    fname<-plotName(paste("patients_table_cluster_based_on_top",N,"_svd1.csv",sep=""))
    colclf<-get.cluster.factor(mplot$colDendrogram,k=3)
    levels(colclf)[1]<-paste("Micro.Pap.new",levels(colclf)[1],sep=".")
    levels(colclf)[-1]<-paste("Conv.TCC.new",levels(colclf)[-1],sep=".")
    df<-data.frame(sample.id=colnames(selDataMatrix),original.label=cancer.label[-to.remove],new.cluster.label=colclf)
    df<-df[hmplot$colInd,]
    rownames(df)<-NULL
    df<-rbind(df,data.frame(sample.id=colnames(dataMatrix[,to.remove]),original.label=cancer.label[to.remove],new.cluster.label="Removed.from.analysis.bad.quality"))
    df<-data.frame(RowNum=1:NROW(df),df)
    #write.csv(df,file=fname,quote=FALSE,row.names=FALSE)
    df[sort(sample(1:NROW(df),10)),] 
}




#Cluster assesment

silhouette.for.hmplot<-function(hmplot,N,k,col=c("violet","lightblue","green","purple")[1:k],...){
    sl<-silhouette(cutree(as.hclust(hmplot$colDendrogram),k=k),hmplot$pi.dist)
    plot(sl,
     col=col,
     main=paste("Silhouette plot for k=",k,"clusters\n based on",N,"top genes"),...)
}

silhouette.for.iter<-function(iter,...){
    hmplot<-iter$hmplot.last
    silhouette.for.hmplot(hmplot,...)
}




### toptable reimplemented



modifTopTableF<-function (fit, number = 10, genelist = fit$genes, adjust.method = "BH", 
    sort.by = "F", p.value = 1, lfc = 0,
    obligatory.genes=NULL) 
{
    if (is.null(fit$coefficients)) 
        stop("Coefficients not found in fit")
    M <- as.matrix(fit$coefficients)
    rn <- rownames(M)
    if (is.null(colnames(M))) 
        colnames(M) <- paste("Coef", 1:ncol(M), sep = "")
    Amean <- fit$Amean
    Fstat <- fit$F
    Fp <- fit$F.p.value
    if (is.null(Fstat)) 
        stop("F-statistics not found in fit")
    if (!is.null(genelist) && is.null(dim(genelist))) 
        genelist <- data.frame(ProbeID = genelist, stringsAsFactors = FALSE)
    if (is.null(rn)) 
        rn <- 1:nrow(M)
    else if (anyDuplicated(rn)) {
        if (is.null(genelist)) 
            genelist <- data.frame(ID = rn, stringsAsFactors = FALSE)
        else if ("ID" %in% names(genelist)) 
            genelist$ID0 <- rn
        else genelist$ID <- rn
        rn <- 1:nrow(M)
    }
    sort.by <- match.arg(sort.by, c("F", "none"))
    adj.P.Value <- p.adjust(Fp, method = adjust.method)
    if (lfc > 0 || p.value < 1) {
        if (lfc > 0) 
            big <- rowSums(abs(M) > lfc, na.rm = TRUE) > 0
        else big <- TRUE
        if (p.value < 1) {
            sig <- adj.P.Value <= p.value
            sig[is.na(sig)] <- FALSE
        }
        else sig <- TRUE
        keep <- big & sig
        if (!is.null(obligatory.genes)){       
            obl.id <- rownames(fit)%in%obligatory.genes
            keep<- keep | obl.id
        }
        if (!all(keep)) {
            M <- M[keep, , drop = FALSE]
            rn <- rn[keep]
            Amean <- Amean[keep]
            Fstat <- Fstat[keep]
            Fp <- Fp[keep]
            genelist <- genelist[keep, , drop = FALSE]
            adj.P.Value <- adj.P.Value[keep]
        }
    }
    if (nrow(M) < number) 
        number <- nrow(M)
    if (number < 1) 
        return(data.frame())
    if (sort.by == "F") 
        o <- order(Fp, decreasing = FALSE)[1:number]
    else o <- 1:number
    if (is.null(genelist)) 
        tab <- data.frame(M[o, , drop = FALSE])
    else tab <- data.frame(genelist[o, , drop = FALSE], M[o,, drop = FALSE])
    tab$AveExpr = fit$Amean[o]
    tab <- data.frame(tab, F = Fstat[o], P.Value = Fp[o], adj.P.Val = adj.P.Value[o])
    rownames(tab) <- rn[o]
    tab
}

