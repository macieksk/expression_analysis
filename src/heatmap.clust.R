library(rpca)

heatmap.clust<-function(pname,genedata,geneSymbol,
                       dend.both=TRUE,
                       pca.type=c("robust.pca","svd","none"),
                       #robust.pca=TRUE,
                       #svd.pca=TRUE,                       
                       plot.only.raw=FALSE,
                       nv=length,
                       sample.labels=NULL,
                       sample.labels.marks=NULL,
                       gene.colors = rep(c("darkred",NA),NROW(genedata)+1)[1:NROW(genedata)],
                       main.title="log(FoldChange) relative to mean gene expression",
                       centralize=TRUE,
                       centralize.for.plot=TRUE,
                       normalize=FALSE,
                       normalize.for.plot=normalize,
                       do.plot=TRUE,
                       breaks=NULL,
                       rect.column=4, rect.row=2*rect.column+2,
                       lhei.colside=0.25,
                       save.sample.dist=TRUE,
                       save.gene.dist=FALSE,
                       clustering.distance=c("angular","L1","oneminus"),
                       absolute.gene.clustering=FALSE,
                       sample.method="ward",
                       density.info='density',
                       genes.method=sample.method,
                       plot.extreme.quantile=0.005,
                       col.panel=colorpanel(50, "lightgreen", "black", "red"),
                       row.cluster.order.by.vars=TRUE,
                       hm.width=1000+10*max(NCOL(genedata)-100,0),
                       hm.height=3400+6*(max(N-300,0)),
                       hm.hprop=c((700+max(800-N,0))/700*350/(800+max(N-800,0)), (700+max(800-N,0))/700*150/(800+max(N-800,0)) ,10),
                       hm.wprop=c(1.5,9),
		       #Predifned clusterings
                       hclust.row = NULL,
                       hclust.col = NULL,
                       #Predefined order
                       row.order = NULL,
                       col.order = NULL,
                       #Separators
                       colsep = NULL,
                       rowsep = NULL,
                       sepwidth = c(0.05, 0.05),
                       create.png=TRUE, #legacy option
                       close.png=TRUE,  #legacy option
                       create.device=ifelse(create.png,"png",""), 
                       close.device=close.png,  
                       eps.pointsize=70,
                       eps.pointtopixel=1400,
                       par.args=list(),
                       legend.fun=NULL, 
                       use.gpu=FALSE,                   
                       ...){

    clustering.distance<-match.arg(clustering.distance)
    if (genes.method=="ward") genes.method<-"ward.D2"
    if (sample.method=="ward") sample.method<-"ward.D2"
    pca.type<-match.arg(pca.type)
    robust.pca<-pca.type=="robust.pca"
    svd.pca<-pca.type=="svd"
    
    #par(mfcol=c(1,2))
    plot.genedata<-genedata
    genedata.centralized<-function(){genedata-rowMedians(genedata)}
    genedata.normalized <-function(){genedata/(rowIQRs(genedata)/1.349)}
    if (centralize) {
        genedata<-genedata.centralized()
        if (centralize.for.plot)
            plot.genedata<-genedata
    } else if (centralize.for.plot) {
            plot.genedata<-genedata.centralized()
    }
    if (normalize){
        genedata<-genedata.normalized()
        if (normalize.for.plot)
            plot.genedata<-genedata
    } else if (normalize.for.plot) {
            plot.genedata<-genedata.normalized()
    }

    if (pca.type=="robust.pca"){
         #sd<-svd(genedata,nu=nu,nv=nv)
         if (!use.gpu)
             tm<-system.time(pres<-rpca(genedata,...))
         else 
            tm<-system.time(pres<-rpca.gpucula(genedata,...))
         print(paste("rpca converged:",pres$convergence$converged,"in:"))
         print(tm)

         sd<-pres$L.svd
         #sd$v<-t(sd$vt)

         d2<-sd$d^2

         #print(cumsum(d2/sum(d2)))
         #print(sd$d)

         if (is.function(nv))
             nv<-nv(sd$d)
         nv<-min(nv,length(sd$d))
         nu<-nv

         print(c(robust.rank=length(d2),rank.used.in.clustering=nv))
         
         #Errorbars
         Sa<-abs(pres$S)
         sSa<-sum(Sa)
         cerr<-colSums(Sa)/sSa
         rerr<-rowSums(Sa)/sSa

         #Genes
         u<-t(t(sd$u[,1:nu])*head(sd$d,nu))
         gi<-u  #/sqrt(rowSums(u^2)) #not needed since we look at the correlation
         #corgi<-cor(t(gi))
         unorm<-u/sqrt(rowSums(u^2))
         corgi<-unorm%*%t(unorm)
         #corgi<-cor(t(gi))
                
         #Patients, 
         v<-t(sd$vt[1:nv,]*head(sd$d,nv))
         pi<-v    #/sqrt(rowSums(v^2)) #we do not normalize patients (columns)
        
         vnorm<-v/sqrt(rowSums(v^2))
         corpi<-vnorm%*%t(vnorm)
         #corpi<-cor(t(pi))

         colnames(pres$L)<-colnames(genedata)
     } else {
        pres<-u<-v<-rerr<-cerr<-NULL
        
        if (pca.type=="svd"){
            sd<-La.svd(genedata)
            d2<-sd$d^2
            if (is.function(nv))
                 nv<-nv(sd$d)
            nv<-min(nv,length(sd$d))
            nu<-nv
            #Genes
            u<-t(t(sd$u[,1:nu])*head(sd$d,nu))
            gi<-u
            #corgi<-cor(t(gi))
            unorm<-u/sqrt(rowSums(u^2))
            corgi<-unorm%*%t(unorm)

            #Patients, 
            v<-t(sd$vt[1:nv,]*head(sd$d,nv))
            pi<-v 
            #corpi<-cor(t(pi))
            vnorm<-v/sqrt(rowSums(v^2))
            corpi<-vnorm%*%t(vnorm)
   
        } else {
           gi<-genedata         
           pi<-t(genedata)
           corgi<-cor(t(genedata))
           corpi<-cor(t(pi))
        }        
        
        colnames(corgi)<-1:NCOL(corgi)
        rownames(corgi)<-1:NCOL(corgi)        
     }
     
     if (absolute.gene.clustering) corgi<-abs(corgi)
     gi.dist<-  if (clustering.distance=="L1") as.dist(sqrt(1-corgi))
                else if (clustering.distance=="angular") as.dist(acos(corgi))
                else                       as.dist(1-corgi)

     #pi.dist<- if (use.L1.for.clustering) dist(pi,method="minkowski",p=1) 
     #          else  dist(pi)^2           
     pi.dist<-  if (clustering.distance=="L1") as.dist(sqrt(1-corpi))
                else if (clustering.distance=="angular") as.dist(acos(corpi))
                else                       as.dist(1-corpi)

     print(c(dim(pi),dim(gi)))
     
     hclustfun<- if(is.null(hclust.row)) function(x)hclust(x,method=genes.method)
                 else function(x)hclust.row
     hclustfun.col<-if(is.null(hclust.col)) function(x)hclust(x,method=sample.method)
                 else function(x)hclust.col
     
     N<-NROW(plot.genedata)
     NC<-NCOL(plot.genedata)
     ffamily<-"sans"
     
     #print(rbind(dim(pres$L),dim(plot.genedata),dim(pres$S)))     
     
     if (do.plot){
         tot.width<-if(!plot.only.raw&&robust.pca) 3*hm.width else hm.width
         if (create.device=="png") png(height=hm.height,width=tot.width,filename=pname,res=130,family=ffamily)
         else if (create.device=="eps") do.call(postscript,list(file=pname,onefile=TRUE,family=ffamily,encoding="ISOLatin1.enc",paper="special",
         horizontal=FALSE,
         pointsize=eps.pointsize,
         height=hm.height/eps.pointtopixel*eps.pointsize,width=tot.width/eps.pointtopixel*eps.pointsize))
         else if (create.device=="pdf") do.call(pdf,list(file=pname,onefile=TRUE,family=ffamily,encoding="ISOLatin1.enc",paper="special",
         pointsize=eps.pointsize,
         height=hm.height/eps.pointtopixel*eps.pointsize,width=tot.width/eps.pointtopixel*eps.pointsize))
                  
         do.call(par,c(list(family=ffamily),par.args))
         
         data.args<-if (!plot.only.raw&&robust.pca) list(pres$L,plot.genedata,pres$S,rerr=rerr,cerr=cerr)
                    else list(plot.genedata,NULL,NULL)
                    
         hmplot<-do.call(myheatmap.3,c(data.args,list(
          col= col.panel,
          distfun=function(x)gi.dist,
          distfun.col=function(x)pi.dist,
          hclustfun=hclustfun,
          hclustfun.col=hclustfun.col,
          dendrogram=c("row","both")[dend.both+1],
          Rowv=if (is.null(row.order)) TRUE else row.order,         
          Colv=if (is.null(col.order)) dend.both else col.order,
          scale="none",#"row",
          #rowsep=41,
          breaks=breaks,
          symbreaks=centralize.for.plot,
          symkey=centralize.for.plot,
          ColSideColors=sample.labels,
          ColSideColorsMarks=sample.labels.marks,
          RowSideColors=gene.colors,
          labRow=geneSymbol,
          key=TRUE,
          density.info=density.info,
          trace="none",#"row",
          keysize=0.8,
          densadj=2,
          #cexCol=1,
          #cexRow=0.8,
          cexRow = 0.27 + 1/log10(N),
          cexCol = 0.7,
          lwid=hm.wprop,#c(1.5,9),
          lhei=hm.hprop[c(1,3)],
          lhei.colside=hm.hprop[2],
          colsep = colsep,
          rowsep = rowsep,
          sepwidth = sepwidth,
          main=main.title,
          subtitle.key=if(robust.pca||svd.pca) paste("Clustered with\n",nv,"singular values")
                       else NULL,
          rect.column=rect.column, rect.row=rect.row,
          plot.extreme.quantile=plot.extreme.quantile,
          row.cluster.order.by.vars=  row.cluster.order.by.vars,
          legend.fun=legend.fun
          )))
          if (close.device) dev.off()
    } else {
          hmplot<-list()
          hcr <- hclustfun(gi.dist)
          ddr <- as.dendrogram(hcr)
          hcc <- hclustfun.col(pi.dist)
          ddc <- as.dendrogram(hcc)

          hmplot$rowDendrogram<-ddr
          hmplot$colDendrogram<-ddc
          hmplot$rowInd<-order.dendrogram(ddr)
          hmplot$colInd<-order.dendrogram(ddc)                  
    }
    if (save.sample.dist){
        hmplot$pi.dist<-pi.dist
    }
    if (save.gene.dist){
        hmplot$gi.dist<-gi.dist
    }    
    invisible(list(svd=pres,u=u,v=v,rerr=rerr,cerr=cerr,hmplot=hmplot))
}

require(gplots)

myheatmap.3<-function (x1,x2,x3,cerr=rep(0,NCOL(x1)),rerr=rep(0,NROW(x1)),
    Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
    distfun = dist, distfun.col = distfun,
    hclustfun = hclust, hclustfun.col = hclustfun,
    dendrogram = c("both",
        "row", "column", "none"), symm = FALSE, scale = c("none",
        "row", "column"), na.rm = TRUE, revC = identical(Colv,
        "Rowv"), add.expr=NULL, breaks=NULL, symbreaks = min(x < 0, na.rm = TRUE) ||
        scale != "none", col = "heat.colors", colsep=NULL, rowsep=NULL,
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote=NULL, notecex = 1,
    notecol = "cyan", na.color = par("bg"), trace = c("column",
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
    vline = median(breaks), linecol = tracecol, margins = c(5,
        5), ColSideColors=NULL, ColSideColorsMarks=NULL, RowSideColors=NULL,
    cexRow = 0.2 + 1/log10(nr),
    cexCol = 0.2 + 1/log10(nc), 
    labRow = NULL, labCol = NULL,
    srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA,
        0), offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 1.5,
    density.info = c("histogram", "density", "none"), denscol = tracecol,
    symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
    subtitle.key=NULL,
     main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL,
    lwid = NULL,lhei.colside=0.25,
    rect.column=4, rect.row=4,
    plot.extreme.quantile=0.005,
    row.cluster.order.by.vars=TRUE,
    legend.fun=NULL,
        ...)
{
    #Is it a triple plot?
    allxs<-list(x1,x2,x3)
    nnull<-which(!sapply(allxs,is.null))
    triple.plot<-length(nnull)==3
    
    x<-allxs[[nnull[1]]]    
    
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (is.null(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Colv)&&!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv)&&!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    apply.color.labels<-function(ddr){
        ## toy example to set colored leaf labels :
        #i <- 0
        #mycols <- grDevices::rainbow(length(labRow))
        colLab <- function(n) {
               if(is.leaf(n)) {
                 a <- attributes(n)
                 #i <<- i+1
                 gnum<-as.numeric(attr(n,"label"))
		 if (!is.na(gnum)){
                 mod.lst<-list(lab.font=1,lab.cex=cexRow,pch='.')
                 if (!is.null(RowSideColors)&&!is.na(RowSideColors[gnum])){
                     mod.lst$pch<-15
                     mod.lst<-c(mod.lst,list(lab.col = RowSideColors[gnum],                                     
                                     col = RowSideColors[gnum],
                                      bg = RowSideColors[gnum]))
                 }
                 attr(n, "nodePar") <- c(a$nodePar, mod.lst)
                 attr(n,"label")<-labRow[gnum]
		 }
               }
               n
        }
        dendrapply(ddr, colLab)                
    }

    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        ddr <- apply.color.labels(ddr)
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- apply.color.labels(ddr)
        ddr <- reorder(ddr, Rowv,agglo.FUN=mean)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        if (row.cluster.order.by.vars)
            Rowv <- -rowVars(x, na.rm = na.rm)
        else
            Rowv <- rowMeans(x, na.rm = na.rm)

        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- apply.color.labels(ddr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }


    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        ddc <- apply.color.labels(ddc)
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun.col(distfun.col(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv,agglo.FUN=mean)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- -colMeans(x, na.rm = na.rm)
        hcc <- hclustfun.col(distfun.col(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    compute.breaks<-function(breaks,x){
        #print("breaks computing")
        if (length(breaks) == 1) {
            if (!symbreaks)
                #breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                breaks <- c(min(x, na.rm = na.rm),seq(quantile(x, plot.extreme.quantile, na.rm = na.rm),
                                                      quantile(x,1-plot.extreme.quantile, na.rm = na.rm),
                           length = breaks-2),max(x,na.rm=na.rm))
            else {
                extreme<-quantile(abs(x), 1-2*plot.extreme.quantile, na.rm = na.rm)
                #extreme <- max(abs(x), na.rm = TRUE)
                breaks <- c(min(x, na.rm = na.rm),seq(-extreme, extreme, length = breaks-2),
                            max(x,na.rm=na.rm))
            }
        }
        #print("breaks computed")
        breaks
    }
    breaks<-compute.breaks(breaks,x)

    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(c(4,2), c(3,1))
        if (!is.null(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) !=
                nc)
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(5, 1), lmat[2, ] + 1)
            lhei <- c(lhei[1], lhei.colside, lhei[2])
        }
        #if (!is.null(RowSideColors)) {
        #    if (!is.character(RowSideColors) || length(RowSideColors) !=
        #        nr)
        #        stop("'RowSideColors' must be a character vector of length nrow(x)")
        #    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
        #        1), 1), lmat[, 2] + 1)
        #    lwid <- c(lwid[1], 0.3, lwid[2])
        #}
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    if (triple.plot){
        #Make space for additional plot
        mlm<-max(lmat)
        lmat<-cbind(lmat,c(mlm+4,mlm+2,mlm+3),c(0,0,mlm+1),c(mlm+7,mlm+5,mlm+6))
        iwid<-tail(lwid,1)
        lwid<-c(lwid,iwid,0.2,iwid)
        #print(list(lmat=lmat, widths = lwid, heights = lhei, respect = FALSE))
    }
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)


    #if (!is.null(RowSideColors)) {
    #    par(mar = c(margins[1], 0, 0, 0.5))
    #    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    #}
   plot.colsidecolors<-expression({
    if (!is.null(ColSideColors)) {        
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        if (!is.null(ColSideColorsMarks)){
                marks<-as.numeric(as.logical(ColSideColorsMarks[colInd]))+1
                image(cbind(1:nc), col = c("transparent","#00000055")[marks], axes = FALSE,add=TRUE)
        }
    }
   })
  eval(plot.colsidecolors) 

  before.par<-par()
  before.env<-new.env()
  for (v in ls(parent.frame())) {assign("v",get(v,envir=parent.frame()),envir=before.env)}
  plot.image.for.x<-expression({
    #print(" plot.image.for.x")
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr

    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!gtools::invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }    
    if (is.null(srtCol)){
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
            padj = 0.7)#adjCol[2])
        axis(3, 1:nc, labels = labCol, las = 2, line = -1 +
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
            padj = 0.7)#adjCol[2])
    }
    else {
        if (is.numeric(srtCol)) {
            if (is.null(adjCol) || is.null(adjCol))
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
                strheight("M"), labels = labCol, adj = adjCol,
                cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2,
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
                srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!is.null(add.expr))
        eval(substitute(add.expr))
    if (!is.null(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5-sepwidth[1]/2, ybottom = 0,
            xright = csep + 0.5 + sepwidth[1]/2, ytop = ncol(x) +
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!is.null(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
            1 - rsep) - 0.5 + sepwidth[2]/2, xright = nrow(x) + 1, ytop = (ncol(x) +
            1 - rsep) - 0.5 - sepwidth[2]/2, lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!is.null(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
   
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
        if (is.numeric(rect.column))
            rect.hclust(as.hclust(ddc),k=rect.column,border="pink")
    }
    else plot.new()        
    if (!is.null(legend.fun)) {
        legend.fun()
    }
   }) #######
    eval(plot.image.for.x)

    if (!is.null(main)) {        
        mmain<-if (triple.plot) paste(main,"(smoothed)")
               else main
        title(mmain, cex.main = 1.5 * op[["cex.main"]])
    }

    par(mar = c(margins[1], 0, 0, 4))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i")#, leaflab = "none")
        if (is.numeric(rect.row))
            #for(k in 2:rect.row) 
            rect.hclust(as.hclust(ddr),k=rect.row,border="pink",horiz=TRUE)
    }
    else plot.new()


    if (triple.plot){
        #Plot key based on original data
        x<-x2[rowInd, colInd]
        breaks<-compute.breaks(length(breaks),x)
    }

    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key and Density",sub=subtitle.key)
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key and Histogram",sub=subtitle.key)
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    
    retval$breaks<-breaks
    
    if (triple.plot){
        #suppressWarnings(do.call(par,before.par))
        par(cex=before.par$cex)

        ncol.err=20
        #Row Errors
        max.err=1/length(rerr)*5
        col.err=colorpanel(ncol.err,"white","cyan","red")
        err.breaks<-c(0,exp(seq(from=log(max.err/ncol.err),to=log(max.err),len=ncol.err-1)),1)

        par(mar = c(margins[1], 0, 0, 0.5))
        image(z=rbind(rerr[rowInd]), col = col.err, axes = FALSE,
             breaks=err.breaks)

        #Column Errors
        max.err=1/length(cerr)*5
        col.err=colorpanel(ncol.err,"white","cyan","red")
        err.breaks<-c(0,exp(seq(from=log(max.err/ncol.err),to=log(max.err),len=ncol.err-1)),1)

        #par(mar = c(0.5, 0, 0, margins[2]))
        #image(z=cbind(cerr[colInd]), col = col.err, axes = FALSE,
        #     breaks=err.breaks)
        eval(plot.colsidecolors)

        x<-x2[rowInd, colInd]
        breaks<-compute.breaks(length(breaks),x)
        retval$breaks<-breaks

        #assign("ColSideColors",cerr,envir=before.env)
        assign("x",x,envir=before.env)
        assign("breaks",breaks,envir=before.env)
        assign("left.gene.labels",FALSE,envir=before.env)
        eval(plot.image.for.x,envir=before.env)

        if (!is.null(main))
             title(paste(main,"(raw data)"), cex.main = 1.5 * op[["cex.main"]])

        #Column Errors
        par(mar = c(0.5, 0, 0, margins[2]))
        image(z=cbind(cerr[colInd]), col = col.err, axes = FALSE,
             breaks=err.breaks)

        x<-x3[rowInd, colInd]
        #breaks<-compute.breaks(length(breaks),x)
        #assign("ColSideColors",cerr,envir=before.env)
        assign("x",x,envir=before.env)
        assign("breaks",breaks,envir=before.env)
        assign("left.gene.labels",FALSE,envir=before.env)
        eval(plot.image.for.x,envir=before.env)

        if (!is.null(main))
            title(paste(main,"(subtracted outliers)"), cex.main = 1.5 * op[["cex.main"]])
    }

    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}





rect.hclust<- function (tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2, 
    cluster = NULL, horiz=FALSE) 
{
    if (length(h) > 1L | length(k) > 1L) 
        stop("'k' and 'h' must be a scalar")
    if (!is.null(h)) {
        if (!is.null(k)) 
            stop("specify exactly one of 'k' and 'h'")
        k <- min(which(rev(tree$height) < h))
        k <- max(k, 2)
    }
    else if (is.null(k)) 
        stop("specify exactly one of 'k' and 'h'")
    if (k < 2 | k > length(tree$height)) 
        stop(gettextf("k must be between 2 and %d", length(tree$height)), 
            domain = NA)
    if (is.null(cluster)) 
        cluster <- cutree(tree, k = k)
    clustab <- table(cluster)[unique(cluster[tree$order])]
    m <- c(0, cumsum(clustab))
    if (!is.null(x)) {
        if (!is.null(which)) 
            stop("specify exactly one of 'which' and 'x'")
        which <- x
        for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
    }
    else if (is.null(which)) 
        which <- 1L:k
    if (any(which > k)) 
        stop(gettextf("all elements of 'which' must be between 1 and %d", 
            k), domain = NA)
    border <- rep(border, length.out = length(which))
    retval <- list()
    for (n in seq_along(which)) {
        if (horiz) {
         ybottom<-m[which[n]] + 0.66 
         ytop<-m[which[n] + 1] + 0.33 

         xright<-0.5 #par("usr")[3L]
         xleft<-(par("usr")[1L]-xright)/2 #mean(rev(tree$height)[(k - 1):k])
         #print(c(xright=xright,xleft=xleft))
        } else {
         xleft<-m[which[n]] + 0.66 
         xright<-m[which[n] + 1] + 0.33 
 
         ybottom<-par("usr")[3L] 
         ytop<-mean(rev(tree$height)[(k - 1):k])
        }                          
        rect(xleft, ybottom, xright, ytop, border = border[n])
        retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
    }
    invisible(retval)
}










### Error estimation of robustPCA

fit.beta.toerr<-function(rerr,do.plot=FALSE,
                         prior=1/1000,
                         signif=0.01,
                         q.check=0.98){
  n<-length(rerr)    
  #Variance for beta(a,b), scales varbeta(ka,kb)~varbeta(a,b)/k
  varbeta<-function(a,b){a*b/((a+b)^2 * (a+b+1))}
    
  #Set a/b so that mean is a/(a+b)=1/n, mean doesn't change with ka,kb scale
  a<-1
  b<-n-1
  vab<-varbeta(a,b)  
  
  q8<-quantile(rerr,q.check)
  tailrerr<-sort(rerr[rerr>q8])
  num.fits=length(tailrerr)
  all.ft<-lapply(tailrerr,function(tr){
   rerr<-rerr[rerr<=tr]
   k<-vab/var(rerr)
   #print(c(a,b,k))
   #print(c(a=a*k,b=b*k))    
   ft<-fitdistr(rerr,"beta",list(shape1=a*k,shape2=b*k),method="L-BFGS")
   ft$n<-length(rerr)
   ft$k<-k
   ft
  })
  llik<-sapply(all.ft,#BIC
               function(ft){
      ft$loglik+dbinom(n-ft$n,n,prior,log=TRUE)
  })
  best<-which.max(llik)
  ft<-all.ft[[best]]
    
  pvals<-pbeta(tailrerr,ft$estimate[1],ft$estimate[2],lower.tail=FALSE)            
  pvals<-p.adjust(pvals,method="BH",n=n)
  outliers<-which(pvals<=signif)
  outlier.regime<-tailrerr[outliers[1]]
  
  #Cut with Bonferroni correction
  #outlier.regime=qbeta(1-signif/length(rerr),ft$estimate[1],ft$estimate[2])
  #outliers<-which(tailrerr<=outlier.regime)
    
  if (do.plot) {
      par(mfcol=c(1,2))
      plot(llik,type="b",ylab="LogLikelihood",main="Beta fits")
      points(best,llik[best],col="red",pch="Z")
      hist(rerr,breaks=40,freq=FALSE)
      curve(dbeta(x,ft$estimate[1],ft$estimate[2]),add=TRUE,col="darkgreen",lw=2)
      abline(v=outlier.regime,lty=2,col="red")
  }    
  list(outlier.regime=outlier.regime,num.outliers=length(outliers),
       beta.fit=ft,num.fits=num.fits)
       #p.values=pvals)
}
   
   
   
###                           
###  Adding a blank leaf to dendrogram                           
###                           
dendrapply.postfix<-function (X, FUN, ...) 
{
    FUN <- match.fun(FUN)
    if (!inherits(X, "dendrogram")) 
        stop("'X' is not a dendrogram")
    Napply <- function(d) {                
        p1<-0
        if (!is.leaf(d)) {
            sub.res <- lapply(d, function(dd){
                nr<-Napply(dd)
                p1<<-p1+nr$p1   
                nr$r
            })                        
        }
        rr <- FUN(d, ...) 
        r<-rr$n
        here<-rr$here       
        if (!is.leaf(d)) {
            if (!is.list(r)) 
                r <- as.list(r)
            if (length(r) < (n <- length(d))) 
                r[seq_len(n)] <- vector("list", n)
            r[]<-sub.res
            nn<-sum(sapply(sub.res,function(s)attr(s,"members")))
            attr(r,"members")<- nn -here+p1
        }
        list(r=r,p1=here)
    }
    rr<-Napply(X)
    rr$r                       
}


add.blank.leaf<-function(ddr,k,before=FALSE,h=0)  {
       i <- 0
       addl <- function(n) {
           if(is.leaf(n)) {
             i <<- i+1             
           }
           here<-FALSE
           if (i==k) {
               here<-TRUE               
               m<-n
               attr(m,"label")<-""
               attr(m,"edgePar")<-list(lwd=0)
               n<- if (before)
                   merge(m,n,height=h)
                   else
                   merge(n,m,height=h)   
               #attr(n,"members")<-attr(n,"members")
               i<<-i+1
           }
           list(n=n,here=as.integer(here))
       }
       
     dendrapply.postfix(ddr, addl)
}

#ddra<-add.blank.leaf(ddr,8,TRUE)
#plot(ddra)





