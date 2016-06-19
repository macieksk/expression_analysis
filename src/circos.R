
library(OmicCircos)

convert.topdf<-function(basefname,density=150,rm.png=FALSE){
   fname.png<-paste(basefname,".png",sep="")    
   system(paste("convert ",fname.png," -fuzz 0 -transparent white -background none -density",density,
             paste(basefname,".pdf",sep="")))
   if (rm.png) system(paste("rm",fname.png))
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


compare.clust.circos<-function(cl1,cl2, 
    labels=NULL,
    clust.names=paste("C",1:2,sep=""),
    do.revert.2nd=TRUE,
    main ="Clustering Comparison",
    chr.lab.diff=-10,
    font.width=5,
    lab.shift=median(nchar(labels)),
    lwd = 45
    ){ 
  stopifnot(length(cl1)==length(cl2))
  cl1<-factor(cl1)
  cl2<-factor(cl2)
  
  N<-length(cl1)
  k<-2 
    
  cldf<-data.frame(cl1=cl1,cl2=cl2)  
  cldf<-cbind(cldf,id=1:length(cl1))   
    
  #Order for the purpose of good element sort  
  #smat<-clust.sim.matrix(cldf)
  #or<-hclust(dist(smat),method="ward.D2")$order #,method="ward")$order  
  #c1df<-cldf[or,]
  
  intsep<-"_.__"    
  inter12<-interaction(cldf[,1],cldf[,2],sep=intsep)
  inter21<-interaction(cldf[,2],cldf[,1],sep=intsep)
  cldf<-cbind(cldf,inter12=as.character(inter12),inter21=as.character(inter21),stringsAsFactors=FALSE)    
  #print(c(inter12=inter12))
  #This sort should be stable to preserve the above order within clusters  
  or1<-order(inter21,1:N)  #21 in reverse as interaction sorts from the last
  or2<-order(inter12,1:N)
  #or1<-order(as.numeric(cldf[,1]),1:N)
  #or1.2<-order(as.numeric(cldf[or1,2]),1:N)
  #or2<-order(as.numeric(cldf[,2]),1:N)
  #or2.1<-order(as.numeric(cldf[or2,1]),1:N)   

  #Interaction-sorted data for both interactions
  cldfs<-list(cldf[or1,],cldf[or2,])
  #cldfs<-list(cldf[or2,][or2.1,],cldf[or1,][or1.2,])
  
  #Counts within clusters for proper seg.f for circos  
  cumsum.cnts<-function(cl1){  
   clcnt1<-tapply(cl1,cl1,FUN=length)
   #clcnt1<-clcnt1[order(names(clcnt1))]
   clcnt1<-cumsum(clcnt1)
   clcnt1<-c(a=0,clcnt1)
   names(clcnt1)<-names(clcnt1)[c(2:length(clcnt1),1)]
   clcnt1
  }
  
  #Seg.f data for circos segAnglePo  
  give.segf<-function(cl,clname,acc.start=NA){
    clcnt<-cumsum.cnts(cl)   
     
    seg.f<-data.frame(seg.name=paste(clname,as.character(cl),sep="."),
                seg.Start=0:(N-1), seg.End=1:N, the.v=acc.start, NO=NA,
                stringsAsFactors=FALSE)
    
   seg.f$seg.Start<-seg.f$seg.Start-clcnt[as.character(cl)]
   seg.f$seg.End<-seg.f$seg.End-clcnt[as.character(cl)]
   seg.f
  }
  
 
    
  seg.fs<-lapply(1:k,function(i)give.segf(cldfs[[i]][,i],clust.names[i]))
  
  dbs<-lapply(seg.fs,function(seg.f1)as.data.frame(segAnglePo(seg.f1,seg=unique(seg.f1$seg.name))
                                                   ,stringsAsFactors=FALSE))
  
  dbs.org<-dbs    
  
  #Rotate and narrow angles so both clusterings fit in one cicrcle    
  narrow.angles<-function(ang,k,offset=0)(ang%%360)/k+offset
      
  dbs<-lapply(1:k,function(i)      
       within(dbs[[i]],{
             ofs.in<-as.numeric(angle.start[1]) 
             #print(cbind(as=angle.start,ae=angle.end))
             angle.start<-(as.numeric(angle.start)-ofs.in)
             angle.end<-(as.numeric(angle.end)-ofs.in)
             rm(ofs.in)
             if (do.revert.2nd){
               if (i%%2==0){                
                 tmp<-angle.start
                 angle.start<- 359-angle.end
                 angle.end<- 359-tmp
                 rm(tmp)
               }
             }
             p<-9/4
             #print(cbind(as=angle.start,ae=angle.end))
             angle.start<-narrow.angles(angle.start,p,(i-1)*360/2+360*((1-2/p)/4))
             angle.end<-narrow.angles(angle.end,p,(i-1)*360/2+360*((1-2/p)/4))
             
             angle.start<-angle.start+90
             angle.end<-angle.end+90
       }))
  db<-do.call(rbind,dbs)
  
  if (do.revert.2nd){
    #We need to revert within cluster order of 2nd clustering,
    #so we have mirror symmetry
    cldfs[[2]]<-with(cldfs[[2]],
      do.call(rbind,lapply(unique(cl2),function(c2c){
         ids<-which(cl2==c2c)
         or2[ids]<<-or2[rev(ids)] #Also revert order vec
         cldfs[[2]][rev(ids),]
      }))
    )
  }
      
              
  link.pg<-lapply(unique(cldfs[[1]]$inter12),function(incl){
    #For each interaction cluster
    spl<-strsplit(incl,intsep)[[1]]   
    segnames<-sapply(1:k,function(i)paste(clust.names[i],spl[i],sep="."))
    
    ids<-which(cldfs[[1]]$inter12==incl)
    n<-length(ids)
    r1<-range(ids)
    incl2<-paste(spl[2],spl[1],sep=intsep) 
        
    r2<-range(which(cldfs[[2]]$inter21==incl2))
    
    clcnt1<-cumsum.cnts(cldfs[[1]]$cl1)     
    clcnt2<-cumsum.cnts(cldfs[[2]]$cl2)     
    #print(spl)    
    #print(c(i1=incl,i2=incl2))    
    #print(c(r1=r1,r2=r2))    
        
    rw<-data.frame(seg1=segnames[1],start1=r1[1]-1-clcnt1[spl[1]],end1=r1[2]-clcnt1[spl[1]],
                   seg2=segnames[2],start2=r2[1]-1-clcnt2[spl[2]],end2=r2[2]-clcnt2[spl[2]],
                   color=NA,n=n)
    #print(rw)
    rw
  })
  
  Ncl1<-length(unique(cl1))        
  Ncl2<-length(unique(cl2))    

  #print(c(nlink=NROW(link.pg)))
  link.pg<-do.call(rbind,link.pg)
  link.pg$color<-rainbow(NROW(link.pg),alpha=0.5)
  link.pg<-link.pg[order(link.pg$n,decreasing=TRUE),]
  
  to.remove<-apply(link.pg,1,function(rw)any((rw==Inf)||(rw==-Inf)))
  link.pg<-link.pg[!to.remove,]

      
  colors = c(rainbow(Ncl1,s=0.8,v=0.8),rainbow(Ncl2,s=0.6,v=0.6))     
      
  par(mar=c(2, 2, 2, 2));
      
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main=main);
  
  abline(v=400,lty=3,col="darkred",lwd=2)    
      
      
  circos(R=337, cir=db, W=10, mapping=link.pg, type="link.pg", lwd=0,
         link.curvature=1)            

  circos(R=350, cir=db, type="chr",  col=colors, print.chr.lab=TRUE, 
        chr.lab.diff=chr.lab.diff,font.width=font.width,
        W=1, scale=F,B=F,lwd=lwd,cex=1.5,rot.ang=0)
         
  
  if (!is.null(labels)){
    labels.df<-do.call(rbind,lapply(db$seg.name,function(sn){
     seg<-db[db$seg.name==sn,]
     #print(seg)
     #print(c(seg$seg.start,seg$seg.end,by=4))
     starts<-seq(as.numeric(seg$seg.start),as.numeric(seg$seg.end)-1,by=1)+0.5
     data.frame(chr=sn,po=starts,Gene=paste("G",starts,sep=""),stringsAsFactors=FALSE)
    }))
    if (labels[1]=="only.tics"){
        labels.df$Gene<-""
    } else{
        labels.df$Gene<-c(labels[or1],labels[or2])
    }
    circos(R=360, cir=db, W=20, mapping=labels.df, type="label", side="out", col="black",
           cex=0.5,lab.shift=lab.shift)
  }
  
  #image(z=matrix(1:100,ncol=10),col=link.pg$color)   
      
  invisible(ret<-list(seg=seg.fs,db.org=dbs.org,db=db,link=link.pg,col=link.pg$color))
}  
              




circos<-function (mapping = mapping, xc = 400, yc = 400, R = 400, W = W, 
    cir = "", type = "n", col.v = 3, B = F, print.chr.lab = T, 
    chr.lab.diff=R/20,
    lab.shift = 12,
    col.bar = F, col.bar.po = "topleft", cluster = F, order = NULL, 
    scale = F, lwd.scale=0.1,
    cutoff = "n", zoom = "", lwd = 1, col = rainbow(10,alpha = 0.5)[7], side = "",
    cex=1,rot.ang=NULL,link.curvature=2,
    font.width=3) 
{
    r <- R
    cir.m <- 0
    if (is.matrix(cir) == T | is.data.frame(cir) == T | class(cir)[1] == 
        "GRanges") {
        if (class(cir)[1] == "GRanges") {
            chr.po <- cbind(as.character(seqnames(cir)), start(cir), 
                names(cir), mcols(cir))
        }
        else {
            chr.po <- cir
        }
        chr.po[, 1] <- gsub("chr", "", chr.po[, 1])
        chr.num <- nrow(chr.po)
        cir.m <- 1
    }
    else if (nchar(cir) >= 1) {
        pofile.s <- paste("UCSC.", cir, ".RData", sep = "")
        pofile <- system.file("data", pofile.s, package = "OmicCircos")
        if (file.exists(pofile)) {
        }
        else {
            stop("cir name?")
        }
        chr.po <- local(get(load(pofile)))
        chr.po[, 1] <- gsub("chr", "", chr.po[, 1])
        chr.num <- nrow(chr.po)
    }
    else {
        stop("cir name ??")
    }
    if (length(zoom) > 1) {
        chr1 <- which(chr.po[, 1] == zoom[1])
        chr2 <- which(chr.po[, 1] == zoom[2])
        chr.num <- length(chr1:chr2)
        chr.po <- zoom.in(cir.in = chr.po, zoom = zoom)
    }
    if (type == "chr") {
        chr.lw = W
        if (cir.m == 0) {
            chrfile.s <- paste("UCSC.", cir, ".chr.RData", sep = "")
            chrfile <- system.file("data", chrfile.s, package = "OmicCircos")
            if (file.exists(chrfile)) {
            }
            else {
                stop("chrfile name?")
            }
            dat.c <- local(get(load(chrfile)))
            dat.c[, 1] <- gsub("chr", "", dat.c[, 1])
            for (chr.i in c(1:chr.num)) {
                chr.s <- chr.po[chr.i, 1]
                v1 <- as.numeric(chr.po[chr.i, 2])
                v2 <- as.numeric(chr.po[chr.i, 3])
                v3 <- as.numeric(chr.po[chr.i, 6])
                v4 <- as.numeric(chr.po[chr.i, 7])
                dat.v <- subset(dat.c, dat.c[, 1] == chr.s)
                dat.v <- dat.v[order(as.numeric(dat.v[, 2])), 
                  ]
                for (i in 1:nrow(dat.v)) {
                  if (dat.v[i, 5] == "gneg") {
                    col <- colors()[351]
                  }
                  else if (dat.v[i, 5] == "acen" | dat.v[i, 5] == 
                    "gvar" | dat.v[i, 5] == "stalk") {
                    col <- colors()[26]
                  }
                  else {
                    col.v <- gsub("gpos", "", dat.v[i, 5])
                    if (col.v >= 50) {
                      col <- colors()[300]
                    }
                    else {
                      col <- colors()[351]
                    }
                  }
                  w1 <- OmicCircos:::scale.v(dat.v[i, 2], v1, v2, v3, v4)
                  w2 <- OmicCircos:::scale.v(dat.v[i, 3], v1, v2, v3, v4)
                  OmicCircos:::draw.arc.s(xc, yc, r, w1, w2, col = col, lwd = chr.lw)
                }
                if (print.chr.lab) {
                  v1 <- as.numeric(chr.po[chr.i, 2])
                  v2 <- as.numeric(chr.po[chr.i, 3])
                  w.m <- (v1 + v2)/2
                  r.j <- chr.lab.diff
                  chr.t <- gsub("chr", "", chr.s)
                  shft = font.width*cex*nchar(chr.t) 
                  draw.text.rt(xc, yc, r + r.j-shft, w.m, chr.t, cex = cex,rot.ang=rot.ang)
                }
                if (scale) {
                  total.num <- as.numeric(chr.po[nrow(chr.po), 5])
                  #print(c(V1 = v1, V2 = v2, V3 = v3, V4 = v4))
                  OmicCircos:::do.scale.cir(xc = xc, yc = yc, the.r = r + 10, 
                  total.num = total.num, col = "blue", 
                    lwd = lwd.scale, V1 = v1, V2 = v2, V3 = v3, V4 = v4)
                }
            }
        }
        else {
            for (chr.i in c(1:chr.num)) {
                w1 <- as.numeric(chr.po[chr.i, 2])
                w2 <- as.numeric(chr.po[chr.i, 3])
                OmicCircos:::draw.arc.s(xc, yc, r, w1, w2, col = col[chr.i], 
                  lwd = lwd)
                if (print.chr.lab) {
                  w.m <- (w1 + w2)/2
                  r.j <- chr.lab.diff
                  chr.t <- gsub("chr", "", chr.po[chr.i, 1])
                  shft <- font.width*cex*nchar(chr.t)
                  draw.text.rt(xc, yc, r + r.j-shft, w.m, chr.t, cex = cex,rot.ang=rot.ang)
                }
                if (scale) {
                  v1 <- as.numeric(chr.po[chr.i, 2])
                  v2 <- as.numeric(chr.po[chr.i, 3])
                  v3 <- as.numeric(chr.po[chr.i, 6])
                  v4 <- as.numeric(chr.po[chr.i, 7])
                  total.num <- as.numeric(chr.po[nrow(chr.po), 5])
                  #print(c(V1 = v1, V2 = v2, V3 = v3, V4 = v4))
                  OmicCircos:::do.scale.cir(xc = xc, yc = yc, the.r = r + 10,
                    total.num = total.num, col = "blue", 
                    lwd = lwd.scale, V1 = v1, V2 = v2, V3 = v3, V4 = v4)
                }
            }
        }
    }
    if (type == "l" | type == "s" | type == "h" | type == "b" | 
        type == "ls" | type == "lh") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        dat.min <- min(as.numeric(dat.in[, col.v]), na.rm = T)
        dat.max <- max(as.numeric(dat.in[, col.v]), na.rm = T)
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
            if (dim(dat)[1] == 0) {
                next
            }
            my.R1 <- R + W/5
            my.R2 <- R + W - W/5
            my.v <- as.numeric(dat[1, col.v])
            v.old <- OmicCircos:::scale.v(my.v, my.R1, my.R2, dat.min, dat.max)
            po <- as.numeric(dat[1, 2])
            w.from <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
            for (i in 2:nrow(dat)) {
                my.v <- as.numeric(dat[i, col.v])
                if (is.na(my.v)) {
                  next
                }
                my.R1 <- R + W/5
                my.R2 <- R + W - W/5
                v <- OmicCircos:::scale.v(my.v, my.R1, my.R2, dat.min, dat.max)
                po <- as.numeric(dat[i, 2])
                w.to <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                if (type == "ls") {
                  if (v.old > 0) {
                    OmicCircos:::draw.line(xc, yc, w.from, v, v.old, col = col, 
                      lwd = lwd, lend = 2)
                  }
                  OmicCircos:::draw.arc.s(xc, yc, v, w.from, w.to, col = col, 
                    lwd = lwd)
                  v.old <- v
                }
                if (type == "l") {
                  if (v.old > 0) {
                    OmicCircos:::draw.line3(xc, yc, w.from, w.to, v.old, v, 
                      col = col, lwd = lwd)
                  }
                  v.old <- v
                }
                if (type == "lh") {
                  if (v.old > 0) {
                    OmicCircos:::draw.line3(xc, yc, w.from, w.to, v, v, col = col, 
                      lwd = lwd)
                  }
                  v.old <- v
                }
                if (type == "s") {
                  OmicCircos:::draw.point.w(xc, yc, v, w.to, col = col, cex = 0.2)
                }
                if (type == "b") {
                  OmicCircos:::draw.line(xc, yc, w.to, r, v, col = col, lwd = lwd)
                }
                if (type == "h") {
                  tmp.v <- v - r - 1
                  OmicCircos:::draw.arc.pg(xc, yc, w.from, w.to, r, r + tmp.v, 
                    col = col, border = col)
                }
                w.from <- w.to
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "label") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        r <- R
        w <- W
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            if (dim(dat)[1] == 0) {
                next
            }
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            dat <- dat[order(dat[, 2]), ]
            label.num <- nrow(dat)
            v.angle <- v2 - v1
            #maxlab <- as.integer((v2 - v1) * 3)
            #if (label.num > maxlab) {
            #    label.num == maxlab
            #}
            v.wide <- v.angle/label.num
            gap.len <- 2
            #if (v.wide > gap.len) {
            #    v.wide <- gap.len
            #    w.po <- v1 + v.angle/2 - (label.num/2) * gap.len
            #}
            #else {
                w.po <- v1
            #}
            for (i in 1:label.num) {
                label.n <- dat[i, 3]
                po <- as.numeric(dat[i, 2])
                w.to <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                if (side == "in") {
                  OmicCircos:::draw.line(xc, yc, w.to, r, r - w/3, col = col, 
                    lwd = lwd)
                  OmicCircos:::draw.line(xc, yc, w.po, r - w + w/3, r - w, 
                    col = col, lwd = lwd)
                  OmicCircos:::draw.line3(xc, yc, w.to, w.po, r - w/3, r - 
                    w + w/3, col = col, lwd = lwd)
                  draw.text.rt(xc, yc, w.po, r = r - w - 40, 
                    n = label.n, col = col, cex = cex, side = "in",shift=lab.shift)
                }
                else if (side == "out") {
                  OmicCircos:::draw.line(xc, yc, w.to, r, r + w/3, col = col, 
                    lwd = lwd)
                  OmicCircos:::draw.line(xc, yc, w.po, r + w - w/3, r + w, 
                    col = col, lwd = lwd)
                  OmicCircos:::draw.line3(xc, yc, w.to, w.po, r + w/3, r + 
                    w - w/3, col = col, lwd = lwd)
                  draw.text.rt(xc, yc, w.po, r = r + w + 5, 
                    n = label.n, col = col, cex = cex, side = "out",shift=lab.shift)
                }
                w.po <- w.po + v.wide
            }
        }
    }
    if (type == "link") {
        chr.po[, 4] <- gsub("chr", "", chr.po[, 4])
        dat.in <- mapping
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[, 4] <- gsub("chr", "", dat.in[, 4])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        dat.in[dat.in[, 4] == 23, 1] <- "X"
        dat.in[dat.in[, 4] == 24, 1] <- "Y"
        dat <- dat.in
        r <- R
        for (i in 1:nrow(dat)) {
            chr1.s <- dat[i, 1]
            chr2.s <- dat[i, 4]
            po1 <- dat[i, 2]
            po2 <- dat[i, 5]
            chr1 <- which(chr.po[, 1] == chr1.s)
            chr2 <- which(chr.po[, 1] == chr2.s)
            v1 <- as.numeric(chr.po[chr1, 2])
            v2 <- as.numeric(chr.po[chr1, 3])
            v3 <- as.numeric(chr.po[chr1, 6])
            v4 <- as.numeric(chr.po[chr1, 7])
            w1 <- OmicCircos:::scale.v(as.numeric(po1), v1, v2, v3, v4)
            v1 <- as.numeric(chr.po[chr2, 2])
            v2 <- as.numeric(chr.po[chr2, 3])
            v3 <- as.numeric(chr.po[chr2, 6])
            v4 <- as.numeric(chr.po[chr2, 7])
            w2 <- OmicCircos:::scale.v(as.numeric(po2), v1, v2, v3, v4)
            if (chr1 == chr2) {
                OmicCircos:::draw.link(xc, yc, r, w1, w2, col = "#FF000088", 
                  lwd = lwd)
            }
            else {
                OmicCircos:::draw.link(xc, yc, r, w1, w2, col = "#00FFFF88", 
                  lwd = lwd)
            }
        }
    }
    if (type == "link.pg") {
        thefile <- system.file("data", "UCSC.chr.colors.RData", 
            package = "OmicCircos")
        if (file.exists(thefile)) {
        }
        else {
            stop("UCSC.chr.colors ?")
        }
        dat.col <- local(get(load(thefile)))
        dat.col[, 1] <- gsub("chr", "", dat.col[, 1])
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[, 4] <- gsub("chr", "", dat.in[, 4])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (i in 1:nrow(dat.in)) {
            the.d <- as.numeric(dat.in[i, 3]) - as.numeric(dat.in[i, 
                2])
            chr1.s <- dat.in[i, 1]
            chr2.s <- dat.in[i, 4]
            cir.i <- sum(grep("mm", cir))
            if (cir.i == 1) {
                chr2.s <- paste("mm.", dat.in[i, 4], sep = "")
            }
            cir.i <- sum(grep("mac", cir))
            if (cir.i == 1) {
                chr2.s <- paste("mac.", dat.in[i, 4], sep = "")
            }
            if ("color"%in%colnames(dat.in)){
                the.color<-dat.in[i,"color"]
            } else {
                the.color <- as.character(dat.col[dat.col[, 1] == 
                    chr1.s, 2])
                the.color <- paste(the.color, "50", sep = "")
            }
            po1.1 <- as.numeric(dat.in[i, 2])
            po1.2 <- as.numeric(dat.in[i, 3])
            po2.1 <- as.numeric(dat.in[i, 5])
            po2.2 <- as.numeric(dat.in[i, 6])
            chr1 <- which(chr.po[, 1] == chr1.s)
            chr2 <- which(chr.po[, 1] == chr2.s)
            v1 <- as.numeric(chr.po[chr1, 2])
            v2 <- as.numeric(chr.po[chr1, 3])
            v3 <- as.numeric(chr.po[chr1, 6])
            v4 <- as.numeric(chr.po[chr1, 7])
            w1.1 <- OmicCircos:::scale.v(po1.1, v1, v2, v3, v4)
            w1.2 <- OmicCircos:::scale.v(po1.2, v1, v2, v3, v4)
            v1 <- as.numeric(chr.po[chr2, 2])
            v2 <- as.numeric(chr.po[chr2, 3])
            v3 <- as.numeric(chr.po[chr2, 6])
            v4 <- as.numeric(chr.po[chr2, 7])
            w2.1 <- OmicCircos:::scale.v(po2.1, v1, v2, v3, v4)
            w2.2 <- OmicCircos:::scale.v(po2.2, v1, v2, v3, v4)
            if (po1.1 > po1.2 & po2.2 < po2.1) {
                draw.link.pg(xc, yc, r, w1.2, w1.1, w2.2, w2.1, 
                  col = the.color, lwd = lwd, s=link.curvature)
            }
            else if (po1.1 < po1.2 & po2.2 > po2.1) {
                draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, 
                  col = the.color, lwd = lwd, s=link.curvature)
            }
            else {
                draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, 
                  col = the.color, lwd = lwd, s=link.curvature)
            }
        }
    }
    if (type == "highlight.link") {
        xc <- mapping[1]
        yc <- mapping[2]
        r <- mapping[3]
        w1.1 <- mapping[4]
        w1.2 <- mapping[5]
        w2.1 <- mapping[6]
        w2.2 <- mapping[7]
        draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col = col, 
            lwd = lwd, s=link.curvature)
    }
    if (type == "heatmap") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        if (cluster) {
            dat.d <- dist(t(dat.m))
            dat.h <- hclust(dat.d)
            dat.m <- dat.m[, dat.h$order]
        }
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        dat.c <- as.numeric(unlist(dat.m))
        dat.min <- min(dat.c, na.rm = T)
        dat.max <- max(dat.c, na.rm = T)
        dat.n <- matrix(dat.c, ncol = num.col, nrow = num.row)
        colnames(dat.n) <- colnames(dat.n)
        rbPal <- colorRampPalette(c("blue", "white", "red"))
        col.i <- rbPal(20)[cut(dat.c, breaks = 20)]
        col.m <- matrix(col.i, ncol = num.col, nrow = num.row)
        w <- W - W/25
        w.i <- w/length(dat.i)
        for (i in 1:nrow(dat.in)) {
            dat.i <- as.numeric(dat.m[i, ])
            chr.i <- which(chr.po[, 1] == dat.in[i, 1])
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            w.to <- OmicCircos:::scale.v(as.numeric(dat.in[i, 2]), v1, v2, 
                v3, v4)
            r.from <- r + w - w/25
            for (j in 1:length(dat.i)) {
                w.to <- OmicCircos:::scale.v(as.numeric(dat.in[i, 2]), v1, 
                  v2, v3, v4)
                r.to <- r.from - w.i
                OmicCircos:::draw.line(xc, yc, w.to, r.from, r.to, col = col.m[i, 
                  j], lwd = lwd)
                r.from <- r.to
            }
        }
        if (col.bar) {
            if (col.bar.po == "topleft") {
                color.bar(40, 700, 50, 780, v.min = dat.min, 
                  v.max = dat.max)
            }
            if (col.bar.po == "bottomright") {
                color.bar(700, 40, 780, 50, v.min = dat.min, 
                  v.max = dat.max)
            }
        }
        if (cluster) {
            heatmap.cluster(0, 20, 150, 70, dat.m)
        }
    }
    if (type == "heatmap2") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        if (cluster) {
            dat.d <- dist(t(dat.m))
            dat.h <- hclust(dat.d)
            dat.m <- dat.m[, dat.h$order]
        }
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        dat.c <- as.numeric(unlist(dat.m))
        dat.min <- min(dat.c, na.rm = T)
        dat.max <- max(dat.c, na.rm = T)
        dat.n <- matrix(dat.c, ncol = num.col, nrow = num.row)
        colnames(dat.n) <- colnames(dat.n)
        rbPal <- colorRampPalette(c("blue", "white", "red"))
        col.i <- rbPal(20)[cut(dat.c, breaks = 20)]
        col.m <- matrix(col.i, ncol = num.col, nrow = num.row)
        col.df <- cbind(dat.in[, c(1:2)], col.m)
        w <- W - W/25
        w <- w - w/5
        w5 <- w/5
        lwd1 <- 360/(num.row + 24) * 1.5
        w.i <- w/length(dat.i)
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            col.c <- subset(col.df, col.df[, 1] == chr.s)
            col.c <- col.c[order(as.numeric(col.c[, 2])), ]
            col.c <- col.c[, -c(1:2)]
            if (dim(dat)[1] == 0) {
                next
            }
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            v.angle <- v2 - v1
            v.wide <- v.angle/nrow(dat)
            w.po <- v1
            for (i in 1:nrow(dat)) {
                r.from <- r + w - w/25
                w.to <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                  v3, v4)
                for (j in 1:length(dat.i)) {
                  r.to <- r.from - w.i
                  if (v.wide > 0.5) {
                    OmicCircos:::draw.arc.pg(xc, yc, w.po, w.po + v.wide, 
                      r.from, r.to, col = col.c[i, j], border = col.c[i, 
                        j])
                  }
                  else {
                    OmicCircos:::draw.line(xc, yc, w.po, r.from, r.to, col = col.c[i, 
                      j], lwd = lwd1)
                  }
                  r.from <- r.to
                }
                OmicCircos:::draw.line(xc, yc, w.to, r + w + w5 * 2/3, r + 
                  w + w5, col = col, lwd = lwd)
                OmicCircos:::draw.line(xc, yc, w.po, r + w, r + w + w5/3, 
                  col = col, lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.to, w.po, r + w + w5 * 2/3, 
                  r + w + w5/3, col = col, lwd = lwd)
                w.po <- w.po + v.wide
            }
        }
        if (col.bar) {
            if (col.bar.po == "topleft") {
                color.bar(40, 700, 50, 780, v.min = dat.min, 
                  v.max = dat.max)
            }
            if (col.bar.po == "bottomright") {
                color.bar(790, 1, 800, 100, v.min = dat.min, 
                  v.max = dat.max)
            }
        }
        if (cluster) {
            heatmap.cluster(0, 20, 150, 70, dat.m)
        }
    }
    if (type == "box") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in <- na.omit(dat.in)
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- data.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        mat.quant <- apply(dat.m, 1, quantile, probs = c(0.05, 
            0.25, 0.5, 0.75, 0.95))
        lwd <- 360/nrow(dat.in)
        if (lwd < 3) {
            lwd <- 3
        }
        color1 <- rainbow(10, alpha = 0.5)[7]
        color2 <- rainbow(10)[7]
        for (i in 1:nrow(dat.in)) {
            chr.i <- which(dat.in[i, 1] == chr.po[, 1])
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            w.to <- OmicCircos:::scale.v(as.numeric(dat.in[i, 2]), v1, v2, 
                v3, v4)
            my.R1 <- r
            my.R2 <- r + W
            mat.num <- length(mat.quant[, i])
            if (mat.num < 5) {
                next
            }
            q1 <- OmicCircos:::scale.v(mat.quant[1, i], my.R1, my.R2, dat.min, 
                dat.max)
            q2 <- OmicCircos:::scale.v(mat.quant[2, i], my.R1, my.R2, dat.min, 
                dat.max)
            qm <- OmicCircos:::scale.v(mat.quant[3, i], my.R1, my.R2, dat.min, 
                dat.max)
            q3 <- OmicCircos:::scale.v(mat.quant[4, i], my.R1, my.R2, dat.min, 
                dat.max)
            q4 <- OmicCircos:::scale.v(mat.quant[5, i], my.R1, my.R2, dat.min, 
                dat.max)
            OmicCircos:::draw.line(xc, yc, w.to, q1, q4, col = color1, lwd = 0.2)
            OmicCircos:::draw.line(xc, yc, w.to, q2, q3, col = color1, lwd = lwd)
            OmicCircos:::draw.line(xc, yc, w.to, qm - 0.1, qm + 0.1, col = color2, 
                lwd = lwd)
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "quant90") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in <- na.omit(dat.in)
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- data.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        mat.quant <- apply(dat.m, 1, quantile, probs = c(0.05, 
            0.25, 0.5, 0.75, 0.95))
        dat.qu <- cbind(dat.in[, c(1, 2)], t(mat.quant))
        colors <- rainbow(10, alpha = 0.5)
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.qu, dat.qu[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            my.R1 <- r
            my.R2 <- r + W
            i <- 1
            q1.from <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                my.R2, dat.min, dat.max)
            q2.from <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                my.R2, dat.min, dat.max)
            qm.from <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, 
                my.R2, dat.min, dat.max)
            q3.from <- OmicCircos:::scale.v(as.numeric(dat[i, 6]), my.R1, 
                my.R2, dat.min, dat.max)
            q4.from <- OmicCircos:::scale.v(as.numeric(dat[i, 7]), my.R1, 
                my.R2, dat.min, dat.max)
            w.from <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                v3, v4)
            for (i in 2:nrow(dat)) {
                q1.to <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                  my.R2, dat.min, dat.max)
                q2.to <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                  my.R2, dat.min, dat.max)
                qm.to <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, 
                  my.R2, dat.min, dat.max)
                q3.to <- OmicCircos:::scale.v(as.numeric(dat[i, 6]), my.R1, 
                  my.R2, dat.min, dat.max)
                q4.to <- OmicCircos:::scale.v(as.numeric(dat[i, 7]), my.R1, 
                  my.R2, dat.min, dat.max)
                w.to <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                  v3, v4)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, q1.from, q1.to, 
                  col = colors[1], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, qm.from, qm.to, 
                  col = colors[4], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, q4.from, q4.to, 
                  col = colors[7], lwd = lwd)
                w.from <- w.to
                q1.from <- q1.to
                q2.from <- q2.to
                qm.from <- qm.to
                q3.from <- q3.to
                q4.from <- q4.to
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "quant75") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in <- na.omit(dat.in)
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- data.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        mat.quant <- apply(dat.m, 1, quantile, probs = c(0.05, 
            0.25, 0.5, 0.75, 0.95))
        dat.qu <- cbind(dat.in[, c(1, 2)], t(mat.quant))
        colors <- rainbow(10, alpha = 0.5)
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.qu, dat.qu[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            my.R1 <- r
            my.R2 <- r + W
            i <- 1
            q1.from <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                my.R2, dat.min, dat.max)
            q2.from <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                my.R2, dat.min, dat.max)
            qm.from <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, 
                my.R2, dat.min, dat.max)
            q3.from <- OmicCircos:::scale.v(as.numeric(dat[i, 6]), my.R1, 
                my.R2, dat.min, dat.max)
            q4.from <- OmicCircos:::scale.v(as.numeric(dat[i, 7]), my.R1, 
                my.R2, dat.min, dat.max)
            w.from <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                v3, v4)
            for (i in 2:nrow(dat)) {
                q1.to <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                  my.R2, dat.min, dat.max)
                q2.to <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                  my.R2, dat.min, dat.max)
                qm.to <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, 
                  my.R2, dat.min, dat.max)
                q3.to <- OmicCircos:::scale.v(as.numeric(dat[i, 6]), my.R1, 
                  my.R2, dat.min, dat.max)
                q4.to <- OmicCircos:::scale.v(as.numeric(dat[i, 7]), my.R1, 
                  my.R2, dat.min, dat.max)
                w.to <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                  v3, v4)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, q2.from, q2.to, 
                  col = colors[1], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, qm.from, qm.to, 
                  col = colors[4], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, q3.from, q3.to, 
                  col = colors[7], lwd = lwd)
                w.from <- w.to
                q1.from <- q1.to
                q2.from <- q2.to
                qm.from <- qm.to
                q3.from <- q3.to
                q4.from <- q4.to
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ci95") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in <- na.omit(dat.in)
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- data.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        get.conf.int <- function(x) t.test(x)$conf.int
        mat.ci <- apply(dat.m, 1, get.conf.int)
        mat.me <- apply(dat.m, 1, mean, rm.na = T)
        dat.ci <- cbind(dat.in[, c(1, 2)], t(mat.ci), mat.me)
        colors <- rainbow(10, alpha = 0.5)
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.ci, dat.ci[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            my.R1 <- r
            my.R2 <- r + W
            i <- 1
            ci1.from <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                my.R2, dat.min, dat.max)
            ci2.from <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                my.R2, dat.min, dat.max)
            m.from <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, my.R2, 
                dat.min, dat.max)
            w.from <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                v3, v4)
            for (i in 2:nrow(dat)) {
                ci1.to <- OmicCircos:::scale.v(as.numeric(dat[i, 3]), my.R1, 
                  my.R2, dat.min, dat.max)
                ci2.to <- OmicCircos:::scale.v(as.numeric(dat[i, 4]), my.R1, 
                  my.R2, dat.min, dat.max)
                m.to <- OmicCircos:::scale.v(as.numeric(dat[i, 5]), my.R1, 
                  my.R2, dat.min, dat.max)
                w.to <- OmicCircos:::scale.v(as.numeric(dat[i, 2]), v1, v2, 
                  v3, v4)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, ci1.from, ci1.to, 
                  col = colors[1], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, m.from, m.to, 
                  col = colors[4], lwd = lwd)
                OmicCircos:::draw.line3(xc, yc, w.from, w.to, ci2.from, ci2.to, 
                  col = colors[7], lwd = lwd)
                w.from <- w.to
                ci1.from <- ci1.to
                ci2.from <- ci2.to
                m.from <- m.to
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "sv") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        dat.me <- mean(as.numeric(dat.m), na.rm = T)
        var1 <- function(x) var(x, na.rm = T)
        all.var <- apply(dat.m, 1, var1)
        var.min <- min(all.var, na.rm = T)
        var.max <- max(all.var, na.rm = T)
        num.col <- ncol(dat.in)
        the.cex <- 360/nrow(dat.in)
        if (the.cex < 3) {
            the.cex <- 3
        }
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            for (i in 1:nrow(dat)) {
                my.v <- as.numeric(dat[i, dat.i])
                d.v <- var(my.v, na.rm = T)
                d.m <- mean(my.v, na.rm = T)
                if (d.m > dat.me) {
                  col <- rainbow(10, alpha = 0.5)[1]
                }
                else {
                  col <- rainbow(10, alpha = 0.5)[7]
                }
                po <- as.numeric(dat[i, 2])
                w.to <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                the.v <- OmicCircos:::scale.v(d.m, my.R1, my.R2, dat.min, 
                  dat.max)
                the.c <- OmicCircos:::scale.v(d.v, 0.1, the.cex, var.min, 
                  var.max)
                OmicCircos:::draw.point.w(xc, yc, the.v, w.to, col = col, 
                  cex = the.c)
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "s.sd") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        dat.me <- mean(as.numeric(dat.m), na.rm = T)
        sd1 <- function(x) sd(x, na.rm = T)
        all.sd <- apply(dat.m, 1, sd1)
        sd.min <- min(all.sd, na.rm = T)
        sd.max <- max(all.sd, na.rm = T)
        num.col <- ncol(dat.in)
        the.cex <- 360/nrow(dat.in)
        if (the.cex < 3) {
            the.cex <- 3
        }
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            for (i in 1:nrow(dat)) {
                my.v <- as.numeric(dat[i, dat.i])
                sd.v <- sd(my.v, na.rm = T)
                d.m <- mean(my.v, na.rm = T)
                if (d.m > dat.me) {
                  col <- rainbow(10, alpha = 0.5)[1]
                }
                else {
                  col <- rainbow(10, alpha = 0.5)[7]
                }
                po <- as.numeric(dat[i, 2])
                w.to <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                the.v <- OmicCircos:::scale.v(d.m, my.R1, my.R2, dat.min, 
                  dat.max)
                the.c <- OmicCircos:::scale.v(sd.v, 0.1, the.cex, sd.min, 
                  sd.max)
                OmicCircos:::draw.point.w(xc, yc, the.v, w.to, col = col, 
                  cex = the.c)
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ss") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        dat.me <- mean(as.numeric(dat.m), na.rm = T)
        num.col <- ncol(dat.in)
        the.cex <- 360/nrow(dat.in)
        if (the.cex < 3) {
            the.cex <- 3
        }
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            for (i in 1:nrow(dat)) {
                my.v <- as.numeric(dat[i, dat.i])
                d.m <- mean(my.v, na.rm = T)
                if (d.m > dat.me) {
                  col <- rainbow(10, alpha = 0.5)[1]
                }
                else {
                  col <- rainbow(10, alpha = 0.5)[7]
                }
                po <- as.numeric(dat[i, 2])
                w.to <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                the.v <- OmicCircos:::scale.v(d.m, my.R1, my.R2, dat.min, 
                  dat.max)
                the.c <- OmicCircos:::scale.v(d.m, 0.1, the.cex, dat.min, 
                  dat.max)
                OmicCircos:::draw.point.w(xc, yc, the.v, w.to, col = col, 
                  cex = the.c)
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ml") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        if (length(col) == num.col) {
            colors <- col
        }
        else {
            colors <- rainbow(num.col, alpha = 0.5)
        }
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            col.i <- 0
            for (j in col.v:ncol(dat)) {
                col.i <- col.i + 1
                col <- colors[col.i]
                my.v <- as.numeric(dat[1, j])
                dat.i.old <- my.v
                v.old <- OmicCircos:::scale.v(my.v, my.R1, my.R2, dat.min, 
                  dat.max)
                po <- as.numeric(dat[1, 2])
                w.from <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                for (k in 1:nrow(dat)) {
                  dat.i <- as.numeric(dat[k, j])
                  if (is.na(dat.i)) {
                    next
                  }
                  v <- OmicCircos:::scale.v(dat.i, my.R1, my.R2, dat.min, 
                    dat.max)
                  w.to <- OmicCircos:::scale.v(as.numeric(dat[k, 2]), v1, 
                    v2, v3, v4)
                  if (w.from > 0) {
                    OmicCircos:::draw.line3(xc, yc, w.from, w.to, v.old, v, 
                      col = col, lwd = lwd)
                  }
                  dat.i.old <- dat.i
                  w.from <- w.to
                  v.old <- v
                }
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ml2") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        if (cutoff == "n") {
            if (length(col) == num.col) {
                colors <- col
            }
            else {
                colors <- rainbow(num.col, alpha = 0.5)
            }
        }
        else {
            if (length(col) == 2) {
                colors <- col
            }
            else {
                colors <- rainbow(10, alpha = 0.5)[c(1, 7)]
            }
        }
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            col.i <- 0
            for (j in col.v:ncol(dat)) {
                my.v <- as.numeric(dat[1, j])
                if (cutoff == "n") {
                  col.i <- col.i + 1
                  col <- colors[col.i]
                }
                else {
                  if (my.v > cutoff) {
                    col <- colors[1]
                  }
                  else {
                    col <- colors[2]
                  }
                }
                dat.i.old <- my.v
                v.old <- OmicCircos:::scale.v(my.v, my.R1, my.R2, dat.min, 
                  dat.max)
                po <- as.numeric(dat[1, 2])
                w.from <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                for (k in 1:nrow(dat)) {
                  dat.i <- as.numeric(dat[k, j])
                  if (is.na(dat.i)) {
                    next
                  }
                  v <- OmicCircos:::scale.v(dat.i, my.R1, my.R2, dat.min, 
                    dat.max)
                  w.to <- OmicCircos:::scale.v(as.numeric(dat[k, 2]), v1, 
                    v2, v3, v4)
                  OmicCircos:::draw.arc.s(xc, yc, v, w.from, w.to, col = col, 
                    lwd = lwd)
                  dat.i.old <- dat.i
                  w.from <- w.to
                  v.old <- v
                }
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ml3") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        if (length(col) == 2) {
            colors <- col
        }
        else {
            colors <- rainbow(10, alpha = 0.5)[c(1, 7)]
        }
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            for (j in col.v:ncol(dat)) {
                my.v <- as.numeric(dat[1, j])
                dat.i.old <- my.v
                if (is.na(my.v) == T) {
                  next
                }
                v.old <- OmicCircos:::scale.v(my.v, my.R1, my.R2, dat.min, 
                  dat.max)
                po <- as.numeric(dat[1, 2])
                w.from <- OmicCircos:::scale.v(po, v1, v2, v3, v4)
                for (k in 1:nrow(dat)) {
                  dat.i <- as.numeric(dat[k, j])
                  if (is.na(dat.i)) {
                    next
                  }
                  if (dat.i > cutoff) {
                    col <- colors[1]
                  }
                  else {
                    col <- colors[2]
                  }
                  v <- OmicCircos:::scale.v(dat.i, my.R1, my.R2, dat.min, 
                    dat.max)
                  w.to <- OmicCircos:::scale.v(as.numeric(dat[k, 2]), v1, 
                    v2, v3, v4)
                  if (v.old > 0) {
                    tmp.v <- dat.i * dat.i.old
                    if (tmp.v < 0) {
                      the0 <- OmicCircos:::scale.v(0, my.R1, my.R2, dat.min, 
                        dat.max)
                      if (dat.i > cutoff) {
                        col <- colors[1]
                        OmicCircos:::draw.line(xc, yc, w.from, v, the0, col = col, 
                          lwd = lwd, lend = 2)
                        col <- colors[2]
                        OmicCircos:::draw.line(xc, yc, w.from, the0, v.old, 
                          col = col, lwd = lwd, lend = 2)
                      }
                      else {
                        col <- colors[2]
                        OmicCircos:::draw.line(xc, yc, w.from, v, the0, col = col, 
                          lwd = lwd, lend = 2)
                        col <- colors[1]
                        OmicCircos:::draw.line(xc, yc, w.from, the0, v.old, 
                          col = col, lwd = lwd, lend = 2)
                      }
                    }
                    else {
                      OmicCircos:::draw.line(xc, yc, w.from, v, v.old, col = col, 
                        lwd = lwd, lend = 2)
                    }
                  }
                  OmicCircos:::draw.arc.s(xc, yc, v, w.from, w.to, col = col, 
                    lwd = lwd)
                  dat.i.old <- dat.i
                  w.from <- w.to
                  v.old <- v
                }
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "ms") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        for (chr.i in 1:chr.num) {
            chr.s <- chr.po[chr.i, 1]
            chr.s <- gsub("chr", "", chr.s)
            dat <- subset(dat.in, dat.in[, 1] == chr.s)
            dat <- dat[order(as.numeric(dat[, 2])), ]
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            if (B) {
                OmicCircos:::draw.arc.pg(xc, yc, v1, v2, r, r + W - 5, col = colors()[245])
            }
            else {
                OmicCircos:::draw.arc.s(xc, yc, r, v1, v2, col = colors()[245], 
                  lwd = lwd)
            }
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- as.matrix(dat.m)
        num.col <- ncol(dat.m)
        num.row <- nrow(dat.m)
        my.R1 <- R + W/5
        my.R2 <- R + W - W/5
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        dat.med <- median(as.numeric(dat.m), na.rm = T)
        colors <- rainbow(10, alpha = 0.3)[c(1, 7)]
        for (j in 1:num.col) {
            for (i in 1:num.row) {
                dat.i <- as.numeric(dat.m[i, j])
                if (is.na(dat.i)) {
                  next
                }
                chr.i <- which(chr.po[, 1] == dat.in[i, 1])
                if (dat.i > dat.med) {
                  col <- colors[1]
                }
                else {
                  col <- colors[2]
                }
                v <- OmicCircos:::scale.v(dat.i, my.R1, my.R2, dat.min, dat.max)
                v1 <- as.numeric(chr.po[chr.i, 2])
                v2 <- as.numeric(chr.po[chr.i, 3])
                v3 <- as.numeric(chr.po[chr.i, 6])
                v4 <- as.numeric(chr.po[chr.i, 7])
                w.to <- OmicCircos:::scale.v(as.numeric(dat.in[i, 2]), v1, 
                  v2, v3, v4)
                OmicCircos:::draw.point.w(xc, yc, v, w.to, col = col, cex = 0.1)
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "hist") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        dat.in[, 1] <- gsub("chr", "", dat.in[, 1])
        dat.in[dat.in[, 1] == 23, 1] <- "X"
        dat.in[dat.in[, 1] == 24, 1] <- "Y"
        dat.i <- c(col.v:ncol(dat.in))
        dat.m <- dat.in[, dat.i]
        dat.m <- data.matrix(dat.m, )
        dat.min <- min(as.numeric(dat.m), na.rm = T)
        dat.max <- max(as.numeric(dat.m), na.rm = T)
        dat.df <- as.data.frame(t(dat.m))
        hist.max <- ncol(dat.m) * 0.75
        cut.num <- 16
        region.num <- cut.num - 1
        breaks <- seq(dat.min, dat.max, length.out = cut.num)
        out <- lapply(dat.df, cut, br = breaks, labels = c(1:region.num))
        color1 <- rainbow(10, alpha = 0.5)[7]
        w.i <- (W - W/10)/region.num
        for (i in 1:nrow(dat.in)) {
            out.i <- table(out[[i]])
            chr.i <- which(dat.in[i, 1] == chr.po[, 1])
            v1 <- as.numeric(chr.po[chr.i, 2])
            v2 <- as.numeric(chr.po[chr.i, 3])
            v3 <- as.numeric(chr.po[chr.i, 6])
            v4 <- as.numeric(chr.po[chr.i, 7])
            w.to <- OmicCircos:::scale.v(as.numeric(dat.in[i, 2]), v1, v2, 
                v3, v4)
            OmicCircos:::draw.line(xc, yc, w.to, r, r + W - 5, col = color1, 
                lwd = 0.01)
            my.R <- r
            for (i in 1:region.num) {
                my.R <- my.R + w.i
                the.n <- out.i[i]
                if (the.n == 0) {
                  next
                }
                w.from <- OmicCircos:::scale.v(the.n, w.to, w.to + 360/nrow(dat.m), 
                  0, hist.max - 1)
                OmicCircos:::draw.arc.s(xc, yc, my.R, w.from, w.to, col = "lightblue", 
                  lwd = w.i, lend = 1)
            }
        }
        if (scale) {
            OmicCircos:::do.scale(xc, yc, dat.min, dat.max, R, W - W/5)
        }
    }
    if (type == "hl") {
        if (class(mapping)[1] == "GRanges") {
            po <- as.integer((start(mapping) + end(mapping))/2 + 
                0.5)
            dat.in <- cbind(as.character(seqnames(mapping)), 
                po, names(mapping), mcols(mapping))
        }
        else {
            dat.in <- mapping
        }
        r1 <- as.numeric(dat.in[1])
        r2 <- as.numeric(dat.in[2])
        seg1 <- dat.in[3]
        seg2 <- dat.in[5]
        po1 <- as.numeric(dat.in[4])
        po2 <- as.numeric(dat.in[6])
        the.col1 <- dat.in[7]
        the.col2 <- dat.in[8]
        seg.i <- which(seg1 == chr.po[, 1])
        v1 <- as.numeric(chr.po[seg.i, 2])
        v2 <- as.numeric(chr.po[seg.i, 3])
        v3 <- as.numeric(chr.po[seg.i, 6])
        v4 <- as.numeric(chr.po[seg.i, 7])
        w1 <- OmicCircos:::scale.v(as.numeric(po1), v1, v2, v3, v4)
        seg.i <- which(seg2 == chr.po[, 1])
        v1 <- as.numeric(chr.po[seg.i, 2])
        v2 <- as.numeric(chr.po[seg.i, 3])
        v3 <- as.numeric(chr.po[seg.i, 6])
        v4 <- as.numeric(chr.po[seg.i, 7])
        w2 <- OmicCircos:::scale.v(as.numeric(po2), v1, v2, v3, v4)
        OmicCircos:::draw.arc.pg(xc, yc, w1, w2, r1, r2, col = the.col1, border = the.col2)
    }
}


draw.link.pg<-function (xc, yc, r, w1.1, w1.2, w2.1, w2.2, col = col, lwd = lwd,s=2) 
{
    w1 <- w1.1
    w2 <- w2.2
    w3 <- (w1 + w2)/2
    w1 <- w1/360 * 2 * pi
    w2 <- w2/360 * 2 * pi
    w3 <- w3/360 * 2 * pi
    x0 <- xc + r * cos(w1)
    y0 <- yc - r * sin(w1)
    x1 <- xc + r * cos(w2)
    y1 <- yc - r * sin(w2)
    x <- c(x0, rep(xc,s), x1)
    y <- c(y0, rep(yc,s), y1)
    bc1 <- bezierCurve(x, y, 100)
    ang.d <- abs(w1.1 - w1.2)
    pix.n <- ang.d * 10
    if (pix.n < 10) {
        pix.n <- 10
    }
    ang.seq <- rev(seq(w1.1, w1.2, length.out = pix.n))
    ang.seq <- ang.seq/360 * 2 * pi
    fan.1.x <- xc + cos(ang.seq) * r
    fan.1.y <- yc - sin(ang.seq) * r
    w1 <- w1.2
    w2 <- w2.1
    w3 <- (w1 + w2)/2
    w1 <- w1/360 * 2 * pi
    w2 <- w2/360 * 2 * pi
    w3 <- w3/360 * 2 * pi
    x0 <- xc + r * cos(w1)
    y0 <- yc - r * sin(w1)
    x1 <- xc + r * cos(w2)
    y1 <- yc - r * sin(w2)
    x <- c(x0, rep(xc,s), x1)
    y <- c(y0, rep(yc,s), y1)
    bc2 <- bezierCurve(x, y, 100)
    ang.d <- abs(w2.1 - w2.2)
    pix.n <- ang.d * 10
    if (pix.n < 10) {
        pix.n <- 10
    }
    ang.seq <- rev(seq(w2.1, w2.2, length.out = pix.n))
    ang.seq <- ang.seq/360 * 2 * pi
    fan.2.x <- xc + cos(ang.seq) * r
    fan.2.y <- yc - sin(ang.seq) * r
    polygon(c(bc1$x, fan.2.x, rev(bc2$x), fan.1.x), c(bc1$y, 
        fan.2.y, rev(bc2$y), fan.1.y), fillOddEven = T, 
        border = col, 
        col = col)
}

draw.arc.s<-function (xc, yc, r, w1, w2, col = "lightblue", lwd = 1, lend = 0) 
{
    ang.d <- abs(w1 - w2)
    pix.n <- ang.d * 5
    if (pix.n < 2) {
        pix.n <- 2
    }
    ang.seq <- rev(seq(w1, w2, length.out = pix.n))
    ang.seq <- ang.seq/360 * 2 * pi
    fan.i.x <- xc + cos(ang.seq) * r
    fan.i.y <- yc - sin(ang.seq) * r
    lines(fan.i.x, fan.i.y, col = col, lwd = lwd, type = "l", 
        lend = lend)
}


draw.text.rt<-function (xc, yc, r, w, n, col = "black", cex = 1, side = "out",
                        rot.ang=NULL,shift=12) 
{
    #r = r-5
    w <- w%%360
#    the.o <- w
    the.w <- 360 - w
    w <- w/360 * 2 * pi
    x <- xc + r * cos(w)
    y <- yc - r * sin(w)
    num2 <- shift
#    if (side == "out") {
#        if (the.w <= 90) {
#            the.pos <- 4
#        }
#        else if (the.w > 90 & the.w <= 180) {
#            the.w <- the.w + 180
#            the.pos <- 2
#        }
#        else if (the.w > 180 & the.w <= 270) {
#            the.w <- the.w%%180
#            the.pos <- 2
#        }
#        else if (the.w > 270 & the.w <= 360) {
#            the.pos <- 4
#        }
#        if (the.pos == 2) {
#            x <- x + num2
#        }
#        if (the.pos == 4) {
#            x <- x - num2
#        }
#    }
#    if (side == "in") {
        if (the.w <= 90) {
            the.pos <- 4
            x <- x - num2 #*cos(w)
#            y <- y - num2*sin(w)
        }
        else if (the.w > 90 & the.w <= 180) {
            the.w <- the.w + 180
            the.pos <- 2
            x <- x + num2 #*sin(w)
#            y <- y + num2*cos(w)
        }
        else if (the.w > 180 & the.w <= 270) {
            the.w <- the.w%%180
            the.pos <- 2
           x <- x + num2#*sin(w)
#            y <- y + num2*cos(w)
        }
        else if (the.w > 270 & the.w <= 360) {
            the.pos <- 4
            x <- x - num2#*cos(w)
#            y <- y - num2*sin(w)
        }
#    }
    if (!is.null(rot.ang)){
        the.w<-the.w+rot.ang
    }
    text(as.double(x), as.double(y), adj = 0, offset = 1, labels = n, srt = as.double(the.w), 
        pos = the.pos, col = col, cex = cex)
}



bezierCurve<-function (x, y, n = 10) 
{
    out<-sapply(seq(0, 1, length.out = n),function(t) bez(x, y, t))
    return(list(x = out["x",], y = out["y",]))
}

bez<-function (x, y, t) 
{
    outx <- 0
    outy <- 0
    n <- length(x) - 1
    for (i in 0:n) {
        outx <- outx + choose(n, i) * ((1 - t)^(n - i)) * t^i * 
            x[i + 1]
        outy <- outy + choose(n, i) * ((1 - t)^(n - i)) * t^i * 
            y[i + 1]
    }
    return(c(x = outx, y = outy))
}


