library(compiler)

F2norm<-cmpfun(function(M)sqrt(sum(M^2)))

thr.op<-cmpfun(function(x,thr){sign(x)*pmax(abs(x)-thr,0)})

La.svd.cmp<-cmpfun(La.svd)

thresh.nuclear<-cmpfun(function(M,thr){
    s<-La.svd.cmp(M)
    dd<-thr.op(s$d,thr)
    id<-which(dd!=0)
    s$d<-dd[id]
    s$u<-s$u[,id,drop=FALSE]
    s$vt<-s$vt[id,,drop=FALSE]
    s$L<-s$u%*%(s$d*s$vt)
    s
})


thresh.l1<-cmpfun(function(M,thr){
    thr.op(M,thr)
})



rpca<-function(M,
             lambda = 1/sqrt(max(dim(M))),
             mu = prod(dim(M))/(4*sum(abs(M))),
             term.delta=10^(-7),
             max.iter=5000,
             trace=FALSE,
             thresh.nuclear.fun=thresh.nuclear){
    dm<-dim(M)
    term.norm<-term.delta*F2norm(M)
    S<-Yimu<-matrix(0,nrow=dm[1],ncol=dm[2])

    imu<-1/mu
    limu<-lambda/mu

    i<-0
    stats<-c()
    converged<-FALSE
    while(TRUE){
        i<-i+1
        L.svd<-thresh.nuclear.fun(M-S+Yimu,imu)
        L<-L.svd$L
        S<-thresh.l1(M-L+Yimu,limu)

        MLS = M-L-S

        resid.norm<-F2norm(MLS)
        stats<-c(stats,resid.norm)
        if (trace) 
            print(c(iter=i,resid.norm=resid.norm))
        converged<-resid.norm<term.norm
        if ((i>max.iter)||converged)
            break;

        Yimu = Yimu + MLS
    }
    final.delta<-resid.norm*term.delta/term.norm
    if (!converged)
        warning(paste("rpca did not converge after",i,"iterations, final.delta=",final.delta))
    list(L=L,S=S,L.svd=L.svd,
         convergence=list(converged=converged,
                          iterations=i,
                          final.delta=final.delta,
                          all.delta=stats))
}

#require(gputools)

thresh.nuclear.gpu<-cmpfun(function(M,thr){
    s<-gpuSvd(M)
    dd<-thr.op(s$d,thr)
    id<-which(dd!=0)
    s$d<-dd[id]
    s$u<-s$u[,id,drop=FALSE]
    s$v<-s$v[,id,drop=FALSE]
    s$vt<-t(s$v)
    s$L<-gpuMatMult(s$u,s$d*s$vt)
    s
})



rpca.gpucula<-function(...,term.delta=10^(-5),gpu.to.choose=NULL){
    stopifnot(require(gputools))
    if (!exists("gpuSvd")) 
        stop("library(gputools) needs to be compiled with --cula, and its version needs to implement gpuSvd")
    if (!is.null(gpu.to.choose))
        chooseGpu(gpu.to.choose)    
    rpca(...,,term.delta=term.delta,thresh.nuclear.fun=thresh.nuclear.gpu)
}





