#RAS, TP53, RB1, FGFR3 (RYC.1A) te jednak olewamy bo to analiza mutacyjna

#KRT1,CD44, CDH3, KRT14, KRT16, KRT5, KRT6A, KRT6B, KRT6C, (RYC1B PART1) - basal_markers_UP

#FGFR3, FOXA1 GPX2, ERBB2, ERBB3, CYP2J2, GATA3, PPARG, KRT19, KRT7, KRT8, FABP4, KRT20, CD24, KRT18, XBP1 (RYC1B PART2) - luminal_markers_UP
#(CD44, KRT5, KRT6, KRT14, CDH3) i (CD24, FOXA1, GATA3, ERBB2, ERBB3, XBP1, KRT20) - też olewamy bo zawierają się w powyższych gr.
# 
#
#KRT14, DSG3, KRT6B, KRT5, KRT6A, KRT6C, LOC653499, LOC728910, P13, S100A7 - basal_top_10_UP_MIBC ACTG2, CNN1, MYH11, MFAP4, PGM 5, FLNC, ACTC1, DES, PCP4,DMN - p53_top_10_UP_MIBC
#
#MAL, FMO9P, BHMT, SNX31, KRT20, SPINK1, DHRS2, UPK2, UPK1A, VSIG2 - luminal_top_10_UP_MIBC.
#
#Dwie pierwsze grupy t niby markery dla basal i lumi, trzy ostatnie to tylko top 10 które wyciągnęli już po podziale na grupy.

marker.genes.list<-list(
basal.top.Choi.et.al. = c('KRT14','DSG3', 'KRT6B','KRT5', 'KRT6A', 'KRT6C', 'LOC653499', 'LOC728910', 'P13', 'S100A7'),
p53.top.Choi.et.al. = c( 'ACTG2', 'CNN1', 'MYH11', 'MFAP4', 'PGM5', 'FLNC', 'ACTC1', 'DES', 'PCP4', 'DMN'),
luminal.top.Choi.et.al.=c( 'MAL', 'FMO9P', 'BHMT', 'SNX31', 'KRT20', 'SPINK1', 'DHRS2', 'UPK2', 'UPK1A', 'VSIG2'),
#basal.suspected.upregulated = c('CD44', 'KRT5', 'KRT6', 'KRT14', 'CDH3'),
#luminal.suspected.upregulated = c('CD24', 'FOXA1', 'GATA3', 'ERBB2', 'ERBB3', 'XBP1', 'KRT20'),
#ras.group = c('RAS', 'TP53', 'RB1', 'FGFR3'),
basal.markers.UP = c('KRT1', 'CD44', 'CDH3', 'KRT14', 'KRT16', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C'),
luminal.markers.UP = c('FGFR3', 'FOXA1', 'GPX2', 'ERBB2', 'ERBB3', 'CYP2J2', 'GATA3', 'PPARG', 'KRT19', 'KRT7', 'KRT8', 'FABP4', 'KRT20', 'CD24', 'KRT18', 'XBP1','UPK1B','UPK3A'),
#PNAS - bez tych do final
Damrauer.luminalvsbasal.markers.base47 = c('SCNN1B','PPARG','TOX3','GATA3','HMGCS2','RAB15',
'AHNAK2','ADIRF','SEMA5A','CHST15','TRAK1','SCNN1G','MT1X','TMPRSS2','VGLL1',
'TBX2','UPK1A','GAREM','BHMT','SPINK1','GPD1L','RNF128','CYP2J2','EMP3','GDPD3','FBP1',
'MSN','MT2A','CDK6','ALOX5AP','PRRX1','SLC27A2','TMEM97','CD14','PLEKHG6','CYP4B1',
'GLIPR1','PDGFC','PRKCDBP','FAP','CAPN5','PALLD','TUBB6','SLC9A2','PPFIBP2','FAM174B'),
#http://onlinelibrary.wiley.com/doi/10.1110/ps.073196308/full
higurashi2008.hub.proteins = c('rnaSA', 'KBP1A', 'KBP1', 'KBP12', 'AN', 'PD1', 'DL235C', '0790', 'D4', 'psF', 'ps6', 'CL2', 'CP1', 'CYA2', 'NASE1', 'NS1', 'PB5', 'PA7', 'PC9', 'BR154C', 'BR1204', 'OU2F1', 'CT1', 'TF1', 'unx1', 'ml1', 'bfa2', 'ebp2ab', 'PO', 'POR', 'ag-pol', 'A', 'mm', 'eiE9', 'ol', 'ei', 'KP2', 'BXL1', 'KP1', 'MC19', 'CP2', 'KP1A', 'CEB1L', 'nf1a', 'nf-1', 'nf-1a', 'cf1', 'NAT1', 'NB1', 'NGT1', 'B1', 'B1', 'MY2', 'AN', 'RA24', 'K/SW-cl.81', 'PNB1', 'TF97', 'WF', '8VWF', 'heY', '1882', 'W1871', 'TNNB1', 'TNNB', 'K/SW-cl.35', 'RO2286', 'CF7L2', 'CF4', 'yp102A1', 'yp102', 'LA-B', 'LAB', '2M', 'DABP0092', 'DCMA22P', 'ya', 'XO1-122', 'XA0141', 'BAA_pXO1_0142', 'MAD2', 'ADH2', 'ADR2', 'la', 'laT-3', 'laT-4', 'laT-5', 'laT-6', 'L2', 'PA1', 'PA', 'EGFA', 'EGF', '', 'amp2', 'yb2', 'tx1a', 'ap', 'nap25', 'nap', 'nap25', 'nap', 'IR2DL1', 'D158A', 'KAT1', 'naN', '3701', 'W3678', 'inB', 'inP', '0231', 'W0221', 'ago', 'gn', 'G9401', 'su', '14', 'G8781', '35', 'P1BA', 'FN1', 'QCRC1', 'QCRC2', 'T-CYB', 'OB', 'YTB', 'TCYB', 'YC1', 'QCRFS1', 'QCRB', 'QCRQ', 'QCRH', 'QCRFS1', 'QCR10', 'QCRC1', 'QCRC2', 'T-CYB', 'OB', 'YTB', 'TCYB', 'YC1', 'QCRFS1', 'QCRB', 'QCRQ', 'QCRH', 'QCRFS1', 'QCR10', 'pr', 'AB5A', 'AB5', 'csA', 'kc1', 'RF1', 'YTH2', 'RNO', 'SCD2', 'SCD2L', 'rkaca', 'kaca', 'GF1', 'GFA', 'cp', 'sfl1c', 'auC', 'mi', 'co', 'ti', '2209', 'W2197', 'rss2', 'ry2', 'lpha-LP', 'nai1', 'nai-1', 'naSA', 'isF', 'M_1036', 'nfsf13', 'pril', 'GFBR1', 'LK5', 'KR4', 'AN', 'RP1', 'AP60', 'NL189W', '1606', 'SE1', 'GL238W', 'RC135', 'ta', 'A1148', 'PARG', 'R1C3', 'r0b2', 'hp', 'QCRC1', 'QCRC2', 'T-CYB', 'OB', 'YTB', 'TCYB', 'YC1', 'QCRFS1', 'QCRB', 'QCRQ', 'QCRH', 'QCRFS1', 'QCR10', 'CTA1', 'CTA', 'NASE1', 'NL1', 'ASF2', 'AVE2', 'P53', '53', 'D', 'S6', 'ufA', 'uf', 'THA1694', 'DK2', 'DKN2', 'CNA2', 'CN1', 'CNA', 'DC6', 'DC18L', 'TP5A1', 'TP5A2', 'TP5B', 'TP5C1', 'TP5C', 'TP5D', 'TP5E', 'L4', 'ASP3', 'PP32', 'ASP3', 'PP32', 'CP1', 'CP', 'PO', 'KR066C', 'us', 'G2944', 'ceb2', 'ceb1', 'CTA1', 'CTA', '3', 'L10', 'uc', 'YZ', 'APS', 'ACWR095', '3R', 'ntB'),
p63.Choi.et.al.=c(
'CCNG1','S100A4','XPC','IGFBP3','MYC','TNFRSF10A',
'ITGA3','SFN','S100A8','KRT14','KRT6C','PI3','KRT6A','S100A7','KRT1',
'CSTA','SERPINB2','CKS2','CCNB1','CCNA2','KIF23','PMAIP1','IGFBP6','PRNP',
'GFI1','RAC2','ADA','CDK6','CD82','TNC','SH2B3','PTPN12','LYN','PLAU','SERPINE1',
'IL1B','IL8','HBEGF','F3','IL1RAP','PTHLH'),
PPARgamma.Choi.et.al.=c(
'PEPD','KRT18','ACADM','MGST1','CYP4B1','ACAA1','CEBPA','ACOX1','ACSL1','KRT20','SCNN1G','ELOVL6','HES1',
'SDC1','IGFBP3','PPARG','GDF15','COL1A1','COL1A2','VIM','ACTA2','TPM2','RARRES2','TCF4','IGFBP5','IL6',
'SOCS3','CTGF','JUN','IGFBP6','CAV1','CES1','PC','CEBPB','MMP9','SERPINA1','TBXAS1')
)

excluded.tccs.cases <- c("8381670042_F", "9287078021_F", "8381670051_H", "9287078094_C", "9287078094_F", "9287078094_I", "9287078102_H", "9287078002_B", "9287078002_H", "9287078002_I", "9287078002_J", "9287078074_I", "9287078074_L", "9287078096_C", "9464921093_H", "9287078074_D", "9464921108_E", "9464921108_I", "9464921093_D"
#)
,"8381670042_H", "9287078096_H")
excluded.tccs.cases.19 <- c("8381670042_F", "9287078021_F", "8381670051_H", "9287078094_C", "9287078094_F", "9287078094_I", "9287078102_H", "9287078002_B", "9287078002_H", "9287078002_I", "9287078002_J", "9287078074_I", "9287078074_L", "9287078096_C", "9464921093_H", "9287078074_D", "9464921108_E", "9464921108_I", "9464921093_D"
)

marker.genes.df<-do.call(rbind,lapply(names(marker.genes.list),function(n){
    data.frame(gene=marker.genes.list[[n]],gene.type=n)
}))
#all.gene.colors<-rainbow(length(marker.genes.list))
all.gene.colors<-c("red","green",rgb(0.2,0.2,1),"darkred","darkblue","violet","purple","magenta","cyan")
marker.genes.df<-transform(marker.genes.df,gene.type=factor(gene.type,levels=names(marker.genes.list)))
marker.genes.df<-transform(marker.genes.df,color=all.gene.colors[gene.type],stringsAsFactors=FALSE)

marker.genes.df<-transform(marker.genes.df,
    bladder.cancer.obligatory=gene.type%in%c(
    "basal.top.Choi.et.al.","p53.top.Choi.et.al.","luminal.top.Choi.et.al.","basal.markers.UP","luminal.markers.UP",
    "Damrauer.luminalvsbasal.markers.base47","Damrauer.luminalvsbasal.markers.base47"
    ),
    choi.et.al. = gene.type%in%c("basal.top.Choi.et.al.","p53.top.Choi.et.al.","luminal.top.Choi.et.al.","basal.markers.UP","luminal.markers.UP")
)#,"higurashi2008.hub.proteins"))

marker.genes.df.uniq<-marker.genes.df[!duplicated(marker.genes.df$gene),]
rownames(marker.genes.df.uniq)<-marker.genes.df.uniq$gene

marker.genes.df.uniq.fromlast<-marker.genes.df[!duplicated(marker.genes.df$gene,fromLast=TRUE),]
rownames(marker.genes.df.uniq.fromlast)<-marker.genes.df.uniq.fromlast$gene


add.present.column<-function(x,name="present",genes.df=marker.genes.df){
    genes.df[,name]<-with(genes.df,gene%in%x)
    genes.df
}

get.color.for.gene<-function(x,genes.df=marker.genes.df.uniq){
#    marger.genes.df$
    genes.df[x,"color"]
}

get.marker.type.for.gene<-function(x,genes.df=marker.genes.df.uniq){
#    marger.genes.df$
    genes.df[x,"gene.type"]
}

######
.__this.dir <- dirname(parent.frame(2)$ofile) 
.__old.dir<-getwd()
setwd(.__this.dir) 
selected.for.validation<-tryCatch(read.csv('List_of_probes_ordered_for_validation.csv',as.is=T),
                            finally=setwd(.__old.dir))

setwd(.__this.dir) 
blca.gnames<-tryCatch(read.csv('BLCA_gnames.txt',sep="|",header=F,as.is=T),
                            finally=setwd(.__old.dir))
colnames(blca.gnames)<-c("gene","count")



