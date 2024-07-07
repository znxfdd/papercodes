
rm(list = ls())
library(AnnoProbe) 
suppressPackageStartupMessages(library(GEOquery)) 
gset=AnnoProbe::geoChina('GSE14407')
gset

eSet=gset[[1]]

probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
probes_expr=log2(probes_expr+1)
boxplot(probes_expr,las=2)

phenoDat <- pData(eSet)
head(phenoDat[,1:4])


(gpl=eSet@annotation)
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)


group_list=factor(c(rep('Normal',12),rep('Tumor',12)))
table(group_list)
library(limma)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)
write.table(DEG,'DEG.txt',row.names = T,sep='\t')


need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,1)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list)
deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('FDX1',genes_expr,group_list)
check_diff_genes('MPP6',genes_expr,group_list)



library(KEGGREST)
cg <- KEGGREST::keggGet("hsa03410")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)



library(ggplot2)
library(ggsignif)
library(ggpubr)
interest<-read.table("gene_2.txt",header = T)
data<-genes_expr[which(rownames(genes_expr) %in% interest$Gene),]
a <- t(data)
t <- as.data.frame(c(rep('Normal',12),rep('Tumor',12)))
colnames(t) <- 'group'
b=cbind(a,t)
c <- data.frame(gene=c(rep('FDX1',24),rep('LIPT1',24),rep('LIAS',24),rep('DLD',24),rep('DBT',24),rep('GCSH',24),rep('DLAT',24),rep('PDHA1',24),rep('PDHB',24),rep('SLC31A1',24),rep('ATP7A',24),rep('ATP7B',24),rep('MTF1',24),rep('GLS',24),rep('CDKN2A',24)),group=rep(b$group,15),expression=c(b$FDX1,b$LIPT1,b$LIAS,b$DLD,b$DBT,b$GCSH,b$DLAT,b$PDHA1,b$PDHB,b$SLC31A1,b$ATP7A,b$ATP7B,b$MTF1,b$GLS,b$CDKN2A))

c$group = factor(c$group)
p=ggplot(c, aes(x = gene, y = expression, fill = group)) + geom_boxplot() +
  scale_fill_manual(values=c("#00C3C6","#FF6C67"))+
  stat_compare_means(aes(group = group),method = "t.test",label = "p.signif") + 
  xlab(NULL) +
theme_classic()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

pdf('geneexpr.pdf',width = 10,height = 6)
print(p)
dev.off()


DEG$lab = as.factor(ifelse(DEG$adj.P.Val>=0.05,"",
                               
                               ifelse(DEG$adj.P.Val>=0.01&DEG$adj.P.Val<0.05,"*",
                                      
                                      ifelse(DEG$adj.P.Val>=0.001&DEG$adj.P.Val<0.01,"**", "***"))))



DEG$new <- paste(rownames(DEG),DEG$lab)
data<-DEG[which(rownames(DEG) %in% interest$Gene),]


library(dplyr)
data <- data %>% 
  arrange(row.names(.))
data <- data[rownames(data),]

rownames(data) <- data$new 


input <- genes_expr[which(rownames(genes_expr) %in% interest$Gene),]
input <- input %>% 
  arrange(row.names(.))
row.names(input) <- row.names(data)

pre_heatdata <- t(scale(t(input)))

pre_heatdata[pre_heatdata> 1] <- 1

pre_heatdata[pre_heatdata< -1] <- -1



annColors <- list()

annColors[['Level']] <- c('Normal'="#00C3C6",'Tumor'="#FF6C67")



library(pheatmap)
annotation=as.data.frame(b$group)
colnames(annotation) <- 'group'
row.names(annotation) <- row.names(b)

pdf(file="pheatmap.pdf",width = 10,height = 12)

pheatmap(pre_heatdata,
         
         color = colorRampPalette(c("blue",'white',"red"))(1000),
         
         annotation_col = annotation,
         
         annotation_colors = annColors,
         
         treeheight_row = 50,
         
      
         
         show_rownames = T,
         
         show_colnames = F,
         
         cluster_rows = T,
         
         cluster_cols = F)

dev.off()


gene <- genes_expr[which(rownames(genes_expr) %in% interest$Gene),]
gene <- gene[,13:ncol(gene)]
gene <- t(gene) 



gene <- log(gene+1)



gene_cor <- cor(gene, method = 'pearson')



diag(gene_cor) <- 0

gene_cor 

gene_cor <- reshape2::melt(gene_cor)

gene_cor <- subset(gene_cor, value != 0) 

head(gene_cor) 

library(circlize)
pdf('cor.pdf')
chordDiagram(gene_cor,
             
             annotationTrack = c('grid', 'name', 'axis'), 
             
             grid.col = c(GABRD = 'green3', PLVAP = 'red', CDKN3 = 'orange', CDC25C = 'purple', UBE2T = 'skyblue', SKA1 = 'blue'), 
             
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), 
             
             annotationTrackHeight = c(0.05, 0.05))
dev.off()

options(stringsAsFactors = F)
library(stringr)
project="TCGA-OV"
if(!dir.exists("clinical"))dir.create("clinical")
if(!dir.exists("expdata"))dir.create("expdata")
dir()


command1 <- "./gdc-client download -m gdc_manifest_clinical.txt -d clinical"
command2 <- "./gdc-client download -m gdc_manifest_expdata.txt -d expdata"

system(command = command1) 
system(command = command2) 

length(dir("./clinical/"))

length(dir("./expdata/"))


library(XML)
xmls = dir("clinical/",pattern = "*.xml$",recursive = T)
cl = list()
for(i in 1:length(xmls)){
  result = xmlParse(paste0("clinical/",xmls[[i]]))
  rootnode = xmlRoot(result)
  cl[[i]] = xmlToDataFrame(rootnode[2])
}
clinical = do.call(rbind,cl)
clinical[1:3,1:3]

count_files = dir("expdata/",pattern = "*.tsv$",recursive = T)

exp = list()
for(i in 1:length(count_files)){
  exp[[i]] = read.table(paste0("expdata/",count_files[[i]]),header=T,sep="\t")
  exp[[i]] = exp[[i]][-(1:4),] 
  exp[[i]] = exp[[i]]$fpkm_unstranded 
}
exp = as.data.frame(do.call(cbind,exp))

dim(exp)

exp[1:4,1:4]



meta = jsonlite::fromJSON("metadata.cart.json")
ID = sapply(meta$associated_entities,
            function(x){x$entity_submitter_id})
file2id = data.frame(file_name = meta$file_name,
                     ID = ID)

count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
table(count_files2 %in% file2id$file_name)


file2id = file2id[match(count_files2,file2id$file_name),]
identical(file2id$file_name,count_files2)


colnames(exp) = file2id$ID

gene_name = data.table::fread(paste0("expdata/",count_files[1]))$gene_name
gene_name = gene_name[-seq(1,4)] 
exp = cbind(gene_name=gene_name,exp)
dim(exp)

exp = exp[!duplicated(exp$gene_name),]
rownames(exp) = exp$gene_name
exp = exp[,-1]
dim(exp)



exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ] 
dim(exp)


library(stringr)
table(str_sub(colnames(exp),14,15))
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)



if(!dir.exists("data"))dir.create("data")
save(exp,clinical,Group,project,file = paste0("data/",project,"_gdc.Rdata"))


write.csv(exp,'mRNAmatrix.',sep = '\t') 
write.table(clinical,'clinical.txt',row.names = F,sep = '\t') 


library(limma)            
expFile="symbol.txt"      
geneFile="gene.txt"      



rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table(geneFile, header=T, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="m6aGeneExp.txt",sep="\t",quote=F,col.names=F)


library(limma)
corFilter=0.2          
pvalueFilter=0.05      



rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=read.table("m6aGeneExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
m6A=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
m6A=avereps(m6A)
m6A=m6A[rowMeans(m6A)>0.1,]


group=sapply(strsplit(colnames(m6A),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
m6A=m6A[,group==0]


outTab=data.frame()
for(i in row.names(lncRNA)){
  
  for(j in row.names(m6A)){
    x=as.numeric(lncRNA[i,])
    y=as.numeric(m6A[j,])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="postive"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="negative"))
      
    }
  }
}


pvalues <- outTab$pvalue
pvalue_adjust <- p.adjust(pvalues, method = "fdr", n = length(pvalues))


write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)

lncNode=data.frame(Node=unique(as.vector(outTab[,"lncRNA"])), Type="lncRNA")
mrnaNode=data.frame(Node=unique(as.vector(outTab[,"m6A"])), Type="m6A")
nodeOut=rbind(lncNode, mrnaNode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)


m6aLncRNA=unique(as.vector(outTab[,"lncRNA"]))
m6aLncRNAexp=data[m6aLncRNA,]
m6aLncRNAexp=rbind(ID=colnames(m6aLncRNAexp), m6aLncRNAexp)
write.table(m6aLncRNAexp,file="m6aLncExp.txt",sep="\t",quote=F,col.names=F)


library(igraph)               
nodefile="net.node.txt"      
edgefile="net.network.txt"   
outfile="network.pdf"        
lncRNAcol="#80CAA6"         
m6Acol="#EE4442"             



node.data=read.table(nodefile, header=T, sep="\t", check.names=F)
edge.data=read.table(edgefile, header=T, sep="\t", check.names=F)
color=ifelse(node.data$Type=="lncRNA", lncRNAcol, m6Acol)
value=ifelse(node.data$Type=="lncRNA", 3, 4)
fontSize=ifelse(node.data$Type=="lncRNA", 0.01, 1.1)
node=data.frame(id=node.data$Node,label=node.data$Node,color=color,shape="dot",value=value,fontSize=fontSize)
edge=data.frame(from=edge.data$lncRNA,to=edge.data$m6A,length=100,arrows="middle",smooth=TRUE,shadow=FALSE,weight=edge.data$cor)


d=data.frame(p1=edge$from, p2=edge$to, weight=abs(edge$weight))
g=graph.data.frame(d,directed = FALSE)
E(g)$color="grey60"
V(g)$size=node$value[match(names(components(g)$membership),node$label)]
V(g)$shape="sphere"
V(g)$lable.cex=node$fontSize[match(names(components(g)$membership),node$label)]
V(g)$color=node$color[match(names(components(g)$membership),node$label)]


pdf(outfile,width=28,height=28)
layout(mat=matrix(c(1,2,1,2),nc=2), height=c(1,7))
par(mar=c(0,5,3,5))
plot(1,type="n",axes=F,xlab="",ylab="")
legend('center',legend=c('LncRNA','Cuproptosis-related gene'),col=c(lncRNAcol,m6Acol),pch=16,bty="n",ncol=2,cex=5)
vertex.frame.color = node$color
edge_col=E(g)$color
plot(g,layout=layout.fruchterman.reingold,vertex.size=V(g)$size,vertex.label=node$label,vertex.label.cex=V(g)$lable.cex,edge.width =0.1,edge.arrow.size=0,vertex.label.color='yellow',vertex.frame.color='black',edge.color=edge_col,vertex.label.font=50)
dev.off()


library(limma)               
lncFile="m6aLncExp.txt"      
cliFile="time.txt"          


rt=read.table(lncFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data=t(data)
data=avereps(data)


cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     


sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)


library(survival)      
pFilter=0.01         
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365   


outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP) )
  }
}


write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)



bioForest=function(coxFile=null, forestFile=null){

  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  

  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'p-value',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="gray80",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, "red", "green")
  points(as.numeric(hr), n:1, pch = 18, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}

bioForest(coxFile="uniCox.txt", forestFile="forest.pdf")



library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
lncFile="uniSigExp.txt"     
expFile="lncRNA.txt"       

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


lncRNA=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[colnames(lncRNA)[3:ncol(lncRNA)],]
exp=data


cl <- read.table('clinical.txt',header=T )
newid <-substring(colnames(data),1,12) 
colnames(data) <- newid
newid <- as.data.frame(newid)
colnames(newid) <- 'id'
newdata <- merge(newid,cl,by.x = 'id',sort = T)
df<- data[,colnames(data) %in% newdata$id] 
early <- newdata[which(newdata$group=='I-II'),]
late <- newdata[which(newdata$group=='III-IV'),]
sampleType=ifelse(colnames(df) %in% early$id,1,2)
conNum=length(which(sampleType=='1'))
treatNum=length(which(sampleType=='2'))
df.early <- df[,early$id]
df.late <- df[,late$id]
data <- cbind(df.early,df.late)


sigVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ sampleType)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec


Type=c(rep("Early",conNum), rep("Late",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type,colnames(data)) 
data=log2(data+1)
pdf("heatmap.pdf", width=7.5, height=4.7)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()


exp=as.data.frame(t(df))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Early", "Late")
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")


p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            
            xlab="",
            legend.title="Type",
            palette = c("blue", "red"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")


pdf(file="boxplot.pdf", width=7.5, height=5)
print(p1)
dev.off()


library(limma)
library(ConsensusClusterPlus)


rt=read.table("uniSigExp.txt", header=T, sep="\t", check.names=F, row.names=1)
data=rt[,(3:ncol(rt))]
data=t(data)


maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


clusterNum=2        
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)


library(survival)
library(survminer)
clusterFile="cluster.txt"    
cliFile="time.txt"           


cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365


sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,], Cluster=cluster[sameSample,])
rt$Cluster=paste0("Cluster", rt$Cluster)


length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line="hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()


library(pheatmap)             
ClusterFile="cluster.txt"       
cliFile="clinical.txt"        
expFile="uniSigExp.txt"        


Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")
Cluster=Cluster[order(Cluster$Cluster),,drop=F] 


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)


samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
Type=cbind(Cluster, cli)
Type$Cluster=paste0("Cluster", Type$Cluster)


sigVec=c("Cluster")
for(clinical in colnames(Type[,2:ncol(Type)])){
  data=Type[c("Cluster", clinical)]
  colnames(data)=c("Cluster", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))

}
colnames(Type)=sigVec

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(exp[,(3:ncol(exp))])
data=data[,row.names(Type)]


colorList=list()

bioCol=c("#0066FF","#FF0000","#FF9900","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
  cliLength=length(levels(factor(Type[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(Type[,cli]))
  if("unknow" %in% levels(factor(Type[,cli]))){
    cliCol["unknow"]="grey75"}
  colorList[[cli]]=cliCol
}


data=log2(data+1)
pdf("heatmap.pdf", width=12, height=10)
pheatmap(data,
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


library(limma)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"            
clusterFile="cluster.txt"      
gene="CDKN2A"                   
showName="CDKN2A"               


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=t(data[gene,,drop=F])


group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])    
Type=c(rep(1,conNum), rep(2,treatNum))


exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")


group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(exp$Type)))]

boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                  xlab="",
                  ylab=paste0(showName, " expression"),
                  legend.title="Type",
                  palette = bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()



group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,,drop=F]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)


cluster=read.table(clusterFile, header=F, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
data=data.frame(gene=data[sameSample,], Cluster=cluster[sameSample,])
data$Cluster=paste0("Cluster", data$Cluster)

group=levels(factor(data$Cluster))
data$Cluster=factor(data$Cluster, levels=group)
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cluster)))]


boxplot=ggboxplot(data, x="Cluster", y="gene", color="Cluster",
                  xlab="",
                  ylab=paste0(showName, " expression"),
                  legend.title="Cluster",
                  palette = bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file=paste0(gene,".Cluster.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

library(limma)
library(corrplot)
expFile="symbol.txt"       
lncFile="uniSigExp.txt"  
gene="CD274"               
showName="PD-L1"           


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]


lncRT=read.table(lncFile, header=T, sep="\t", check.names=F, row.names=1)
lncRNA=c(gene, colnames(lncRT)[3:ncol(lncRT)])
sameGene=intersect(lncRNA, row.names(data))
data=data[sameGene,]
row.names(data)[1]=showName


data=t(data)
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)


pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()



library("limma")         
expFile="symbol.txt"    



rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


v <-voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        


source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

library(limma)
library(estimate)
inputFile="symbol.txt"  


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
out=data[rowMeans(data)>0,]


out=rbind(ID=colnames(out), out)
write.table(out, file="uniq.symbol.txt", sep="\t", quote=F, col.names=F)


filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds="commonGenes.gct",
              output.ds="estimateScore.gct")


scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
scores=scores[,1:3]
out=rbind(ID=colnames(scores), scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)


library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"     
cluFile="cluster.txt"                  
pFilter=0.05         


immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)


cluster=read.table(cluFile, header=F, sep="\t", row.names=1, check.names=F)


lowName=row.names(cluster)[cluster[,1]==1]
highName=row.names(cluster)[cluster[,1]==2]


lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)


outTab=data.frame()
pdf("vioplot.pdf", width=13, height=9)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Cluster1", "Cluster2"),
       lwd=4.5,bty="n",cex=1.5,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()


write.table(outTab,file="diff.result.txt",sep="\t",row.names=F,quote=F)


library(limma)
library(ggpubr)
immuneFile="CIBERSORT-Results.txt"     
cluFile="cluster.txt"                
pFilter=0.05       

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)


cluster=read.table(cluFile, header=F, sep="\t", row.names=1, check.names=F)


sameSample=intersect(row.names(immune), row.names(cluster))
immune1=immune[sameSample,,drop=F]
cluster1=cluster[sameSample,,drop=F]
colnames(cluster1)=c("Cluster")
data=cbind(immune1, cluster1)
data$Cluster=paste0("Cluster", data$Cluster)


type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cluster)))]



colnames(data)[9] <- 'T cells regulatory'


for(i in colnames(data)[1:(ncol(data)-1)]){

  boxplot=ggboxplot(data, x="Cluster", y=i, fill="Cluster",
                    xlab="",
                    ylab=i,
                    legend.title="Cluster",
                    palette=bioCol
  )+ 
    stat_compare_means(comparisons=my_comparisons)

  pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}


library(limma)
library(ggpubr)
scoreFile="scores.txt"       
cluFile="cluster.txt"      


score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

score$rownames=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score <- score[-which(duplicated(score$rownames)==T),]
row.names(score) <- score$rownames
score <- subset(score,select=-rownames)
score=avereps(score)


cluster=read.table(cluFile, header=F, sep="\t", row.names=1, check.names=F)


sameSample=intersect(row.names(score), row.names(cluster))
score1=score[sameSample,,drop=F]
cluster1=cluster[sameSample,,drop=F]
colnames(cluster1)=c("Cluster")
data=cbind(score1, cluster1)
data$Cluster=paste0("Cluster", data$Cluster)


type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$Cluster)))]


for(i in colnames(data)[1:(ncol(data)-1)]){

  boxplot=ggboxplot(data, x="Cluster", y=i, fill="Cluster",
                    xlab="",
                    ylab=i,
                    legend.title="Cluster",
                    palette=bioCol
  )+ 
    stat_compare_means(comparisons=my_comparisons)

  pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}

exp <- read.table('symbol.txt',header = F) 
colnames(exp) <- exp[1,]
exp=exp[-1,]
colnames(exp) [2:ncol(exp)]<- substring(colnames(exp)[2:ncol(exp)],1,12)

group <- ifelse(duplicated(colnames(exp)),1,2)
exp <- exp[,-which(group==1)]

da <- t(read.table('cluster.txt',header = F,row.names = 1))
exp=exp[,-which(ifelse(colnames(exp)[2:ncol(exp)] %in% colnames(da),1,2)=='2')]

colnames(exp)[1] <- 'Name'
exp$Description <- rep('na',nrow(exp))
exp <- cbind(exp[1],exp['Description'],exp[2:(ncol(exp)-1)])
n_gene=nrow(exp)
n_sample=ncol(exp)-2
data2 <- rbind(c('#1.2',rep('',ncol(exp)-1)),c(paste(n_gene), paste(n_sample), rep('',ncol(exp)-2) ),colnames(exp), exp[1:nrow(exp),])
write.table(data2,file = 'Cluster.gct',col.names = F,quote = F, row.names = F,sep = '\t')     

a <- read.table('cluster.txt',header = F)
a=t(a)[-1,]
a=ifelse(a=='1','C1','C2')
a=as.data.frame(t(a))
data3 <- rbind(c(ncol(exp)-2,'2','1',rep('',ncol(exp)-5)),c('#','C1',"C2",rep('',ncol(exp)-5)),a)
write.table(data3,file = 'Cluster.cls',col.names = F,quote = F,row.names = F,sep = '\t')

library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

rt=read.table("uniSigExp.txt", header=T, sep="\t", check.names=F, row.names=1)    

inTrain<-createDataPartition(y=rt[,2], p=0.5, list=F) 
train<-rt[inTrain,]
test<-rt[-inTrain,]
trainOut=cbind(id=row.names(train), train)
testOut=cbind(id=row.names(test), test)
  
x=as.matrix(train[,c(3:ncol(train))])
y=data.matrix(Surv(train$futime,train$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
 
trainFinalGeneExp=train[,lassoGene]

trainScore=predict(cvfit, newx=as.matrix(train[,c(3:ncol(train))]), s="lambda.min", type="response")
outCol=c("futime", "fustat", lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore), "high", "low"))
train=cbind(train[,outCol], riskScore=as.vector(trainScore), risk)
trainRiskOut=cbind(id=rownames(train), train)
  
 
testFinalGeneExp=test[,lassoGene]
 
testScore=predict(cvfit, newx=as.matrix(test[,c(3:ncol(test))]), s="lambda.min", type="response")
outCol=c("futime", "fustat", lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore), "high", "low"))
test=cbind(test[,outCol], riskScore=as.vector(testScore), risk)
testRiskOut=cbind(id=rownames(test), test)
  
 
diff=survdiff(Surv(futime, fustat) ~risk, data=train)
pValue=1-pchisq(diff$chisq, df=1)
diffTest=survdiff(Surv(futime, fustat) ~risk, data=test)
pValueTest=1-pchisq(diffTest$chisq, df=1)
  

predictTime=1   
roc=timeROC(T=train$futime, delta=train$fustat,
              marker=trainScore, cause=1,
              weighting='aalen',
              times=c(predictTime), ROC=TRUE)
rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=testScore, cause=1,
                  weighting='aalen',
                  times=c(predictTime), ROC=TRUE)	
  
  
write.table(trainOut,file="train.data.txt",sep="\t",quote=F,row.names=F)
write.table(testOut,file="test.data.txt",sep="\t",quote=F,row.names=F)
   
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
  
write.table(geneCoef, file="geneCoef.txt", sep="\t", quote=F, row.names=F)
write.table(trainRiskOut,file="trainRisk.txt",sep="\t",quote=F,row.names=F)
write.table(testRiskOut,file="testRisk.txt",sep="\t",quote=F,row.names=F)
  
allRiskOut=rbind(trainRiskOut, testRiskOut)
write.table(allRiskOut,file="allRisk.txt",sep="\t",quote=F,row.names=F)
    
library(survival)
library(survminer)


bioSurvival=function(inputFile=null,outFile=null){
 
  rt=read.table(inputFile, header=T, sep="\t")

  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  

  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("red", "blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")

library(survival)
library(survminer)
library(timeROC)


bioROC=function(inputFile=null, rocFile=null){
  predictTime=1   

  rt=read.table(inputFile, header=T, sep="\t")

  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(predictTime), ROC=TRUE)
  pdf(file=rocFile, width=5, height=5)
  plot(ROC_rt, time=predictTime, col='red', title=FALSE, lwd=2)
  
  legend('bottomright', cex=1.3,
         paste0('AUC=',sprintf("%.03f",ROC_rt$AUC[2])),
         col="white", lwd=1, bty = 'n')
  dev.off()
}

bioROC(inputFile="trainRisk.txt",rocFile="train.ROC.pdf")
bioROC(inputFile="testRisk.txt",rocFile="test.ROC.pdf")


library('survivalROC')

sRocFuction=function(inputFile=null){
  rt=read.table(inputFile, header=T, sep="\t")
  par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)
  sROC=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
  plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
       xlab="1-Specificity", ylab="Sensitivity",
       lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
  sROC3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM")
  lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
  aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") 
  sROC5=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")
  lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
  aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") 
  legend("bottomright", c(aucText,aucText3,aucText5),
         lwd=2,bty="n",col=c("red","green","blue"))
}
pdf(file="train.ROC.pdf",width=6,height=5)
sRocFuction(inputFile="trainRisk.txt") 
dev.off() 

pdf(file="test.ROC.pdf",width=6,height=5) 
sRocFuction(inputFile="testRisk.txt") 
dev.off() 

library(pheatmap)        

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)  
  rt=rt[order(rt$riskScore),]     
  
 
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile, width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  

  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile, width=7, height=4)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  

  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  ann_colors = list(
    type = c(high="red", low="blue")
  )
  pdf(file=heatmapFile, width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation, 
           annotation_colors = ann_colors,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("green",3.5), "white", rep("red",3.5)))(50),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

bioRiskPlot(inputFile="trainRisk.txt",riskScoreFile="train.riskScore.pdf",survStatFile="train.survStat.pdf",heatmapFile="train.heatmap.pdf")

bioRiskPlot(inputFile="testRisk.txt",riskScoreFile="test.riskScore.pdf",survStatFile="test.survStat.pdf",heatmapFile="test.heatmap.pdf")

library(survival)        

bioForest=function(coxFile=null,forestFile=null,forestCol=null){

  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  

  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
 
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'p-value',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="gray80",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 18, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}



indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)  
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     
  

  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
 
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")
  
  

  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}



indep(riskFile="trainRisk.txt",
      cliFile="clinical.txt",
      uniOutFile="train.uniCox.txt",
      multiOutFile="train.multiCox.txt",
      uniForest="train.uniForest.pdf",
      multiForest="train.multiForest.pdf")



indep(riskFile="testRisk.txt", 
      cliFile="clinical.txt",
      uniOutFile="test.uniCox.txt",
      multiOutFile="test.multiCox.txt",
      uniForest="test.uniForest.pdf",
      multiForest="test.multiForest.pdf")


library(survival)
library(survminer)
riskFile="allRisk.txt"     
cliFile="clinical.txt"     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)   
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)    


sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=risk[,1],fustat=risk[,2],cli,risk=risk[,"risk"])

for(i in colnames(data[,3:(ncol(data)-1)])){
  rt=data[,c("futime","fustat",i,"risk")]
  rt=rt[(rt[,i]!="unknow"),]
  colnames(rt)=c("futime","fustat","clinical","risk")
  tab=table(rt[,"clinical"])
  tab=tab[tab!=0]
 
  for(j in names(tab)){
    rt1=rt[(rt[,"clinical"]==j),]
    tab1=table(rt1[,"risk"])
    tab1=tab1[tab1!=0]
    labels=names(tab1)
    if(length(labels)==2){
      titleName=j
      if((i=="age") | (i=="Age") | (i=="AGE")){
        titleName=paste0("age",j)
      }
      diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
      pValue=1-pchisq(diff$chisq,df=1)
      if(pValue<0.001){
        pValue="p<0.001"
      }else{
        pValue=paste0("p=",sprintf("%.03f",pValue))
      }
      fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
     
      surPlot=ggsurvplot(fit, 
                         data=rt1,
                         conf.int=F,
                         pval=pValue,
                         pval.size=6,
                         title=paste0("Patients with ",titleName),
                         legend.title="Risk",
                         legend.labs=labels,
                         font.legend=12,
                         xlab="Time(years)",
                         break.time.by = 1,
                         palette=c("red", "blue"),
                         risk.table=TRUE,
                         risk.table.title="",
                         risk.table.col = "strata",
                         risk.table.height=.25)
    
      j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
      pdf(file=paste0("survival.",i,"_",j,".pdf"),onefile = FALSE,
          width = 6,      
          height =5)        
      print(surPlot)
      dev.off()
    }
  }
}

library(limma)
library(pheatmap)
ClusterFile="cluster.txt"   
cliFile="clinical.txt"       
riskFile="allRisk.txt"       
scoreFile="scores.txt"       


Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)


name <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score <- score[-which(duplicated(gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score)))),]
name1 <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
row.names(score)=name1
score=avereps(score)


samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
risk=risk[samSample,,drop=F]
score=score[samSample,,drop=F]
score[,"ImmuneScore"]=ifelse(score[,"ImmuneScore"]>median(score[,"ImmuneScore"]), "High", "Low")
data=cbind(risk, Cluster, score[,"ImmuneScore",drop=F], cli)
data=data[order(data$riskScore),,drop=F]    
Type=data[,(ncol(risk):ncol(data))]     
exp=data[,(3:(ncol(risk)-2))]      
Type$Cluster=paste0("Cluster", Type$Cluster)


sigVec=c("risk")
for(clinical in colnames(Type[,2:ncol(Type)])){
  data=Type[c("risk", clinical)]
  colnames(data)=c("risk", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))

}
colnames(Type)=sigVec


colorList=list()

bioCol=c("#FF0000","#0066FF","#0066FF","#FF0000","#FF9900","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
  cliLength=length(levels(factor(Type[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(Type[,cli]))
  if("unknow" %in% levels(factor(Type[,cli]))){
    cliCol["unknow"]="grey75"}
  colorList[[cli]]=cliCol
}

pdf("heatmap.pdf", height=10, width=9)
pheatmap(t(exp),
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

library(limma)
library(ggpubr)
ClusterFile="cluster.txt"    
cliFile="clinical.txt"      
riskFile="allRisk.txt"        
scoreFile="scores.txt"       


Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)


name <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score <- score[-which(duplicated(gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score)))),]
name1 <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
row.names(score)=name1
score=avereps(score)


samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
risk=risk[samSample,,drop=F]
score=score[samSample,,drop=F]
score[,"ImmuneScore"]=ifelse(score[,"ImmuneScore"]>median(score[,"ImmuneScore"]), "High", "Low")
data=cbind(risk, Cluster, score[,"ImmuneScore",drop=F], cli)
rt=data[order(data$riskScore),,drop=F] 
rt=rt[,((ncol(risk)-1):ncol(rt))]
rt=rt[,-2]
rt$Cluster=paste0("Cluster", rt$Cluster)


for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
 
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
 
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab=clinical,
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
 
  pdf(file=paste0(clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}

library(limma)
library(ggpubr)


gene=showName='GLS'


riskGene=function(riskFile=null, expFile=null, boxFile=null){

  
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0,]
  

  group=sapply(strsplit(colnames(data),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group=gsub("2", "1", group)
  data=data[,group==0]
  

  data=rbind(data, gene=data[gene,])
  exp=t(data[c("gene",gene),])
  row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
  exp=avereps(exp)

 
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  
  
  sameSample=intersect(row.names(exp), row.names(risk))
  exp=exp[sameSample,]
  
  risk=risk[sameSample,]
  data=cbind(as.data.frame(exp), as.data.frame(risk))
  
  
  data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
  group=levels(factor(data$risk))
  data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
 
  
  boxplot=ggboxplot(data, x="risk", y="gene", color="risk",
                    xlab="",
                    ylab=paste0(showName, " expression"),
                    legend.title="",
                    palette = c("blue", "red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  

  pdf(file=boxFile, width=5, height=4.5)
  print(boxplot)
  dev.off()
}


riskGene(riskFile="allRisk.txt", expFile="symbol.txt", boxFile=paste0(gene,".boxplot.pdf"))


library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
immFile="CIBERSORT-Results.txt"       
riskFile="allRisk.txt"             
pFilter=0.05     


immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(immune),row.names(risk))
immune1=immune[sameSample,]
risk1=risk[sameSample,]


outTab=data.frame()
x=as.numeric(risk1[,"riskScore"])

for(j in colnames(immune1)){
  y=as.numeric(immune1[,j])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    p1=ggplot(df1, aes(x, y)) + 
      xlab("Risk score")+ ylab(j)+
      geom_point(color='black',alpha=0.7)+ geom_smooth(method="lm",formula=y~x) +
      stat_cor(method = 'spearman', aes(x =x, y =y),color='red')+ theme_bw()
    if(pValue<pFilter){
      pdf(file=paste0(j,".pdf"), width=5, height=4.6)
      print(p1)
      dev.off()
      outTab=rbind(outTab,cbind(Cell=j, pValue))
    }
  }
}
write.table(outTab,file="immuneCor.result.txt",sep="\t",row.names=F,quote=F)

