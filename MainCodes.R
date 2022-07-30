library("ggplot2",quietly = TRUE)
library("gplots",quietly = TRUE)
library("pheatmap",quietly = TRUE)
library("RColorBrewer",quietly = TRUE)
library("reshape2",quietly = TRUE)
library("ggpubr",quietly = TRUE)
library("cowplot",quietly = TRUE)
library(ggsci,quietly = TRUE)



#############################################
##  Fig2

load("diversity_MM.RData")
source("functions.friedman.R")
friedman.shannon <- friedman(data=diversity_MM, "shannon")
friedman.simpson <- friedman(data=diversity_MM, "simpson")

options(repr.plot.width = 5, repr.plot.height = 4)
ggboxplot(diversity_MM, x = "Stage2", y = "shannon", add = "jitter",color = "Stage2", fill = "Stage2",
          palette =pal_nejm()(8)[c(6,4,3,2,1)],alpha=0.1, lwd=1, outlier.shape = NA,
          font.label = list(size = 16, color = "black"),) +
  stat_pvalue_manual(friedman.shannon$pwc, hide.ns = TRUE,tip.length = 0.02) +
  labs(
    subtitle = rstatix::get_test_label(friedman.shannon$res.fried,  detailed = TRUE),
    caption = rstatix::get_pwc_label(friedman.shannon$pwc))+
  theme(legend.position = "none")+
  scale_y_continuous("Shannon")+
  scale_x_discrete("")+
  my_theme2+
  theme(title = element_text(size=8))


load("otu_MM.RData")
library("GGally")
data3 <- dcast(melt(otu_MM[,-c(3:4)], ids=c("subject","Stage2")), subject+variable ~ Stage2, median)

options(repr.plot.width = 5, repr.plot.height = 4)
ggcorr(data3[,-c(1,2)], nbreaks = 10, method = c("pairwise", "spearman"),
            label_round = 2,limits = c(0,1),
            geom = "tile", #text, blank, tile, circle
            label = TRUE, label_size = 4, label_color = "white", 
            hjust = 0.75, size = 4.5, color = "black", 
            low = pal_lancet("lanonc")(9)[2],
            high = pal_lancet("lanonc")(9)[1],
            midpoint = 0)



load("relat.abun_phylum_MM.RData")
tmp1 <- melt(dcast(melt(relat.abun_MM[,-c(1,3:4)], ids="Stage2"),Stage2~variable, mean))
tmp1 <- tmp1[which(tmp1$variable!="k__Archaea.p__Euryarchaeota"),]
tmp1$variable <- as.character(tmp1$variable)
tmp1$variable <- unlist(lapply(tmp1$variable,function(x){strsplit(x,"p__")[[1]][2]}))
tmp1$variable <- factor(tmp1$variable, levels = c("Firmicutes","Proteobacteria","Bacteroidetes","Actinobacteria","Fusobacteria","TM7",".Thermi.","Synergistetes","Verrucomicrobia","Tenericutes","Cyanobacteria"))

options(repr.plot.width = 7, repr.plot.height = 4)
ggplot(tmp1, 
       aes(x=Stage2, 
           y=value/10000, 
           stratum = variable, 
           alluvium = variable,
           fill=variable))+
  ggalluvial::geom_flow(width = 0.5, alpha=0.4) +
  ggalluvial::geom_stratum(width = 0.5, cex=0.2, alpha=0.8)+
  scale_fill_manual(values = c(colorRampPalette(pal_lancet()(8))(27),pal_simpsons()(7)[c(1:4)]))+
  theme_cowplot()+
  scale_x_discrete("")+
  scale_y_continuous("Relative abundance (%)")


load("qiime2_net.avg.change.RData")
options(repr.plot.width = 7, repr.plot.height = 4)
ggplot(qiime2_net.avg.change,aes(x=log10(Net.Avg.Change+1), y=importance))+
  geom_point(aes(color = id3), size=2)+
  geom_vline(xintercept = 0, color="black", cex=0.2)+
  theme_test()+
  scale_color_lancet()+
  scale_y_log10()+
  theme(legend.position = "right", 
        panel.border  = element_rect(colour = "black", size=0.6),
        axis.ticks = element_line(size=0.8, ),
        text = element_text(color = "black", size =15),
        axis.text.x = element_text(vjust = -1),
        axis.ticks.length = unit(1.5, "mm"))


# Process of taxomomy files

load("meta.RData")
load("meta2.RData")
load("gg_level1-7.RData")


# calculate relative abundance 
sumdata <- apply(data,1,sum)
data_relabun <- data/sumdata*1000000
write.table(data.frame(ID=rownames(data_relabun), data_relabun), file = "gg_level1-7-relabun.txt", quote = F, sep="\t", row.names = F)

#remove samples with sum < 2000 
freq <-  read.table("sample-frequency-detail.csv", header = F, stringsAsFactors = F, colClasses = c("character","numeric"), row.names = 1, sep=",")
excldeID <- paste("X",rownames(freq)[which(freq$V2<2000)], sep="")
data_relabun2 <- data_relabun[!rownames(data_relabun) %in% excldeID,]

# calculate average value for repeated measured samples 
data0 <- merge(x=meta[,-1], y=data_relabun2, by="row.names")
rownames(data0) <- data0$Row.names
data0$Row.names <- NULL
data_averaged <- dcast(melt(data0, id.vars = c(colnames(data0)[c(1:8)])), subject+Stage2+disease1+disease2+CRSgrade+Outcome1+Outcome2+Outcome3 ~ variable, mean)
rownames(data_averaged) <- paste(data_averaged$subject,data_averaged$Stage2,data_averaged$disease1, sep=".")

write.table(data.frame(ID=rownames(data_averaged), data_averaged), file = "gg_level1-7-relabun_averaged_filter.freq.2000.txt", quote = F, sep="\t", row.names = F)




data <- read.table("gg_level1-7-relabun_averaged_filter.freq.2000-filter.unclassified2.txt", header = T, row.names=1, stringsAsFactors = F)
id = names(which(table(diversity_MM$subject) == 5))
re.MM.g <- friedman.taxon(data, "MM", "g")
re1 <- re.MM.g$re.friedman
re2 <- re.MM.g$re.paired

re1.sig <- subset(re1, p.fdr<0.05) 
re1.sig$id <- exid(rownames(re1.sig), "g")
re1.sig[,c(1:4)] <- apply(re1.sig[,c(1:4)] , 2, as.numeric)
re1.sig <- re1.sig[order(re1.sig$p.fdr),]
re1.sig$id <- factor(re1.sig$id, levels = rev(re1.sig$id))
 ggplot(data=re1.sig, aes(x=id,y=effect.size))+
  geom_col(fill=pal_lancet()(3)[1], alpha=0.9) +
  coord_flip()+
  theme_classic()+
  my_theme2+
  scale_x_discrete("")+
  scale_y_continuous("Effect size")



## calculate log2fc
re3 <- re2[rownames(re1.sig),]
fc.genus <-  fc(data, "MM", "g")
fc.genus.sig <- fc.genus[rownames(re3),]
fc.genus.sig[re3>0.05] <- 0
#re3[re3>0.05] <- 0
rownames(fc.genus.sig) = rownames(re3) <- exid(rownames(fc.genus.sig),"g")

## heatmap
re4 <- re3
re4[re4>0.05] <- ""
pheatmap(as.matrix(fc.genus.sig), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         # border_color ="white",
         color = unique(c(colorRampPalette(colors = c(pal_lancet()(9)[1],"white"))(length(seq(-5,-0.1,by=0.1))),"white",colorRampPalette(colors = c("white",pal_lancet()(9)[7]))(length(seq(0,10,by=1))))),
         breaks = unique(c(seq(-5,-0.1,by=0.1),seq(0,10,by=1))),
         show_rownames = T,
         show_colnames = T,
         legend = T, 
         annotation_legend=F,
         annotation_names_col = T,
         display_numbers = re4,
         fontsize_number = 10,
         number_color = "black" )



#############################################
## fig 3

load("MM_diversity.RData")
ggplot(data=data_averaged[which(data_averaged$Outcome1!="VGPR"),], aes(y=shannon , x=Stage2))+
  geom_boxplot(aes(y=shannon, x=Stage2, col=Outcome1), alpha=0.1, cex=0.4 ,outlier.color = NA, width=0.7)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), cex=0.3, col="darkgray")+
  geom_point(aes(y=shannon, x=Stage2, col=Outcome1), size=1.4, shape=16, alpha=0.4, position=position_jitterdodge(0.3))+
  geom_smooth(aes(y=shannon, x=Stage2, group=Outcome1, color=Outcome1, fill=Outcome1), cex=0.8 , alpha=0.1,se=F)+
  stat_compare_means(method="wilcox",  aes(label= ..p.format.., group=Outcome1,comparisons = list(c("CR","PR"))), label.y = max(data_averaged$shannon)+0.2)+
  scale_x_discrete("")+
  scale_y_continuous("shannon index", limits = c(min(data_averaged$shannon)-0.05, max(data_averaged$shannon)+0.4))+
  theme_bw(base_rect_size = 0.8)+
  scale_color_lancet()+
  theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank(), text=element_text(size=14), axis.ticks.x = element_blank())+
  ggtitle("wilcox.test")



### OTU-based PCoA
data <-read.table("feature-table-relabun_averaged_filter.freq.2000.txt", header = T, stringsAsFactors = F, row.names = 1)
dt <- subset(data, disease1=="MM"& Stage2=="CRSb"& Outcome1 %in% c("CR","PR"))

mt <- dt[,-c(1:8)]
distance <- vegan::vegdist(mt, method = 'canberra')
meta_sub_MM <- meta2[rownames(mt),]
pcoa <- cmdscale(distance, k = (nrow(mt) - 1), eig = TRUE)
plot_data <- data.frame({pcoa$point})[1:2]
plot_data$Sample_name <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig = pcoa$eig
plot_data$Outcome1 <- meta2[rownames(plot_data),"Outcome1"]

options(repr.plot.width = 6, repr.plot.height = 4)
ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, color=Outcome1)) +
  geom_point( size=2) +
  #stat_chull(geom="polygon",aes(fill=Outcome1), alpha=0.3) +
  stat_ellipse(level = 0.9, size=0.7)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme_classic()+
  scale_color_lancet()+
  my_theme2

# adonis test

vegan::adonis(distance ~ Outcome1, meta_sub_MM, permutations = 9999) 


######## maSigPro
library("maSigPro")
library("MASS")
load("otu_table_0.97.RData")
source("masigpro.R")
# extract MM
data_MM <- data[which(data$disease1=="MM"),]
data_MM1 <- data_MM[which(data_MM$Outcome1 %in% c("PR","CR")),]
#remove 90% samples = 0
temp <- apply(data_MM1[,c(9:ncol(data_MM1))], 2 , function (x) {sum(x==0)})
temp2 <- temp[which(temp > nrow(data_MM1)*0.9)]
data_MM11 <- data_MM1[,-which(colnames(data_MM1) %in% names(temp2))]
data_averaged <- data_MM11
data_averaged$Stage2 <- factor(data_averaged$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc"))
data_averaged <- data_averaged[order(data_averaged$Stage2),]

# maSigPro PR vs CR
data_averaged_sig2 <- data_averaged[which(data_averaged$Outcome1 %in% c("PR","CR")),]
meta <- subset(data_averaged_sig2, select = "Stage2")
meta$Stage2 <- as.numeric(factor(meta$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc")))
meta$Replicate <- meta$Stage2
meta$group <- data_averaged_sig2$Outcome1
meta$CR <- data_averaged_sig2$Outcome1
meta$PR <- data_averaged_sig2$Outcome1
meta$PR[which(meta$PR=="PR")] <- 1
meta$PR[which(meta$PR=="CR")] <- 0
meta$CR[which(meta$CR=="PR")] <- 0
meta$CR[which(meta$CR=="CR")] <- 1
meta$group <- NULL
meta$Replicate <- paste(meta$Stage2,meta$CR, sep="")
meta$Replicate <- paste(meta$Replicate,meta$PR, sep="")
meta$Replicate <- as.numeric(as.factor(meta$Replicate))
meta <- as.matrix(apply(meta,2,as.numeric))
rownames(meta) <- rownames(data_averaged_sig2)
matrix <- data.frame(t(data_averaged_sig2[,-c(1:8)]))

design<- make.design.matrix(meta, degree=4)
NBp <- p.vector(matrix, design, counts = T,family=negative.binomial(theta = 10)) 
NBt <- T.fit2(NBp,  nvar.correction=T)
get<- get.siggenes(NBt, rsq = 0, vars="groups")
ID <- get$summary$PRvsCR

see.re <-  see.genes2(get$sig.genes$PRvsCR, 
                      edesign = get$sig.genes$PRvsCR$edesign,
          k=3,
          show.fit = T,
          dis = design$dis, 
          groups.vector = design$groups.vector,
          cluster.method="hclust",
          # distance = "euclidean",
           #distance = "euclidean",
          agglo.method="ward.D2",
          # summary.mode="representative",
          cluster.data = 1,
           #show.lines = T,
          newX11 = F,
          legend = F,
          ylim=c(-0.3,1.5))

cl <- data.frame(see.re$cut)
tt <- data.frame(get$sig.genes$PRvsCR$sig.pvalues)
tt2 <- data.frame(get$sig.genes$PRvsCR$coefficients)

# heatmap
cl <- data.frame(see.re$cut)
cl.1 <-  data.frame(t(log10(data_averaged[,rownames(cl)[which(cl$see.re.cut=="1")]]+1)))
cl.2 <-  data.frame(t(log10(data_averaged[,rownames(cl)[which(cl$see.re.cut=="2")]]+1)))
cl.3 <-  data.frame(t(log10(data_averaged[,rownames(cl)[which(cl$see.re.cut=="3")]]+1)))
# order
temp1 <- apply(cl.1,1,mean)
cl.12 <- cl.1[order(-temp1),]
temp2<- apply(cl.12,2,mean)
cl.13 <- cl.12[,order(-temp2)]
cl.1_meta <- data_averaged[colnames(cl.13),c("Stage2","Outcome1")]
cl.1_meta <- cl.1_meta[order(cl.1_meta$Outcome1,cl.1_meta$Stage2),]
cl.14 <- cl.13[,rownames(cl.1_meta)]

temp1 <- apply(cl.2,1,mean)
cl.22 <- cl.2[order(-temp1),]
cl.22 <- cl.22[,colnames(cl.14)]

temp1 <- apply(cl.3,1,mean)
cl.32 <- cl.3[order(-temp1),]
cl.32 <- cl.32[,colnames(cl.14)]
# plot
library("pheatmap")
library("RColorBrewer")
data.p <- rbind(cl.14, cl.22, cl.32)
 pheatmap(data.p, 
         cluster_rows = F, 
         cluster_cols = F,
         # border_color = "white",
         border=F,
         color = colorRampPalette(c(brewer.pal(8,"PuBu")[6], brewer.pal(8,"Set3")[2],brewer.pal(8,"Set1")[1]))(50),
         gaps_row = c(67,119,125), 
         gaps_col = c(99), 
         show_rownames = F,
         show_colnames = F,
         annotation_col = cl.1_meta,
         annotation_colors = list(Outcome1=c(CR=pal_lancet()(7)[1], PR=pal_lancet()(7)[7]), 
                                  Stage2=c(FCa=colorRampPalette(brewer.pal(8, "YlGn"))(8)[4],
                                           FCb=colorRampPalette(brewer.pal(8, "YlGn"))(8)[5],
                                           CRSa=colorRampPalette(brewer.pal(8, "YlGn"))(8)[6],
                                           CRSb=colorRampPalette(brewer.pal(8, "YlGn"))(8)[7],
                                           CRSc=colorRampPalette(brewer.pal(8, "YlGn"))(8)[8])))




############################################
#OTU: ID transform，composition

IDs <- read.table("gg13_OTU_annotation.tsv", sep="\t", header = T, stringsAsFactors = F)
rownames(IDs) <- IDs$ID0

cl$ID <- IDs[rownames(cl),"taxonomy"]
cl$ID1 <- unlist(lapply(cl$ID,function(x){strsplit(x,split = "; ")[[1]][2]}))
cl$ID2 <- unlist(lapply(cl$ID,function(x){strsplit(x,split = "; ")[[1]][3]}))
cl$ID3 <- unlist(lapply(cl$ID,function(x){strsplit(x,split = "; ")[[1]][4]}))
cl$ID4 <- unlist(lapply(cl$ID,function(x){strsplit(x,split = "; ")[[1]][5]}))
cl$ID5 <- unlist(lapply(cl$ID,function(x){strsplit(x,split = "; ")[[1]][6]}))

##### pie
dd2 <- melt(dcast(data=cl, ID1+ID3 ~ see.re.cut, value.var= "ID3",length, margins = F), ids="ID3")
dd2 <- dd2[order(dd2$variable ,dd2$ID1, decreasing = T),]
dd2$ID1 <- factor(dd2$ID1, levels = unique(dd2$ID1))
dd2$ID3 <- factor(dd2$ID3, levels = unique(dd2$ID3))
dd2$variable <- as.numeric(dd2$variable)
dd3 =dcast(dd2, ID1~variable, sum)
colnames(dd3) <- c("ID1","v1","v2","v3")

t0 <- colorRampPalette(pal_rickandmorty()(9))(19)
names(t0) <- as.character(dd2[which(dd2$variable==1),"ID3"])
t1 <- colorRampPalette(pal_lancet()(8))(9)
names(t1) <- c("p__Firmicutes","p__Proteobacteria","p__Bacteroidetes","p__Actinobacteria","p__Fusobacteria","p__Euryarchaeota", "p__Cyanobacteria", "p__TM7","p__[Thermi]")

  
p1 <- ggplot(data=dd2[which(dd2$variable==1),])+
  geom_bar(data=dd3 ,aes(x = 1, fill = ID1, y = v1), stat = "identity",width=1.3, color="white")+
  geom_bar(aes(x = 2, fill = ID3, y = value), stat = "identity", width=0.5, color="white", alpha=0.7)+
  # geom_text(data=dd3[rev(which(dd3$v3!=0)),], aes(x=1, y = v1, label=unlist(lapply(ID1, function(x){gsub("p__","",x)}))), position = position_stack(vjust=0.5),size=3)+
  # geom_text(data=dd2[rev(which(dd2$variable==1 & dd2$value!=0)),], aes(x=2, y = value, label=unlist(lapply(ID3, function(x){gsub("c__","",x)}))),position = position_stack(reverse = TRUE,vjust=0.5),size=5, vjust=1,hjust=0)+
  coord_polar(theta = "y")+
  scale_fill_manual(values = c(t0,t1))+
  theme_void() + 
  labs(x = "", y = "", title = "") + 
  theme(legend.position = "none",
        panel.grid = element_blank()) 


p2 <- ggplot(data=dd2[which(dd2$variable==2),])+
  geom_bar(data=dd3 ,aes(x = 1, fill = ID1, y = v2), stat = "identity",width=1.3, color="white")+
  geom_bar(aes(x = 2, fill = ID3, y = value), stat = "identity", width=0.5, color="white", alpha=0.7)+
  # geom_text(data=dd3[rev(which(dd3$v3!=0)),], aes(x=1, y = v2, label=unlist(lapply(ID1, function(x){gsub("p__","",x)}))), position = position_stack(vjust=0.5),size=3)+
  # geom_text(data=dd2[rev(which(dd2$variable==2 & dd2$value!=0)),], aes(x=2, y = value, label=unlist(lapply(ID3, function(x){gsub("c__","",x)}))),position = position_stack(reverse = TRUE,vjust=0.5),size=5, vjust=1,hjust=0)+
  coord_polar(theta = "y")+
  scale_fill_manual(values = c(t0,t1))+
  theme_void() + 
  labs(x = "", y = "", title = "") + 
  theme(legend.position = "none",
        panel.grid = element_blank()) 

p3 <- ggplot(data=dd2[which(dd2$variable==3),])+
  geom_bar(data=dd3 ,aes(x = 1, fill = ID1, y = v3), stat = "identity",width=1.3, color="white")+
  geom_bar(aes(x = 2, fill = ID3, y = value), stat = "identity", width=0.5,color="white", alpha=0.7)+
  # geom_text(data=dd3[rev(which(dd3$v3!=0)),], aes(x=1, y = v3, label=unlist(lapply(ID1, function(x){gsub("p__","",x)}))), position = position_stack(vjust=0.5),size=3)+
  # geom_text(data=dd2[rev(which(dd2$variable==3 & dd2$value!=0)),], aes(x=2, y = value, label=unlist(lapply(ID3, function(x){gsub("c__","",x)}))),position = position_stack(reverse = TRUE,vjust=0.5),size=5, vjust=1,hjust=0)+
  scale_fill_manual(values = c(t0,t1))+
  theme_void() + 
  labs(x = "", y = "", title = "") + 
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  coord_polar(theta = "y")

p <- plot_grid(p1,p2,p3, ncol=1, align="h", axis="l", labels = c("Cluster1","Cluster2","Cluster3"))

# legend
legend <- get_legend(p1 + scale_shape(guide = FALSE) + theme(legend.position = "right"))
p <- plot_grid(p, legend, ncol=2)
p




#############################################
## fig 4

# maSigPro
data <- read.table("gg_level1-7-relabun_averaged_filter.freq.2000-filter.unclassified2.txt", header = T, row.names=1, stringsAsFactors = F)
data <- data[,c(1:8,grep("g__",colnames(data)))]
data <- data[,c(grep("s__",colnames(data), invert = T))]
# extract MM
data_MM <- data[which(data$disease1=="MM"),]
data_MM1 <- data_MM[which(data_MM$Outcome1 %in% c("PR","CR")),]
#remove 90% samples = 0
temp <- apply(data_MM1[,c(9:ncol(data_MM1))], 2 , function (x) {sum(x==0)})
temp2 <- temp[which(temp > nrow(data_MM1)*0.9)]
data_MM11 <- data_MM1[,-which(colnames(data_MM1) %in% names(temp2))]
data_averaged <- data_MM11
data_averaged$Stage2 <- factor(data_averaged$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc"))
data_averaged <- data_averaged[order(data_averaged$Stage2),]

## maSigPro PR vs CR
data_averaged_sig2 <- data_averaged
meta <- subset(data_averaged_sig2, select = "Stage2")
meta$Stage2 <- as.numeric(factor(meta$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc")))
meta$Replicate <- meta$Stage2
meta$group <- data_averaged_sig2$Outcome1
meta$CR <- data_averaged_sig2$Outcome1
meta$PR <- data_averaged_sig2$Outcome1
meta$PR[which(meta$PR=="PR")] <- 1
meta$PR[which(meta$PR=="CR")] <- 0
meta$CR[which(meta$CR=="PR")] <- 0
meta$CR[which(meta$CR=="CR")] <- 1
meta$group <- NULL
meta$Replicate <- paste(meta$Stage2,meta$CR, sep="")
meta$Replicate <- paste(meta$Replicate,meta$PR, sep="")
meta$Replicate <- as.numeric(as.factor(meta$Replicate))
meta <- as.matrix(apply(meta,2,as.numeric))
rownames(meta) <- rownames(data_averaged_sig2)
matrix <- data.frame(t(data_averaged_sig2[,-c(1:8)]))

design<- make.design.matrix(meta, degree=4)
NBp <- p.vector(matrix, design, counts = T,family=negative.binomial(theta = 10)) 
NBt <- T.fit(NBp,nvar.correction=T)
get<- get.siggenes(NBt, rsq = 0, vars="groups")
ID <- get$summary$PRvsCR

re.pvalues <- data.frame(get$sig.genes$PRvsCR$sig.pvalues)
re.coefficients <- data.frame(get$sig.genes$PRvsCR$coefficients)



### glmm
ma2 <- re.pvalues[-grep(".s_",rownames(re.pvalues)),]
ma2$fdr <- p.adjust(ma2$p.value, method = "fdr", n = 96)

source("glme.taxon.functions.R")
data1 <-  sub1(data, "MM", "g")
data11 <- cbind(data1[,c(1:8)],data1[,rownames(ma2)])
data2 <- subset(data11, Outcome1 %in% c("CR","PR")) 

#  generalized linear-mixed models, GLMM（before， after）
data.Bf <- subset(data2, Outcome1 %in% c("CR","PR")) %>% subset(.,Stage2 %in% c("FCa","FCb"))
genus.lmer.Bf <- taxa.lmer(data.Bf)
genus.lmer.Bf2 <- lapply(genus.lmer.Bf, function(x){x$coefficients[,4]}) %>% unlist %>% matrix(.,ncol=2,byrow = T) %>% as.data.frame()
rownames(genus.lmer.Bf2) <- names(genus.lmer.Bf)
colnames(genus.lmer.Bf2) <- rownames(genus.lmer.Bf[[1]]$coefficients)
genus.lmer.Bf2$fdr <- p.adjust(genus.lmer.Bf2$Outcome1PR, method = "fdr")
genus.lmer.Bf2$coef <- lapply(genus.lmer.Bf, function(x){x$coefficients[2,1]}) %>% unlist 

data.Aft <- subset(data2, Outcome1 %in% c("CR","PR")) %>% subset(.,Stage2 %in% c("CRSa","CRSb","CRSc"))
genus.lmer.Aft <- taxa.lmer(data.Aft)
genus.lmer.Aft2 <- lapply(genus.lmer.Aft, function(x){x$coefficients[,4]}) %>% unlist %>% matrix(.,ncol=2,byrow = T) %>% as.data.frame()
rownames(genus.lmer.Aft2) <- names(genus.lmer.Aft)
colnames(genus.lmer.Aft2) <- rownames(genus.lmer.Aft[[1]]$coefficients)
genus.lmer.Aft2$fdr <- p.adjust(genus.lmer.Aft2$Outcome1PR, method = "fdr")
genus.lmer.Aft2$coef <- lapply(genus.lmer.Aft, function(x){x$coefficients[2,1]}) %>% unlist 

## merge tables
mt <- merge(ma2[,c("p.value","fdr")], genus.lmer.Bf2[,-1], by="row.names", all = T )
rownames(mt) <- mt$Row.names
mt <- mt[,-1]
mt2 <- merge(mt, genus.lmer.Aft2[,-1], by="row.names", all = T )
rownames(mt2) <- exid(mt2$Row.names,"g")
mt2 <- mt2[,-1]
mt3 <- merge(mt2, lefse2, by="row.names", all = T )
colnames(mt3) <- c("ID","MaSigPro p", "MaSigPro fdr","Before CART p","Before CART fdr","Before CART coefficient","After CART p","After CART fdr","After CART coefficient","LDA Before CART","LDA Before CART p","LDA After CART","LDA After CART p")



## random forest
load("data_RF.RData")
FCa <- data.frame(t(FCa.0), stringsAsFactors = F)
FCb <- data.frame(t(FCb.0), stringsAsFactors = F)
data <- merge(FCa,FCb,by="subject")
data$Outcome1.x <- data$Outcome1.y
data$Outcome1.x[9] <- "PR"
data$Outcome1.y <- NULL
data$subject <- NULL

colnames(data)[1] <- c("DiseaseStat2")
g <- colnames(data)[-1][grep("g__",colnames(data)[-1])][grep("s__",colnames(data)[-1][grep("g__",colnames(data)[-1])],invert = T)]
data <- data[,c("DiseaseStat2",g)]
data$DiseaseStat2 <- as.factor(data$DiseaseStat2)
set.seed(12345)

# setseed
set.seed(123456)
subsetSizes <- c(seq(1,ncol(data)-1,by=1))
seeds <- vector(mode = "list", length = 50)
for(i in 1:50) seeds[[i]] <- sample.int(10000, length(subsetSizes) + 1)
seeds[[51]] <- sample.int(100, 1)

#5-CV
ctrl.cv= rfeControl(functions = rfFuncs, method = "repeatedcv", verbose = FALSE, number=5,repeats = 5, seeds = seeds)
# feature selection
fs_nb13 <- rfe(x = data[,-1],
             y = data[,1],
             sizes = c(seq(1,ncol(data)-1,by=1)),
             rfeControl = ctrl.cv,
             rerank=T)

# res
cv.result <- fs_nb13$results
ggplot(cv.result, aes(x=Variables, y=Accuracy), color="blue")+
  geom_line(col="grey", linetype=1, size=1)+
  geom_point(size=2, shape=21, fill="grey", col="black" ) + 
  xlab("Number of variables") +
  ylab("Mean CV Error") +
  ggtitle("5-fold rfcv,10times") +
  theme_bw() +
  theme(panel.grid=element_blank(), panel.border = element_rect(size=0.8),axis.text = element_text(size=12, color="black"))        

# important variables
imp.var <- fs_nb13$optVariables
imp <- fs_nb12$fit$importance

# re-construct
control <- trainControl(method="loocv")
data_se <- data[,c("DiseaseStat2",imp.var[1:3])]
model <- train(x = data_se[,-1], y = data_se[,1], method="rf",trControl=control)

# loocv
data_se <- data[,c("DiseaseStat2",imp.var[1:8])]
data2 <- data_se
predictions <- 1:22
for (k in 1:22){
  set.seed(1234) 
  predictions[k] <- predict(randomForest(DiseaseStat2 ~. , 
                                         data = data2[-k,], ntree=500),type="prob", 
                            newdata = data2[k,])[2] 
}
auc(data2$DiseaseStat2,predictions,direction="<", 
    levels = levels(data2$DiseaseStat2))
pred=prediction(predictions,data2$DiseaseStat2)
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]
perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve

##  AUC 
rocobj <- roc(data2$DiseaseStat2,predictions, ci=T)
cis <- ci(rocobj, of="auc")
##  ROC 
ggplot(data=data.frame("Specificity"=perf_ROC@x.values[[1]], "Sensitivity"=perf_ROC@y.values[[1]]), aes(x=Specificity, y=Sensitivity))+
  geom_line(color=pal_lancet("lanonc")(9)[1],size=0.8)+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.border = element_rect(size=0.8),axis.text = element_text(size=12, color="black"), axis.title = element_text(size=12, color="black"))+
  geom_text(x=0.8,y=0.15,col=pal_lancet("lanonc")(9)[1],label=paste("AUC:",format(AUC, digits=2, scientific=FALSE)))+
  geom_text(x=0.8,y=0.05,col="darkgray",label=paste("[",format(cis[1], digits=2, scientific=FALSE),"-",format(cis[3], digits=2, scientific=FALSE),"]"))+
  ggtitle("ROC for FCa top importance (LOOCV)")


## survival
load("survival.RData")
library("survival")
library("survminer")
g <- "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Alcaligenaceae.g__Sutterella"
temp <- data_MM1[,c("subject","Outcome1","Stage2",g)]
colnames(temp)[4] <- "feature"
temp$Stage2 <- factor(temp$Stage2, levels = c("FCa","FCb","CRSa","CRSb","CRSc"))
temp$PFS <- meta3[rownames(temp),"PFS"]
temp$PFS_status <- meta3[rownames(temp),"PFS_status"]

# Mean
tmp2 <- melt(temp[,-3], id=c("subject","Outcome1","PFS","PFS_status"))
temp2 <- dcast(tmp2,subject+Outcome1+PFS+PFS_status~variable, mean)  
temp2$f_group <- temp2$feature
temp2$f_group[which(temp2$feature<quantile(temp2$feature,0.33))] <- "Low"
temp2$f_group[which(temp2$feature>=quantile(temp2$feature,0.33) & temp2$feature<quantile(temp2$feature,0.66))] <- "Median"
temp2$f_group[which(temp2$feature>=quantile(temp2$feature,0.66))] <- "High"
fit <- survfit(Surv(PFS, PFS_status) ~ f_group, data = temp2)
ggsurvplot(fit, temp2,
           pval = TRUE,
           pval.method=T,
           palette  = c(pal_lancet("lanonc")(9)[c(1:3)]),
           legend="right",
           title="Mean Abundance")




#############################################
##  Fig 5

# maSigPro
data_MM1 <- data_MM
#remove 90% samples = 0
temp <- apply(data_MM1[,c(9:ncol(data_MM1))], 2 , function (x) {sum(x==0)})
temp2 <- temp[which(temp > nrow(data_MM1)*0.9)]
data_MM11 <- data_MM1[,-which(colnames(data_MM1) %in% names(temp2))]
data_averaged <- data_MM11
data_averaged$Stage2 <- factor(data_averaged$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc"))
data_averaged <- data_averaged[order(data_averaged$Stage2),]

# maSigPro，CRSgrade 
data_averaged_sig2 <- data_averaged
meta <- subset(data_averaged_sig2, select = "Stage2")
meta$Stage2 <- as.numeric(factor(meta$Stage2, levels = c("FCa","FCb", "CRSa","CRSb","CRSc")))
meta$Replicate <- meta$Stage2
meta$group <- data_averaged_sig2$CRSgrade
meta$CRS1 <- data_averaged_sig2$CRSgrade
meta$CRS2 <- data_averaged_sig2$CRSgrade
meta$CRS3 <- data_averaged_sig2$CRSgrade
meta$CRS1[which(meta$group==1)] <- 1
meta$CRS1[which(meta$group!=1)] <- 0
meta$CRS2[which(meta$group==2)] <- 1
meta$CRS2[which(meta$group!=2)] <- 0
meta$CRS3[which(meta$group==3)] <- 1
meta$CRS3[which(meta$group!=3)] <- 0
meta$group <- NULL
meta$Replicate <- paste(meta$Stage2,meta$CRS1, sep="")
meta$Replicate <- paste(meta$Replicate,meta$CRS2, sep="")
meta$Replicate <- paste(meta$Replicate,meta$CRS3, sep="")
meta$Replicate <- as.numeric(as.factor(meta$Replicate))
meta <- as.matrix(apply(meta,2,as.numeric))
rownames(meta) <- rownames(data_averaged_sig2)
matrix <- data.frame(t(data_averaged_sig2[,-c(1:8)]))

design<- make.design.matrix(meta, degree=4)
NBp <- p.vector(matrix, design, counts = T,family=negative.binomial(theta = 10)) 
NBt <- T.fit2(NBp,nvar.correction=T)
get<- get.siggenes(NBt, rsq = 0, vars="groups")
ID <- get$summary$CRS3vsCRS1

### res
re.p <- get$sig.genes$CRS3vsCRS1$sig.pvalues
re.coef <- get$sig.genes$CRS3vsCRS1$group.coeffs



# glmm
ma2 <- re.pvalues[-grep(".s_",rownames(re.pvalues)),]
ma2$fdr <- p.adjust(ma2$p.value, method = "fdr", n = 96)
data2 <- subset(data11, CRSgrade %in% c("1","2","3")) 
# Poisson generalized linear-mixed models, GLMM（before， after）
data.Bf <- subset(data2,CRSgrade %in% c("1","2","3") ) %>% subset(.,Stage2 %in% c("FCa","FCb","CRSa"))
genus.lmer.Bf <- taxa.lmer2(data.Bf)
genus.lmer.Bf2 <- lapply(genus.lmer.Bf, function(x){x$coefficients[,4]}) %>% unlist %>% matrix(.,ncol=2,byrow = T) %>% as.data.frame()
rownames(genus.lmer.Bf2) <- names(genus.lmer.Bf)
colnames(genus.lmer.Bf2) <- rownames(genus.lmer.Bf[[1]]$coefficients)
genus.lmer.Bf2$fdr <- p.adjust(genus.lmer.Bf2$CRSgrade, method = "fdr")
genus.lmer.Bf2$coef <- lapply(genus.lmer.Bf, function(x){x$coefficients[2,1]}) %>% unlist 
  
data.Aft <- subset(data2, CRSgrade %in% c("1","2","3") ) %>% subset(.,Stage2 %in% c("CRSb"))
genus.lmer.Aft <- taxa.lmer3(data.Aft)
genus.lmer.Aft2 <- lapply(genus.lmer.Aft, function(x){x$coefficients[,4]}) %>% unlist %>% matrix(.,ncol=2,byrow = T) %>% as.data.frame()
rownames(genus.lmer.Aft2) <- names(genus.lmer.Aft)
colnames(genus.lmer.Aft2) <- rownames(genus.lmer.Aft[[1]]$coefficients)
genus.lmer.Aft2$fdr <- p.adjust(genus.lmer.Aft2$CRSgrade, method = "fdr")
genus.lmer.Aft2$coef <- lapply(genus.lmer.Aft, function(x){x$coefficients[2,1]}) %>% unlist 

## merge tables
mt <- merge(ma2[,c("p.value","fdr")], genus.lmer.Bf2[,-1], by="row.names", all = T )
rownames(mt) <- mt$Row.names
mt <- mt[,-1]
mt2 <- merge(mt, genus.lmer.Aft2[,-1], by="row.names", all = T )
rownames(mt2) <- exid(mt2$Row.names,"g")
mt2 <- mt2[,-1]
mt3 <- merge(mt2, lefse2, by="row.names", all = T )
colnames(mt3) <- c("ID","MaSigPro p", "MaSigPro fdr","Before CRS p","Before CRS fdr","Before CRS coefficient","After CRS p","After CRS fdr","After CRS coefficient","LDA Before CRS","LDA Before CRS p","LDA After CRS","LDA After CRS p")


# correlation net, rmcorr
source("rmcorr_functions.R")
load("cyto.RData")
taxa2 <- data_MM1
taxa2[,c(3:8)] <- NULL
taxa2$ID <- paste(taxa2$subject, taxa2$Stage2, sep=".")
all <- merge(data[,-c(2,3)], taxa2[,-c(1,2)], by="ID", all=T)
data3 <- all[,c(1:ncol(data[,-c(2,3)]))]
taxa3 <- all[,c(1,(ncol(data[,-c(2,3)])+1):ncol(all))]
id <- colnames(taxa3) %>% grep("g__",., value = T) %>% grep("s__",.,invert = T, value = T)
taxa3 <- taxa3[,c("ID",id)]
colnames(taxa3) <- lapply(colnames(taxa3), function(x){strsplit(x,"g__")[[1]][2]}) %>% unlist() 
colnames(taxa3)[1] <- "ID"
id3 <- c("Sutterella","Collinsella","Paraprevotella","Bifidobacterium","Weissella","Anaerotruncus","Prevotella","Oscillospira","Megasphaera","Megamonas","Methanobrevibacter","Faecalibacterium","Enterococcus","Gemmiger","Bulleidia","Clostridium","Leuconostoc","Pseudoramibacter_Eubacterium","Odoribacter","Roseburia","Dialister","Staphylococcus","Enhydrobacter","Ruminococcus","Lactobacillus","Dorea","Atopobium","Bacillus","Meiothermus","Sphingomonas","Bradyrhizobium","Leuconostoc","Bifidobacterium","Anaerotruncus","Lactococcus","Morganella","Butyricicoccus","Lachnoanaerobaculum","Fusobacterium","Staphylococcus","Veillonella","Sphingomonas","Akkermansia","Peptostreptococcus","Holdemania","Oscillospira","Alistipes","Finegoldia","Enhydrobacter","Novosphingobium",".Prevotella.","Enterococcus","Megamonas","Thermus","Curvibacter","Parabacteroides","Methanobrevibacter","Butyricimonas") %>% unique()
taxa3 <- taxa3[,c("ID",id3)] 
subject <- unlist(lapply(data3$ID,function(x){strsplit(x,"\\.")[[1]][1]}))

# rmcorr
library(rmcorr)
data4 <- apply(data3[,-1],2,as.numeric) %>% as.data.frame()
taxa4 <- apply(taxa3[,-1],2,as.numeric) %>% as.data.frame()
res_rmcorr <- sapply(log10(data4+1),
                     function(x) Map(function(a,b) my_rmcorr(a,b,subject),
                                     log10(taxa4+1), as.data.frame(x)))

res_rmcorr2 <- unlist(res_rmcorr) %>% matrix(.,ncol = 2, byrow = T) %>% as.data.frame()
res_rmcorr2$cytokine <- rep(colnames(data3)[-1], each=(ncol(taxa3)-1))
res_rmcorr2$taxa <- rep(colnames(taxa3)[-1], n=(ncol(data3)-1))
colnames(res_rmcorr2)[c(1,2)] <- c("r","p")
res_rmcorr2$FDR <- p.adjust(res_rmcorr2$p, method = "fdr")

res_rmcorr3 <- res_rmcorr2[which(res_rmcorr2$FDR<0.05),]

pp <- cbind(subset(data4, select = "CCL2.MCP.1"),subset(taxa4, select = "Lactobacillus"))
pp <- log10(pp+1)
pp$subject <- subject
ggplot(data=pp, aes(x=CCL2.MCP.1, y=Lactobacillus))+
  geom_point(aes(color=subject))+
  geom_smooth(aes(color=subject), method = "lm", se=F, lwd=0.7, lty=18)+
  geom_smooth( method = "lm", lwd=1, color="darkred", alpha=0.35)+
  theme_classic()+
  my_theme2+
  theme(legend.position = "none")


pp <- cbind(subset(data4, select = "leukocyte"),subset(taxa4, select = "Veillonella"))
pp <- log10(pp+1)
pp$subject <- subject
 ggplot(data=pp, aes(x=	leukocyte	, y=Veillonella))+
  geom_point(aes(color=subject))+
  geom_smooth(aes(color=subject), method = "lm", se=F, lwd=0.7, lty=18)+
  geom_smooth( method = "lm", lwd=1, color="darkred", alpha=0.35)+
  theme_classic()+
  my_theme2+
  theme(legend.position = "none")



















