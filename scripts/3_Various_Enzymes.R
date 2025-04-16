# Nano-FIRSTseq 
# R Analysis
# Oguzhan Begik
# Part 4. Various Enzymes

library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggcorrplot)



#Color scheme
protoscript_col <- "#aa6f73"
ss2_col <- "#eea990"
ss3_col <- "#f6e0b5"
ss4_col <- "#D193C1"
maxima_col <- "#9C8CC3"
tgirt_col <- "#5294a3"
induro_col <- "#a3d9d9"
marathon_col <- "#60bfae"
#Palette vector
enzyme_palette <- c(protoscript_col,ss2_col,ss3_col,ss4_col,maxima_col,tgirt_col,induro_col,marathon_col)



# 1.COUNTS ANALYSIS
process_count <- function(input, enzyme,buffer) {
	data <- read.delim(input, header=FALSE)
	data2 <- data[,c("V1", "V3")]
	colnames(data2) <- c("Chr", "Count")
	data2$Count <- ((data2$Count)/sum(data2$Count))*1000000
	data2$Enzyme <-  enzyme
	data2$Buffer <-  buffer
	return(data2)
}

#Mg
protoscript_mg.cnt <- process_count("data/chrcount/FIRST_Mg_PS2.chrcount.tsv", "PS2","Mg")
ss2_mg.cnt <- process_count("data/chrcount/FIRST_Mg_SS2.chrcount.tsv", "SS2","Mg")
ss3_mg.cnt <- process_count("data/chrcount/FIRST_Mg_SS3.chrcount.tsv", "SS3","Mg")
ss4_mg.cnt <- process_count("data/chrcount/FIRST_Mg_SS4.chrcount.tsv", "SS4","Mg")
maxima_mg.cnt <- process_count("data/chrcount/FIRST_Mg_Maxima.chrcount.tsv", "Maxima","Mg")
tgirt_mg.cnt <- process_count("data/chrcount/FIRST_Mg_TGIRT.chrcount.tsv", "TGIRT","Mg")
induro_mg.cnt <- process_count("data/chrcount/FIRST_Mg_Induro.chrcount.tsv", "Induro","Mg")
marathon_mg.cnt <- process_count("data/chrcount/FIRST_Mg_Marathon.chrcount.tsv", "Marathon","Mg")

#Mn
protoscript_mn.cnt <- process_count("data/chrcount/FIRST_Mn_PS2.chrcount.tsv", "PS2","Mn")
ss2_mn.cnt <- process_count("data/chrcount/FIRST_Mn_SS2.chrcount.tsv", "SS2","Mn")
ss3_mn.cnt <- process_count("data/chrcount/FIRST_Mn_SS3.chrcount.tsv", "SS3","Mn")
ss4_mn.cnt <- process_count("data/chrcount/FIRST_Mn_SS4.chrcount.tsv", "SS4","Mn")
maxima_mn.cnt <- process_count("data/chrcount/FIRST_Mn_Maxima.chrcount.tsv", "Maxima","Mn")
tgirt_mn.cnt <- process_count("data/chrcount/FIRST_Mn_TGIRT.chrcount.tsv", "TGIRT","Mn")
induro_mn.cnt <- process_count("data/chrcount/FIRST_Mn_Induro.chrcount.tsv", "Induro","Mn")
marathon_mn.cnt <- process_count("data/chrcount/FIRST_Mn_Marathon.chrcount.tsv", "Marathon","Mn")



enzymes.cnt <- rbind(protoscript_mg.cnt,ss2_mg.cnt,ss3_mg.cnt,ss4_mg.cnt,maxima_mg.cnt,tgirt_mg.cnt,induro_mg.cnt,marathon_mg.cnt,protoscript_mn.cnt,ss2_mn.cnt,ss3_mn.cnt,ss4_mn.cnt,maxima_mn.cnt,tgirt_mn.cnt,induro_mn.cnt,marathon_mn.cnt) 
enzymes.cnt2 <-  dcast(enzymes.cnt, Chr ~ Enzyme+Buffer, value.var = "Count")
enzymes.cntfilt <- enzymes.cnt2[apply(enzymes.cnt2[,-1], 1, function(x) !all(x==0)),]
enzymes.cntfilt$Biotype <- gsub(".*_", "", enzymes.cntfilt$Chr )
enzymes.cntfilt <- subset(enzymes.cntfilt, Biotype!="oligo3" & Biotype!="oligo5")
enzymes.cntfilt <-  enzymes.cntfilt %>% select("Biotype", everything())

columns_mg <-c("PS2_Mg", "SS2_Mg","SS3_Mg", "SS4_Mg","Maxima_Mg", "TGIRT_Mg","Induro_Mg", "Marathon_Mg")
columns_mn <-c("PS2_Mn", "SS2_Mn","SS3_Mn", "SS4_Mn","Maxima_Mn", "TGIRT_Mn","Induro_Mn", "Marathon_Mn")


#Now PCA
pca_plot <- function(data, label, biotype, levels) {
	data2 <- data[,c("Chr","Biotype", levels)]
	data3 <- subset(data2, Biotype %in% biotype)
	rownames(data3) <- data3$Chr
	data4 <- data3[,c(3:10)]
	data5 <- t(data4)
	data5 <- data5[, colSums(data5 != 0) > 0]
	pca_data <- prcomp(data5, center = TRUE, scale=TRUE)
	pca_loadings <- data.frame(RNAs = rownames(pca_data$rotation), pca_data$rotation)
	summ <- summary(pca_data)
	variance <- summ$importance[2,]
	df_out <- as.data.frame(pca_data$x)
	df_out$Enzymes <- rownames(df_out)
	df_out$Enzymes <- factor(df_out$Enzymes, levels = levels)
	pdf(file= paste(label, "Enzymes_ChrCounts_PCA.pdf", sep="_"),height=4,width=5.5,onefile=FALSE)
    print(ggplot(df_out,aes(x=PC1,y=PC2, color=Enzymes))+
    	#geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = (PC1*50),yend = (PC2*50)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
    	#annotate("text", x = (pca_loadings$PC1*50000), y = (pca_loadings$PC2*50000),label = pca_loadings$RNAs)+
        geom_point(size=3)+
        theme_bw()+
        ggtitle(label)+
        scale_color_manual(values=enzyme_palette)+
        xlab(paste("PC1", as.character(round(variance[1]*100)),"%"))+
        ylab(paste("PC2", as.character(round(variance[2]*100)),"%"))+
        geom_text_repel(data=df_out, aes(label=rownames(df_out),x=PC1, y=PC2), colour="black",segment.size  = 0.4,segment.color = "grey50",size=2))
    dev.off()
}

pca_plot(enzymes.cntfilt, "tRNA_counts_Mg", "tRNA",columns_mg)
pca_plot(enzymes.cntfilt, "tRNA_counts_Mn", "tRNA",columns_mn)

pca_plot(enzymes.cntfilt, "norRNA_counts_Mg",  c("tRNA","IVT", "ncRNA", "snRNA","snoRNA", "Sequins"),columns_mg)
pca_plot(enzymes.cntfilt, "norRNA_counts_Mn",  c("tRNA","IVT", "ncRNA", "snRNA","snoRNA", "Sequins"),columns_mn)

pca_plot(enzymes.cntfilt, "All_counts_Mg", unique(enzymes.cntfilt$Biotype),columns_mg)
pca_plot(enzymes.cntfilt, "All_counts_Mn", unique(enzymes.cntfilt$Biotype),columns_mn)


bar_plot <- function(data, label, levels) {
	data <- data[,c("Chr","Biotype", levels)]
	data2 <- melt(data)
	data3 <- aggregate(data2$value, by=list(variable=data2$variable,Biotype=data2$Biotype), FUN=sum)
	data3$Biotype <- factor(data3$Biotype, levels = c("ncRNA","snRNA", "snoRNA","mtrRNA","tRNA", "Sequins", "IVT", "rRNA"))
	colnames(data3) <- c("Enzymes", "Biotype", "Count")
	pdf(file= paste(label,"BiotypeCount_Barplot.pdf", sep="_"),height=5,width=10,onefile=FALSE)
		print(ggplot(data3,aes(x = Enzymes, y = Count, fill = Biotype)) +
			geom_col(colour = "black",position = "fill")+
			#scale_y_continuous(labels = scales::percent) +
			scale_fill_brewer(palette = "Pastel1")+
			#coord_cartesian(ylim = c(0.9,1))+
			#scale_y_break(c(0, 0.8)) + 
			theme_classic())
			dev.off()
}

bar_plot(enzymes.cntfilt, "Mg",columns_mg)
bar_plot(enzymes.cntfilt, "Mn",columns_mn)



### Correlation Matrix of the enzymes
#library(viridis)

enzymes.cntfilt.norRNA <- subset(enzymes.cntfilt, Biotype!="rRNA")
enzymes.cntfilt.norRNA.Mg<- enzymes.cntfilt.norRNA[,columns_mg]
enzymes.cntfilt.norRNA.Mg.log <- log(enzymes.cntfilt.norRNA.Mg+1)
corr.enzymes.cntfilt.norRNA.Mg.scaled <- round(cor(enzymes.cntfilt.norRNA.Mg.log),2)

soft_blue_gradient2 <- c("#EAE4E0",  "#51618D")
pdf(file="Correlation_Matrix_Enzymes_Mg.pdf",height=5,width=5,onefile=FALSE)
	print(ggcorrplot(corr.enzymes.cntfilt.norRNA.Mg.scaled, hc.order = TRUE,lab = TRUE,lab_size = 3,type = "lower",outline.col = "white")+
	scale_fill_gradientn(colors = soft_blue_gradient2, na.value = "#EAE4E0",breaks=c(0.9,0.95, 1), limit=c(0.9, 1)))
	dev.off()

soft_blue_gradient2 <- c("#EAE4E0",  "#51618D")
pdf(file="Correlation_Matrix_Enzymes_Mg_var2.pdf",height=5,width=5,onefile=FALSE)
	print(ggcorrplot(corr.enzymes.cntfilt.norRNA.Mg.scaled, hc.order = TRUE,outline.col = "white")+
	scale_fill_gradientn(colors = soft_blue_gradient2, na.value = "#EAE4E0",breaks=c(0.9,0.95, 1), limit=c(0.9, 1)))
	dev.off()


enzymes.cntfilt.norRNA <- subset(enzymes.cntfilt, Biotype!="rRNA")
enzymes.cntfilt.norRNA.Mn<- enzymes.cntfilt.norRNA[,columns_mn]
enzymes.cntfilt.norRNA.Mn.log <- log(enzymes.cntfilt.norRNA.Mn+1)
corr.enzymes.cntfilt.norRNA.Mn.scaled <- round(cor(enzymes.cntfilt.norRNA.Mn.log),2)


soft_blue_gradient2 <- c("#EAE4E0",  "#51618D")
pdf(file="Correlation_Matrix_Enzymes_Mn.pdf",height=5,width=5,onefile=FALSE)
	print(ggcorrplot(corr.enzymes.cntfilt.norRNA.Mn.scaled, hc.order = TRUE,outline.col = "white")+
	scale_fill_gradientn(colors = soft_blue_gradient2, na.value = "#EAE4E0",breaks=c(0.9,0.95, 1), limit=c(0.9, 1)))
	dev.off()





# 2.Modifications
mod_positions <- read.delim("data/references/yeast_mod_positions.tsv")

# 1. Read data and clean up
########################################
process_func <- function(stats,readends, enzyme,buffer) {
	#Stats file input
	stats_file <- read.delim(stats, sep="")
	#Readends file input
	readends_file <- read.delim(readends, sep="")
	# Merge two files
	merged <- merge(stats_file,readends_file[c("Ref", "Pos", "Ends")], by.x=c("chr","pos"), by.y=c("Ref", "Pos"))
	#Take out unnecessary columns
	clean <- merged[,c("chr", "pos", "ref_nuc", "coverage","rtstop","Ends", "ins", "del", "A", "T","C","G")]
	#Calculate mismatch per position
	bases <- c("A", "T","C","G")
	mis_all <- vector()
	for (base in unique(clean$ref_nuc)){
		subs <- subset(clean, ref_nuc==base)
		mis_bases <-  bases[ !bases == base]
		subs$mismatch <- subs[, mis_bases[1]] + subs[, mis_bases[2]]+ subs[, mis_bases[3]]
		mis_all <- rbind(mis_all, subs)
	}
	#Remove first 20 nts 
	mis_all2 <- subset(mis_all, pos > 20)
	#Create another column with Biotype
	mis_all2$Biotype <- gsub(".*_","",mis_all2$chr)
	#Create another column with condition
	mis_all2$Enzyme <- enzyme
	#Create another column with Ion
	mis_all2$Buffer <- buffer	
	#Create a unique coordinate column
	mis_all2$uniq_coord <- paste(mis_all2$chr, mis_all2$pos, sep="_")
	#Normalize coverage by each chromosome
	merged<-vector() 
	for (ref in (unique(mis_all2$chr))){
		subs<- subset(mis_all2, chr==ref)
		subs_3end<- subset(subs, pos> nrow(subs)/2) #Normalize by 3end coverage
		subs$norm_cov<- round(subs$coverage/max(subs_3end$coverage), 3)
		merged<- rbind(merged, subs)
	}
	merged$rtstop <- as.numeric(merged$rtstop)
	merged$norm_rt<- round(merged$rtstop/merged$coverage, 3) #to normalize RTdrop by the coverage at that position
	merged$norm_rt_end<- round(merged$Ends/merged$coverage, 3) #to normalize RTdrop by the coverage at that position
	merged$mis_freq <- round(merged$mismatch/merged$coverage,3)
	merged2 <- merge(merged,mod_positions, by.x=c("chr","pos"), by.y=c("Reference","Position"),all=TRUE )
	merged2$Unique <- merged2$Unique %>% replace_na('Unm')
	merged2$Mod <- merged2$Mod %>% replace_na('Unm')
	#Relabel unmodified bases 
	unm <- subset(merged2 , Mod=="Unm")
	unm$Mod <- unm$ref_nuc
	mod <- subset(merged2, Mod!="Unm")
	final <- rbind(unm, mod)
	return(final)

}

#Process Mg runs
protoscript_mg.processed <- process_func("data/stats/FIRST_Mg_PS2.STATS", "data/readends/FIRST_Mg_PS2.readends.tsv","PS2", "Mg")
ss2_mg.processed <- process_func("data/stats/FIRST_Mg_SS2.STATS", "data/readends/FIRST_Mg_SS2.readends.tsv","SS2", "Mg")
ss3_mg.processed <- process_func("data/stats/FIRST_Mg_SS3.STATS", "data/readends/FIRST_Mg_SS3.readends.tsv","SS3", "Mg")
ss4_mg.processed <- process_func("data/stats/FIRST_Mg_SS4.STATS", "data/readends/FIRST_Mg_SS4.readends.tsv","SS4", "Mg")
maxima_mg.processed <- process_func("data/stats/FIRST_Mg_Maxima.STATS", "data/readends/FIRST_Mg_Maxima.readends.tsv","Maxima", "Mg")
tgirt_mg.processed <- process_func("data/stats/FIRST_Mg_TGIRT.STATS", "data/readends/FIRST_Mg_TGIRT.readends.tsv","TGIRT", "Mg")
induro_mg.processed <- process_func("data/stats/FIRST_Mg_Induro.STATS", "data/readends/FIRST_Mg_Induro.readends.tsv","Induro", "Mg")
marathon_mg.processed <- process_func("data/stats/FIRST_Mg_Marathon.STATS", "data/readends/FIRST_Mg_Marathon.readends.tsv","Marathon", "Mg")


mg.processed <- rbind(protoscript_mg.processed,ss2_mg.processed,ss3_mg.processed,ss4_mg.processed,maxima_mg.processed,tgirt_mg.processed,induro_mg.processed,marathon_mg.processed)


protoscript_mn.processed <- process_func("data/stats/FIRST_Mn_PS2.STATS", "data/readends/FIRST_Mn_PS2.readends.tsv","PS2", "Mn")
ss2_mn.processed <- process_func("data/stats/FIRST_Mn_SS2.STATS", "data/readends/FIRST_Mn_SS2.readends.tsv","SS2", "Mn")
ss3_mn.processed <- process_func("data/stats/FIRST_Mn_SS3.STATS", "data/readends/FIRST_Mn_SS3.readends.tsv","SS3", "Mn")
ss4_mn.processed <- process_func("data/stats/FIRST_Mn_SS4.STATS", "data/readends/FIRST_Mn_SS4.readends.tsv","SS4", "Mn")
maxima_mn.processed <- process_func("data/stats/FIRST_Mn_Maxima.STATS", "data/readends/FIRST_Mn_Maxima.readends.tsv","Maxima", "Mn")
tgirt_mn.processed <- process_func("data/stats/FIRST_Mn_TGIRT.STATS", "data/readends/FIRST_Mn_TGIRT.readends.tsv","TGIRT", "Mn")
induro_mn.processed <- process_func("data/stats/FIRST_Mn_Induro.STATS", "data/readends/FIRST_Mn_Induro.readends.tsv","Induro", "Mn")
marathon_mn.processed <- process_func("data/stats/FIRST_Mn_Marathon.STATS", "data/readends/FIRST_Mn_Marathon.readends.tsv","Marathon", "Mn")

mn.processed <- rbind(protoscript_mn.processed,ss2_mn.processed,ss3_mn.processed,ss4_mn.processed,maxima_mn.processed,tgirt_mn.processed,induro_mn.processed,marathon_mn.processed)


# 2.1. Lets first try to cluster based on their error signature at the modified sites
mg.onlymod <- subset(mg.processed, Mod %in% c("m1acp3Y", "m1A", "I", "m3C", "m1G", "m22G", "m3U"))
mn.onlymod <- subset(mn.processed, Mod %in% c("m1acp3Y", "m1A", "I", "m3C", "m1G", "m22G", "m3U"))

mg.onlymodmis <- na.omit(mg.onlymod[,c("Enzyme","Unique", "Mod", "mis_freq")])
mn.onlymodmis <- na.omit(mn.onlymod[,c("Enzyme","Unique", "Mod", "mis_freq")])

mg.onlymodmistable <- dcast(mg.onlymodmis, Mod+Unique ~ Enzyme, value.var="mis_freq" )
mn.onlymodmistable <- dcast(mn.onlymodmis, Mod+Unique ~ Enzyme, value.var="mis_freq" )


levels_mg <-c("PS2", "SS2","SS3", "SS4","Maxima", "TGIRT","Induro", "Marathon")
levels_mn <-c("PS2", "SS2","SS3", "SS4","Maxima", "TGIRT","Induro", "Marathon")


#Now PCA
pca_plot <- function(data, label,levels) {
	rownames(data) <- data$Unique
	data2 <- na.omit(data[,c(3:10)])
	data3 <- t(data2)
	pca_data <- prcomp(data3, center = TRUE)
	pca_loadings <- data.frame(RNAs = rownames(pca_data$rotation), pca_data$rotation)
	summ <- summary(pca_data)
	variance <- summ$importance[2,]
	df_out <- as.data.frame(pca_data$x)
	df_out$Enzymes <- rownames(df_out)
	df_out$Enzymes <- factor(df_out$Enzymes, levels = levels)
	pdf(file= paste(label, "Enzymes_ModOnly_Mismatch_PCA.pdf", sep="_"),height=4,width=5.5,onefile=FALSE)
    print(ggplot(df_out,aes(x=PC1,y=PC2, color=Enzymes))+
    	#geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = (PC1*50),yend = (PC2*50)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
    	#annotate("text", x = (pca_loadings$PC1*50000), y = (pca_loadings$PC2*50000),label = pca_loadings$RNAs)+
        geom_point(size=3)+
        theme_bw()+
        ggtitle(label)+
        scale_color_manual(values=c(protoscript_col,ss2_col,ss3_col,ss4_col,maxima_col,tgirt_col,induro_col,marathon_col))+
        xlab(paste("PC1", as.character(round(variance[1]*100)),"%"))+
        ylab(paste("PC2", as.character(round(variance[2]*100)),"%"))+
        geom_text_repel(data=df_out, aes(label=rownames(df_out),x=PC1, y=PC2), colour="black",segment.size  = 0.4,segment.color = "grey50",size=2))
    dev.off()
}
pca_plot(mg.onlymodmistable, "Mg",levels_mg )
pca_plot(mn.onlymodmistable, "Mn",levels_mn)



## Clean mismatch Aggregated
clean_mis_input_agg<-function(data) {
data <- subset(data, coverage> 20)
data<- data[,c("ref_nuc","Mod", "Enzyme","Unique","A", "T","C","G")]
data<- subset(data, ref_nuc %in% c("A","T","C","G"))
data[,c("A", "T", "C","G")] <- round(data[,c("A", "T", "C","G")]/(data$A+data$T+data$C+data$G),3)
data_count <- data %>% count(c("Mod","Enzyme"))
data2<- aggregate(data[5:8],by=list(data$Mod,data$ref_nuc, data$Enzyme),FUN=mean, na.rm=TRUE)
colnames(data2) <- c("Mod", "ref_nuc", "Enzyme", "A", "T","C", "G")
data3 <- melt(data2)
  ref_base_data<- vector()
  for (enzyme in unique(data3$Enzyme)) {
  	subs <- subset(data3, Enzyme==enzyme)
  for (mod in unique(subs$Mod)) {
      subs2<- subset(subs, Mod==mod)
      subs2$variable<- gsub(as.character(unique(subs2$ref_nuc)),"Ref", subs2$variable)
      ref_base_data<- rbind(ref_base_data, subs2)
    }
}
data_final <- merge(ref_base_data,data_count, by.x=c("Enzyme", "Mod"),  by.y=c("Enzyme", "Mod"))
return(data_final)
}

mg.mis <-clean_mis_input_agg(mg.processed)
mn.mis <-clean_mis_input_agg(mn.processed)




# Create a PCA plot based on these profiles

mg.mis$Mod_variable <- paste(mg.mis$Mod, mg.mis$variable, sep="_")
mn.mis$Mod_variable <- paste(mn.mis$Mod, mn.mis$variable, sep="_")

mg.mis.onlymod <- subset(mg.mis, Mod %in% c("m1acp3Y", "m1A", "I", "m3C", "m1G", "m22G", "m3U"))
mn.mis.onlymod <- subset(mn.mis, Mod %in% c("m1acp3Y", "m1A", "I", "m3C", "m1G", "m22G", "m3U"))

mg.mis.table <- dcast(mg.mis.onlymod[,c("Enzyme","Mod_variable", "value")], Mod_variable~Enzyme,value.var = "value" )
mn.mis.table <- dcast(mn.mis.onlymod[,c("Enzyme","Mod_variable", "value")], Mod_variable~Enzyme,value.var = "value" )

pca_plot_basedonbasefreq <- function(data, label) {
	rownames(data) <- data$Mod_variable
	data2 <- data[,c(2:9)]
	data3 <- t(data2)
	pca_data <- prcomp(data3, center = TRUE)
	pca_loadings <- data.frame(RNAs = rownames(pca_data$rotation), pca_data$rotation)
	summ <- summary(pca_data)
	variance <- summ$importance[2,]
	df_out <- as.data.frame(pca_data$x)
	df_out$Enzymes <- rownames(df_out)
	pdf(file= paste(label, "Enzymes_ModOnly_BaseFreq_PCA.pdf", sep="_"),height=4,width=5.5,onefile=FALSE)
    print(ggplot(df_out,aes(x=PC1,y=PC2, color=Enzymes))+
    	#geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = (PC1*50),yend = (PC2*50)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
    	#annotate("text", x = (pca_loadings$PC1*50000), y = (pca_loadings$PC2*50000),label = pca_loadings$RNAs)+
        geom_point(size=3)+
        theme_bw()+
        ggtitle(label)+
        scale_color_manual(values=c(protoscript_col,ss2_col,ss3_col,ss4_col,maxima_col,tgirt_col,induro_col,marathon_col))+
        xlab(paste("PC1", as.character(round(variance[1]*100)),"%"))+
        ylab(paste("PC2", as.character(round(variance[2]*100)),"%"))+
        geom_text_repel(data=df_out, aes(label=rownames(df_out),x=PC1, y=PC2), colour="black",segment.size  = 0.4,segment.color = "grey50",size=2))
    dev.off()
}
pca_plot_basedonbasefreq(mg.mis.table, "Mg")
pca_plot_basedonbasefreq(mn.mis.table, "Mn")








#COmparison pairwise Mg vs Mn
mn.processed.allmod <- subset(mn.processed[,c("Enzyme","Unique", "Mod", "mis_freq")], Unique!="Unm")
mg.processed.allmod <- subset(mg.processed[,c("Enzyme","Unique", "Mod", "mis_freq")], Unique!="Unm")

both.processed.allmod <- merge(mg.processed.allmod,mn.processed.allmod,by=c("Enzyme", "Unique","Mod"))


scatter_plot <- function(data,enzyme){
	subs <- subset(data, Enzyme==enzyme)
	colnames(subs) <- c("Enzyme", "Unique","Mod", "x", "y")
	subs$Diff <- abs(subs$x-subs$y)
    pdf(file=paste(enzyme, "Mg_Mn_Mismatch_Scatter.pdf", sep="_"),height=5,width=5,onefile=FALSE)
    print(ggplot(subs, aes(x=x, y=y)) +
     geom_point(size=1)+
     geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
     geom_point(data=subset(subs, Diff> 0.1), color="red")+
     geom_text_repel(data=subset(subs, Diff> 0.2), aes(label=Mod), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
     xlim(0,1)+
     ylim(0,1)+
     xlab("Mg Mismatch Frequency")+
     ylab("Mn Mismatch Frequency") +
     theme_classic()+
     theme(axis.text.x = element_text(face="bold", color="black",size=11),
             axis.text.y = element_text(face="bold", color="black", size=11),
     plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
     axis.title.x = element_text(color="black", size=15, face="bold"),
     axis.title.y = element_text(color="black", size=15, face="bold"),
     panel.background = element_blank(),
     axis.line = element_line(colour = "black", size=0.5),
     legend.title = element_text(color = "black", size = 20,face="bold"),
     legend.text = element_text(color = "black", size=20)))
dev.off()
}

scatter_plot(both.processed.allmod,"PS2")




library(ggpubr)

#Second approach : Boxplots
bound_both_long <- rbind(mg.processed[,c("Enzyme","Buffer", "Unique", "Mod", "mis_freq")], mn.processed[,c("Enzyme","Buffer", "Unique", "Mod", "mis_freq")])
bound_both_long_allmod <- subset(bound_both_long, Unique!="Unm" & Enzyme!="NA")

pdf(file="Mg_Mn_Mismatch_Boxplot.pdf",height=10,width=10,onefile=FALSE)
            print(ggplot(bound_both_long_allmod, aes(x=Mod, y=mis_freq)) + 
              geom_boxplot(aes(colour=Buffer), position = position_dodge(0.9),) +
              stat_compare_means(aes(group = Buffer),label = "p.signif", method = "kruskal.test")+
              #stat_n_text(y.pos=0)+
              #scale_fill_manual(values = alpha(palette_enzymes, 0.5))+
              #scale_colour_manual(values = alpha(palette_enzymes, 1))+
              facet_wrap(~Enzyme,nrow=4, ncol=2)+
              theme_classic())
dev.off()




levels <-c("PS2", "SS2","SS3", "SS4","Maxima", "TGIRT","Induro", "Marathon")


plot_area_coverage_facet<- function(merged, name,palette) {
	for (chro in c("18s_rRNA", "25s_rRNA")){
	subs<- subset(merged, chr==chro)
	subs$pos <- subs$pos-14
	subs$position <- paste(subs$chr, subs$pos, sep="_")
	subs$Enzyme <- factor(subs$Enzyme, levels = levels_mg)
	pdf(file= paste(chro,name,"_coverage_facet.pdf",sep=""),height=10,width=7 ,onefile=FALSE)
	print(ggplot(subs, aes(x=pos, y=norm_cov,colour=Enzyme)) +
    geom_line() +
    geom_area(aes(fill = Enzyme, group = Enzyme),alpha = 1/3, position = 'identity')
	+theme_classic()
	+scale_fill_manual(values=palette)
	+scale_colour_manual(values=palette)
	+geom_vline(data=subset(subs, position=="18s_rRNA_1191"| position=="25s_rRNA_645" | 
			position=="25s_rRNA_2142"| position=="25s_rRNA_2634"|position=="25s_rRNA_2843")
			, aes(xintercept=pos),linetype="dashed")
	+facet_wrap(~Enzyme,scales="free", nrow=length(unique(subs$Enzyme))))
	dev.off()
	}
}

plot_area_coverage_facet(mg.processed,"Various_Enzymes_Mg",enzyme_palette)
plot_area_coverage_facet(mn.processed,"Various_Enzymes_Mn",enzyme_palette)







plot_area_rtdrop_facet<- function(data, name,palette) {
	data2 <- subset(data, Biotype=="rRNA")
	data3 <- subset(data2, coverage>25)
	mg <- subset(data3, Buffer=="Mg")
	mn <- subset(data3, Buffer=="Mn")
	mn$norm_rt <- mn$norm_rt*-1 
	merged <- rbind(mg, mn)
	for (chro in c("18s_rRNA", "25s_rRNA")){
	subs<- subset(merged, chr==chro)
	subs$pos <- subs$pos-14
	subs$position <- paste(subs$chr, subs$pos, sep="_")
	subs$Enzyme <- factor(subs$Enzyme, levels = levels)
	pdf(file= paste(chro,name,"_rtdrop_facet.pdf",sep=""),height=10,width=7,onefile=FALSE)
	print(ggplot(subs, aes(x=pos, y=norm_rt,colour=Enzyme)) +
    geom_line(stat="identity")
	+scale_fill_manual(values=palette)
	+scale_colour_manual(values=palette)
    +geom_area(aes(fill = Enzyme, group = Enzyme),alpha = 1/2, position = 'identity')
	+theme_classic()
	+geom_vline(data=subset(subs, position=="18s_rRNA_1191"| position=="25s_rRNA_645" | 
				position=="25s_rRNA_2142"| position=="25s_rRNA_2634"|position=="25s_rRNA_2843")
				, aes(xintercept=pos),linetype="dashed",size=0.1)
	+geom_hline(aes(yintercept=0), linetype="dashed",size=0.1)
	+facet_wrap(~Enzyme, nrow=length(unique(subs$Enzyme)))
	)
	dev.off()
	}
}

plot_area_rtdrop_facet(mg.processed,"Various_Enzymes_Mg",enzyme_palette)
plot_area_rtdrop_facet(mn.processed,"Various_Enzymes_Mn",enzyme_palette)






#Create a table for modifications on rRNAs
#Name of modifications
mod_name<- c("m1acp3Y", "m1A", "m1A", "m3U", "m3U")
#Position of mods by chromosome
mod_position<-c("18s_rRNA_1191", "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843")
#Chromosome
mod_chr<- c("18s_rRNA","25s_rRNA", "25s_rRNA","25s_rRNA","25s_rRNA")
#Numeric Position
pos_num<-c( 1191, 645,2142,2634,2843)
#Merge the data
mod<-as.data.frame(cbind(mod_name, mod_position,mod_chr))
#Rename the columns
colnames(mod)<- c("ModName", "position","Chr")
#Add the numeric positions
mod<- cbind(mod, pos_num)

# Only focus on the WC positions with a 10 nt shift (with a 6 nt window)
rtdrop_window<- function(enzymes) {
	collapse_window<-vector()
	for (modification in 1:nrow(mod)){
		mod_pos <- mod[modification,c("pos_num")]
		neigh<- vector()
		neigh<- c(mod_pos+8, mod_pos+9, mod_pos+10,mod_pos+11,mod_pos+12, mod_pos+13, mod_pos+14, mod_pos+15, mod_pos+16,mod_pos+17,mod_pos+18)
		ref<- as.data.frame(c("8", "9", "10","11", "12", "13", "14", "15","16", "17", "18"))
		neigh_ref<- cbind(neigh,ref)
		colnames(neigh_ref) <- c("pos", "ref")
		only_mod_chr<- subset(enzymes,chr == as.character(mod[modification,c("Chr")]))
		only_mod_chr_pos<- join(only_mod_chr, neigh_ref, by="pos")
		only_mod_chr_pos$Unique <- NULL
		only_mod_chr_pos$Mod <- NULL
		only_mod_chr_pos<- na.omit(only_mod_chr_pos)
		only_mod_chr_pos$pos_num<- rep(mod_pos,nrow(only_mod_chr_pos))
		only_mod_chr_pos<- only_mod_chr_pos[,c("Enzyme", "pos_num", "norm_rt", "norm_cov")]
		aggdata <-aggregate(only_mod_chr_pos, by=list(Enzyme=only_mod_chr_pos$Enzyme), FUN=max)
		collapse_window<- rbind(collapse_window,aggdata)
		collapse_window2 <-join(collapse_window,mod,by="pos_num")
		collapse_window2$unique<- paste(collapse_window2$ModName, collapse_window2$position)
		collapse_window2$Enzyme <- NULL
		print(modification)
	}
	return(collapse_window2)
}

mg.rtdrop <- rtdrop_window(mg.processed)
mg.rtdrop.table <- dcast(mg.rtdrop[,c("Enzyme", "unique", "norm_rt")], unique ~Enzyme,value.var = "norm_rt")

mn.rtdrop <- rtdrop_window(mn.processed)
mn.rtdrop.table <- dcast(mn.rtdrop[,c("Enzyme", "unique", "norm_rt")], unique ~Enzyme,value.var = "norm_rt")


pca_plot_basedonrtdrop <- function(data, label,levels) {
	rownames(data) <- data$unique
	data2 <- data[,c(2:9)]
	data3 <- t(data2)
	pca_data <- prcomp(data3, center = TRUE)
	pca_loadings <- data.frame(RNAs = rownames(pca_data$rotation), pca_data$rotation)
	summ <- summary(pca_data)
	variance <- summ$importance[2,]
	df_out <- as.data.frame(pca_data$x)
	df_out$Enzymes <- rownames(df_out)
	df_out$Enzymes <- factor(df_out$Enzymes, levels = levels)
	pdf(file= paste(label, "Enzymes_ModOnly_RTDrop_PCA.pdf", sep="_"),height=4,width=5.5,onefile=FALSE)
    print(ggplot(df_out,aes(x=PC1,y=PC2, color=Enzymes))+
    	#geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = (PC1*50),yend = (PC2*50)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
    	#annotate("text", x = (pca_loadings$PC1*50000), y = (pca_loadings$PC2*50000),label = pca_loadings$RNAs)+
        geom_point(size=3)+
        theme_bw()+
        ggtitle(label)+
        scale_color_manual(values=c(protoscript_col,ss2_col,ss3_col,ss4_col,maxima_col,tgirt_col,induro_col,marathon_col))+
        xlab(paste("PC1", as.character(round(variance[1]*100)),"%"))+
        ylab(paste("PC2", as.character(round(variance[2]*100)),"%"))+
        geom_text_repel(data=df_out, aes(label=rownames(df_out),x=PC1, y=PC2), colour="black",segment.size  = 0.4,segment.color = "grey50",size=2))
    dev.off()
}
pca_plot_basedonrtdrop(mg.rtdrop.table, "Mg",levels_mg)



### Mismatch Barplots
order_mods <- c("Am", "i6A","m26A", "t6A", "ac4C", "Cm", "m5C", "Gm","m7G", "m5U","Um","Y","m1A","I","m3C","m1G","m22G","m1acp3Y","m3U")
mismatch_profile_barplot_facet<- function(data,label) {
data <- subset(data, Mod %in% order_mods)
data$variable <- factor(data$variable, levels = c("A", "T", "C", "G", "Ref"))
data <- transform(data, Mod = factor(Mod, levels = order_mods))
pdf(file= paste(label,"_mismatch_profile_barplot.pdf",sep=""),height=8,width=6, onefile=FALSE)
  print(ggplot(data, aes(x=Mod, y=value, fill=variable), width=1) +
      scale_fill_manual(values=c("#1fab89","#eb4d55","#1e56a0", "#f0cf85","#888888"))+
        geom_bar(stat='identity', colour="black")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45,hjust = 1))+
        geom_text(data=data[!duplicated(data$Mod), ],aes(y=1.2, label = freq), size=3, check_overlap = FALSE)+
        facet_wrap(~Enzyme,nrow=length(unique(data$Enzyme))))
        #theme_bw()
  dev.off()
}

mismatch_profile_barplot_facet(mg.mis,"Enzymes_Mg")
mismatch_profile_barplot_facet(mn.mis,"Enzymes_Mn")



#### Dotplot of the modificaitons
library(ggbeeswarm)
tail_dotplot_facet <- function(data, label){
	data <- subset(data, (!is.na(data[,"Enzyme"])))
	data2 <- subset(data, Mod=="m1A" | Mod=="m3C")
  pdf(file=paste(label, "Dotplot.pdf", sep="_"),height=3,width=5,onefile=FALSE)
  print(ggplot(data2, aes(x = Enzyme, y=mis_freq , color = Enzyme)) +
    labs(x = "Method", y = "Mismatch Frequency") +
    geom_quasirandom(varwidth = TRUE, aes(), alpha=0.6)+
    geom_boxplot(width=0.3, fill = "gray80", size = 0.5, alpha = .8) +
    scale_color_manual(values=enzyme_palette)+
    facet_wrap(~ Mod,nrow=1)+
    theme_classic())
  dev.off()
}
tail_dotplot_facet(mn.processed,  "AllEnzymes_Modifications_Mismatch_Mn")



rtdrop_barplot<- function(data,label,palette) {
	data$Enzyme <- factor(data$Enzyme, levels = c("PS2","SS2","SS3","SS4","Maxima","TGIRT","Induro","Marathon"))
	data <- transform(data, Mod = factor(Enzyme, levels = c("PS2","SS2","SS3","SS4","Maxima","TGIRT","Induro","Marathon")))
	pdf(file= paste(label,"_rtdrop_onlymod_barplot.pdf",sep=""),height=3,width=10,onefile=FALSE)
	print(ggplot(data, aes(x=unique, y=norm_rt, fill=Enzyme, width=0.8))
			+scale_fill_manual(values=palette)
			+ylim(c(0,0.5))
			+scale_colour_manual(values=palette)+
			theme_classic()+
	  		geom_bar(stat='identity', position=position_dodge(0.8), colour="black", alpha=0.8))
	dev.off()
}

rtdrop_barplot(mg.rtdrop,"Enzymes_Mg", enzyme_palette)
rtdrop_barplot(mn.rtdrop,"Enzymes_Mn", enzyme_palette)



#### WHICH ENZYME IS BEST FOR WHAT? 

mg.processed.m1a.m3C <- subset(mg.processed[,c("Enzyme", "Unique", "mis_freq", "Mod")], Mod=="m1A" |  Mod=="m3C")
mn.processed.m1a.m3C <- subset(mn.processed[,c("Enzyme", "Unique", "mis_freq", "Mod")], Mod=="m1A" |  Mod=="m3C")

mg.rtdrop.m1Aonly <- subset(mg.rtdrop, ModName=="m1A")
mn.rtdrop.m1Aonly <- subset(mn.rtdrop, ModName=="m1A")


mg.rtdrop.mis.m1a.m3C <- merge(mg.rtdrop.m1Aonly[,c("Enzyme", "position", "norm_rt")],mg.processed.m1a.m3C[,c("Enzyme", "Unique", "mis_freq", "Mod")], by.y=c("Enzyme", "Unique"),by.x=c("Enzyme", "position"), all.y=TRUE)
mn.rtdrop.mis.m1a.m3C <- merge(mn.rtdrop.m1Aonly[,c("Enzyme", "position", "norm_rt")],mn.processed.m1a.m3C[,c("Enzyme", "Unique", "mis_freq", "Mod")], by.y=c("Enzyme", "Unique"),by.x=c("Enzyme", "position"), all.y=TRUE)


mg.processed.m1a.m3C.average <- mg.processed.m1a.m3C %>%
  group_by(Enzyme) %>%
  dplyr::summarise(mean_mis_freq = mean(mis_freq, na.rm = TRUE))


mn.processed.m1a.m3C.average <- mn.processed.m1a.m3C %>%
  group_by(Enzyme) %>%
  dplyr::summarise(mean_mis_freq = mean(mis_freq, na.rm = TRUE))


mg.rtdrop.mis.m1a.m3C.average <- mg.rtdrop.mis.m1a.m3C %>%
  group_by(Enzyme) %>%
  dplyr::summarise(mean_rt = mean(norm_rt, na.rm = TRUE))


mn.rtdrop.mis.m1a.m3C.average <- mn.rtdrop.mis.m1a.m3C %>%
  group_by(Enzyme) %>%
  dplyr::summarise(mean_rt = mean(norm_rt, na.rm = TRUE))


mg.rtdrop.mis.m1a.m3C.average <- mg.rtdrop.mis.m1a.m3C %>%
  group_by(Enzyme) %>%
  dplyr::summarise(mean_rt = mean(norm_rt, na.rm = TRUE),mean_mis_freq = mean(mis_freq, na.rm = TRUE) )

mg.rtdrop.mis.m1a.m3C.average$Sum_Error <- mg.rtdrop.mis.m1a.m3C.average$mean_rt+mg.rtdrop.mis.m1a.m3C.average$mean_mis_freq

all_merged_final_decision <- merge(mg.rtdrop.mis.m1a.m3C.average,mn.processed.m1a.m3C.average, by="Enzyme")

all_merged_final_decision_final <- merge(all_merged_final_decision,mn.rtdrop.mis.m1a.m3C.average,by="Enzyme" )


colnames(all_merged_final_decision_final) <- c("Enzyme", "Mg_RT", "Mg_Mismatch", "Mg_Sum", "Mn_Mismatch", "Mn_RT")

all_merged_final_decision_final$Mn_Sum <- all_merged_final_decision_final$Mn_Mismatch+all_merged_final_decision_final$Mn_RT




barplot_plot <- function(data, feature,ylim) {
	data2 <- data[,c("Enzyme", feature)]
	colnames(data2) <- c("Enzyme", "Feature")
	data2 <- data2[order(-data2$Feature),]
	data2$Enzyme <- factor(x = data2$Enzyme, levels = data2$Enzyme) 

	pdf(file= paste(feature,"mean_barplot.pdf",sep=""),height=3,width=5,onefile=FALSE)
	print(ggplot(data2, aes(x=Enzyme, y=Feature, fill=Enzyme, width=0.8))
		+scale_fill_manual(values=c("PS2"="#aa6f73","SS2"="#eea990","SS3"="#f6e0b5","SS4"="#D193C1","Maxima"="#9C8CC3","TGIRT"="#5294a3","Induro"="#a3d9d9","Marathon"="#60bfae"))
		+scale_colour_manual(values=c("PS2"="#aa6f73","SS2"="#eea990","SS3"="#f6e0b5","SS4"="#D193C1","Maxima"="#9C8CC3","TGIRT"="#5294a3","Induro"="#a3d9d9","Marathon"="#60bfae"))+
		ylab(feature)+
		ylim(ylim)+
		theme_classic()+
			geom_bar(stat='identity', position=position_dodge(0.8), colour="black", alpha=0.8))
	dev.off()
}
barplot_plot(all_merged_final_decision_final,"Mg_RT", c(0,0.3))
barplot_plot(all_merged_final_decision_final,"Mg_Mismatch", c(0,0.7))
barplot_plot(all_merged_final_decision_final,"Mg_Sum",c(0,0.9))


barplot_plot(all_merged_final_decision_final,"Mn_RT", c(0,0.3) )
barplot_plot(all_merged_final_decision_final,"Mn_Mismatch", c(0,0.7) )
barplot_plot(all_merged_final_decision_final,"Mn_Sum",c(0,0.9) )






### GROUPED BOXPLOT MG and MN on m1A and m3C Seperately

# 3.1. Lets first try to cluster based on their error signature at the modified sites
mg.onlydms <- subset(mg.processed, Mod %in% c( "m1A", "m3C"))
mn.onlydms <- subset(mn.processed, Mod %in% c( "m1A", "m3C"))

mg.onlydmsmis <- na.omit(mg.onlydms[,c("Enzyme","Unique", "Mod", "mis_freq")])
mn.onlydmsmis <- na.omit(mn.onlydms[,c("Enzyme","Unique", "Mod", "mis_freq")])


mg.onlydmsmis$Buffer <- "Mg"
mn.onlydmsmis$Buffer <- "Mn"

both.onlydmsmis <- rbind(mg.onlydmsmis,mn.onlydmsmis)

both.onlydmsmis <-c("PS2", "SS2","SS3", "SS4","Maxima", "TGIRT","Induro", "Marathon")


#Color scheme
protoscript_col <- "#aa6f73"
ss2_col <- "#eea990"
ss3_col <- "#f6e0b5"
ss4_col <- "#D193C1"
maxima_col <- "#9C8CC3"
tgirt_col <- "#5294a3"
induro_col <- "#a3d9d9"
marathon_col <- "#60bfae"
#Palette vector
enzyme_palette <- c(protoscript_col,ss2_col,ss3_col,ss4_col,maxima_col,tgirt_col,induro_col,marathon_col)

levels <-c("PS2", "SS2","SS3", "SS4","Maxima", "TGIRT","Induro", "Marathon")
both.onlydmsmis$Enzyme <- factor(both.onlydmsmis$Enzyme, levels = levels)


both.onlydmsmis.m1A <- subset(both.onlydmsmis, Mod=="m1A")
both.onlydmsmis.m3C <- subset(both.onlydmsmis, Mod=="m3C")




#### Dotplot of the modificaitons
library(ggpubr)


  pdf(file="m1A_Mg_Mn_Dotplot.pdf",height=3,width=9,onefile=FALSE)
  print(ggplot(both.onlydmsmis.m1A, aes(x = Enzyme, y=mis_freq, color = Buffer, fill=Buffer)) +
    labs(x = "Method", y = "Mismatch Frequency") +
    geom_jitter(position=position_dodge(width = 0.8),  alpha = 0.6) +   # Using jitter for dodged points
    geom_boxplot(position=position_dodge(width = 0.8), size = 0.5, alpha = .8) +
    scale_fill_manual(values=c("white", "grey"))+
    scale_color_manual(values=c("black", "black"))+
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..))+
    theme_classic())
  dev.off()


  pdf(file="m3C_Mg_Mn_Dotplot.pdf",height=3,width=9,onefile=FALSE)
  print(ggplot(both.onlydmsmis.m3C, aes(x = Enzyme, y=mis_freq , color = Buffer, fill=Buffer)) +
    labs(x = "Method", y = "Mismatch Frequency") +
    geom_jitter(position=position_dodge(width = 0.8),  alpha = 0.6) +   # Using jitter for dodged points
    geom_boxplot(position=position_dodge(width = 0.8), size = 0.5, alpha = .8) +
    scale_fill_manual(values=c("white", "grey"))+
    scale_color_manual(values=c("black", "black"))+
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..))+
    theme_classic())
  dev.off()