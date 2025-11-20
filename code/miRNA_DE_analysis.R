# set working environment - change to your local one #

install.packages("extrafont")
install.packages("showtext")
setwd("D:/yahel/phd/Collaboration/Lenveberg/Hagit")

# load the libraries with required functions #

library(reshape2)
library(DESeq2)
library(sva)
library(gamlss)
library(readxl)
library(rolr)
library(ggplot2)
library(GGally)
library(plotly)
library(memisc)
library(rms)
library(pec)
library(dplyr)
library(glmnet)
library(viridis)
library(ggrepel)
library(gridExtra)
library(factoextra)
library(NbClust)
library(psych)  
library(cowplot)
library(egg)
library(grid)
library(ggpubr)
library(ggbreak)
library(stringr)
library(DescTools)
library(enrichR)
library(pheatmap)
library(dplyr)
library(fgsea)
library(stats)
library(edgeR)
library(showtext)
library(grDevices)
library(ggbeeswarm)

set.seed(111)

# Add Calibri font by specifying the path
font_add("Calibri", 
         regular = "C:/Windows/Fonts/calibri.ttf", 
         bold = "C:/Windows/Fonts/calibrib.ttf",
         italic = "C:/Windows/Fonts/calibrii.ttf", 
         bolditalic = "C:/Windows/Fonts/calibriz.ttf")
# Automatically use showtext for new plots
showtext_auto()
###############################
#     Built in functions      #
###############################

## Function that runs Gene set enrichment analysis and plots it based on ggplot & fgsea ##
# gene_list - named vector: numeric; ordered vector by chosen parameter with names corresponding to names of gmt files
# GO_file - gmt list: a gene set file in the gmt format
# pval - numeric: the adusjted p-value threshold to filter significant results by
# min_set - numeric: the minimal set size to include in analysis
# max_set - numeric: the maximal set size to include in analysis
# nperm - numeric: num of permutation to run in the fgsea function
# top_sets - numeric: top down/unregulated sets to show in plot

GSEA = function(gene_list, GO_file, pval,
                min_set = 15, max_set = 500,
                top_sets= 10,color_up = 'firebrick2', color_down = 'dodgerblue2') {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=min_set, ## minimum gene set size
                        maxSize=max_set, ## maximum gene set size
                        scoreType = 'pos') %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  if (nrow(fgRes) == 0) {
    warning("No GO terms enriched after first filtration.")
    return(NULL)
  }
  
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  
  if (nrow(fgRes) == 0) {
    warning("No GO terms enriched after collapsing pathways.")
    return(NULL)
  }
  
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Shared-targets", "non-target")
  filtRes = rbind(head(fgRes, n = top_sets),
                  tail(fgRes, n = top_sets ))
  
  total_up = sum(fgRes$Enrichment == "Shared-targets")
  total_down = sum(fgRes$Enrichment == "non-target")
  header = paste0("Top 10 (Total pathways: Shared=", total_up,", non-target=",    total_down, ")")
  
  colos = setNames(c(color_up, color_down),
                   c("Shared-targets", "non-target"))
  
 g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) + theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}



## Function the ranks miRNAs targets based on how they are shared ##
# target_score_df - dataframe: 
#                     columns: values of mRNA targets, types of experimental supported , miRNA name
#                     rows: mRNA target per miRNA per evidence
# MOI: character vector: if not null (default) subset from target_score_df relevant miRNAs names must match the ones in target_score_df
# support: character - name of column with types of experimental supported in target_score_df (default 'Support')
# target:  character - name of column with mRNA targets in target_score_df (default 'Target...Target.pathway')
# miR_name: character - name of column with miRNAs in target_score_df (default 'miRNA')
# predicted_sign: character - indicator of predicted  types in "support", all other values are experimentally supported (default "V")
# strong_evidence: character - indicator of strong experimentally supported targets in "support" (default "Functional MTI")
# plot: boolean - if TRUE (default) than plots  abar plot of ranked targets

rank_mir_target <- function(target_score_df, MOI= NULL,
                            support = 'Support',target = "Target...Target.pathway", miR_name = "miRNA",
                            predicted_sign = 'V',strong_evidence = 'Functional MTI',
                            plot = T){
  
  # filter by miRs of interest
  if (!is.null(MOI)) {
    target_score_df <- target_score_df[target_score_df[[miR_name]] %in% MOI,]
  }
  
  
  # separate between targets with validation and targets with prediction #
  targetlink_perdicted <- target_score_df[target_score_df[[support]] == predicted_sign,]
  targetlink_validated <- target_score_df[target_score_df[[support]] != predicted_sign,]
  
  
  # count number of evidance 
  targetlink_perdicted_frq <- as.data.frame(table(targetlink_perdicted[[miR_name]],targetlink_perdicted[[target]]))
  targetlink_validated_frq <- as.data.frame(table(targetlink_validated[[miR_name]],targetlink_validated[[target]]))
  
  
  # remove empty values which are results of the table # 
  targetlink_perdicted_frq <- targetlink_perdicted_frq[targetlink_perdicted_frq$Var1 != '',]
  targetlink_validated_frq <- targetlink_validated_frq[targetlink_validated_frq$Var2 != '',]
  
  # sum evidence and penalize for predication or reward for strong validation #
  # penalize by division of the validated counts mean
  # reward by validated counts mean * counts per proteins if strong evidence
  # validated weak evidence stayed the same 
  # this is to use as a score for gsea #
  
  predicted_targetlink_sum <- tapply(targetlink_perdicted_frq$Freq, targetlink_perdicted_frq$Var2, sum)
  
  if (is.null(MOI)) {
    validated_targetlink_sum <- tapply(targetlink_validated_frq$Freq, targetlink_validated_frq$Var2, sum)[-1]
  }else{
    validated_targetlink_sum <- tapply(targetlink_validated_frq$Freq, targetlink_validated_frq$Var2, sum)
  }
  
  
  # penalize / reward #
  predicted_targetlink_sum <- predicted_targetlink_sum/mean(validated_targetlink_sum) #penalize
  
  
  validated_targetlink_sum <- ifelse(names(validated_targetlink_sum) %in%  targetlink_validated[[target]][targetlink_validated[[support]] == strong_evidence],
                                     validated_targetlink_sum + validated_targetlink_sum* mean(validated_targetlink_sum),
                                     validated_targetlink_sum)
  names(validated_targetlink_sum) <- unique(targetlink_validated_frq$Var2)
  
  
  # join scores together #
  
  joint_targetlink_score <- c()
  
  for (name in names(validated_targetlink_sum)) {
    if (name %in% names(predicted_targetlink_sum)) {
      joint_targetlink_score[name] <- validated_targetlink_sum[name] + predicted_targetlink_sum[name]
    }else{
      joint_targetlink_score[name] <- validated_targetlink_sum[name]
    }
  }
  
  joint_targetlink_score <- sort(c(joint_targetlink_score,
                                   predicted_targetlink_sum[!names(predicted_targetlink_sum) %in% names(joint_targetlink_score)]),
                                 decreasing = T)
  if (plot) {
    # plot the targetlink score #
    hist_df <- cbind.data.frame('score' = joint_targetlink_score , 'proteins' = names(joint_targetlink_score))
    hist_df$proteins <- factor(hist_df$proteins, levels = hist_df$proteins)
    hist_df$likelihood <- cut(hist_df$score,breaks = c(0,10,20,40,100),labels = c('Low','Medium','High','Very High'))
    
    hist_tar <- ggplot(data = hist_df[1:50,], aes(y = score,x = proteins,fill = likelihood)) +geom_col() +
      labs(y = expression('miRNA-target interaction score'), x = expression('Targets')) + guides(fill =  guide_legend(title = expression("Interaction likelihood"))) + 
      theme(panel.background = element_blank(),panel.grid = element_blank(),
            axis.line = element_line(colour = 'black',linewidth = 2),
            axis.text.x = element_blank(), axis.ticks.x =  element_blank(),
            axis.text.y = element_text(size = 20),
            axis.title = element_text(face = 'bold',size = 34),
            legend.position = 'top',
            legend.title = element_text(face = 'bold',size = 28),
            legend.text = element_text(size = 24),
            text = element_text(family = "Calibri")) +
      scale_fill_manual(values = c("Low" = "#C6DBEF", "Medium" = "#6BAED6", "High" = "#2171B5", "Very High" = "#08306B")) +
      geom_segment(aes(x = 2, y = hist_df$score[1] + 2.7, xend = 1, yend = hist_df$score[1] + 0.5),
                   arrow = arrow(length = unit(0.25, "cm")),linewidth = 2) +
      geom_segment(aes(x = 3, y = hist_df$score[2] + 3, xend = 2, yend = hist_df$score[2] + 0.5),
                   arrow = arrow(length = unit(0.25, "cm")),linewidth = 2) + 
      geom_segment(aes(x = 5.5, y = hist_df$score[3] + 2, xend = 3, yend = hist_df$score[3] + 0.5),
                   arrow = arrow(length = unit(0.25, "cm")),linewidth = 2) +
      annotate(geom = 'label',label = hist_df$proteins[1:3], x = c(5.3,5.5,8) , y =hist_df$score[1:3] + c(3.75,4,2.75),size = 10,family = "Calibri")
    return(list('score' = joint_targetlink_score, 'plot' = hist_tar))
  }else{
    return(joint_targetlink_score)
  }

}

###############################
#      load & clean data      #
###############################


# read data #
miR_data <- as.data.frame(read_excel('data/expression/GLUT4_EV_miRNA_data.xlsx', sheet = 'miRNA_piRNA',n_max = 1169))
# change the rownames to miRNAs and omit he columns #
rownames(miR_data) <- miR_data$miRNA ; miR_data <- miR_data[,-1]

# filter only miRNA names 
miR_data <- miR_data[,str_detect(colnames(miR_data),pattern = 'UMI|miRNA')]

# define groups #
treatment <- sapply(colnames(miR_data), function(x)(paste(str_split(x,'-')[[1]][1:2],collapse = '_')))
names(treatment) <- NULL

# remove treatments with less than 3 repeats

miR_data_WT2D <- miR_data[,treatment == names(table(treatment)[table(treatment)<3])]
miR_data <- miR_data[,treatment != names(table(treatment)[table(treatment)<3])]

treatment <- treatment[treatment != names(table(treatment)[table(treatment)<3])]

# ---------------------------------------- #
# filter by expression in at least 1 group #
# ---------------------------------------- #

# calculate mean per group per miRNA
miRs_mean_per_group <- apply(miR_data, 1, function(x)(tapply(x, treatment, mean)))


# filter by a threshold of mean 50 UMI #

species <- c()
for (thresh in 1:length(seq(0,100,5))) {
  species[thresh] <- sum(apply(miRs_mean_per_group, 2, function(x)(sum(x >= seq(0,100,5)[thresh]) > 0)))
}

barplot(species~seq(0,100,5), col = 'maroon',
        xlab = 'UMI threshold',ylab = 'miRNA species passing filtration')

umi_filter <- apply(miRs_mean_per_group, 2, function(x)(sum(x >= 25) > 0)) # change the 50 to change threshold 
miR_data_filtered <- miR_data[umi_filter,]


miR_data_WT2D <- miR_data_WT2D[rownames(miR_data_WT2D) %in% rownames(miR_data_filtered),]

# reorder columns based on treatment #
miR_data_filtered <- miR_data_filtered[,match(sort(colnames(miR_data_filtered)),colnames(miR_data_filtered))]

treatment <- sort(treatment)

# add WT_2D to the treatment vector and count matrix #

miR_data_filtered <- cbind(miR_data_filtered,miR_data_WT2D)

treatment <- c(treatment,rep('WT_2D',2))


apply(miR_data_filtered[rownames(miR_data_filtered) %in% c('hsa-miR-122-5p','hsa-miR-16-5p','hsa-miR-486-5p'),],
      1, function(x)(tapply(x, rep(c('G4','WT'),c(7,6)), mean)))

####################################
#       Differential analysis      #
####################################

miR_data_filtered <- miR_data_filtered +1

# prepare data for analysis #

coldata <- data.frame('OE' = sapply(treatment, function(x)(str_split(x,'_')[[1]][1])),
           'formation' = sapply(treatment, function(x)(str_split(x,'_')[[1]][2])),
           'per_group' = paste(sapply(treatment, function(x)(str_split(x,'_')[[1]][1])),sapply(treatment, function(x)(str_split(x,'_')[[1]][2])),sep = '_'),
           row.names = colnames(miR_data_filtered))

# build a DESeq object #
# notice the desgin taks into account both factors and the??r interaction #
dds <- DESeqDataSetFromMatrix(countData = miR_data_filtered,
                              colData = coldata,
                              design = ~ OE + formation + OE:formation)
dds$formation <- factor(dds$formation, levels = c("3D","2D"))

# run DE analysis
DE_analysis <- DESeq(dds)

# data exploration prior to DE analysis #

# extract normalized counts # 
norm_miRs <- counts(DE_analysis,normalized=TRUE)

###############
# QC plotting #
###############

# plot cooks distance for outliers #

par(mar=c(10,5,2,2))
boxplot(log10(assays(DE_analysis)[['cooks']]),range = 0, las =2, ylab = "Cook's Distance")

#plot PCA#
vst <- varianceStabilizingTransformation(dds, fitType = 'local')
pcaData <- plotPCA(vst, intgroup=c("formation", "OE"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=OE, shape=formation)) +
  geom_point(size=4) + theme_minimal_grid() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#plot dispersion plots #
# idealy it will be as flat as you can on y axis
plotDispEsts(DE_analysis,cex = 1.5)

# plot heat map for remaining miRs #

vsd <- varianceStabilizingTransformation(as.matrix(miR_data_filtered),blind = F)

counts(vsd)

# plot the whole data # 
# notice one sample has higher expression in 2 miRs 
df <- as.data.frame(colData(DE_analysis)[,c("OE","formation")])
pheatmap(scale(vsd), cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)

# lets find out who... #
hetmap_outlier <- rownames(vsd)[apply(vsd >= 11, 1, function(x)(sum(x) != 0))]
hetmap_outlier # "hsa-miR-30d-5p" "hsa-miR-3135b"
# these are the ones corresponding to the dispersion plots

# now lets plot the heatmap without them #
pheatmap(scale(vsd[!rownames(vsd) %in% hetmap_outlier,]), cluster_rows=T, show_rownames=TRUE,
         cluster_cols=T, annotation_col=df,cellheight = 10, fontsize = 10)

###########################
# results for DE analysis #
###########################
res <- results(DE_analysis,contrast=c("OE","GLU4","WT"),cooksCutoff = TRUE)
#res <- results(DE_analysis,contrast=c("formation","3D","2D"),cooksCutoff = TRUE)
summary(res)

# save and filter outliers in analysis #
outlier_miR <- rownames(res)[is.na(res$padj)]
res <- res[!is.na(res$padj),]


res.sig <- res[res$padj <= 0.05,]
res.sig <- res.sig[order(res.sig$log2FoldChange,decreasing = T),]

# top unregulated miRs - Not significant
upregulated_miRs <- res[res$log2FoldChange > 0,]
upregulated_miRs <- upregulated_miRs[order(upregulated_miRs$log2FoldChange,decreasing = T),]

head(upregulated_miRs,n =10)

# plot volcano plot plot #

volcan_df <- cbind.data.frame(res,'color' = ifelse(res$padj <= 0.05 & res$log2FoldChange < 0,yes = 'Downregulated',
                                         ifelse(res$padj <= 0.05 & res$log2FoldChange > 0,yes = 'upregulated','Steady')))

volcan <- ggplot(data = volcan_df,aes(x = log2FoldChange,y = -log10(padj),color = color, fill = color)) + 
  geom_jitter(size = 4,alpha = 0.7) + ylab(expression(paste(-log[10],'(FDR)'))) + 
  xlab(expression(paste(log[2],(FC[GLUT4 - WT]))))

volcan <- volcan + theme(panel.background = element_blank(), axis.title = element_text(face = 'bold', size = 34),
                         axis.text = element_text(size = 20),
                         legend.position = 'top',legend.title = element_text(face = 'bold',size = 28),
                         legend.text = element_text(size = 24),panel.grid.major = element_blank(),
                         axis.line = element_line(colour = 'black'),
                         text = element_text(family = "Calibri")) 

volcan <- volcan + scale_color_manual(values = (c('blue','black','red'))) + guides(col = guide_legend('DE miRNAs'),fill="none") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red",linewidth = 2.5)

volcan + geom_label_repel(label = ifelse(!volcan_df$color %in% 'Steady',
                                         sapply(rownames(volcan_df),function(x)(paste(strsplit(x,'-')[[1]][-1],collapse = '-'))),''),
                          size = 8,show.legend = F,max.overlaps = 1000,alpha = 0.8,seed = 42,family = "Calibri") + 
  scale_fill_manual(values = setNames(c("cadetblue1", "NA","mistyrose"), levels(as.factor(volcan_df$color))))

# plot boxplots of chosen miRs #
#I chose only mu miRs of interest  #

miR_OI <- c('hsa-miR-122-5p','hsa-miR-486-5p','hsa-miR-16-5p')

miR_OI_df <- res.sig[rownames(res.sig) %in% miR_OI,]

box_df = cbind.data.frame('miR_exp' = c(t(norm_miRs[rownames(norm_miRs) %in% rownames(miR_OI_df),])),
                 'miR' = rep(sapply(rev(rownames(miR_OI_df)), function(x)(paste(strsplit(x,'-')[[1]][-1],collapse = '-')))
                             ,each = ncol(norm_miRs)),
                 'OE' = rep(sapply(treatment, function(x)(strsplit(x,'_')[[1]][1])),3),
                 'formation' = as.character(rep(sapply(treatment, function(x)(strsplit(x,'_')[[1]][2])),3)))

box_df$OE <- factor(box_df$OE, levels = c('WT','GLU4'))
box_df$formation <- factor(box_df$formation, levels = c('2D','3D'))

# filter out the 2D #
box_df <- box_df[box_df$formation == '3D',]

num_miRNAs <- length(unique(box_df$miR))
dodge_width <- 2.6 / num_miRNAs

# add this line to geom_point if 2D was not filtered out #
# aes(shape = formation, group = interaction(OE, miR))

box_plt <- ggplot(data = box_df, aes(x = OE, y = miR_exp, color = miR)) + 
  geom_boxplot(aes(group = interaction(OE, miR)), 
               position = position_dodge(width = dodge_width), size = 1, alpha = 0) + 
  geom_point(position = position_dodge(width = dodge_width), size = 3) +
  ylab(expression('miRNA expression (UMI)')) +
  theme(panel.background = element_blank(),
        axis.title = element_text(face = 'bold', size = 34),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(face = 'bold',color = 'black', size = 34),
        axis.title.x = element_blank(),
        legend.position = 'top', 
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = 'black'),
        text = element_text(family = "Calibri")) +
  scale_x_discrete(labels = c(expression('WT'),expression('GLUT4')))


########################
# miR targets analysis #
########################

# --------------------- #
# from Target-link 2.0  #
# --------------------- #

# reads targets of miRs downloaded from:
# https://ccb-compute.cs.uni-saarland.de/mirtargetlink2/ #

# targets of chosen miRs 
path_targetlink <- paste(getwd(),'data','miRNAs_targets','miRTargetLink 2.0.csv',sep = '/')
targetlink_targets <- read.csv(path_targetlink)

# experimantlly validated targets
path_targetlink_validated <- paste(getwd(),'data','miRNAs_targets','miRTargetLink 2.0 all miRs  not predicted.csv',sep = '/')
targetlink_targets_validated <- read.csv(path_targetlink_validated)

# predicted targets
path_targetlink_predicted <- paste(getwd(),'data','miRNAs_targets','miRTargetLink 2.0 all miRs  predicted.csv',sep = '/')
targetlink_targets_predicted <- read.csv(path_targetlink_predicted)

#combine both datasets
targetlink_targets_all <- rbind(targetlink_targets_predicted,targetlink_targets_validated)

original_IS <- rank_mir_target(targetlink_targets,plot = T)

#rank_mir_target(targetlink_targets_downregulated,
#                MOI = unique(targetlink_targets_downregulated$miRNA))

rank_mir_target(targetlink_targets_all,
                MOI = unique(targetlink_targets_all$miRNA))

# get all possible trios from the 10 downregulated miRs
mir_combi <- combn(unique(targetlink_targets_all$miRNA), 3)

# Run interaction score analysis for all trios and combine to a matrix
# 657359
# permutation_targets_df <- NULL
# for (trio in 1:ncol(mir_combi)) {
#   if (is.null(permutation_targets_df)) {
#     permutation_targets_df <-  cbind(rank_mir_target(targetlink_targets_all,
#                                                  MOI = mir_combi[,trio],plot = F))
#   }else{
#     add <-  cbind(rank_mir_target(targetlink_targets_all,
#                                                 MOI = mir_combi[,trio],plot = F))
#     permutation_targets_df <- merge(permutation_targets_df,add,by = 0, all = T)
#     colnames(permutation_targets_df)[2:ncol(permutation_targets_df)] <- as.character(1:(ncol(permutation_targets_df)-1))
#     rownames(permutation_targets_df) <- permutation_targets_df[,1]
#     permutation_targets_df <- permutation_targets_df[,-1]
#   }
#   if ((trio %% 1000) == 0 ) {
#     print(trio)
#   }
#   
# }
# 
# save(mir_combi,targetlink_targets_all,rank_mir_target,file = 'for_permutation_all_calc.RData')

# load sequential run 
load('data/miRNAs_targets/permutation_targets_df_sequential.RData')
permutation_targets_df_sequntial <- permutation_targets_df

# load random run 
load('data/miRNAs_targets/permutation_targets_results.RData')

# comapre the location of thge random df and remove overlap
# convert the values to string to use match for faster computation
combi_str <- apply(mir_combi, 2, paste, collapse = "_")
#get location - if location is smaller than the column number of permutation_targets_df_sequntial there's an overlap, than remove
permutation_targets_df <- permutation_targets_df[,match(apply(chosen_miRs, 2, paste, collapse = "_"),
                                                        combi_str) > ncol(permutation_targets_df_sequntial)]

permutation_targets_all <- cbind(permutation_targets_df_sequntial,permutation_targets_df)


# get missing value percentage for each gene target
# low values means that it appeared at least once in a trio and is shared across different trios
missing_from_trio <- sort(apply(permutation_targets_df, 1, function(x)(sum(is.na(x)))),decreasing = T) / ncol(permutation_targets_df)
1 - missing_from_trio[c('IGF1R','AKT3', 'VEGFA','BCL2','TARBP2')]

# count the number of occurrences chosen genes were in the top 10 most likely interactors across the trios 
apply(apply(permutation_targets_df, 2, function(x)(names(sort(x,decreasing = T))[1:10]) %in% c('IGF1R','AKT3', 'VEGFA')), 2, sum)

# count the number of occurrences ALL genes were in the top 10 most likely interactors across the trios
top_ten <- c()
for (gene in rownames(permutation_targets_df)) {
  top_ten[gene] <- sum(apply(apply(permutation_targets_df, 2, function(x)(names(sort(x,decreasing = T))[1:10]) %in% gene), 2, sum))
} # VEGFA was second most abundant while IGF1R was only at 12 combination & AKT3 at 6. BCL2 was the most abundant with 64 occurrences from 120  


# scale interaction score across all permutations by min-max within each trio
scaled_permutation_df <- apply(permutation_targets_df, 2, function(x)( (x - min(x,na.rm = T)) / (max(x,na.rm = T) - min(x,na.rm = T)) ))
rownames(scaled_permutation_df) <- rownames(permutation_targets_df)


# get a permutation p-value based on scoring higher scaled value in other trios -  this is similar to top_ten
original_scaled <- (original_IS$score - min(original_IS$score)) / (max(original_IS$score) - min(original_IS$score))
for (i in c('IGF1R','VEGFA','AKT3',)) {
  print(sum(scaled_permutation_df[rownames(scaled_permutation_df) %in% i,] >= original_scaled[i],na.rm = T))
}


# calculate the density kernal for each trio across all existant targets
density_df <- NULL
for (gene in rownames(scaled_permutation_df)) {
  # remove missing values
  tmp <- unlist(scaled_permutation_df[match(gene,rownames(scaled_permutation_df)),])
  tmp <- tmp[!is.na(tmp)]
  if (is.null(density_df)) {
    density_df <- cbind.data.frame(density(tmp)$y) # get kernel density
    colnames(density_df) <- gene
  }else{
    density_df <- cbind(density_df,density(tmp)$y) # get kernel density
    colnames(density_df)[ncol(density_df)] <- gene
  }
}

# remove targets that have average score across all trios #
# these are likley not a shared target of the trio and mostly predicted or unsupported

mean_target_IS <- apply(scaled_permutation_df, 1, function(x)(mean(x,na.rm = T)))

par(family = "Calibri", cex.lab = 2, cex.axis = 1.5,mgp = c(2.3, 0.8, 0))
hist(mean_target_IS,breaks = 50,
     xlab = expression('Average scaled intreaction score'), ylab = expression('Number of mRNA tragets'),
     main = '')

abline(v = 0.02, col = 'red' ,lty = 2, lwd = 2)

# targets with low mean interaction score across all trios
eligable_targets <- rownames(scaled_permutation_df)[apply(scaled_permutation_df, 1, function(x)(mean(x,na.rm = T))) > 0.02]

# ------------------- #
# cluster with Kmeans #
# ------------------- #

# scale all densities by normalizing to max value to cluster 
scaled_density <- apply(density_df, 2, function(x)(x/max(x)))

# remove targets with low mean IS 
subset_scaled_density <- scaled_density[,eligable_targets]

kmeans_result <- kmeans(t(subset_scaled_density),iter.max = 20,
                        centers = 4, nstart = 25)
table(kmeans_result$cluster)

kmeans_result$cluster[kmeans_result$cluster == 1]
missing_from_trio[names(kmeans_result$cluster[kmeans_result$cluster == 1])]

colors = c('red','orange','blue','green')
for (i in 1:nrow(kmeans_result$centers)) {
  if (i == 1) {
    print(plot(kmeans_result$centers[i,],col = colors[i],type="l",
               lwd=2,ylim = c(0,1), ylab ='Normalized density', xlab = 'Normalized interaction score'))
  }else{
    lines(kmeans_result$centers[i,], col = colors[i], lwd = 2)
  }
}

# the first,red, cluster is a just noise, the tragets mean IS is identical to the general one 
# and all targets have the highest missing value rate of 0.7
kmeans_result$cluster[c('BCL2','VEGFA','AKT3','IGF1R','TARBP2')]


plot(unlist(scaled_permutation_df[match('POGZ',rownames(scaled_permutation_df)),]),x = 1:120,ylim = c(0,1),main = 'POGZ',)
hist(unlist(scaled_permutation_df[match('POGZ',rownames(scaled_permutation_df)),]))
which(diff(sign(diff(density_val$y))) == -2)

# GSEA analysis #
# downloaded the gene sets from https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

enrichmet_results <- list()
for (set in dir(paste(getwd(),'data/Gene_sets',sep = '/'))) {
  print(set)
  pathway_path <- paste(getwd(),'data/Gene_sets/',set,sep = '/')
  pathways <- gmtPathways(pathway_path)
  
  enrichmet_results[[set]] <- GSEA(joint_targetlink_score, GO_file = pathway_path, pval = 0.05)
}


enr_react <- enrichmet_results$reactome.gmt$Plot
new_y_labels <- sapply(enrichmet_results$reactome.gmt$Results$pathway, function(x)(paste(strsplit(x,'_')[[1]][-1],collapse = '_')))
new_y_labels[4] <- 'DISEASES_OF_SIGNAL_TRANSDUCTION_BY\nGROWTH_FACTOR_RECEPTORS'
names(new_y_labels) <- NULL

enr_react + theme(panel.grid = element_blank(),
                  plot.title = element_text(face = 'bold',size = 30,hjust = 1),
                  axis.title.x = element_text(face = 'bold',size = 28), axis.title.y =  element_blank(),
                  axis.text = element_text(size = 22),
                  axis.line = element_line(colour = 'black'),
                  legend.position = 'top',legend.justification = c(1,0),
                  legend.text = element_text(size = 22,hjust = 2),
                  legend.title = element_text(face = 'bold',size = 24,hjust = 2),
                  text = element_text(family = "Calibri")) + guides(fill = 'none') + 
  scale_x_discrete(labels = new_y_labels) +labs(size = 'Set size')


enrichmet_results$wikipathwaysgmt.gmt$Plot
enrichmet_results$GO_terms_all.gmt


enrichmet_results_cds$KEGG.gmt$Results
# plot a single set enrich,emt 
enr.plot <- plotEnrichment(pathway = pathways[['KEGG_MEDICUS_REFERENCE_GF_RTK_PI3K_SIGNALING_PATHWAY']],
               stats = joint_targetlink_score,ticksSize = 0.2)
enr.plot <- enr.plot + labs(title = 'KEGG_MEDICUS_REFERENCE_GF_RTK_PI3K_SIGNALING_PATHWAY') + 
  theme(plot.title = element_text(size = 26,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 20,face = 'bold',hjust = 0.5),
        axis.text = element_text(size = 16,face = 'bold'))

enr.plot$layers[[5]] <- geom_line(color = 'blue', size = 2)
enr.plot$layers[[2]] <- geom_hline(yintercept = enr.plot$layers[[2]]$data$yintercept, color = 'red', linetype = 'dashed', size = 1.4)
enr.plot$layers[[3]] <- geom_hline(yintercept = enr.plot$layers[[3]]$data$yintercept, color = 'red', linetype = 'dashed', size = 1.4)
enr.plot$layers[[4]] <- geom_hline(yintercept = enr.plot$layers[[4]]$data$yintercept, size = 1)

# --------------------- #
# from MicroT-CDS 2.0   #
# --------------------- #

# based on data downloaded from:
# https://dianalab.e-ce.uth.gr/microt_webserver/#/interactions #

path_microT_CDS <- paste(getwd(),'data','miRNAs_targets','microT_Interactions_miRs.tsv',sep = '/')
microT_CDS <- read.table(path_microT_CDS,header = T)[,-1]

# order data by gene name
microT_CDS <- microT_CDS[order(microT_CDS$gene_symbol),]

# filter by genes appearing at least twice - i.e. targeted by several miRs 
duplicated_names <- names(table(microT_CDS$gene_symbol))[table(microT_CDS$gene_symbol) > 1]
microT_CDS_dup <- microT_CDS[microT_CDS$gene_symbol %in% duplicated_names,]

# num of miRs targeting  a gene
apperance_CDS <- colSums(table(microT_CDS_dup$mirna_name,microT_CDS_dup$gene_symbol))  

# mean Interaction score per gene across all miRs targeting it
mean_CDS_score <- tapply(microT_CDS_dup$interaction_score, microT_CDS_dup$gene_symbol, mean) 

# reward factor based on the prodcut of the above
CDS_factor <- apperance_CDS * mean_CDS_score 

# filter based on gene found in both databases 
CDS_factor_filtered <- CDS_factor[names(CDS_factor) %in% names(joint_targetlink_score)]

updated_targetlink_score <- joint_targetlink_score[names(CDS_factor_filtered)] * CDS_factor_filtered 
sort(updated_targetlink_score,decreasing = T)

# run fgsea on all CDS scores and updated ones #

cds_score <- microT_CDS$interaction_score[!microT_CDS$gene_symbol %in% names(CDS_factor)]
names(cds_score) <- microT_CDS$gene_symbol[!microT_CDS$gene_symbol %in% names(CDS_factor)]

cds_score <- c(cds_score,CDS_factor) ; cds_score <- sort(cds_score,decreasing = T)

enrichmet_results_cds <- list()
for (set in dir(paste(getwd(),'data/Gene_sets',sep = '/'))) {
  pathway_path <- paste(getwd(),'data/Gene_sets/',set,sep = '/')
  pathways <- gmtPathways(pathway_path)
  
  enrichmet_results_cds[[set]] <- GSEA(cds_score, GO_file = pathway_path, pval = 0.05)
}

enrichmet_results_cds$go_mf.gmt
enrichmet_results_cds$wikipathwaysgmt.gmt
enrichmet_results_cds$GO_terms_all.gmt



pathways <- gmtPathways(paste(getwd(),'data/Gene_sets/KEGG.gmt',sep = '/'))
# plot a single set enrich,emt 
enr.plot <- plotEnrichment(pathway = pathways[['KEGG_MEDICUS_REFERENCE_GF_RTK_PI3K_SIGNALING_PATHWAY']],
                           stats = joint_targetlink_score,ticksSize = 0.2)
enr.plot <- enr.plot + labs(title = 'KEGG_MEDICUS_REFERENCE_GF_RTK_PI3K_SIGNALING_PATHWAY') + 
  theme(plot.title = element_text(size = 26,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 20,face = 'bold',hjust = 0.5),
        axis.text = element_text(size = 16,face = 'bold'))

enr.plot$layers[[5]] <- geom_line(color = 'blue', size = 2)
enr.plot$layers[[2]] <- geom_hline(yintercept = enr.plot$layers[[2]]$data$yintercept, color = 'red', linetype = 'dashed', size = 1.4)
enr.plot$layers[[3]] <- geom_hline(yintercept = enr.plot$layers[[3]]$data$yintercept, color = 'red', linetype = 'dashed', size = 1.4)
enr.plot$layers[[4]] <- geom_hline(yintercept = enr.plot$layers[[4]]$data$yintercept, size = 1)


##############################
#  analysis of Karus's data  #
##############################

# data from https://www.frontiersin.org/articles/10.3389/fphys.2022.937899/full

Karus_data <- as.data.frame(read_excel('data/Bursak et al/Burask_miRNA_data.XLSX',sheet = 'RawData'))
rownames(Karus_data) <- Karus_data$mature ; Karus_data <- Karus_data[,-1]

Karus_data <- Karus_data[,str_which(colnames(Karus_data),'EV')]

#remove outlier found by Cook's distance 
#Karus_data <- Karus_data[,-match('S116_EV_CLFS',table = colnames(Karus_data))]

Karus_groups <- sapply(colnames(Karus_data), function(x)(strsplit(x,'_')[[1]][3]))

means_miRs_Karus <- apply(Karus_data, 1, function(x)(tapply(x, Karus_groups, mean)))



species <- c()
for (thresh in 1:length(seq(0,100,5))) {
  species[thresh] <- sum(apply(means_miRs_Karus, 2, function(x)(sum(x > seq(0,100,5)[thresh]) > 0)))
}

barplot(species~seq(0,100,5), col = 'navy',
        xlab = 'UMI threshold',ylab = 'miRNA species passing filtration')

Karus_data_filtered <- Karus_data[apply(means_miRs_Karus, 2, function(x)(sum(x > 25) > 0)),]

# build a DESeq object #

dds_karus <- DESeqDataSetFromMatrix(countData = Karus_data_filtered,
                              colData = as.data.frame(Karus_groups),
                              design = ~ Karus_groups)


# run DE analysis
DE_karus_analysis <- DESeq(dds_karus)

# plot cooks distance for outliers #

par(mar=c(10,5,2,2))
boxplot(log10(assays(DE_karus_analysis)[['cooks']]),range = 0, las =2, ylab = "Cook's Distance")

#plot PCA#
vst_kraus <- varianceStabilizingTransformation(dds_karus, fitType = 'local')
pcaData <- plotPCA(vst_kraus, intgroup=c("Karus_groups"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Karus_groups)) +
  geom_point(size=4) + theme_minimal_grid() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#plot dispersion plots #
# idealy it will be as flat as you can on y axis
plotDispEsts(DE_karus_analysis,cex = 1.5)

vst_counts_karus <- assay(vst_kraus)


pheatmap(scale(vst_counts_karus), cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=as.data.frame(Karus_groups),cellheight = 10)

# only 3 comparison so it's ok to save it as so

volvano_plots <- list()
for (first in 1:length(unique(Karus_groups))) {
  mid_list <- list()
  for (second in 1:length(unique(Karus_groups))) {
    if (first >= second) {
      next
    }else{
      res <- results(DE_karus_analysis,
                     contrast=c("Karus_groups",unique(Karus_groups)[first],
                                unique(Karus_groups)[second]),cooksCutoff = TRUE,alpha = 0.1)
      res['color'] <- ifelse(res$log2FoldChange > 0 & res$pvalue <= 0.05,'red',
                             ifelse(res$log2FoldChange < 0 & res$pvalue <= 0.05,'blue','grey'))
      
      p_title <- paste(unique(Karus_groups)[first], unique(Karus_groups)[second],sep = ' vs ')
      
      p <- local(ggplot(data = as.data.frame(res),aes(x = log2FoldChange, y = -log10(pvalue),color = color)) + geom_point(size = 5, alpha = 0.5) + 
        theme(panel.background = element_blank(),axis.line = element_line(colour = 'black'), panel.grid = element_line(colour = 'lightgrey'),
              axis.title = element_text(face = 'bold',size = 20),axis.text = element_text(size = 16),
              plot.title = element_text(face = 'bold',size = 24,hjust = 0.5)) + 
        geom_hline(yintercept = -log10(0.05),color = 'red',linewidth = 2, linetype = 'dashed') + labs(title =p_title) +
        scale_color_manual(values = c('blue' = 'blue','grey' = 'black','red' = 'red'),
                           labels = c( 'blue' = 'Downregulated','grey' = 'Stable','red' = 'Upregulated')) +
        geom_label_repel(label = ifelse(res$pvalue <= 0.05,rownames(res),''),
                         show.legend = F,
                         color = ifelse(res$log2FoldChange <= 0 & res$pvalue < 0.05,'blue',
                                        ifelse(res$log2FoldChange >= 0 & res$pvalue < 0.05,'red','grey')), size = 7))
      
      p_built <- ggplotGrob(p)
      mid_list[[second]] <- local(p_built)
    }

  }
  volvano_plots[[first]] <- mid_list
}

sapply(volvano_plots[[1]][-1], plot)
plot(volvano_plots[[2]][[3]])

summary(res)


##############################
#     proteomic FC plot      #
##############################

EV_FC_proteom <- read.csv('data/FC_proteomics_EVs_mat.csv')
rownames(EV_FC_proteom) <- EV_FC_proteom$Gene.Symbol ; EV_FC_proteom <- EV_FC_proteom[,-c(1,2)]

#reverse the FC since now negative is higher in GLU4 and I want positive value to indicate that
EV_FC_proteom$FC <- -EV_FC_proteom$FC

EV_FC_proteom$color <- ifelse(EV_FC_proteom$FDR <= 0.1 & EV_FC_proteom$FC < 0,yes = 'Depleted',
                              ifelse(EV_FC_proteom$FDR <= 0.1 & EV_FC_proteom$FC > 0,
                                     yes = 'Enriched','Steady'))

volcan_protein <- ggplot(data = EV_FC_proteom,aes(x = FC,y = -log10(FDR),color = color, fill = color)) + 
  geom_jitter(size = 4,alpha = 0.7) + ylab(expression(paste(-log[10],'(FDR)'))) + 
  xlab(expression(paste(log[2],(FC[GLUT4 - WT]))))

volcan_protein <- volcan_protein + theme(panel.background = element_blank(), axis.title = element_text(face = 'bold', size = 34),
                         axis.text = element_text(size = 20),
                         legend.position = 'top',legend.title = element_text(face = 'bold',size = 28),
                         legend.text = element_text(size = 24),panel.grid.major = element_blank(),
                         axis.line = element_line(colour = 'black'),
                         text = element_text(family = "Calibri")) 

volcan_protein <- volcan_protein + scale_color_manual(values = (c('blue','red','black'))) + guides(col = guide_legend('DE miRNAs'),fill="none") + 
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color = "red",size = 2.5)

volcan_protein + geom_label_repel(label = ifelse(abs(EV_FC_proteom$FC) >= 1 & EV_FC_proteom$FDR <= 0.1,
                                         rownames(EV_FC_proteom),''),
                          size = 8,show.legend = F,max.overlaps = 5000,alpha = 0.8,
                          nudge_x = ifelse(EV_FC_proteom$FC > 1,1,0),
                          nudge_y = ifelse(-log10(EV_FC_proteom$FDR) > 1.25 & EV_FC_proteom$FC > 1,0.1,0),family = "Calibri") + 
  scale_fill_manual(values = setNames(c("cadetblue1", "NA","mistyrose"), levels(as.factor(EV_FC_proteom$color))))

# ------------------------------------------------------------- #
# GO analysis results by enrichR for enriched proteins in GLUT4 #
# ------------------------------------------------------------- #

GLUT4_enriched <- read.csv('data/GLUT4_enriched_Terms.csv')[,-1]
#GLUT4_enriched <- read.csv('data/GLUT4_depleted_Terms.csv')[,-1]

GLUT4_enriched <- GLUT4_enriched[GLUT4_enriched$Gene_set %in% c("GO_Biological_Process_2023", "GO_Molecular_Function_2023"),]

# function to choose top k names from a df based in a given rank
# decreasing - to sort in a decreasing order or not
get_top_k <- function(df,rank,k = 10,value, decreasing =F){
  df_ = df[order(df[[rank]],decreasing = decreasing),]
  if (nrow(df_) < k) {
    return(df_[[value]])
  }else{
  return(df_[[value]][1:k])
  }
}

top_GOs <- c()
for (p_set in unique(GLUT4_enriched$Gene_set)) {
  sub_df <- GLUT4_enriched[GLUT4_enriched$Gene_set %in% p_set,]
  top_GOs <- append(top_GOs,get_top_k(sub_df,rank = 'Adjusted.P.value',value = 'Term', k = 10))
  
}

subset_enr <- GLUT4_enriched[GLUT4_enriched$Term %in% top_GOs,]

# remove the suffix of some og the GO-terms 
subset_enr$Term <- sapply(subset_enr$Term,function(x)(paste(strsplit(x,' ')[[1]][-length(strsplit(x,' ')[[1]])],collapse = ' ')))


subset_enr$Term[c(6,10,13,14)] <- c("Regulation Of Vascular Endothelial\nGrowth Factor Signaling Pathway",
                                    "Antigen Processing And Presentation Of\nExogenous Peptide Antigen Via MHC Class II",
                                    "Cysteine-Type Endopeptidase Activator\nActivity Involved In Apoptotic Process",
                                    "Peptidase Activator Activity\nInvolved In Apoptotic Process")


# subset_enr$Term[c(24)] <- c("Reactome 2022: Extracellular Matrix Organization")
# Calculate -log10(Adjusted.P.value)
subset_enr <- subset_enr %>%
  mutate(neg_log10_pval = -log10(Adjusted.P.value))

subset_enr$Gene_set <- factor(subset_enr$Gene_set,
                              levels = c("GO_Molecular_Function_2023", 'GO_Biological_Process_2023'))

# Create a combined factor of Gene_set and Term
subset_enr <- subset_enr %>%
  arrange(Gene_set, neg_log10_pval)

# Set the levels of Term to ensure correct ordering in the plot
subset_enr$Term <- factor(subset_enr$Term, levels = subset_enr$Term)




ggplot(data = subset_enr, aes(y = Term, x = neg_log10_pval, fill = Gene_set)) +
  geom_col() + xlab(expression(paste(-log[10],'(FDR)'))) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold',size = 30,hjust = 1),
        axis.title.x = element_text(face = 'bold',size = 28), axis.title.y =  element_blank(),
        axis.text = element_text(size = 22),
        axis.line = element_line(colour = 'black'),
        legend.position = 'top',legend.justification = c(1,0),
        legend.text = element_text(size = 22),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.y = element_text(lineheight = 0.6,
                                   margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5)))

# ----------------------------------- #
# EVs vs CM QC MISEV2018 protein lvls #
# ----------------------------------- #

color1 <- rgb(0.8662821991541715, 0.2901191849288735, 0.2978085351787774)
color2 <- rgb(0.3686274509803922, 0.30980392156862746, 0.6352941176470588)

EVs_qc_boxplot_df <- read.csv('data/EVs_qc_boxplot_df_GLUT_log2.csv')[,-1]


EVs_qc_boxplot_df$affiliation <- factor(EVs_qc_boxplot_df$affiliation,
                                        levels = levels(as.factor(EVs_qc_boxplot_df$affiliation))[c(2,1,3,5,4)],ordered = T)


anova_qc_protein <- aov(as.formula("proteins_level~Fraction + affiliation + Fraction:affiliation"),data = EVs_qc_boxplot_df)
summary(anova_qc_protein)

for (affiliation in unique(EVs_qc_boxplot_df$affiliation)) {
  subset_qc_dat <- EVs_qc_boxplot_df[EVs_qc_boxplot_df$affiliation == affiliation,]
  sub_t_test <- t.test(as.formula(proteins_level~Fraction),subset_qc_dat,var.equal = T)
  print(paste('for', affiliation,'the p-value is',sub_t_test$p.value,sep = ' '))
}


ggplot(data = EVs_qc_boxplot_df, aes(y = proteins_level, x = affiliation, fill = Fraction)) +
  geom_boxplot(alpha = 0.75) + theme_minimal()+ labs(y = expression(Log[2]~"(intensities)")) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold',size = 30,hjust = 1),
        axis.title.y = element_text(face = 'bold',size = 28), axis.title.x =  element_blank(),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 65,hjust = 1, color = 'black',size = 28),
        axis.line = element_line(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        text = element_text(family = "Calibri")) +
  scale_fill_manual(values = c(color2,color1))  +
  annotate(geom = 'segment', x= , 1.85, xend = 2.15, y = 38, yend = 38, linewidth = 2) + 
  annotate(geom = 'segment', x= , 1.85, xend = 1.85, y = 38.05, yend = 37, linewidth = 2) +
  annotate(geom = 'segment', x= , 2.15, xend = 2.15, y = 38.05, yend = 37, linewidth = 2) +
  annotate(geom = 'segment', x= , 2.85, xend = 3.15, y = 38, yend = 38, linewidth = 2) + 
  annotate(geom = 'segment', x= , 2.85, xend = 2.85, y = 38.05, yend = 37, linewidth = 2) +
  annotate(geom = 'segment', x= , 3.15, xend = 3.15, y = 38.05, yend = 37, linewidth = 2) +
  annotate(geom = 'text',x = c(2,3),y = c(38.2,38.2),label = c('*','***'),size = c(12,12))
  
########################
#       AFM plots      #
########################
AFM_path <- "D:/yahel/phd/Collaboration/Lenveberg/Hagit/data/AFM _size_quantification.xlsx"
wt_afm <- as.data.frame(read_excel(AFM_path,sheet = 'WT',range = "I1:M130")) # B1:F135 for 4um ; I1:M130 for 10um
g4_afm <- as.data.frame(read_excel(AFM_path,sheet = 'Glut1',range = "I1:M174")) # B1:F131 for 4um ; I1:M174 for 10um

# remove top X EVS looks like outliers# 
# wt_afm <- wt_afm[order(wt_afm$`eqR (m)`,decreasing = T),]
# wt_afm <- wt_afm[-c(1:5),]

all_afm <- cbind(rbind(wt_afm,g4_afm), "Treatment" = rep(c('WT','G4OE'),c(nrow(wt_afm),nrow(g4_afm))))
all_afm$Treatment <- factor(all_afm$Treatment, levels = c('WT','G4OE'))
t.test(`eqR (m)` ~ Treatment,data = all_afm,var.equal = T)

EVS_size_means <- tapply(all_afm$`eqR (m)`, all_afm$Treatment, mean) * 1e9
EVS_size_means
# pooled std for Cohen's D
pooled_std <- sqrt(((nrow(wt_afm)-1)*var(wt_afm$`eqR (m)`) + (nrow(g4_afm)-1)*var(g4_afm$`eqR (m)`)) / (nrow(wt_afm) + nrow(g4_afm) - 2)) * 1e9

# Cohen's D
CD <- (EVS_size_means[1] - EVS_size_means[2]) / pooled_std

# beeswarm plot #
ggplot(all_afm, aes(x = Treatment, y = `eqR (m)` * 1e9, color = Treatment)) +
  geom_beeswarm(cex = 2,size=2.5) +  theme_minimal() +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 5, color = "black",fill = 'orange', stroke = 2.5) +
  annotate("segment", x = 0.75, xend = 1.25, y = EVS_size_means[1], yend = EVS_size_means[1], linewidth = 1.5, linetype = 'dashed') +
  annotate("segment", x = 1.75, xend = 2.25, y = EVS_size_means[2], yend = EVS_size_means[2], linewidth = 1.5, linetype = 'dashed') +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold',size = 30,hjust = 1),
        axis.title.y = element_text(face = 'bold',size = 28), axis.title.x =  element_blank(),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(color = 'black',size = 28),
        axis.line = element_line(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        text = element_text(family = "Calibri")) + 
  scale_color_manual(values = c("red", "green")) + guides(color = 'none') + 
  labs(y = 'Radius (nm)') +
  annotate("segment", x = 1, xend = 2, y = 330, yend = 330, linewidth = 2) + 
  annotate("segment", x = 1, xend = 1, y = 331, yend = 310, linewidth = 2) + 
  annotate("segment", x = 2, xend = 2, y = 331, yend = 310, linewidth = 2) + 
  annotate('text', x = 1.5, y = 335, label = '**', size = 20)

# boxplot #
ggplot(all_afm, aes(x = Treatment, y = `eqR (m)` * 1e9, color = Treatment)) +
  geom_boxplot(size = 1.5, outlier.size = 3) +  theme_minimal() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, color = "black",fill = 'orange', stroke = 2.5) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold',size = 30,hjust = 1),
        axis.title.y = element_text(face = 'bold',size = 28), axis.title.x =  element_blank(),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(color = 'black',size = 28),
        axis.line = element_line(colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        text = element_text(family = "Calibri")) + 
  scale_color_manual(values = c("red", "green")) + guides(color = 'none') + 
  labs(y = 'Radius (nm)') +
  annotate("segment", x = 1, xend = 2, y = 330, yend = 330, linewidth = 2) + 
  annotate("segment", x = 1, xend = 1, y = 331, yend = 310, linewidth = 2) + 
  annotate("segment", x = 2, xend = 2, y = 331, yend = 310, linewidth = 2) + 
  annotate('text', x = 1.5, y = 335, label = '**', size = 20)


