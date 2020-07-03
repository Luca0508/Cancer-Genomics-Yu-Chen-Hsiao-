library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dplyr)
data=fread("genie_data_mutations_extended.txt", stringsAsFactors = TRUE)
sample_size=nrow(data)

###############################################
## select SNP and missense mutations
missnse=data[which(data$Variant_Classification=="Missense_Mutation" & data$Variant_Type=="SNP") ]


n_tumor=length(unique(missnse$Tumor_Sample_Barcode))
n_tumor

#############################################################
## create a new column that combine Hugo_Symbol and HGVSp_Short

library(stringr)
missnse$HGVS=str_sub(as.character(missnse$HGVSp_Short), 3, -2)
missnse$Hugo_HGVS=paste(as.character(missnse$Hugo_Symbol), missnse$HGVS, sep="_")

## since there is one row that has missing value in HGVSp_Short, we delete it
missnse=missnse%>%
  filter(Hugo_HGVS!="PRKDC_")


#####################################################
## find the frequency by Hugo_HGVS

count=missnse %>%
  select(Hugo_HGVS)%>%
  group_by(Hugo_HGVS)%>%
  summarize(N=n())

#############################################################
## read pickle file 
## expected probability to get a mutation in each codon 

library(data.table)
## transfer the pickle file in python and export to CSV
mutability=read.csv("mutability.csv")
mutability=as.data.table(mutability)
codon_chr=as.character(mutability$codon)
mutability=mutability[,"codon":=codon_chr]

## merge the count and mutability dataset
mut_merge=merge(mutability,count, by.x = "codon", by.y="Hugo_HGVS", all.x = TRUE, all.y = FALSE)

######################################
## binomial test
get_binom_p=function(x, p, n){
  binom.test(x=x,n=n, p=p, alternative = "greater")$p.value
}

binom_p_value=mapply( get_binom_p, mut_merge$N, mut_merge$mutability, n_tumor)
mut_merge=mut_merge[, "binom_p_value":=binom_p_value]


## accept or reject H0
compare_p_value=function(x){
  if(x>0.01){
    return("not significant")
  }else if (x<0.01){
    return("significant")
  }
}

significant=sapply(mut_merge$binom_p_value, compare_p_value)
mut_merge=mut_merge[, "significant":=significant]


nrow(mut_merge[which(mut_merge$significant=="significant")])
## 24061 codon_id are significant from binomial test with alpha=0.01

## Benjamini & Yekutieli multiple correction method 
library(stats)

by_adjusted_p = p.adjust(p=mut_merge$binom_p_value, method = "BY", n = length(mut_merge$binom_p_value))
significant_by=sapply(by_adjusted_p, compare_p_value)
mut_merge=mut_merge[, "significant_by":=significant_by]
mut_merge=mut_merge[, "by_ajusted_p":=by_adjusted_p]
nrow(mut_merge[which(mut_merge$significant_by=="significant")])
## 5011 codons are significant from by adjusted with q-value=0.01

sum(by_adjusted_p==1)
## adjust 151118 p-value to 1
sum(by_adjusted_p==1)/(length(by_adjusted_p))
## 0.8580936

hist_by=ggplot(data = mut_merge, aes(x=by_adjusted_p))+
  geom_histogram(breaks=seq(0,1,by=0.02))+
  geom_density()+
  labs(x="Benjamini & Yekutieli p-value")

## Benjamini & Hochberg  multiple correction method 
library(stats)

bh_adjusted_p = p.adjust(p=mut_merge$binom_p_value, method = "BH", n = length(mut_merge$binom_p_value))
significant_bh=sapply(bh_adjusted_p, compare_p_value)
mut_merge=mut_merge[, "significant_bh":=significant_bh]
mut_merge=mut_merge[, "bh_adjusted_p":=bh_adjusted_p]

nrow(mut_merge[which(mut_merge$significant_bh=="significant")])
## 9573 codons are significant from bh adjusted with q-value=0.01

sum(bh_adjusted_p==1)
## 0
sum(bh_adjusted_p==1)/(length(bh_adjusted_p))
## 0

hist_bh=ggplot(data = mut_merge, aes(x=bh_adjusted_p))+
  geom_histogram(breaks=seq(0,1,by=0.02))+
  geom_density()+
  labs(x="Benjamini & Hochberg p-value")



## Bonferroni  multiple correction method 
library(stats)

bf_adjusted_p = p.adjust(p=mut_merge$binom_p_value, method = "bonferroni", n = length(mut_merge$binom_p_value))
significant_bf=sapply(bf_adjusted_p, compare_p_value)
mut_merge=mut_merge[, "significant_bf":=significant_bf]
mut_merge=mut_merge[, "bf_adjusted_p":=bf_adjusted_p]


nrow(mut_merge[which(mut_merge$significant_bf=="significant")])
## 2599 codons are significant from bh adjusted with q-value=0.01

sum(bf_adjusted_p==1)
## adjust 171861 p-value to 1
sum(bf_adjusted_p==1)/(length(bf_adjusted_p))
## 0.9758786 

hist_bf=ggplot(data = mut_merge, aes(x=bf_adjusted_p))+
  geom_histogram(breaks=seq(0,1,by=0.02))+
  geom_density()+
  labs(x="Bonferroni p-value")


## combine three histgram from three procedures
library(gridExtra)
grid.arrange(hist_bf, hist_bh, hist_by, ncol=3, nrow=1)


#########################################################
## find gene and hgsv
gene=sapply(strsplit(mut_merge$codon, "_"), "[", 1)
HGVS=sapply(strsplit(mut_merge$codon, "_"), "[", 2)

mut_merge=mut_merge[, "gene":=gene]
mut_merge=mut_merge[, "HGVS":=HGVS]

## rerange the column order
setcolorder(mut_merge, c("codon", "gene", "HGVS",  "N", "mutability", 
                         "binom_p_value", "significant", 
                         "by_ajusted_p", "significant_by" ,
                         "bh_adjusted_p", "significant_bh",
                         "bf_adjusted_p", "significant_bf" ))


## create a dataset for hotspot by BY correction
hot_by=mut_merge[which(mut_merge$significant_by=="significant")]


## create a dataset for hotspot by BH correction
hot_bh=mut_merge[which(mut_merge$significant_bh=="significant")]


library(dplyr)

## count gene in hotspot
ct_gene_bh=
  hot_bh %>%
  select(gene, N) %>%
  group_by(gene) %>%
  summarise(sum_n_mutation=sum(N),count=n())%>%
  arrange(desc(count))
ct_gene_bh

## 858 gene in 9573 hotspot

######################################################
## Top 50 number of significant mutation in genes

ct_gene_bh_50=ct_gene_bh %>%
  slice(1:50)

library(ggplot2)
top50gene=ggplot(ct_gene_bh_50, aes(x=reorder(gene, count), count))+
  geom_bar(stat="identity")+
  xlab("Gene")+
  ylab("Number of Driver Mutations in Genes")+
  coord_flip()+
  ggtitle("Top 50 Number of Driver Mutation \n across Genes")
top50gene



## number of driver mutations across genes
hot_across_gene=ct_gene_bh%>%
  select(gene, count)%>%
  group_by(count)%>%
  summarize(count_count=n())
hot_across_gene %>%
  arrange(desc(count_count))

ggplot(ct_gene_bh, aes(x=count))+
  geom_histogram(binwidth = 1)+
  xlab("Number of Driver Mutation")+
  ylab("Number of Genes")


## there are 275 gene with 1 occurence in hotspot, which is the most category
## 32.05%
275/858
## there are 91 gene with 2 occurence in hotspot, which is the second most
## there are 57 gene with 3 occurence in hotspot, which is the third most
(275+91+57)/858

## top 50 frequently mutated gene
fq_gene_bh_50=hot_bh %>%
  select(gene, N) %>%
  group_by(gene) %>%
  summarise(sum_n_mutation=sum(N),count=n())%>%
  arrange(desc(sum_n_mutation))%>%
  slice(1:50)
fq_gene_bh_50

top50fqgene=ggplot(fq_gene_bh_50, aes(x=reorder(gene, sum_n_mutation), sum_n_mutation))+
  geom_bar(stat="identity")+
  xlab("Gene")+
  ylab("Total Number of Mutations in Genes")+
  coord_flip()+
  ggtitle("Top 50 Freqeuntly Mutated Genes")

top50fqgene

19475/sum(ct_gene_bh$sum_n_mutation)


grid.arrange(top50gene, top50fqgene,
             ncol=2, nrow=1)


## find most significant hotspot

library(ggplot2)
ggplot(ct_gene_bh, aes(x=reorder(gene, -sum_n_mutation), sum_n_mutation))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_blank())+
  xlab("gene")+
  ylab("total number of mutation per gene")

##################################################
## distribution of hotspot

fq_hotspot=hot_bh%>%
  select(codon, N)%>%
  arrange(desc(N))


ggplot(fq_hotspot, aes(x=reorder(codon, -N), N))+
  geom_bar(stat="identity", fill="red")+
  theme(axis.text.x = element_blank())+
  xlab("Driver Mutations")+
  ylab("Total Number of Mutations ")

nrow(fq_hotspot[fq_hotspot$N<=3,])
##2205
nrow(fq_hotspot[fq_hotspot$N==3,])/nrow(fq_hotspot)
## 0.2303

nrow(fq_hotspot[fq_hotspot$N<=10,])
## 8139
nrow(fq_hotspot[fq_hotspot$N<=10,])/nrow(fq_hotspot)
## 0.8502


##################################################
## most frequently codon id 
ct_codon_bh=hot_bh %>%
  select(codon, N)%>%
  arrange(desc(N))
  
ct_codon_bh_50=ct_codon_bh %>%
  slice(1:50)

ggplot(ct_codon_bh_50, aes(x=reorder(codon, N), N))+
  geom_bar(stat="identity")+
  theme(axis.text.y = element_text(size=7))+
  xlab("Driver Mutations")+
  ylab("Total Number of Mutation")+
  coord_flip()

8088/sum(ct_codon_bh$N)

################################################################

nrow(hot_bh[hot_bh$bh_adjusted_p==0])
## there are 93 driver mutation with p-value =0

sign_hotspot=hot_bh[hot_bh$bh_adjusted_p==0]

length(intersect(sign_hotspot$codon, ct_codon_bh_50$codon))
## all top 50 frequently driver mutation are with p-value=0




################################################################
################################################################
## test file
library(data.table)
test=fread("pre_tcga_mutations_data.txt", stringsAsFactors = TRUE)

sample_size_test=nrow(test)
test_mis=test[which(test$Variant_Classification=="Missense_Mutation" & test$Variant_Type =="SNP") ]

n_tumor_test=length(unique(test_mis$Tumor_Sample_Barcode))

#############################################################
## create a new column that combine Hugo_Symbol and HGVSp_Short

library(stringr)
test_mis$HGVS=str_sub(as.character(test_mis$HGVSp_Short ), 3, -2)
test_mis$Hugo_HGVS=paste(as.character(test_mis$Hugo_Symbol), test_mis$HGVS, sep="_")


test_count=test_mis[, .N, by=.(Hugo_HGVS)]

test_merge=merge(mutability,test_count, by.x = "codon", by.y="Hugo_HGVS", all.x = FALSE, all.y = FALSE)


get_binom_p_test=function(x, p){
  binom.test(x=x,n=n_tumor_test, p=p, alternative = "greater")$p.value
}

binom_p_value_test=mapply( get_binom_p_test, test_merge$N, test_merge$mutability)
test_merge=test_merge[, "binom_p_value_test":=binom_p_value_test]

#### compare_p_value function are in read_gene file
significant_test=sapply(test_merge$binom_p_value_test, compare_p_value)
test_merge=test_merge[, "significant_test":=significant_test]

nrow(test_merge[which(significant_test=="significant")])
## 3459 codons are significant from binomial test with alpha=0.01


## Benjamini & Hochberg  multiple correction method 
library(stats)

bh_adjusted_p_test = p.adjust(p=test_merge$binom_p_value_test, method = "BH", n = length(test_merge$binom_p_value_test))
significant_bh_test=sapply(bh_adjusted_p_test, compare_p_value)
test_merge=test_merge[, "significant_bh_test":=significant_bh_test]
nrow(test_merge[which(test_merge$significant_bh_test=="significant")])
## 788 codons are significant from bh adjusted with q-value=0.01

##############################################
## find gene and hgsv
gene=sapply(strsplit(test_merge$codon, "_"), "[", 1)
HGVS=sapply(strsplit(test_merge$codon, "_"), "[", 2)

test_merge=test_merge[, "gene":=gene]
test_merge=test_merge[, "HGVS":=HGVS]

#####################################################################
## create a dataset for hotspot by BH correction
hot_bh_test=test_merge[which(test_merge$significant_bh_test=="significant")]

inter_bh=intersect(hot_bh$codon, hot_bh_test$codon)
length(inter_bh)
## only 306 hotspot are intersect with the original data
## only less than half of codon are intersect

#####################################################################
library(dplyr)

## count gene in hotspot
ct_gene_bh_test=
  hot_bh_test %>%
  select(gene, N) %>%
  group_by(gene) %>%
  summarise(sum_n_mutation=sum(N),count=n())%>%
  arrange(desc(count))
ct_gene_bh_test
## 364 gene in 788 hotspot
sum(ct_gene_bh_test$sum_n_mutation)

ct_gene_bh_50_test=ct_gene_bh_test %>%
  slice(1:50)

## Top 50 number of significant mutation in genes

inter_top50gene_hot=intersect(ct_gene_bh_50_test$gene, ct_gene_bh_50$gene ) 
length(inter_top50gene_hot)
inter_top50gene_hot

compare_include_gene_hot=function(x){
  if(x %in% inter_top50gene_hot){
    return("include")
  }else{
    return("not include")
  }
}


ct_gene_bh_50_test$include=mapply(compare_include_gene_hot, ct_gene_bh_50_test$gene)

nrow(ct_gene_bh_50_test[ct_gene_bh_50_test$include=="include",])



top50gene_test=ggplot(ct_gene_bh_50_test, aes(x=reorder(gene, count), count, 
                                              fill=ifelse(include=="include", 
                                                          "Included in TOP50 \nNumber of Driver Mutations \n across Gene \nin GENIE dataset", 
                                                          "Not Included")))+
  geom_bar(stat="identity")+
  xlab("Genes")+
  ylab("Number of Driver Mutations in Genes")+
  coord_flip()+
  scale_fill_manual(values=c("red", "gray40"))+
  theme(legend.title = element_blank())+
  ggtitle("Top 50 Number of Driver Mutations \nacross Genes in TCGA dataset")

top50gene_test

###############################################
## Top 50 frequently mutated genes

fq_gene_bh_test=
  hot_bh_test %>%
  select(gene, N) %>%
  group_by(gene) %>%
  summarise(sum_n_mutation=sum(N),count=n())%>%
  arrange(desc(sum_n_mutation))
fq_gene_bh_test

fq_gene_bh_test50=fq_gene_bh_test%>%
  slice(1:50)

inter_top50gene=intersect(fq_gene_bh_test50$gene, fq_gene_bh_50$gene ) 
length(inter_top50gene)
inter_top50gene

compare_include_gene=function(x){
  if(x %in% inter_top50gene){
    return("include")
  }else{
    return("not include")
  }
}

fq_gene_bh_test50$include=mapply(compare_include_gene, fq_gene_bh_test50$gene)

nrow(fq_gene_bh_test50[fq_gene_bh_test50$include=="include",])



top50fqgene_test=ggplot(fq_gene_bh_test50, aes(x=reorder(gene, sum_n_mutation), sum_n_mutation, 
                                               fill=ifelse(include=="include", 
                                                           "Included in TOP50 \nfrequently mutated Gene \nof GENIE dataset", 
                                                           "Not Included")))+
  geom_bar(stat="identity")+
  xlab("Genes")+
  ylab("Total Number of Mutations")+
  coord_flip()+
  scale_fill_manual(values=c("red", "gray40"))+
  theme(legend.title = element_blank())+
  ggtitle("Top 50 Freqeuntly Mutated Genes \nin TCGA dataset")

#########################################
grid.arrange(top50gene_test, top50fqgene_test,
             ncol=2, nrow=1)

#################################################
##########################################
## merge the primary site with the hot_bh
test_primary=test_mis %>%
  select(Hugo_HGVS, `Primary site`)%>%
  distinct()


hot_bh_p=merge(hot_bh,test_primary, by.x = "codon", by.y="Hugo_HGVS", all.x = TRUE, all.y = FALSE)
sum(is.na(hot_bh_p$`Primary site`))

length(unique(test_primary$Hugo_HGVS))


#######################################
## distribution of hotspot
fq_hotspot_test=hot_bh_test%>%
  select(codon, N)%>%
  arrange(desc(N))
fq_hotspot_test
ggplot(fq_hotspot_test, aes(x=reorder(codon, -N), N))+
  geom_bar(stat="identity", fill="red")+
  theme(axis.text.x = element_blank())+
  xlab("Driver Mutations")+
  ylab("Total Number of Mutations ")


fq_hotspot_test50=fq_hotspot_test%>%
  slice(1:50)
fq_hotspot_test50


inter_top50hotspot=intersect(ct_codon_bh_50$codon, fq_hotspot_test50$codon ) 
length(inter_top50hotspot)
inter_top50hotspot

compare_include_hot=function(x){
  if(x %in% inter_top50hotspot){
    return("include")
  }else{
    return("not include")
  }
}


fq_hotspot_test50

fq_hotspot_test50$include=mapply(compare_include_hot, fq_hotspot_test50$codon)




ggplot(fq_hotspot_test50, aes(x=reorder(codon, N), N, 
                              fill=ifelse(include=="include", 
                                          "Included in TOP50 Driver \nMutations of GENIE dataset", 
                                          "Not Included")))+
  geom_bar(stat="identity")+
  xlab("Driver Mutations")+
  ylab("Total Number of Mutations")+
  coord_flip()+
  scale_fill_manual(values=c("red", "gray40"))+
  theme(legend.title = element_blank())



###########################################################
## lineage diversity
tumor_test=test_mis%>%
  select(Hugo_HGVS, Hugo_Symbol, HGVS,`Primary site`)


### TP53, KRAS, BRAF


tp53_bh=hot_bh_test %>%
  select(gene, HGVS) %>%
  filter(gene=="TP53")%>%
  distinct()

TP53_tumor=tumor_test %>%
  filter(Hugo_Symbol=="TP53")%>%
  filter_at(vars(HGVS), any_vars(. %in% tp53_bh$HGVS))

TP53_tumor%>%
  group_by(HGVS)%>%
  summarize(hotspot_count=n())%>%
  arrange(desc(hotspot_count))

## R273 R248 R175


TP53_tumor%>%
  filter_at(vars(HGVS), any_vars(. %in% c("R273", "R248", "R175")))%>%
  select(`Primary site`)%>%
  distinct()

TP53_tumor%>%
  group_by(`Primary site`)%>%
  summarize(tumor_count=n())%>%
  arrange(desc(tumor_count))

## "lusc", "ov", "hnsc", "lgg"

tp=TP53_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("lusc", "ov", "hnsc", "lgg")))%>%
  group_by(`Primary site`, HGVS)%>%
  summarize(count=n())%>%
  arrange(count)

### ov: R273 R248 C176 Y220

tp53_ov=TP53_tumor%>%
  filter(`Primary site`=="ov")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("R273", "R248", "C176", "Y220")))
tp53_ov

tp53_ov=tp53_ov%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

tp53_ov=tp53_ov%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()

### hnsc: R273 R248 R175 H179

tp53_hnsc=TP53_tumor%>%
  filter(`Primary site`=="hnsc")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("R273", "R248", "R175", "H179")))


tp53_hnsc=tp53_hnsc%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

tp53_hnsc=tp53_hnsc%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()


### lusc: R158 V157 R273 G245

tp53_lusc=TP53_tumor%>%
  filter(`Primary site`=="lusc")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("R273", "R158", "V157", "G245")))


tp53_lusc=tp53_lusc%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

tp53_lusc=tp53_lusc%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
tp53_lusc


### lgg: R273 H193 P151 R248 	V173 Y220

tp53_lgg=TP53_tumor%>%
  filter(`Primary site`=="lgg")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("R273", "R248", "V173", "H193", "P151", "Y220")))


tp53_lgg=tp53_lgg%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

tp53_lgg=tp53_lgg%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
tp53_lgg

tp53_comb=rbind(tp53_ov, tp53_hnsc, tp53_lgg, tp53_lusc)

tp53_draw=ggplot(tp53_comb, aes(fill=HGVS, x=`Primary site`, y=percentage))+ 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  ggtitle("TP53")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("blue", "red", "dark green", 
                               "wheat", "plum", "#d8b365",
                               "yellowgreen", "#5ab4ac", "pink",
                               "purple", "brown", "tomato"))+
  guides(fill=guide_legend(ncol=2)) +
  xlab("")+
  scale_x_discrete(labels = c("Head and Neck \nsquamous cell carcinoma\n",
                              "Brain Lower \nGrade Glioma\n",
                              "Lung squamous \ncell carcinoma\n",
                              "Ovarian serous \ncystadenocarcinoma\n"))

tp53_draw

#######################################
## KRAS
kras_bh=hot_bh_test %>%
  select(gene, HGVS) %>%
  filter(gene=="KRAS")%>%
  distinct()

kras_tumor=tumor_test %>%
  filter(Hugo_Symbol=="KRAS")%>%
  filter_at(vars(HGVS), any_vars(. %in% kras_bh$HGVS))

kras_tumor%>%
  group_by(HGVS)%>%
  summarize(hotspot_count=n())%>%
  arrange(desc(hotspot_count))

kras_tumor%>%
  group_by(`Primary site`)%>%
  summarize(tumor_count=n())%>%
  arrange(desc(tumor_count))

## coadread luad paad 

kras=kras_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("coadread", "luad",  "paad" )))%>%
  group_by(`Primary site`, HGVS)%>%
  summarize(count=n())%>%
  arrange(count)


kras_comb=kras_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("coadread", "luad",  "paad" )))%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

kras_comb=kras_comb%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
kras_comb


Kras_draw=ggplot(kras_comb, aes(fill=HGVS, x=`Primary site`, y=percentage))+ 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  ggtitle("KRAS")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("purple", "red", "skyblue", 
                               "wheat", "plum", "#d8b365"))+
  guides(fill=guide_legend(ncol=2)) +
  xlab("")+scale_x_discrete(labels = c( "Colorectal adenocarcinoma",
                                        "Lung adenocarcinoma",
                                        "Pancreatic adenocarcinoma"))

Kras_draw



##########################################################
## BRAF
braf_bh=hot_bh_test %>%
  select(gene, HGVS) %>%
  filter(gene=="BRAF")%>%
  distinct()

BRAF_tumor=tumor_test %>%
  filter(Hugo_Symbol=="BRAF")%>%
  filter_at(vars(HGVS), any_vars(. %in% braf_bh$HGVS))

BRAF_tumor%>%
  group_by(HGVS)%>%
  summarize(hotspot_count=n())%>%
  arrange(desc(hotspot_count))

BRAF_tumor%>%
  group_by(`Primary site`)%>%
  summarize(tumor_count=n())%>%
  arrange(desc(tumor_count))

## thca skcm coadread luad

braf=BRAF_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("thca", "skcm", "coadread", "luad" )))%>%
  group_by(`Primary site`, HGVS)%>%
  summarize(count=n())%>%
  arrange(count)

BRAF_tumor%>%
  group_by(`Primary site`, HGVS)%>%
  summarize(count=n())%>%
  arrange(`Primary site`)


BRAF_comb=BRAF_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("thca", "skcm", "coadread", "luad")))%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

BRAF_comb=BRAF_comb%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
BRAF_comb


BRAF_draw=ggplot(BRAF_comb, aes(fill=HGVS, x=`Primary site`, y=percentage))+ 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  ggtitle("BRAF")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c( "red", "skyblue", 
                                "pink", "plum", "wheat"))+
  guides(fill=guide_legend(ncol=2)) +
  xlab("")+scale_x_discrete(labels = c("Colorectal adenocarcinoma",
                                       "Lung adenocarcinoma",
                                       "Skin Cutaneous Melanoma",
                                       "Thyroid carcinoma"))



BRAF_draw


###################################################
## PIK3CA

pik_bh=hot_bh_test %>%
  select(gene, HGVS) %>%
  filter(gene=="PIK3CA")%>%
  distinct()

pik_tumor=tumor_test %>%
  filter(Hugo_Symbol=="PIK3CA")%>%
  filter_at(vars(HGVS), any_vars(. %in% pik_bh$HGVS))

pik_tumor%>%
  group_by(HGVS)%>%
  summarize(hotspot_count=n())%>%
  arrange(desc(hotspot_count))

pik_tumor%>%
  group_by(`Primary site`)%>%
  summarize(tumor_count=n())%>%
  arrange(desc(tumor_count))

## coadread hnsc lusc stad blca
pik=pik_tumor%>%
  filter_at(vars(`Primary site`), any_vars(. %in% c("coadread",  "hnsc",  "lusc", "stad", "blca" )))%>%
  group_by(`Primary site`, HGVS)%>%
  summarize(count=n())%>%
  arrange(count)

### blca: 	E545 H1047 E542
pik_blca=pik_tumor%>%
  filter(`Primary site`=="blca")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("E545", "H1047", "E542")))


pik_blca=pik_blca%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

pik_blca=pik_blca%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
pik_blca

### 	coadread: E545 E542	H1047 	Q546 L239 	N345

pik_coad=pik_tumor%>%
  filter(`Primary site`=="coadread")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("E545", "H1047", "E542",	"Q546", "L239", "N345")))


pik_coad=pik_coad%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

pik_coad=pik_coad%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
pik_coad


### 	hnsc: E545 E542	H1047 	E365 G1007 	M1043

pik_hnsc=pik_tumor%>%
  filter(`Primary site`=="hnsc")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("E545", "H1047", "E542",	"E365", "G1007", 	"M1043")))


pik_hnsc=pik_hnsc%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

pik_hnsc=pik_hnsc%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
pik_hnsc


### 	lusc: E545 E542	H1047 

pik_lusc=pik_tumor%>%
  filter(`Primary site`=="lusc")%>%
  filter_at(vars(HGVS), any_vars(. %in% c("E545", "H1047", "E542")))


pik_lusc=pik_lusc%>%
  select( HGVS, `Primary site`)%>%
  group_by(`Primary site`)%>%
  mutate(tumor_n=n())

pik_lusc=pik_lusc%>%
  select( HGVS, `Primary site`, tumor_n)%>%
  group_by(`Primary site`, HGVS)%>%
  mutate(hgvs_n=n())%>%
  mutate(percentage=hgvs_n/tumor_n)%>%
  distinct()
pik_lusc





pik_comb=rbind(pik_lusc, pik_hnsc, pik_coad)




pik_draw=ggplot(pik_comb, aes(fill=HGVS, x=`Primary site`, y=percentage))+ 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  ggtitle("PIK3CA")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("wheat", "plum", "#d8b365",
                               "yellowgreen", "#5ab4ac", "dark green",
                               "purple", "pink", "tomato"))+
  xlab("")+
  guides(fill=guide_legend(ncol=2)) +
  scale_x_discrete(labels = c("Colorectal adenocarcinoma",
                              "Head and Neck \n squamous cell carcinoma",
                              "Lung squamous \ncell carcinoma"))

pik_draw


grid.arrange(tp53_draw, Kras_draw, BRAF_draw, pik_draw, ncol=1, nrow=4)



