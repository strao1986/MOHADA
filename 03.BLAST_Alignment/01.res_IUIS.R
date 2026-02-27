library(data.table)

library(dplyr)
library(stringr)

library(tximport)

rm(list = ls())
gc()

path <- 'D:/all/004.Graduation/250428_continue/result'
setwd(path)

bp_path <- file.path(path, "01.BLAST","Result.blastp.tsv")

allo_path <- "../Database/AllergenOnline.clean.317.tsv"
iuis_path <- "../Database/IUIS.clean.1510.tsv"

allo <- fread(allo_path) %>% as.data.frame
iuis <- fread(iuis_path) %>% as.data.frame


btype <- "blastp"
blast_path <- bp_path
eval_cutoff <- 1e-7
identity_cutoff <- 60

db <- iuis
db_ind <- "accession2"  # accs without version
###########################
blast1 <- fread(blast_path) %>% as.data.frame()
colnames(blast1) <- c('Query_id', 'Subject_id', 'identity', 'alignment_length', 'mismatchs', 'gap_openings',
                      'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score', 'query_length', 'subject_length',
                      'query_cov_per_subject', 'query_cov_per_hsp')

### rm useless column
blast1 <- blast1[, -15:-16]

### filter blast result
blast2 <- subset(blast1, evalue < eval_cutoff & identity > identity_cutoff)

nrow(blast2) 
length(unique(blast2$Query_id)) 
length(unique(blast2$Subject_id)) 

### rm version number for consistency...
blast2$subid_clean <- sub("\\..*", "", blast2$Subject_id)
d1 <- merge(db, blast2, by.x = db_ind, by.y = 'subid_clean')
#d1 <- d1 %>% select(db, everything())

if(!file.exists('iuis1.align.tsv')){
  fwrite(d1, 'iuis1.align.tsv', sep = '\t')
}


d2 <- d1

##### load salmon data ----
files <- list.files(pattern = 'quant.sf', recursive = T, path = '02.Salmon', full.names = T)
list2 <- sapply(strsplit(files, '\\/'), function(x) x[length(x) - 1])
list2 <- gsub('_quant', '', list2)

txi <- tximport(files, type="salmon", countsFromAbundance="lengthScaledTPM", txOut = T)

quant <- txi$abundance %>% as.data.frame()
colnames(quant) <- list2

quant$Name <- rownames(quant)

d3 <- merge(d2, quant, by.x = 'Query_id', by.y = 'Name')

d3$avgTPM <- apply(d3[, (ncol(d3) - length(list2) + 1):ncol(d3)], 1, mean)

if(!file.exists('iuis3.tpm.tsv')){
  fwrite(d3, 'iuis3.tpm.tsv', sep = '\t')
}

d4<-d3[which(d3$avgTPM>10),]
nrow(d4) # 425
length(unique(d4$Query_id)) # 82
length(unique(d4$accession2)) # 162

if(!file.exists('iuis4.tpmfilt.tsv')){
  fwrite(d4, 'iuis4.tpmfilt.tsv', sep = '\t')
}

#### adding cleaned & consistent biochemical name for each record 
### 1. all lower case
d4$Biochemical_name <- tolower(d4$Biochemical_name)
tt <- table(tolower(d4$Biochemical_name)) %>% as.data.frame %>% arrange(Var1) %>% mutate(ptgroup = tt$Var1)
## manually fix the biochemical name...
if(!file.exists('iuis.bioname.clean.tb.tsv')){
  fwrite(tt, "iuis.bioname.clean.tb.tsv", sep = "\t")
}
#### should be manually fixed...----
### adjusting ptgroup...
rm(tt)
tt <- fread("iuis.bioname.clean.tb.tsv") %>% as.data.frame()
####

tt_group <- tt %>% select(Var1, ptgroup)

d5 <- merge(d4, tt_group, by.x = "Biochemical_name", by.y = "Var1")
if(!file.exists('iuis5.protGroup.tsv')){
  fwrite(d5, 'iuis5.protGroup.tsv', sep = '\t')
}

tt2 <- table(d5$ptgroup) %>% as.data.frame %>% arrange(desc(Freq)) %>% rename(Allergen = Var1, FreqTPM = Freq)

tb <- tt2

#### summary
list3 <- tb$Allergen[!is.na(tb$FreqTPM)] %>% as.character()

for (i in 1: length(list3)) {
  
  tmp_aller <- list3[i]
  
  if(is.na(tmp_aller)){
    tmp_ind3 <- which(is.na(d5$ptgroup))
  } else {
    tmp_ind3 <- which(d5$ptgroup == tmp_aller)
  }
  
  tb[i, "N_prots"] <- length(unique(d5$Query_id[tmp_ind3]))
  tb[i, "N_allers"] <- length(unique(d5$Subject_id[tmp_ind3]))
  
}

fwrite(tb, "iuis.summary.tsv", sep = "\t")


