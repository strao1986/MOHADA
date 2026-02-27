library(data.table)
library(openxlsx)

library(dplyr)
library(stringr)

library(tximport)

rm(list = ls())
gc()

path <- 'D://all/004.Graduation/250428_continue/result'
setwd(path)

bp_path <- file.path(path, "01.BLAST","Result.blastp.tsv")

allo_path <- "../Database/AllergenOnline.clean.317.0729.tsv"
iuis_path <- "../Database/IUIS.clean.1510.tsv"

allo <- fread(allo_path) %>% as.data.frame
iuis <- fread(iuis_path) %>% as.data.frame

btype <- "blastp"
blast_path <- bp_path
eval_cutoff <- 1e-7
identity_cutoff <- 60

db <- allo
db_ind <- "Acc_clean"
###########################
blast1 <- fread(blast_path) %>% as.data.frame()
colnames(blast1) <- c('Query_id', 'Subject_id', 'identity', 'alignment_length', 'mismatchs', 'gap_openings',
                  'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit_score', 'query_length', 'subject_length',
                  'query_cov_per_subject', 'query_cov_per_hsp')
### other
#nrow(blast1)
#length(unique(blast1[,1]))
#length(unique(blast1[,2]))


# i suddenly forget what these columns do... so remove them for now
blast1 <- blast1[, -15:-16]

### filter blast(np) result
blast2 <- subset(blast1, evalue < eval_cutoff & identity > identity_cutoff)

nrow(blast2)
length(unique(blast2$Query_id))
length(unique(blast2$Subject_id))

### removing version number for consistency...
blast2$subid_clean <- sub("\\..*", "", blast2$Subject_id)
d1 <- merge(db, blast2, by.x = db_ind, by.y = 'subid_clean')
d1 <- d1 %>% select(db, everything())


if(!file.exists('allo1.align.tsv')){
  fwrite(d1, 'allo1.align.tsv', sep = '\t')
}


### remove iuis ids that were overlapped with IUIS database
d2 <- d1 %>% select(Group, iuis_id, Biochemical_name, everything())

  tmpind2 <- d2$iuis_id[!is.na(d2$iuis_id)] %>% unique
  d2 <- d2[which(!d2$iuis_id %in% intersect(tmpind2, iuis$Allergen_name)),]

#nrow(d2)
#length(unique(d2$Query_id))
#length(unique(d2$Accession))


if(!file.exists('allo2.clean.tsv')){
  fwrite(d2, 'allo2.clean.tsv', sep = '\t')
}

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

if(!file.exists('allo3.tpm.tsv')){
  fwrite(d3, 'allo3.tpm.tsv', sep = '\t')
}

length(which(d3$avgTPM>10))

d4<-d3[which(d3$avgTPM>10),]
nrow(d4)
length(unique(d4$Query_id))
length(unique(d4$Accession))

if(!file.exists('allo4.tpmfilt.tsv')){
  fwrite(d4, 'allo4.tpmfilt.tsv', sep = '\t')
}

### summary
tb1 <- table(d2$Biochemical_name, useNA = "always") %>% as.data.frame %>% arrange(desc(Freq)) %>% rename(Allergen = Var1, Freq1 = Freq)
tb2 <- table(d4$Biochemical_name, useNA = "always") %>% as.data.frame %>% arrange(desc(Freq)) %>% rename(Allergen = Var1, FreqTPM = Freq)

tb <- merge(tb1, tb2, by = "Allergen", all.x = T, sort = F)

list3 <- tb$Allergen[!is.na(tb$FreqTPM)] %>% as.character()

for (i in 1: length(list3)) {
  
  tmp_aller <- list3[i]
  
  if(is.na(tmp_aller)){
    tmp_ind3 <- which(is.na(d4$Biochemical_name))
  } else {
    tmp_ind3 <- which(d4$Biochemical_name == tmp_aller)
    }
  
  tb[i, "N_trans"] <- length(unique(d4$Query_id[tmp_ind3]))
  tb[i, "N_allers"] <- length(unique(d4$Subject_id[tmp_ind3]))
  
}

fwrite(tb, "allergenonline.summary.tsv", sep = "\t")


