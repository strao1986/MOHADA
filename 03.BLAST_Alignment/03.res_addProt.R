library(data.table)
#library(openxlsx)

library(dplyr)
library(stringr)


rm(list = ls())
gc()

path <- "D:/all/004.Graduation/250428_continue/result"
setwd(path)
getwd()

allo_res <- fread("allo4.tpmfilt.tsv") %>% as.data.frame
iuis_res <- fread("iuis5.protGroup.tsv") %>% as.data.frame

### adjust illustration from allergenonline
allo_res[which(allo_res$Biochemical_name == "heat shock 70 kda protein"), "Biochemical_name"] <- "heat shock protein 70 kda (hsp70)"
allo_res[which(allo_res$Biochemical_name == "cyclophilin"), "Biochemical_name"] <- "cyclophilin, peptidyl-prolyl cis-trans isomerase"
allo_res[which(allo_res$Biochemical_name == "Euphausia"), "Biochemical_name"] <- "tropomyosin"

### adjust... IUIS
### 25.12 no need as i already formatted the db..
#iuis_res[which(iuis_res$ptgroup == "tropomyosin alpha"), "ptgroup"] <- "tropomyosin"
###############################################

prot <- fread("03.Protein/report.pg_matrix.tsv") %>% as.data.frame()

colnames(prot)[7:15] <- paste0("group", 1:9)
prot <- prot %>% mutate(ptgroup = str_split_i(Protein.Group, ";", 1))


### dont need other information
allo1 <- allo_res %>% select(Query_id, Biochemical_name, Acc_clean, avgTPM) %>% mutate(db = "allergenonline") %>% rename(protein = Biochemical_name, id = Acc_clean)
iuis1 <- iuis_res %>% select(Query_id, ptgroup, accession2, avgTPM) %>% mutate(db = "iuis") %>% rename(protein = ptgroup, id = accession2)

total <- rbind(allo1, iuis1)
fwrite(total, "CombinedDB_BLAST.tsv", sep = "\t")

tb <- table(total$protein) %>% sort(decreasing = T) %>% as.data.frame() %>% rename(Allergen = Var1, FreqTPM = Freq)

#### same protocol as allo...
list3 <- tb$Allergen

for (i in 1: length(list3)) {
  
  tmp_aller <- list3[i]
  
  if(is.na(tmp_aller)){
    tmp_ind3 <- which(is.na(total$protein))
  } else {
    tmp_ind3 <- which(total$protein == tmp_aller)
  }
  
  tb[i, "N_Prots"] <- length(unique(total$Query_id[tmp_ind3]))
  tb[i, "N_Allers"] <- length(unique(total$id[tmp_ind3]))
  
}
tb <- tb %>% arrange(desc(FreqTPM))

fwrite(tb, "CombinedDB_summary.tsv", sep = "\t")


ff <- function(x){ 
  y <- sum(!is.na(x))
  
  if(y >= 5){
    res <- "yes"
    } else {res <- "no"}
    
  return(res)
  }

prot$included <- apply(prot[, 7:15], 1, ff)
prot <- prot %>% filter(included == "yes")

iuis_map <- iuis_res %>% select(Query_id, ptgroup) %>% rename(group = ptgroup)
iuis_map <- iuis_map[!duplicated(iuis_map),]
head(iuis_map)

allo_map <- allo_res %>% select(Query_id, Biochemical_name) %>% rename(group = Biochemical_name)
allo_map <- allo_map[!duplicated(allo_map),]
head(allo_map)

total_map <- rbind(iuis_map, allo_map)
total_map <- total_map[!duplicated(total_map),]


prot2 <- merge(prot, total_map, by.x = "ptgroup", by.y = "Query_id")
prot2[is.na(prot2)] <- 0
prot2[, 8:16] <- log2(prot2[, 8:16] + 1) %>% round(., 2)

tt <- table(prot2$group) %>% sort(decreasing = T) %>% as.data.frame()

fwrite(prot2, "03.Protein/test.protRES.tsv", sep = "\t")
fwrite(tt, "03.Protein/test.protRES_tb.tsv", sep = "\t")


######### 0930. creating a table for illustration
iuis_map2 <- iuis_res %>% select(Query_id, ptgroup, Allergen_name, avgTPM, evalue, identity, Species, accession2) %>% rename(group = ptgroup)
iuis_map2 <- iuis_map2[!duplicated(iuis_map2),]
head(iuis_map2)

allo_map2 <- allo_res %>% select(Query_id, Biochemical_name, IUIS_Allergen, avgTPM, evalue, identity, Species, Acc_clean) %>% 
  rename(group = Biochemical_name, Allergen_name = IUIS_Allergen, accession2 = Acc_clean)
allo_map2 <- allo_map2[!duplicated(allo_map2),]
head(allo_map2)

total_map2 <- rbind(iuis_map2, allo_map2)
total_map2 <- total_map2[!duplicated(total_map2),]


prot3 <- merge(prot, total_map2, by.x = "ptgroup", by.y = "Query_id")
prot3[is.na(prot3)] <- 0
prot3[, 8:16] <- log2(prot3[, 8:16] + 1) %>% round(., 2)

tt2 <- table(prot3$group) %>% sort(decreasing = T) %>% as.data.frame()
prot3$avgProtExp <- apply(prot3[, 8:16], 1, mean)

fwrite(prot3, "total_blast_tb.prot.tsv", sep = "\t")


