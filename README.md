# Multi-Omics integration framework for High-abundance Allergen Discovery in Abalone (MOHADA)  
  MOHADA is an integrated analytical workflow that combines transcriptome assembly, genome annotation, alignment of allergen databases, and post-alignment analysis. This workflow not only involves the systematic analysis and identification of the coding and amino acid sequences for each protein but also determines their actual expression levels by integrating transcriptional and translational data.  
  
  By cross-referencing all identified proteins against authoritative allergen databases, we have first constructed a high-abundance allergen profile for Haliotis discus hannai (HDH). This pipeline could be readily applied to the identification of allergen repertoires in other molluscan and crustacean species. The allergen profile provides specific allergen data to support clinical testing for patients with aquatic food allergies.  

# Overview  
## 1. Transcriptome assembly and Quality assessment  
   1.1. Data preprocessing  
   1.2. Denovo transcriptome assembly  
   1.3. Assemble quality assessment  

## 2. Genome Annotation  
   2.1. Genome preparation  
   2.2. Funannotate pipeline  

## 3. Alignment to representative allergen database  

## 4. Post-alignment analysis 
   4.1. pairwise alignment  
         tools: EMBOSS-Needle  
           
   4.2. Multiple sequence alignment & Phlyogenetic tree construction  
         tools: Clustal Omega, MEGA12  
           
   4.3. B cell epitope prediction  
         Tools for linear B cell epitope prediction: ABCPred, BcePred and Bepipred;  
               for conformational B cell epitope: DiscoTope, SEPPA and CBTOPE.  
         
   
   
Upon completion of the entire analytical workflow, not only can a reliable allergen repertoire be identified for a given species, but also high-confidence, high-abundance allergens can be determined based on transcriptomic and proteomic expression levels.
