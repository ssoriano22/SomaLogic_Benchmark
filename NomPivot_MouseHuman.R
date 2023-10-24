library(dplyr)
library(tidyr)
library(rlang)
library(tidyverse)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(KEGGREST)

#Example gene sets for in-script testing ("simple" and "complex")
musGenes = c("Hmmr", "Tlx3", "Cpeb4")
musGenes2 = c("Alb","Amy2a5","Lyz1","Reg1","Reg3g")
musGenes3 = c("Alb","Amy1","Lyz1","Reg2","Reg3b")
musGenes4 = c("Hmmr", "Tlx3", "ONYXAMBER")
humGenes = c("HMMR", "TLX3", "CPEB4")
humGenes2 = c("ALB","AMY2A","LYZ","REG1A","REG3G")

#FUNCTION: Human --> Mouse Nomenclature Pivot
nomConvert_Human_Mouse <- function(genes){
  print("Starting mapping...")
  #Get EntrezIDs of submitted human genes (key for orthologyDB search)
  egs = mapIds(org.Hs.eg.db, genes, "ENTREZID","SYMBOL")
  #Check to make sure valid genes are submitted
  if(any(is.na(egs))){
    #stop(simpleError("Pivot failed: All genes submitted must be valid in either NCBI (EntrezID) or KEGG DB."))
    print("Not all genes submitted are valid in NCBI (EntrezID) DB. Proceeding with pivot for valid genes only.")
    egs = na.omit(egs)
  }
  print("Mouse genes mapped to EntrezID...")
  #Perform orthology selection with human EntrezIDs
  mapped_genes = select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  print("Mouse --> Human Entrez ID Orthology Match performed...")
  #Add columns for each annotation of interest (human and mouse)
  mapped_genes$HOM_Gene = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "SYMBOL", "ENTREZID")
  mapped_genes$HOM_EnsemblID = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ENSEMBL", "ENTREZID")
  mapped_genes$HOM_EnsemblProtID = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ENSEMBLPROT", "ENTREZID")
  mapped_genes$HOM_GeneDescription = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "GENENAME", "ENTREZID")
  mapped_genes$HOM_PATH = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "PATH", "ENTREZID")
  mapped_genes$HOM_GO = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "GO", "ENTREZID")
  mapped_genes$HOM_Ontology = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ONTOLOGY", "ENTREZID")
  mapped_genes$HOM_UniProt = mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "UNIPROT", "ENTREZID")
  mapped_genes$MUS_Gene = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "SYMBOL", "ENTREZID"))
  mapped_genes$MUS_EnsemblID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ENSEMBL", "ENTREZID"))
  mapped_genes$MUS_EnsemblProtID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ENSEMBLPROT", "ENTREZID"))
  mapped_genes$MUS_MGI = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "MGI", "ENTREZID"))
  mapped_genes$MUS_GeneDescription = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "GENENAME", "ENTREZID"))
  mapped_genes$MUS_GO = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "GO", "ENTREZID"))
  mapped_genes$MUS_Ontology = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ONTOLOGY", "ENTREZID"))
  mapped_genes$MUS_UniProt = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "UNIPROT", "ENTREZID"))
  print("All annotation columns successfully added to NCBI-Mapped DF...")
  #Reorganize/rename columns for clarity
  mapped_genes = mapped_genes %>% dplyr::select(Homo.sapiens,Mus.musculus,HOM_Gene,everything())
  colnames(mapped_genes)[1] = "HOM_EntrezID"
  colnames(mapped_genes)[2] = "MUS_EntrezID"
  #Separate out genes wit no mouse orthology found using entrez ID/orthology.eg.db
  EnzID_mapped_genes = mapped_genes %>% filter(!is.na(MUS_EntrezID)) %>% mutate(MUS_EntrezID = as.character(MUS_EntrezID))
  unmapped_genes = mapped_genes %>% filter(is.na(MUS_EntrezID))
  if(length(unmapped_genes$MUS_EntrezID) == 0){
    print("Gene Nomenclature Conversion: Human --> Mouse Complete. NCBI DB used for all gene pivots.")
    return(EnzID_mapped_genes)
  }
  else{
  #Check for orthologs not detected by entrez search using KEGG data
  KEGG_mapped_genes = unmapped_genes %>% mutate(HOM_NCBIID = paste("ncbi-geneid",HOM_EntrezID,sep = ":")) %>%
                                          mutate(MUS_EntrezID = as.character(MUS_EntrezID)) %>%
                                          #Generate NCBI and KEGGIDs for searching KEGG
                                          rowwise() %>%
                                          mutate(HOM_KEGGID = list(head(keggConv("genes",HOM_NCBIID)))) %>%
                                          #Retrieve human gene results (as list) from human KEGGID
                                          mutate(HOM_KEGG_Res = ifelse(is_empty(HOM_KEGGID), NA, keggGet(HOM_KEGGID))) %>%
                                          drop_na(HOM_KEGG_Res) %>%
                                          #Save KEGG Ontology (KO) ID
                                          rowwise() %>%
                                          mutate(HOM_KO = list(names(HOM_KEGG_Res[[4]]))) %>%
                                          mutate(HOM_KO = ifelse(grepl("[0-9]",HOM_KO), HOM_KO, NA)) %>%
                                          filter(!is.na(HOM_KO)) %>%
                                          #Use KEGG db link to find mouse (mmu) version of gene by KO number (selects 1st result)
                                          rowwise() %>%
                                          mutate(MUS_KEGGID = head(tryCatch(keggLink("mmu", HOM_KO), error=function(e) NA))[1]) %>%
                                          filter(!is.na(MUS_KEGGID)) %>%
                                          #Convert mouse KEGGID to NCBI ID
                                          #rowwise() %>% 
                                          mutate(MUS_NCBIID = sub("mmu:","ncbi-geneid:",MUS_KEGGID)) %>%
                                          # #Fill in MUS_EntrezID column with NCBI ID
                                          mutate(MUS_EntrezID = as.character(sub(".*:", "", MUS_NCBIID))) %>%
                                          #Use Orthology.eg.db to fill in other mouse data columns
                                          mutate(MUS_Gene = mapIds(org.Mm.eg.db, MUS_EntrezID, "SYMBOL", "ENTREZID"),
                                                 MUS_EnsemblID = mapIds(org.Mm.eg.db, MUS_EntrezID, "ENSEMBL", "ENTREZID"),
                                                 MUS_EnsemblProtID = mapIds(org.Mm.eg.db, MUS_EntrezID, "ENSEMBLPROT", "ENTREZID"),
                                                 MUS_MGI = mapIds(org.Mm.eg.db, MUS_EntrezID, "MGI", "ENTREZID"),
                                                 MUS_GeneDescription = mapIds(org.Mm.eg.db, MUS_EntrezID, "GENENAME", "ENTREZID"),
                                                 MUS_GO = mapIds(org.Mm.eg.db, MUS_EntrezID, "GO", "ENTREZID"),
                                                 MUS_Ontology = mapIds(org.Mm.eg.db, MUS_EntrezID, "ONTOLOGY", "ENTREZID"),
                                                 MUS_UniProt = mapIds(org.Mm.eg.db, MUS_EntrezID, "UNIPROT", "ENTREZID")
                                                )
  # TESTING FOR LOOP
  # for (gene in KEGG_mapped_genes$HOM_KO){
  #   #HOM_KEGGID = head(keggConv("genes",gene))
  #   #MUS_KEGGID = head(keggLink("mmu",gene))[1]
  #   paste("KO:",as.character(gene))
  #   #paste("MUS_KEGGID:",as.character(MUS_KEGGID))
  # }
  #Combine mapped and KEGG-mapped genes
  combo_mapped_genes = bind_rows(EnzID_mapped_genes,KEGG_mapped_genes)
  #Print final statement confirming completed execution of the above
  print("Gene Nomenclature Conversion: Human --> Mouse Complete. KEGG Orthology DB used to address gaps in NCBI DB.")
  return(combo_mapped_genes)
  }
}

#orthDB_HtM_mapped_res_df = nomConvert_Human_Mouse(humGenes2)

#FUNCTION: Mouse --> Human Nomenclature Pivot
nomConvert_Mouse_Human <- function(genes){
  print("Starting mapping...")
  #Get EntrezIDs of submitted mouse genes (key for orthologyDB search)
  egs = mapIds(org.Mm.eg.db, genes, "ENTREZID","SYMBOL")
  #Check to make sure valid genes are submitted
  if(any(is.na(egs))){
    #stop(simpleError("Pivot failed: All genes submitted must be valid in either NCBI (EntrezID) or KEGG DB."))
    print("Not all genes submitted are valid in NCBI (EntrezID) DB. Proceeding with pivot for valid genes only.")
    egs = na.omit(egs)
  }
  print("Mouse genes mapped to EntrezID...")
  #Perform orthology selection with mouse EntrezIDs
  #mmug = keys(Orthology.eg.db,"Mus.musculus")
  mapped_genes = select(Orthology.eg.db, egs, "Homo.sapiens", "Mus.musculus")
  print("Mouse --> Human Entrez ID Orthology Match performed...")
  #Add columns for each annotation of interest (human and mouse)
  mapped_genes$MUS_Gene = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "SYMBOL", "ENTREZID"))
  mapped_genes$MUS_EnsemblID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ENSEMBL", "ENTREZID"))
  mapped_genes$MUS_EnsemblProtID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ENSEMBLPROT", "ENTREZID"))
  mapped_genes$MUS_MGI = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "MGI", "ENTREZID"))
  mapped_genes$MUS_GeneDescription = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "GENENAME", "ENTREZID"))
  mapped_genes$MUS_GO = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "GO", "ENTREZID"))
  mapped_genes$MUS_Ontology = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "ONTOLOGY", "ENTREZID"))
  mapped_genes$MUS_UniProt = as.character(mapIds(org.Mm.eg.db, as.character(mapped_genes$Mus.musculus), "UNIPROT", "ENTREZID"))
  mapped_genes$HOM_Gene = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "SYMBOL", "ENTREZID"))
  mapped_genes$HOM_EnsemblID = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ENSEMBL", "ENTREZID"))
  mapped_genes$HOM_EnsemblProtID = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ENSEMBLPROT", "ENTREZID"))
  mapped_genes$HOM_GeneDescription = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "GENENAME", "ENTREZID"))
  mapped_genes$HOM_PATH = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "PATH", "ENTREZID"))
  mapped_genes$HOM_GO = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "GO", "ENTREZID"))
  mapped_genes$HOM_Ontology = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "ONTOLOGY", "ENTREZID"))
  mapped_genes$HOM_UniProt = as.character(mapIds(org.Hs.eg.db, as.character(mapped_genes$Homo.sapiens), "UNIPROT", "ENTREZID"))
  print("All annotation columns successfully added to NCBI-Mapped DF...")
  #Reorganize/rename columns for clarity
  mapped_genes = mapped_genes %>% dplyr::select(Mus.musculus,Homo.sapiens,everything())
  colnames(mapped_genes)[1] = "MUS_EntrezID"
  colnames(mapped_genes)[2] = "HOM_EntrezID"
  #Separate out genes wit no human orthology found using entrez ID/orthology.eg.db
  EnzID_mapped_genes = mapped_genes %>% filter(!is.na(HOM_EntrezID)) %>% mutate(HOM_EntrezID = as.character(HOM_EntrezID))
  unmapped_genes = mapped_genes %>% filter(is.na(HOM_EntrezID))
  if(length(unmapped_genes$MUS_EntrezID) == 0){
      print("Gene Nomenclature Conversion: Mouse --> Human Complete. NCBI DB used for all gene pivots.")
      return(EnzID_mapped_genes)
  }
  else{
      #Check for orthologs not detected by entrez search using KEGG data
      KEGG_mapped_genes = unmapped_genes %>% mutate(MUS_NCBIID = paste("ncbi-geneid",MUS_EntrezID,sep = ":")) %>%
                                              mutate(MUS_EntrezID = as.character(MUS_EntrezID)) %>%
                                              mutate(HOM_EntrezID = as.character(HOM_EntrezID)) %>%
                                              #Generate NCBI and KEGGIDs for searching KEGG
                                              rowwise() %>%
                                              ##################
                                              #NOTE: IF FUNCTION ERRORS WITH HTTP 403 CODE - switch comment/uncomment status for 
                                              #      the two mutate lines directly below and re-run.
                                              mutate(MUS_KEGGID = paste("mmu",MUS_EntrezID,sep = ":")) %>%
                                              # mutate(MUS_KEGGID = keggConv("genes",MUS_NCBIID)[1]) %>%
                                              ##################
                                              filter(!is.na(MUS_KEGGID)) %>%
                                              #Retrieve mouse gene results (as list) from mouse KEGGID
                                              mutate(MUS_KEGG_Res = keggGet(MUS_KEGGID)) %>%
                                              drop_na(MUS_KEGG_Res) %>%
                                              #Save KEGG Ontology (KO) ID
                                              rowwise() %>%
                                              mutate(MUS_KO = ifelse(grepl("K.*",unlist(list(names(MUS_KEGG_Res[[4]]))[1])),
                                                                           unlist(list(names(MUS_KEGG_Res[[4]]))[1]),
                                                                           NA)[1]) %>%
                                              filter(!is.na(MUS_KO)) %>%
                                              #Use KEGG db link to find human (hsa) version of gene by KO number (selects 1st result)
                                              rowwise() %>%
                                              mutate(HOM_KEGGID = head(tryCatch(keggLink("hsa", MUS_KO), error=function(e) NA))[1]) %>%
                                              filter(!is.na(HOM_KEGGID)) %>%
                                              #Convert human KEGGID to NCBI ID
                                              #rowwise() %>%
                                              mutate(HOM_NCBIID = sub("hsa:","ncbi-geneid:",HOM_KEGGID)) %>%
                                              #mutate(HOM_NCBIID = tryCatch(keggConv("ncbi-geneid",HOM_KEGGID), error=function(e) NA)) %>% Not working for some reason, manually creating NCBIID in line above
                                              #Fill in HOM_EntrezID column with NCBI ID
                                              mutate(HOM_EntrezID = as.character(sub(".*:", "", HOM_NCBIID))) %>%
                                              #Use Orthology.eg.db to fill in other mouse data columns
                                              mutate(HOM_Gene = mapIds(org.Hs.eg.db, HOM_EntrezID, "SYMBOL", "ENTREZID"),
                                                     HOM_EnsemblID = mapIds(org.Hs.eg.db, HOM_EntrezID, "ENSEMBL", "ENTREZID"),
                                                     HOM_EnsemblProtID = mapIds(org.Hs.eg.db, HOM_EntrezID, "ENSEMBLPROT", "ENTREZID"),
                                                     HOM_GeneDescription = mapIds(org.Hs.eg.db, HOM_EntrezID, "GENENAME", "ENTREZID"),
                                                     HOM_GO = mapIds(org.Hs.eg.db, HOM_EntrezID, "GO", "ENTREZID"),
                                                     HOM_Ontology = mapIds(org.Hs.eg.db, HOM_EntrezID, "ONTOLOGY", "ENTREZID"),
                                                     HOM_UniProt = mapIds(org.Hs.eg.db, HOM_EntrezID, "UNIPROT", "ENTREZID")
                                              )
      #FOR TESTING KEGG
      # for (gene in KEGG_mapped_genes$MUS_KEGGID){
      #   print(gene)
      #   #HOM_NCBIID = tryCatch(keggConv("ncbi-geneid",gene), error=function(e) NA)
      #   MUS_KEGG_Res = ifelse(is_empty(gene), NA, keggGet(gene))
      #   print(MUS_KEGG_Res)
      # }
      #Combine mapped and KEGG-mapped genes
      combo_mapped_genes = bind_rows(EnzID_mapped_genes,KEGG_mapped_genes)
      #Print final statement confirming completed execution of the above
      print("Gene Nomenclature Conversion: Mouse --> Human Complete. KEGG Orthology DB used to address gaps in NCBI DB.")
      return(combo_mapped_genes)
  }
}

#orthDB_MtH_mapped_res_df = nomConvert_Mouse_Human(musGenes2)
