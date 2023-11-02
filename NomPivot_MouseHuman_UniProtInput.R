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

#Example uniprot sets for in-script testing ("simple")
musUniProt = c("Q00547","O55144","Q7TN98")
humUniProt = c("O75330","O43711","Q17RY0")

#FUNCTION: Human --> Mouse Nomenclature Pivot
nomConvert_Human_Mouse <- function(uniprots){
  print("Starting mapping...")
  #Get EntrezIDs of submitted human uniprots (key for orthologyDB search)
  egs = unlist(unname(mapIds(org.Hs.eg.db, uniprots, "ENTREZID","UNIPROT")))
  #Check to make sure valid uniprots are submitted
  if(anyNA(egs)){
    #stop(simpleError("Pivot failed: All uniprots submitted must be valid in either NCBI (EntrezID) or KEGG DB."))
    print("Not all UniProts submitted are valid in NCBI (EntrezID) DB. Proceeding with pivot for valid UniProts only.")
    egs = na.omit(egs)
  }
  print("Human uniprots mapped to EntrezID...")
  #Perform orthology selection with human EntrezIDs
  mapped_uniprots = select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  print("Mouse --> Human Entrez ID Orthology Match performed...")
  #Add columns for each annotation of interest (human and mouse)
  mapped_uniprots$HOM_Gene = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "SYMBOL", "ENTREZID")
  mapped_uniprots$HOM_EnsemblID = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ENSEMBL", "ENTREZID")
  mapped_uniprots$HOM_EnsemblProtID = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ENSEMBLPROT", "ENTREZID")
  mapped_uniprots$HOM_GeneDescription = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "GENENAME", "ENTREZID")
  mapped_uniprots$HOM_PATH = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "PATH", "ENTREZID")
  mapped_uniprots$HOM_GO = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "GO", "ENTREZID")
  mapped_uniprots$HOM_Ontology = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ONTOLOGY", "ENTREZID")
  mapped_uniprots$HOM_UniProt = mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "UNIPROT", "ENTREZID")
  mapped_uniprots$MUS_Gene = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "SYMBOL", "ENTREZID"))
  mapped_uniprots$MUS_EnsemblID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ENSEMBL", "ENTREZID"))
  mapped_uniprots$MUS_EnsemblProtID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ENSEMBLPROT", "ENTREZID"))
  mapped_uniprots$MUS_MGI = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "MGI", "ENTREZID"))
  mapped_uniprots$MUS_GeneDescription = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "GENENAME", "ENTREZID"))
  mapped_uniprots$MUS_GO = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "GO", "ENTREZID"))
  mapped_uniprots$MUS_Ontology = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ONTOLOGY", "ENTREZID"))
  mapped_uniprots$MUS_UniProt = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "UNIPROT", "ENTREZID"))
  print("All annotation columns successfully added to NCBI-Mapped DF...")
  #Reorganize/rename columns for clarity
  mapped_uniprots = mapped_uniprots %>% dplyr::select(Homo.sapiens,Mus.musculus,HOM_Gene,everything())
  colnames(mapped_uniprots)[1] = "HOM_EntrezID"
  colnames(mapped_uniprots)[2] = "MUS_EntrezID"
  #Separate out uniprots with no mouse orthology found using entrez ID/orthology.eg.db
  EnzID_mapped_uniprots = mapped_uniprots %>% filter(!is.na(MUS_EntrezID)) %>% mutate(MUS_EntrezID = as.character(MUS_EntrezID))
  KEGG_mapped_uniprots = mapped_uniprots %>% filter(is.na(MUS_EntrezID)) %>% mutate(HOM_KEGGID = paste("hsa",HOM_EntrezID,sep = ":"))
  if(length(KEGG_mapped_uniprots$MUS_EntrezID) == 0){
    print("Uniprots Nomenclature Conversion: Human --> Mouse Complete. NCBI DB used for all uniprots pivots.")
    return(EnzID_mapped_uniprots)
  }
  else{
    #Check for orthologs not detected by entrez search using KEGG data - split data into subsets of 100 for KEGG map request
    #Create keggLink urls to get KO#
    sub_size = 100
    total_length = nrow(KEGG_mapped_uniprots)
    r  = rep(1:ceiling(total_length/sub_size),each=sub_size)[1:total_length]
    data_frac = split(KEGG_mapped_uniprots,r)
    link_urls = vector("list", length = length(data_frac))
    link_urls = lapply(link_urls,function(x) paste(x,"https://rest.kegg.jp/link/ko/",sep = ""))
    j = 1 #Initialize data fraction counter
    for (sub_df in data_frac){
      i = 1 #Initialize KEGG ID counter - restart per sub_df/url
      current_link = link_urls[j]
      for (ezid in sub_df$HOM_KEGGID){
          if (i == 1){current_link = paste(current_link,ezid,sep = "")}
          else{current_link = paste(current_link,ezid,sep = "+")}
          i = i + 1
      }
      link_urls[j] = current_link
      j = j + 1
    }
    #Submit urls for KO# results
    HOM_KO_reslist = lapply(link_urls,function(x) read_tsv(x,col_names = c("HOM_KEGG_ID","KO_ID")))
    #Create keggLink urls to get mmu KEGGIDs
    mmu_link_urls = vector("list", length = length(data_frac))
    mmu_link_urls = lapply(mmu_link_urls,function(x) paste(x,"https://rest.kegg.jp/link/mmu/",sep = ""))
    j = 1 #Initialize data fraction counter
    for (ko_sub_df in HOM_KO_reslist){
      i = 1 #Initialize KEGG ID counter - restart per sub_df/url
      current_mmu_link = mmu_link_urls[j]
      for (koid in ko_sub_df$KO_ID){
        if (i == 1){current_mmu_link = paste(current_mmu_link,koid,sep = "")}
        else{current_mmu_link = paste(current_mmu_link,koid,sep = "+")}
        i = i + 1
      }
      mmu_link_urls[j] = current_mmu_link
      j = j + 1
    }
    #Submit mmu urls for mmu KEGGID# results
    MUS_KEGGID_res = lapply(mmu_link_urls,function(x) read_tsv(x,col_names = c("KO_ID","MUS_KEGG_ID")))
    #Turn all mmu KEGGID# to entrezIDs
    new_MUS_KEGGID_res = lapply(MUS_KEGGID_res, function(x) cbind(x, MUS_EntrezID = sub("mmu:","",x$MUS_KEGG_ID)))
    #Combine w/ human EntrezID
    MUS_res = do.call("rbind",new_MUS_KEGGID_res)
    HOM_res = do.call("rbind",HOM_KO_reslist)
    HOM_res = HOM_res %>% mutate(HOM_EntrezID = sub("hsa:","",HOM_KEGG_ID))
    HOM_MUS_res = inner_join(HOM_res,MUS_res, by = "KO_ID")
    HOM_MUS_EzID_res = HOM_MUS_res %>% dplyr::select(c("HOM_EntrezID","MUS_EntrezID"))
    #Combine with main df and fill in Orthology.eg.db MUS columns for filled in MUS EntrezIDs
    KEGG_mapped_uniprots = inner_join(KEGG_mapped_uniprots,HOM_MUS_EzID_res,by = "HOM_EntrezID")
    KEGG_mapped_uniprots = KEGG_mapped_uniprots %>% mutate(MUS_EntrezID.x = ifelse((is.na(MUS_EntrezID.x) & !is.na(MUS_EntrezID.y)),
                                                                                   MUS_EntrezID.y,
                                                                                   MUS_EntrezID.x)) %>%
                                                    dplyr::select(-c("MUS_EntrezID.y"))
    colnames(KEGG_mapped_uniprots)[colnames(KEGG_mapped_uniprots) == "MUS_EntrezID.x"] =  "MUS_EntrezID"
    KEGG_mapped_uniprots = KEGG_mapped_uniprots %>% mutate(MUS_Gene = mapIds(org.Mm.eg.db, MUS_EntrezID, "SYMBOL", "ENTREZID"),
                                                            MUS_EnsemblID = mapIds(org.Mm.eg.db, MUS_EntrezID, "ENSEMBL", "ENTREZID"),
                                                            MUS_EnsemblProtID = mapIds(org.Mm.eg.db, MUS_EntrezID, "ENSEMBLPROT", "ENTREZID"),
                                                            MUS_MGI = mapIds(org.Mm.eg.db, MUS_EntrezID, "MGI", "ENTREZID"),
                                                            MUS_GeneDescription = mapIds(org.Mm.eg.db, MUS_EntrezID, "GENENAME", "ENTREZID"),
                                                            MUS_GO = mapIds(org.Mm.eg.db, MUS_EntrezID, "GO", "ENTREZID"),
                                                            MUS_Ontology = mapIds(org.Mm.eg.db, MUS_EntrezID, "ONTOLOGY", "ENTREZID"),
                                                            MUS_UniProt = mapIds(org.Mm.eg.db, MUS_EntrezID, "UNIPROT", "ENTREZID")
                                                          )
  #Combine mapped and KEGG-mapped uniprots
  combo_mapped_uniprots = bind_rows(EnzID_mapped_uniprots,KEGG_mapped_uniprots)
  combo_mapped_uniprots = combo_mapped_uniprots %>% unique()
  #Print final statement confirming completed execution of the above
  print("UniProt Nomenclature Conversion: Human --> Mouse Complete. KEGG Orthology DB used to address gaps in NCBI DB.")
  return(combo_mapped_uniprots)
  }
}

# orthDB_HtM_mapped_res_df = nomConvert_Human_Mouse(humUniProt)

#FUNCTION: Mouse --> Human Nomenclature Pivot
nomConvert_Mouse_Human <- function(uniprots){
  print("Starting mapping...")
  #Get EntrezIDs of submitted mouse uniprots (key for orthologyDB search)
  egs = mapIds(org.Mm.eg.db, uniprots, "ENTREZID","UNIPROT")
  #Check to make sure valid uniprots are submitted
  if(any(is.na(egs))){
    #stop(simpleError("Pivot failed: All uniprots submitted must be valid in either NCBI (EntrezID) or KEGG DB."))
    print("Not all uniprots submitted are valid in NCBI (EntrezID) DB. Proceeding with pivot for valid uniprots only.")
    egs = na.omit(egs)
  }
  print("Mouse uniprots mapped to EntrezID...")
  #Perform orthology selection with mouse EntrezIDs
  #mmug = keys(Orthology.eg.db,"Mus.musculus")
  mapped_uniprots = select(Orthology.eg.db, egs, "Homo.sapiens", "Mus.musculus")
  print("Mouse --> Human Entrez ID Orthology Match performed...")
  #Add columns for each annotation of interest (human and mouse)
  mapped_uniprots$MUS_Gene = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "SYMBOL", "ENTREZID"))
  mapped_uniprots$MUS_EnsemblID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ENSEMBL", "ENTREZID"))
  mapped_uniprots$MUS_EnsemblProtID = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ENSEMBLPROT", "ENTREZID"))
  mapped_uniprots$MUS_MGI = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "MGI", "ENTREZID"))
  mapped_uniprots$MUS_GeneDescription = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "GENENAME", "ENTREZID"))
  mapped_uniprots$MUS_GO = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "GO", "ENTREZID"))
  mapped_uniprots$MUS_Ontology = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "ONTOLOGY", "ENTREZID"))
  mapped_uniprots$MUS_UniProt = as.character(mapIds(org.Mm.eg.db, as.character(mapped_uniprots$Mus.musculus), "UNIPROT", "ENTREZID"))
  mapped_uniprots$HOM_Gene = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "SYMBOL", "ENTREZID"))
  mapped_uniprots$HOM_EnsemblID = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ENSEMBL", "ENTREZID"))
  mapped_uniprots$HOM_EnsemblProtID = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ENSEMBLPROT", "ENTREZID"))
  mapped_uniprots$HOM_GeneDescription = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "GENENAME", "ENTREZID"))
  mapped_uniprots$HOM_PATH = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "PATH", "ENTREZID"))
  mapped_uniprots$HOM_GO = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "GO", "ENTREZID"))
  mapped_uniprots$HOM_Ontology = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "ONTOLOGY", "ENTREZID"))
  mapped_uniprots$HOM_UniProt = as.character(mapIds(org.Hs.eg.db, as.character(mapped_uniprots$Homo.sapiens), "UNIPROT", "ENTREZID"))
  print("All annotation columns successfully added to NCBI-Mapped DF...")
  #Reorganize/rename columns for clarity
  mapped_uniprots = mapped_uniprots %>% dplyr::select(Mus.musculus,Homo.sapiens,everything())
  colnames(mapped_uniprots)[1] = "MUS_EntrezID"
  colnames(mapped_uniprots)[2] = "HOM_EntrezID"
  #Separate out uniprots with no human orthology found using entrez ID/orthology.eg.db
  EnzID_mapped_uniprots = mapped_uniprots %>% filter(!is.na(HOM_EntrezID)) %>% mutate(HOM_EntrezID = as.character(HOM_EntrezID))
  KEGG_mapped_uniprots = mapped_uniprots %>% filter(is.na(HOM_EntrezID)) %>% mutate(MUS_KEGGID = paste("mmu",MUS_EntrezID,sep = ":"))
  if(length(KEGG_mapped_uniprots$MUS_EntrezID) == 0){
      print("UniProt Nomenclature Conversion: Mouse --> Human Complete. NCBI DB used for all uniprots pivots.")
      return(EnzID_mapped_uniprots)
  }
  else{
      #Check for orthologs not detected by entrez search using KEGG data - split data into subsets of 100 for KEGG map request
      #Create keggLink urls to get KO#
      sub_size = 100
      total_length = nrow(KEGG_mapped_uniprots)
      r  = rep(1:ceiling(total_length/sub_size),each=sub_size)[1:total_length]
      data_frac = split(KEGG_mapped_uniprots,r)
      link_urls = vector("list", length = length(data_frac))
      link_urls = lapply(link_urls,function(x) paste(x,"https://rest.kegg.jp/link/ko/",sep = ""))
      j = 1 #Initialize data fraction counter
      for (sub_df in data_frac){
        i = 1 #Initialize KEGG ID counter - restart per sub_df/url
        current_link = link_urls[j]
        for (ezid in sub_df$MUS_KEGGID){
          if (i == 1){current_link = paste(current_link,ezid,sep = "")}
          else{current_link = paste(current_link,ezid,sep = "+")}
          i = i + 1
        }
        link_urls[j] = current_link
        j = j + 1
      }
      #Submit urls for KO# results
      MUS_KO_reslist = lapply(link_urls,function(x) read_tsv(x,col_names = c("MUS_KEGG_ID","KO_ID")))
      #Create keggLink urls to get mmu KEGGIDs
      hsa_link_urls = vector("list", length = length(data_frac))
      hsa_link_urls = lapply(hsa_link_urls,function(x) paste(x,"https://rest.kegg.jp/link/hsa/",sep = ""))
      j = 1 #Initialize data fraction counter
      for (ko_sub_df in MUS_KO_reslist){
        i = 1 #Initialize KEGG ID counter - restart per sub_df/url
        current_hsa_link = hsa_link_urls[j]
        for (koid in ko_sub_df$KO_ID){
          if (i == 1){current_hsa_link = paste(current_hsa_link,koid,sep = "")}
          else{current_hsa_link = paste(current_hsa_link,koid,sep = "+")}
          i = i + 1
        }
        hsa_link_urls[j] = current_hsa_link
        j = j + 1
      }
      #Submit hsa urls for hsa KEGGID# results
      HOM_KEGGID_res = lapply(hsa_link_urls,function(x) read_tsv(x,col_names = c("KO_ID","HOM_KEGG_ID")))
      #Turn all hsa KEGGID# to entrezIDs
      new_HOM_KEGGID_res = lapply(HOM_KEGGID_res, function(x) cbind(x, HOM_EntrezID = sub("hsa:","",x$HOM_KEGG_ID)))
      #Combine w/ mouse EntrezID
      HOM_res = do.call("rbind",new_HOM_KEGGID_res)
      MUS_res = do.call("rbind",MUS_KO_reslist)
      MUS_res = MUS_res %>% mutate(MUS_EntrezID = sub("mmu:","",MUS_KEGG_ID))
      MUS_HOM_res = inner_join(MUS_res,HOM_res, by = "KO_ID")
      MUS_HOM_EzID_res = MUS_HOM_res %>% dplyr::select(c("MUS_EntrezID","HOM_EntrezID"))
      #Combine with main df and fill in Orthology.eg.db HOM columns for filled in HOM EntrezIDs
      KEGG_mapped_uniprots = inner_join(KEGG_mapped_uniprots,MUS_HOM_EzID_res,by = "MUS_EntrezID")
      KEGG_mapped_uniprots = KEGG_mapped_uniprots %>% mutate(HOM_EntrezID.x = ifelse((is.na(HOM_EntrezID.x) & !is.na(HOM_EntrezID.y)),
                                                                                     HOM_EntrezID.y,
                                                                                     HOM_EntrezID.x)) %>%
                                                      dplyr::select(-c("HOM_EntrezID.y"))
      colnames(KEGG_mapped_uniprots)[colnames(KEGG_mapped_uniprots) == "HOM_EntrezID.x"] =  "HOM_EntrezID"
      KEGG_mapped_uniprots = KEGG_mapped_uniprots %>% mutate(HOM_Gene = mapIds(org.Hs.eg.db, HOM_EntrezID, "SYMBOL", "ENTREZID"),
                                                             HOM_EnsemblID = mapIds(org.Hs.eg.db, HOM_EntrezID, "ENSEMBL", "ENTREZID"),
                                                             HOM_EnsemblProtID = mapIds(org.Hs.eg.db, HOM_EntrezID, "ENSEMBLPROT", "ENTREZID"),
                                                             HOM_GeneDescription = mapIds(org.Hs.eg.db, HOM_EntrezID, "GENENAME", "ENTREZID"),
                                                             HOM_GO = mapIds(org.Hs.eg.db, HOM_EntrezID, "GO", "ENTREZID"),
                                                             HOM_Ontology = mapIds(org.Hs.eg.db, HOM_EntrezID, "ONTOLOGY", "ENTREZID"),
                                                             HOM_UniProt = mapIds(org.Hs.eg.db, HOM_EntrezID, "UNIPROT", "ENTREZID")
      )
      #Combine mapped and KEGG-mapped uniprots
      combo_mapped_uniprots = bind_rows(EnzID_mapped_uniprots,KEGG_mapped_uniprots)
      combo_mapped_uniprots = combo_mapped_uniprots %>% unique()
      #Print final statement confirming completed execution of the above
      print("UniProt Nomenclature Conversion: Mouse --> Human Complete. KEGG Orthology DB used to address gaps in NCBI DB.")
      return(combo_mapped_uniprots)
  }
}

# orthDB_MtH_mapped_res_df = nomConvert_Mouse_Human(musUniProt)
