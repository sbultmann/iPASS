library(shiny)
library(seqinr)
library(tidyverse)
library(cowplot)
library(DT)
library(Biostrings)
library(gtools)
library(fst)
library(data.table)

df <- fst::read_fst(path = "df_small.fst")

#df <- df[c(-2,-3,-4,-6)]

#fst::write.fst(df,"df_small.fst" )
dt <- data.table(df)
setkey(dt, sequence)

#setkey(dt, aa_seq)


#system.time(df%>% 
         #     filter(aa_seq=="RIGE"))


codonList <- as.list(GENETIC_CODE)

  #df <- data.table::as.data.table(df)

#data.table::setkey(dt, aa_seq)
#system.time(dt %>% 
          #    filter(aa_seq=="RIGE"))


check_fasta <- function(x){
  tryCatch({
    read.fasta(x, seqtype = "DNA")
    return(TRUE)
    #read.fasta(x, seqtype = "DNA")
  }, warning = function(war) {
    return(FALSE)
  }, error = function(err) {
    
    return(FALSE)
  }
  
  )}



parseInput <- function(x, InputType){
  if (InputType == "fasta"){
    Seq <- read.fasta(x, seqtype = "DNA")
    Seq <- toupper(as.vector(Seq[[1]]))
    Seq <- Seq[grepl("[ACGT]", Seq)]
  }
  if (InputType == "text"){
    Seq <-  toupper(s2c(x))
    Seq <- Seq[grepl("[ACGT]", Seq)]
  }
  return(Seq)
}

#Seq = "GCAGTATGCAGTATTTGCAGCA"
#Seq <- parseInput(Seq, "text")
CreateDataFrame <- function(Seq){
  sequence <- map_df(as.data.frame(embed(rev(Seq), 15)), rev)
  sequence <- sequence %>% 
    mutate(position = seq(1,nrow(sequence)),
           triplet = paste(V7,V8,V9, sep=""),
           aa = unlist(lapply(triplet, function(x) codonList[[x]]))
           ) %>% 
    dplyr::select(-V7,-V8,-V9) %>% 
    mutate(sequence = apply(sequence[1:15],1, function(x) paste0(paste(x[1:6], collapse = ""),paste(x[10:15], collapse = "")))) %>% 
    select(position, sequence, aa, triplet)
  df_seq <- sequence[seq(1,nrow(sequence), 3),]
  df_seq <-  df_seq %>% 
      left_join(., df, by="sequence") %>% 
    dplyr::rename(position="position.x") %>% 
    mutate(pos = position, 
           aa_position=seq(3, (nrow(df_seq)+2)),
           position = paste0((pos+6), "\n",triplet,"\n",aa )
           
           )
  
  
  df_subset<- df[df$aa_seq %in% c(df_seq$aa_seq),]
  df_seq <- df_subset %>%  
    #filter(aa_seq %in% c(df_seq$aa_seq), !grepl("*", aa_seq,fixed = T)) %>% 
    group_by(aa_seq) %>% 
    summarise(alt_sequence=sequence[which.max(score)], alt_score = score[which.max(score)]) %>% 
    left_join(.,df_seq, by="aa_seq") %>% 
    mutate(sequence = paste0(substr(sequence,1,6), "<font color='red'>TAG</font>", substr(sequence,7,12)),
           alt_sequence = paste0(substr(alt_sequence,1,6), "<font color='red'>TAG</font>", substr(alt_sequence,7,12))) %>% 
    select(score, position, aa_position, pos, sequence, alt_sequence, alt_score)
  return(df_seq)
}

CreateDataFrameDT <- function(Seq){
  sequence <- map_df(as.data.frame(embed(rev(Seq), 15)), rev)
  sequence <- sequence %>% 
    mutate(position = seq(1,nrow(sequence)),
           triplet = paste(V7,V8,V9, sep=""),
           aa = unlist(lapply(triplet, function(x) codonList[[x]]))
    ) %>% 
    dplyr::select(-V7,-V8,-V9) %>% 
    mutate(sequence = apply(sequence[1:15],1, function(x) paste0(paste(x[1:6], collapse = ""),paste(x[10:15], collapse = "")))) %>% 
    select(position, sequence, aa, triplet)
  df_seq <- sequence[seq(1,nrow(sequence), 3),]
  df_seq <- data.table(df_seq)
  setkey(df_seq, sequence)
  
  
  df_seq <-  dt[df_seq] %>% 
    #dplyr::rename(position="i.position") %>% 
    mutate(pos = position, 
           aa_position=seq(3, (nrow(df_seq)+2)),
           position = paste0((pos+6), "\n",triplet,"\n",aa )
           
    )
  
  
  df_subset<- dt[aa_seq %in% c(df_seq$aa_seq),
                 .(alt_sequence=sequence[which.max(score)], alt_score = score[which.max(score)])
                 ,by=.(aa_seq)]
  df_seq <- df_subset %>%  
    left_join(.,df_seq, by="aa_seq") %>% 
    mutate(sequence = paste0(substr(sequence,1,6), "<font color='red'>TAG</font>", substr(sequence,7,12)),
           alt_sequence = paste0(substr(alt_sequence,1,6), "<font color='red'>TAG</font>", substr(alt_sequence,7,12))) %>% 
    select(score, position, aa_position, pos, sequence, alt_sequence, alt_score)
  return(df_seq)
}

#CreateDataFrameDT(Seq)
#microbenchmark(
 # DT = CreateDataFrameDT(Seq),
#  DF = CreateDataFrame(Seq),
#  times = 100
#)


iPASS <-function(x, InputType){
  Seq <- parseInput(x, InputType)
  result <- CreateDataFrameDT(Seq)
  return(result)
}
