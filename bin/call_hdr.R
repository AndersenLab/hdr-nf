#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(valr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# This function intersects g2g long-read alignments with genome bins to extract bin-level coverage and identity 
# a miriad of bin metrics are calculated (see mutate cluster)
alignmentPartition <- function(df) {
  a2a_bed <- df %>% dplyr::select(REF,S1,E1,IDY,STRAIN) %>%
    dplyr::rename(chrom=REF,start=S1,end=E1,Identity=IDY)
  
  df_LRbins <- valr::bed_intersect(bins_1kb_CB_stripped,a2a_bed) %>%
    dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                  start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y) %>%
    dplyr::group_by(CHROM,START_BIN,END_BIN)  %>% 
    dplyr::filter(coverage==max(coverage)) %>% 
    dplyr::filter(Identity==max(Identity)) %>% 
    dplyr::mutate(group_id=cur_group_id()) %>%
    dplyr::distinct(group_id,.keep_all = T) %>%
    dplyr::ungroup()
  
  df_LRbins_unclass <- bins_1kb_CB_stripped %>% 
    dplyr::rename(CHROM=chrom, START_BIN=start, END_BIN=end) %>%
    dplyr::left_join(df_LRbins,by=c("CHROM"="CHROM","START_BIN"="START_BIN","END_BIN"="END_BIN")) %>%
    dplyr::arrange(CHROM,START_BIN) %>% 
    dplyr::mutate(coverage=ifelse(is.na(coverage),0,coverage)) %>%
    dplyr::mutate(STRAIN.y=ifelse(is.na(STRAIN.y),unique(a2a_bed$STRAIN),STRAIN.y)) %>%
    dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
    dplyr::rename(bin_IDY=Identity,bin_COV=coverage,STRAIN=STRAIN.y) %>%
    dplyr::select(CHROM,START_BIN,END_BIN,bin_IDY,bin_COV,STRAIN) 
  
  return(df_LRbins_unclass)
}

# given a set of coverage and identity thresholds, this function classifies bins as divergent or non-divergent
# divergent bins have sub-classifications based on the variant determinants 
classifyPartition <- function(df,cf,idy) {
  df_bins_clasif_LR <- df %>%
    dplyr::mutate(div=ifelse(bin_COV<(cf*10),"C",
                             ifelse(!(is.na(bin_IDY)) & bin_IDY<=idy,"I","nondiv"))) %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(div_class = ifelse((lag(div) == "C" | lag(div) == "I") & (lead(div) == "C" | lead(div) == "I") & div=="nondiv","G", div)) %>%
    dplyr::mutate(div_class = ifelse(div_class=="C" & bin_COV == 0,"Z",div_class)) %>%
    dplyr::mutate(div_gf=ifelse(div_class=="nondiv","nondiv","div")) %>%
    dplyr::ungroup()
  
  return(df_bins_clasif_LR)
}

## this function estimates the frequency at which each bin is classified as hyperdivergent in the sampled population
## low frequency bins can be re discarded to avoid super-clustering of divergent regions
frequencyEstimates <- function(df) {
  denom<-length(unique(df$STRAIN))
  df2 <- df %>% 
    dplyr::group_by(CHROM,START_BIN) %>%
    dplyr::mutate(binCt=stringr::str_count(div_gf,"nondiv")) %>%
    dplyr::mutate(binCt=ifelse(is.na(binCt),1,binCt)) %>%
    dplyr::mutate(binCt=ifelse(binCt==1,0,1)) %>%
    dplyr::mutate(binCtSum=sum(binCt)) %>%
    dplyr::mutate(binFreq=(binCtSum/denom)*100)
  return(df2)
}

## This function clusters contiguous div bins into div regions, 
## and assembles a divergent 'footprint' from the string of bin classifications
clusterBins <- function(df,mode) {
  
  if(mode=="SRF") {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) %>%
      dplyr::mutate(div_gf=ifelse(div_gf =="div" & binFreq < 5,"nondiv",div_gf)) %>%
      dplyr::select(-binCt,-binCtSum,-binFreq)
  } else {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) 
  }
  
  temp$enum <- sequence(rle(as.character(temp$div_gf))$lengths)
  
  div_bins <- temp %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM,data.table::rleid(div_gf)) %>%
    dplyr::mutate(gid=cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(clusTips=ifelse(div_gf == "div" & lead(enum) == 2 & lead(div_gf)=="div","clust_start", 
                                  ifelse(div_gf == "div" & lead(enum) > 1 & lead(div_gf)=="div","clust_center",
                                         ifelse(div_gf == "div" & enum > 1 & lead(div_gf)=="nondiv","clust_end","unclust")))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!(clusTips=="unclust"))
  
  if (mode == "LR") {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(prop_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(prop_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(prop_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(prop_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      dplyr::mutate(
        meanIDY = if (all(is.na(bin_IDY))) NA_real_ else
          sum(bin_IDY * bin_COV, na.rm = TRUE) / sum(bin_COV[!is.na(bin_IDY)])) %>%
      dplyr::mutate(meanCov=mean(bin_COV),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart,group_size=n()) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanIDY,meanCov,prop_covz,prop_lowcov,prop_gf,prop_idy,bin_foot,group_size)
  } else {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(prop_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(prop_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(prop_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(prop_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      dplyr::mutate(meanVC=mean(COUNT),meanCF=mean(pc1X),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart,group_size=n()) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,prop_covz,prop_lowcov,prop_gf,prop_idy,bin_foot,group_size)
  }
  
  out <- list(div_bins,div_regions)
  return(out)
}


coverage_file <- args[[1]]
varct_file    <- args[[2]]
windows_file  <- args[[3]]
REF           <- args[[4]]
GROUP         <- args[[5]]
COV_thresh    <- as.numeric(args[[6]])
VC_thresh     <- as.numeric(args[[7]])
RPO_thresh    <- as.numeric(args[[8]])

bins_1kb_CB_stripped <- readr::read_tsv(
  windows_file,
  col_names = FALSE
) %>%
  dplyr::rename(
    chrom = X1,
    start = X2,
    end = X3
  )

coverage_df <- readr::read_table(
  coverage_file,
  col_names = c(
    "GROUP",
    "STRAIN",
    "CHROM",
    "START_BIN",
    "END_BIN",
    "NAME",
    "c1X",
    "c2X",
    "c5X"
  )
) %>% 
dplyr::filter(GROUP == .env$GROUP) %>% 
dplyr::select(-GROUP)

varct_df <- readr::read_table(
  varct_file,
  col_names = c(
    "GROUP",
    "STRAIN",
    "CHROM",
    "START_BIN",
    "END_BIN",
    "COUNT"
  )
) %>% 
dplyr::filter(GROUP == .env$GROUP) %>% 
dplyr::select(-GROUP)

strainL <- varct_df %>% dplyr::distinct(STRAIN) %>% dplyr::pull(STRAIN)

#merge bin-based variant counts and coverage data
SR_stats_WI <- varct_df %>%
  dplyr::left_join(coverage_df,by=c("STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(STRAIN_count=sum(COUNT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(COUNT_ADJ = COUNT / STRAIN_count) 

#classify bins
all_stats <- SR_stats_WI %>%
  dplyr::mutate(div_class=ifelse(COUNT >= VC_thresh,"I", #over variant count thresh
                                 ifelse(pc1X < COV_thresh,"C","R"))) %>%  #under coverage thresh
  dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>% #fill 1kb gaps
  dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>% #assign a classification to gaps
  dplyr::ungroup()

#cluster bins (no frequency filter)
regList <- list()
binList <- list()
for (i in 1:length(strainL)) {
  #print((i/length(strainL))*100)
  temp <- all_stats %>% dplyr::filter(STRAIN==strainL[i])
  div_call <- clusterBins(temp,"SR")
  #div_call[[1]]$STRAIN <- strainL[[i]]
  div_call[[2]]$STRAIN <- strainL[i]
  #binList[[i]] <- div_call[[1]]
  regList[[i]] <- div_call[[2]]
}

all_calls_SR <- ldply(regList, data.frame)

#QX calls are spurious (likely repeats) and should be removed from other strain calls
ref_calls <- all_calls_SR %>% dplyr::filter(STRAIN==REF) %>% dplyr::select(CHROM,minStart,maxEnd)

#find overlaps and remove
dt_all <- data.table::as.data.table(
  all_calls_SR[all_calls_SR$STRAIN != REF, ]
) # all other strain calls

dt_qx <- data.table::as.data.table(ref_calls) # REF calls

data.table::setnames(
  dt_qx,
  c("CHROM", "minStart", "maxEnd"),
  c("CHROM", "start", "end")
) # rename keys

data.table::setnames(
  dt_all,
  c("minStart", "maxEnd"),
  c("start", "end")
) # rename keys

data.table::setkey(dt_qx, CHROM, start, end) # set the keys for matching
data.table::setkey(dt_all, CHROM, start, end) # set the keys for matching

overlaps <- data.table::foverlaps(
  dt_all,
  dt_qx,
  nomatch = 0L,
  type = "any"
) # find any overlap

overlaps_pct <- data.table::copy(overlaps)

overlaps_pct[, overlap_pct :=
               100 * (pmin(end, i.end) - pmax(start, i.start)) /
               (i.end - i.start)
]

to_drop <- unique(
  overlaps_pct[
    overlap_pct > RPO_thresh,
    .(CHROM, start = i.start, end = i.end, STRAIN)
  ]
) # grab all other strain calls that overlap

key_cols <- c("CHROM", "start", "end", "STRAIN")

dt_filtered <- dt_all[
  !to_drop,
  on = key_cols
] # remove them

data.table::setnames(
  dt_filtered,
  c("start", "end"),
  c("minStart", "maxEnd")
)

all_calls_SR_noREF <- as.data.frame(dt_filtered) # convert back to df

#flag adjacent regions that are 5kb apart
gap_clust <- all_calls_SR_noREF %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

#get flagged and merge them
joinClust <- gap_clust %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gapFoot=paste0(rep("R",forGapSize/1e3),collapse = "")) %>%
  dplyr::mutate(new_foot=ifelse(flag3g=="clust",paste0(bin_foot,gapFoot),bin_foot)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(clust_foot=paste0(new_foot,collapse = "")) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(newMeanVC=mean(meanVC)) %>%
  dplyr::mutate(newMeanCF=mean(meanCF)) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize,-meanVC,-meanCF,-bin_foot) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize,meanVC=newMeanVC,meanCF=newMeanCF,bin_foot=clust_foot) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN,nclust)

#keep unflagged 
nojoin <- gap_clust %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN) %>%
  dplyr::mutate(nclust=1)

#bind unflagged and merged calls, arrange strains by number of HDRs
all_calls_SR_clustered <- rbind(joinClust,nojoin) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 

#filter by size, rename columns for print
all_calls_SR_clustered_sfilt <- all_calls_SR_clustered %>%
  dplyr::filter(divSize >= 5e3) %>%
  dplyr::select(-nclust,-ncalls,-sorter,-rleID,-ystrain) %>%
  dplyr::mutate(group=GROUP) %>%
  dplyr::select(CHROM,start=minStart,end=maxEnd,group,STRAIN,size=divSize,meanVC,meanCF,bin_footprint=bin_foot)

#write
write.table(all_calls_SR_clustered_sfilt, "hdrs.tsv",row.names = F,quote = F,sep = '\t')