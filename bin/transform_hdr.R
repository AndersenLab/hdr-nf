#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(valr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

hdr_file     <- args[1]
coords_file  <- args[2]
group_id     <- args[3]
output_file  <- args[4]

all_calls_SR_clustered_NR <- readr::read_tsv(
  hdr_file,
  show_col_types = FALSE
)

strain <- dplyr::case_match(
  group_id,
  "AD"        ~ "ECA2670",
  "KD"        ~ "JU1348",
  "Temperate" ~ "JU2536",
  "TH"        ~ "NIC1660",
  "TD2"       ~ "BRC20492",
  "TD1"       ~ "BRC20530",
  .default    = NA_character_
)

if (is.na(strain)) {
  stop("Unknown group_id: ", group_id)
}

if (!all(all_calls_SR_clustered_NR$group == group_id)) {
  stop(
    "HDR file contains group(s) inconsistent with group_id = ",
    group_id
  )
}

nuc <- readr::read_tsv(
  coords_file,
  col_names = c(
    "GROUP", "S1", "E1", "S2", "E2", "L1", "L2",
    "IDY", "LENR", "LENQ", "REF", "HIFI"
  ),
  show_col_types = FALSE
) %>%
  dplyr::filter(GROUP == group_id) %>%
  dplyr::filter(!grepl("ptg", HIFI)) %>%
  dplyr::filter(REF == HIFI) %>%
  dplyr::filter(L2 > 3e3) %>%
  dplyr::mutate(STRAIN = strain) %>%
  dplyr::select(-GROUP) %>%
  dplyr::group_by(STRAIN, HIFI) %>%
  dplyr::arrange(STRAIN, HIFI, S2) %>%
  dplyr::mutate(
    leadS1 = lead(S1),
    leadE1 = lead(E1),
    leadS2 = lead(S2),
    leadE2 = lead(E2),
    lagS1  = lag(S1),
    lagE1  = lag(E1),
    lagS2  = lag(S2),
    lagE2  = lag(E2)
  ) %>%
  dplyr::ungroup()

if (nrow(nuc) == 0) {
  stop(
    "No NUCmer alignments remain for group ", group_id,
    " after filtering"
  )
}

if (nrow(all_calls_SR_clustered_NR) == 0) {
  stop("HDR file contains no calls for group ", group_id)
}

calls <- as.data.table(all_calls_SR_clustered_NR)[
  , .(
    rowid = .I,
    region_source = STRAIN,
    genome = strain,
    contig = CHROM,
    start = as.integer(start),
    end = as.integer(end)
  )
]

nuc_dt <- as.data.table(nuc)[, idx := .I]
nuc_dt <- nuc_dt[, .(idx,
                     genome = STRAIN,
                     contig = HIFI,
                     start = pmin(S2, E2),
                     end = pmax(S2, E2),
                     orig_start=S2,
                     orig_end=E2,
                     L1=L1,
                     L2=L2,
                     leadS1=leadS1,
                     leadE1=leadE1,
                     leadS2=leadS2,
                     leadE2=leadE2,
                     lagS1=lagS1,
                     lagE1=lagE1,
                     lagS2=lagS2,
                     lagE2=lagE2,
                     refchrom=REF,
                     refstart=S1,
                     refend=E1
)]

nuc_dt[, `:=`(start = as.integer(start), end = as.integer(end))]

# Set proper keys
setkey(calls, genome, contig, start, end)
setkey(nuc_dt, genome, contig, start, end)

# Run overlap
matched <- foverlaps(
  x = calls,
  y = nuc_dt,
  type = "any",
  nomatch = 0
)

matched_df <- as.data.frame(matched) %>%
  dplyr::mutate(hdr_chrom=contig) %>%
  dplyr::mutate(INV=ifelse(orig_start==start,F,T)) %>%
  dplyr::select(refstart,
                refend,
                orig_start,
                orig_end,
                L1,
                L2,
                refchrom,
                contig,
                genome,
                INV,
                start,
                end,
                rowid,
                region_source,
                hdr_chrom,
                i.start,
                i.end,
                leadS1,
                leadE1,
                leadS2,
                leadE2,
                lagS1,
                lagE1,
                lagS2,
                lagE2) %>%
  dplyr::rename(S1=refstart,
                E1=refend,
                S2=orig_start,
                E2=orig_end,
                REF=refchrom,
                HIFI=contig,
                HIFI_strain=genome,
                St2=start,
                Et2=end,
                group_id=rowid,
                hdr_strain=region_source,
                hdr_start=i.start,
                hdr_end=i.end) %>% 
  dplyr::arrange(group_id,S2) %>%
  dplyr::mutate(HDRid = paste0(hdr_strain,hdr_chrom,hdr_start,hdr_end))


tigFilt2 <- matched_df %>%
  dplyr::arrange(group_id,S1) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(leadDiff=lead(S1)-E1) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 5E4,1,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(E1-S1)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group_id,S1) 


trim_spacer = 1e3
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- tigFilt2 %>%
  dplyr::group_by(group_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(scale_distortion = ((L2 - L1)/L1)) %>%
  dplyr::mutate(rboundDist=hdr_start-min(S2,E2)) %>%
  #dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==F,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(S1=ifelse(rboundDist>trim_spacer & INV==F ,(S1+(rboundDist-trim_spacer+(rboundDist*scale_distortion))),S1)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==T ,(E1-(rboundDist-trim_spacer+(rboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(S2=ifelse(rboundDist>trim_spacer & INV==F,(S2+(rboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & INV==T,(E2+(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(lboundDist=max(S2,E2)-hdr_end) %>%
  dplyr::mutate(E1=ifelse(lboundDist>trim_spacer & INV==F ,(E1-(lboundDist-trim_spacer+(lboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer & INV==T ,(S1+(lboundDist-trim_spacer+(lboundDist*scale_distortion))),S1)) %>%
  dplyr::mutate(E2=ifelse(lboundDist>trim_spacer & INV==F,(E2-(lboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & INV==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::ungroup()

tigMarkExtend <- tigTrim %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(E2_extend=ifelse(INV==F & E2 == max(E2) & E2 < hdr_end, T, F),
                S2_extend=ifelse(INV==F & S2 == min(S2) & S2 > hdr_start, T,F),
                iE2_extend=ifelse(INV==T & E2 == min(E2) & E2 > hdr_start,T,F),
                iS2_extend=ifelse(INV==T & S2 == max(S2) & S2 < hdr_end, T,F)) %>%
  dplyr::mutate(any_extend=ifelse(E2_extend == T | S2_extend == T | iE2_extend==T | iS2_extend ==T,T,F)) %>%
  dplyr::ungroup()
 

tigToExtend <- tigMarkExtend %>% 
  dplyr::filter(any_extend==T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(extend_length_WI_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-E2,
                                             ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-S2,NA)),
                extend_length_WI_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, S2-max(lagS2,lagE2),
                                               ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, E2-max(lagS2,lagE2),NA)),
                extend_length_REF_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, 
                                              ifelse(leadS1 >= E1,leadS1-E1,ifelse(leadE1>=E1,leadE1-E1,NA)),
                                              ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, ifelse(S1>=leadE1,S1-leadE1,ifelse(S1>=leadS1,S1-leadS1,NA)),NA)),
                extend_length_REF_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, 
                                            ifelse(lagE1<=S1,S1-lagE1,ifelse(lagS1<=S1,S1-lagS1,NA)),
                                            ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, ifelse(lagS1>=E1,lagS1-E1,ifelse(lagE1>=E1,lagE1-E1,NA)),NA))) %>%
  dplyr::ungroup()


extendDat <- rbind(tigToExtend %>% 
                     dplyr::select(extend_length_REF_lead,extend_length_WI_lead) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lead,extend_length_REF=extend_length_REF_lead),
                   tigToExtend %>% 
                     dplyr::select(extend_length_REF_lag,extend_length_WI_lag) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lag,extend_length_REF=extend_length_REF_lag)) %>%
  dplyr::filter(!is.na(extend_length_WI) & !is.na(extend_length_REF))

tigExtensions <- rbind(tigToExtend,tigMarkExtend %>% dplyr::filter(any_extend==F) %>% dplyr::mutate(extend_length_WI_lead=NA,extend_length_REF_lead=NA,extend_length_WI_lag=NA,extend_length_REF_lag=NA))

tigExtended_50kb <- tigExtensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(E1=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(leadS1 >= E1,leadS1,ifelse(leadE1>=E1,leadE1,E1)),E1)) %>%
  dplyr::mutate(E1=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagS1>=E1,lagS1,ifelse(lagE1>=E1,lagE1,E1)),E1)) %>%
  dplyr::mutate(S1=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagE1<=S1,lagE1, ifelse(lagS1<=S1,lagS1,S1)),S1)) %>%
  dplyr::mutate(S1=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(S1>=leadE1,leadE1,ifelse(S1>=leadS1,leadS1,S1)),S1)) %>%
  dplyr::mutate(E2=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),E2)) %>%
  dplyr::mutate(E2=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),E2)) %>%
  dplyr::mutate(S2=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),S2)) %>%
  dplyr::mutate(S2=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),S2)) %>% 
  dplyr::ungroup()

counts50kb <- tigExtended_50kb %>%
  dplyr::filter(
    any_extend == TRUE,
    (!is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag)) | (!is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead)),
    (extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4) | (extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4)
  ) %>%
  dplyr::group_by(hdr_strain, HDRid) %>%
  dplyr::summarise(has_extension = any(any_extend), .groups = "drop") %>%  # optional, since filter already ensures this
  dplyr::group_by(hdr_strain) %>%
  dplyr::summarise(count_true = n(), .groups = "drop") %>%
  dplyr::mutate(window_size = "50kb")

hdr_counts <- tigTrim %>%
  dplyr::group_by(hdr_strain) %>%
  dplyr::summarise(num_unique_HDRid = n_distinct(HDRid), .groups = "drop")

hdr_transformed_50ext <- tigExtended_50kb %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    newS1 = min(c(S1, E1), na.rm = TRUE),
    newE1 = max(c(S1, E1), na.rm = TRUE),
    across(
      .cols = -c(S1, E1, S2, E2, St2, Et2),
      .fns = dplyr::first
    ),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    S1 = newS1,
    E1 = newE1
  )

gap_clust_NR_TR <- hdr_transformed_50ext %>%
  dplyr::select(REF,S1,E1,hdr_strain,HIFI_strain) %>%
  dplyr::rename(CHROM=REF,minStart=S1,maxEnd=E1,STRAIN=hdr_strain,REF=HIFI_strain) %>%
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::arrange(STRAIN,CHROM,minStart) %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()
 
joinClust_NR_TR <- gap_clust_NR_TR %>% 
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
  dplyr::group_by(conID) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,nclust,REF)

nojoin_NR_TR <- gap_clust_NR_TR %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,REF) %>%
  dplyr::mutate(nclust=1)

all_calls_SR_clustered_NR_TR <- rbind(joinClust_NR_TR,nojoin_NR_TR) %>%
  dplyr::filter(divSize/1e3 >= 5) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(REF) %>%
  dplyr::arrange(REF,desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("ptg",CHROM))

hdr_trans <- all_calls_SR_clustered_NR_TR %>%
  dplyr::mutate(group=group_id) %>%
  dplyr::select(CHROM,start=minStart,end=maxEnd,group,STRAIN,REF) %>%
  dplyr::rename(source=REF) %>%
  dplyr::mutate(start=round(start),end=round(end),size=end-start) %>%
  dplyr::relocate(source, .after = dplyr::last_col())

write.table(hdr_trans, output_file, row.names = F,quote = F,sep = '\t')

