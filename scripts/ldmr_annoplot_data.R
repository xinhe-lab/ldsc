anno_fun <- function(feature,gr,lev){
  tff <- feature
  t_feat_name <- dplyr::filter(all_feat_df,feature==tff) %>% dplyr::pull(category)
  feat_name=paste0("ra_",make.names(feature))
  plyranges::mutate(readd(feat_name,character_only=TRUE,cache=cc),feature=t_feat_name) %>% 
    plyranges::join_overlap_inner(x = .,gr) %>% plyranges::set_genome_info(.data = .,genome="hg19")
}

grp_anno_fun <- function(df){
  feat_fact <- factor(df$feature)
  pmap(df,anno_fun,lev=levels(feat_fact)) %>% set_names(make.names(df$feature)) %>% 
     plyranges::bind_ranges() %>% mutate(feature_fact=as.integer(factor(feature)))
}

best_fdr_df$stop[2] <- best_fdr_df$stop[2]+20000
#best_fdr_df <- mutate(best_fdr_df,stop=if_else(region_id==356,))
ld_dfrl <- best_fdr_df %>% 
  transmute(seqnames=paste0("chr",chrom),start=start,width=stop-start,region_id=region_id,p_rank=p_rank) %>%
  rowwise() %>% 
  do(.,tibble(
    gr=list(plyranges::set_genome_info(plyranges::as_granges(.,keep_mcols=FALSE),"hg19")),
    region_id=.$region_id,p_rank=.$p_rank)) %>% ungroup()


bait_hic <- transmute(hic,seqnames=bait_chr,start=bait_start,width=bait_end-start,N_reads=N_reads,map_id=map_id,target_bait="bait") %>%
  plyranges::as_granges()
target_hic <- transmute(hic,seqnames=otherEnd_chr,start=otherEnd_start,width=otherEnd_end-start,N_reads=N_reads,map_id=map_id,target_bait="target") %>%
  plyranges::as_granges()


midpoint_df <- function(gr){
  tb <- gr$target_bait[1]
  midpoint_name <- paste0("midpoint_",tb)
  if(tb=="target")
  return(gr %>% as_tibble() %>% 
    dplyr::transmute(seqnames=seqnames,midpoint_target=start+width/2,map_id=map_id,credible_set_target=credible_set,reads=N_reads))
  return(gr %>% as_tibble() %>% 
    dplyr::transmute(seqnames=seqnames,midpoint_bait=start+width/2,map_id=map_id,credible_set_bait=credible_set))
}


join_gr_left <- function(gr,cred_df){
    plyranges::join_overlap_left(gr,cred_df) %>% 
    mutate(credible_set=if_else(is.na(credible_set),FALSE,credible_set))
}


hic_l <- map2(ld_dfrl$gr,susie_cred,function(gr,sgr){
  n_bait <- plyranges::filter_by_overlaps(bait_hic,gr)
  
  n_target <- plyranges::filter_by_overlaps(target_hic,gr) 
  hic_gr <- plyranges::bind_ranges(n_bait,n_target) %>% 
    mutate(group="Hi-C") %>% join_gr_left(sgr)
  
  
  hic_df <- midpoint_df(join_gr_left(n_bait,sgr)) %>% 
    inner_join(midpoint_df(join_gr_left(n_target,sgr))) %>% mutate(any_credible_set=credible_set_target|credible_set_bait,end_height=1) %>% dplyr::filter(midpoint_bait!=midpoint_target)
  
  with_cred <- filter(hic_df,any_credible_set)
  no_cred <-  filter(hic_df,!any_credible_set)
  
  gpr <-   ggplot() +
    geom_rect(data = hic_gr,aes(group=group,alpha=N_reads,col=credible_set),rect.height=0.1,group.selfish=FALSE)
  if(nrow(no_cred)>0){
    gpr <- gpr + geom_curve(
      data = no_cred,
      aes(x=midpoint_bait,y=end_height,xend=midpoint_target,yend=end_height,alpha=0.4),
      curvature=0.5,
      arrow=arrow(length=unit(0.025,"npc"))
      )
  }
  if(nrow(with_cred)>0){
     gpr <- gpr + geom_curve(
      data = with_cred,
      aes(x=midpoint_bait,y=end_height,xend=midpoint_target,yend=end_height),
      curvature=0.9,col="red",
      arrow=arrow(length=unit(0.025,"npc"))
      )
  }
  gpr <- gpr+xlab("")+theme(legend.position = "none",
                            axis.text.y=element_blank(),
                            axis.title.y = element_text(angle=90),
                            axis.ticks.y=element_blank()) + 
    scale_color_manual(values=c("FALSE"="black","TRUE"="red")) + 
    ylim(0,2)+
    ylab("Hi-C\nInteractions")
  return(gpr)
})



anno_rl <- unnest(best_model_df,cols=c(features)) %>% 
  dplyr::select(feature=features) %>% dplyr::filter(str_detect(feature,"hic_all",negate=TRUE)) %>% 
  dplyr::mutate(ina=NA_integer_) %>%
  dplyr::inner_join(mutate(ld_dfrl,ina=NA_integer_)) %>% 
  dplyr::select(-ina) %>% 
  tidyr::nest(data=c(gr,feature)) %>% 
  mutate(anno_l=map(data, ~grp_anno_fun(.x)))

anno_l <- map2(anno_rl$anno_l,susie_cred,function(x,y){
  nx <- join_gr_left(x,y) 
  ggplot(nx) + 
    geom_rect(aes(group=feature,col=credible_set,fill=credible_set)) + 
    scale_color_manual(values=c("FALSE"="black","TRUE"="red")) +
    scale_fill_manual(values=c("FALSE"="black","TRUE"="red")) + 
    theme(legend.position = "none")
})
