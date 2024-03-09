# RR_visualization
#SeanKinard
#2023-07-01
#-----------------------------------------------------------------------------
#Setup
#-----------------------------------------------------------------------------
source('exploration/toolkit.R')

RR_biomass_com<-read_csv("exploration/output/RR_biomass_com.csv")
RR_density_com<-read_csv("exploration/output/RR_density_com.csv")
RR_shan_com<-read_csv("exploration/output/RR_shan_com.csv")
RR_rich_com<-read_csv("exploration/output/RR_rich_com.csv")
RR_cdis_com<-read_csv("exploration/output/RR_cdis_com.csv")
RR_density_fam<-read_csv("exploration/output/RR_density_fam.csv")
RR_density_spe<-read_csv("exploration/output/RR_density_spe.csv")
#-----------------------------------------------------------------------------
#Visualization:ResponseVsBaseline(Community)
#-----------------------------------------------------------------------------
p_RR_RvB<-function(my_data,my_variable){
  
  colnames(my_data)<-str_replace_all(
    colnames(my_data),my_variable,'X')
  
  categorical_colors<-c('cyan','red','grey90')
  
  my_data%>%
    fix_site_order()%>%
    mutate(upper_lim=X_b_qtr_fill+X_sd20,
           lower_lim=X_b_qtr_fill-X_sd20,
           Comparison=case_when(
             X>upper_lim~'Above',
             X<lower_lim~'Below',
             X<=upper_lim&X>=lower_lim~'Inside'))%>%
    ggplot()+
    facet_wrap(~site_code,scales='free_y')+
    geom_point(aes(x=collection_period,y=X,fill=Comparison),
               shape=21,size=3,color='black',alpha=.3)+
    geom_point(aes(x=collection_period,y=X,color=Comparison),
               shape=21,size=3)+
    geom_line(aes(x=collection_period,y=upper_lim),
              linetype=2,linewidth=.4,color='grey40')+
    geom_line(aes(x=collection_period,y=lower_lim),
              linetype=2,linewidth=.4,color='grey40')+
    dark_theme_gray(base_size=12)+
    scale_fill_manual(values=categorical_colors)+
    scale_color_manual(values=categorical_colors)+
    scale_x_date(date_labels="%m")+
    ylab(my_variable)+
    xlab('Month')+
    theme(legend.position='none')
}

RR_fig_RvB_com_biomass<-p_RR_RvB(RR_biomass_com,'AFDM')

RR_fig_RvB_com_density<-p_RR_RvB(RR_density_com,'dens')

RR_fig_RvB_com_shan<-p_RR_RvB(RR_shan_com,'shan')

RR_fig_RvB_com_rich<-p_RR_RvB(RR_rich_com,'rich')

RR_fig_RvB_com_cdis<-p_RR_RvB(RR_cdis_com,'cdis')

#-----------------------------------------------------------------------------
#Visualization:ResponseVsBaseline(Family)
#-----------------------------------------------------------------------------
common_families<-c("Leuciscidae","Poeciliidae","Ictaluridae","Centrarchidae")

RR_fig_RvB_leuc_density<-p_RR_RvB(
  my_data=RR_density_fam%>%filter(family=='Leuciscidae'),
  my_variable='dens')

RR_fig_RvB_poec_density<-p_RR_RvB(
  my_data=RR_density_fam%>%filter(family=='Poeciliidae'),
  my_variable='dens')

RR_fig_RvB_icta_density<-p_RR_RvB(
  my_data=RR_density_fam%>%filter(family=='Ictaluridae'),
  my_variable='dens')

RR_fig_RvB_cent_density<-p_RR_RvB(
  my_data=RR_density_fam%>%filter(family=='Centrarchidae'),
  my_variable='dens')

#-----------------------------------------------------------------------------
#Visualization:variableVsBaseline(Species)
#-----------------------------------------------------------------------------
RR_fig_RvB_Gaff_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="G. affinis"),
  my_variable='dens')

RR_fig_RvB_Plat_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="P. latipinna"),
  my_variable='dens')

RR_fig_RvB_Lcya_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. cyanellus"),
  my_variable='dens')

RR_fig_RR_fig_RvB_Lmac_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. macrochirus"),
  my_variable='dens')

RR_fig_RR_fig_RvB_Lmeg_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. megalotis"),
  my_variable='dens')

RR_fig_RR_fig_RvB_Laur_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. auritus"),
  my_variable='dens')

RR_fig_RvB_Lgul_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. gulosus"),
  my_variable='dens')

RR_fig_RvB_Msal_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="M. salmoides"),
  my_variable='dens')

RR_fig_RvB_Locu_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. oculatus"),
  my_variable='dens')

RR_fig_RvB_Anat_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="A. natalis"),
  my_variable='dens')

RR_fig_RvB_Ngyr_density<-p_RR_RvB(
  my_data=RR_density_spe%>%filter(lowest_taxon=="N. gyrinus"),
  my_variable='dens')
 
#-----------------------------------------------------------------------------
#Visualization:ResponseRatio
#-----------------------------------------------------------------------------

p_RR<-function(my_data,my_variable){
  
  colnames(my_data)<-str_replace_all(
    colnames(my_data),my_variable,'X')
  
  categorical_colors<-c('cyan','red','grey90')
  
  my_data%>%
    fix_site_order()%>%
    mutate(upper_lim=X_b_qtr_fill+X_sd20,
           lower_lim=X_b_qtr_fill-X_sd20,
           Comparison=case_when(
             X>upper_lim~'Above',
             X<lower_lim~'Below',
             X<=upper_lim&X>=lower_lim~'Inside'))%>%
    ggplot()+
    facet_wrap(~site_code,scale='fixed')+
    geom_smooth(aes(x=collection_period,y=X_RR_b_qtr),
                method='loess',span=.8,
                se=F,linetype=2,color='grey45',linewidth=.5)+
    geom_point(aes(x=collection_period,y=X_RR_b_qtr,fill=Comparison),
               shape=21,size=3,color='black',alpha=.3)+
    geom_point(aes(x=collection_period,y=X_RR_b_qtr,color=Comparison),
               shape=21,size=3)+
    dark_theme_gray(base_size=12)+
    scale_fill_manual(values=categorical_colors)+
    scale_color_manual(values=categorical_colors)+
    scale_x_date(date_labels="%m")+
    ylab(my_variable)+
    xlab('Time')+
    theme(legend.position='none')
}

RR_fig_com_biomass<-p_RR(my_data=RR_biomass_com,
                         my_variable='AFDM')

RR_fig_com_density<-p_RR(my_data=RR_density_com,
                         my_variable='dens')

RR_fig_com_rich<-p_RR(my_data=RR_rich_com,
                      my_variable='rich')

RR_fig_com_shan<-p_RR(my_data=RR_shan_com,
                      my_variable='shan')

RR_fig_com_cdist<-p_RR(my_data=RR_cdis_com,
                       my_variable='cdis')

RR_fig_leuc<-p_RR(
  my_data=RR_density_fam%>%filter(family=='Leuciscidae'),
  my_variable='dens')+
  facet_wrap(~site_code,scale='free_y')

RR_fig_poec<-p_RR(
  my_data=RR_density_fam%>%filter(family=='Poeciliidae'),
  my_variable='dens')+
  facet_wrap(~site_code,scale='free_y')

RR_fig_icta<-p_RR(
  my_data=RR_density_fam%>%filter(family=='Ictaluridae'),
  my_variable='dens')+
  facet_wrap(~site_code,scale='free_y')

RR_fig_cent<-p_RR(
  my_data=RR_density_fam%>%filter(family=='Centrarchidae'),
  my_variable='dens')+
  facet_wrap(~site_code,scale='free_y')

RR_fig_Gaff<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="G. affinis"),
  my_variable='dens')

RR_fig_Plat<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="P. latipinna"),
  my_variable='dens')

RR_fig_Lcya<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. cyanellus"),
  my_variable='dens')

RR_fig_Lmac<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. macrochirus"),
  my_variable='dens')

RR_fig_Lmeg<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. megalotis"),
  my_variable='dens')+
  facet_wrap(~site_code,scale='free_y')

RR_fig_Laur<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. auritus"),
  my_variable='dens')+
  ggtitle('L. auritus Density Response Ratio')

RR_fig_Lgul<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. gulosus"),
  my_variable='dens')

RR_fig_Msal<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="M. salmoides"),
  my_variable='dens')

RR_fig_Locu<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="L. oculatus"),
  my_variable='dens')

RR_fig_Anat<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="A. natalis"),
  my_variable='dens')

RR_fig_Ngyr<-p_RR(
  my_data=RR_density_spe%>%filter(lowest_taxon=="N. gyrinus"),
  my_variable='dens')

#------------------------------------------------------------------------------
# Export figures
#------------------------------------------------------------------------------
# my_objects<-ls()
# my_figure_names <- my_objects[str_detect(my_objects,'RR_fig')]
# 
# my_figures <- list(
#   RR_fig_cent, RR_fig_Gaff, RR_fig_icta,
#   RR_fig_Lcya, RR_fig_leuc, RR_fig_Lmac,
#   RR_fig_poec, RR_fig_Anat, RR_fig_com_cdist,
#   RR_fig_com_density, RR_fig_com_rich, RR_fig_Laur,
#   RR_fig_Lgul, RR_fig_Lmeg, RR_fig_Locu,
#   RR_fig_Msal, RR_fig_Ngyr, RR_fig_RR_fig_RvB_Laur_density,
#   RR_fig_RR_fig_RvB_Lmac_density, RR_fig_RR_fig_RvB_Lmeg_density, 
#   RR_fig_RvB_Anat_density,RR_fig_RvB_cent_density, 
#   RR_fig_RvB_com_biomass, RR_fig_RvB_com_cdis,
#   RR_fig_RvB_com_density, RR_fig_RvB_com_rich, 
#   RR_fig_RvB_com_shan, RR_fig_RvB_Gaff_density, 
#   RR_fig_RvB_icta_density, RR_fig_RvB_Lcya_density, 
#   RR_fig_RvB_leuc_density, RR_fig_RvB_Lgul_density, 
#   RR_fig_RvB_Locu_density, RR_fig_RvB_Msal_density, 
#   RR_fig_RvB_Ngyr_density, RR_fig_RvB_poec_density )
# 
# names(my_figures) <- my_figure_names
# 
# for (i in 1:length(my_figures)) {
#   my_place <- paste('exploration/visualization/', names(my_figures[i]), ".png", sep='')
#   my_object <- my_figures[[i]]
#   ggsave(my_place,
#          plot = my_object,
#          width = 10,
#          height = 10,
#          units = c("in")) }

#------------------------------------------------------------------------------
# End RR_visualization