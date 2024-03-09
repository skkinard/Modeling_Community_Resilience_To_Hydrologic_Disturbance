# 06_environment_PCA
# Sean Kinard
# 2023-06-23

# use PCA identify patterns in environmental variables among sample sites
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R')
library(factoextra) # extract prcomp info

d_l <- read_csv("exploration/output/environment_longterm_with_ln.csv")
d_s <- read_csv("exploration/output/environment_shortterm_with_ln.csv") %>%
  select(-wk4_max_30mu)

# fill missing environmental predictors with linear interpolation
d_s <- d_s %>%
  rowwise() %>%
  as_tsibble(key=site_code,
             index=collection_period) %>%  
  na_ma(weighting="linear") %>%
  ungroup() %>%
  as_tibble() %>%
  pivot_longer(cols=-any_of(c('site_code', 'collection_period')),
               names_to='xname',
               values_to = 'x') %>%
  group_by(site_code, collection_period, xname) %>%
  dplyr::summarize(x=mean(x ,na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from = x) %>%
  ungroup()

d_s20 <- d_s %>%
  filter(collection_period >= ymd('2020-01-01')) %>%
  select(-collection_period)

d_s19 <- d_s %>%
  filter(collection_period >= ymd('2019-01-01') &
           collection_period < ymd('2020-01-01')) %>%
  select(-collection_period)

d_s18 <- d_s %>%
  filter(collection_period >= ymd('2018-01-01') &
           collection_period < ymd('2019-01-01')) %>%
  select(-collection_period)

d_s17 <- d_s %>% # outlier flow data at WMC 2017-09-30
  filter(collection_period < ymd('2018-01-01') &
           collection_period >= ymd('2017-01-01')) %>%
  select(-collection_period)

d_s <- d_s %>%
  select(-collection_period)

#------------------------------------------------------------------------------
# PCA functions
#------------------------------------------------------------------------------
my_pca <- function(my_data) {
  
  d_x <- my_data
  
  d_x <- d_x %>%
    fix_site_order() %>%
    arrange(site_code)
  
  # PCA and plots
  PCA <- prcomp(d_x[,-1], scale = TRUE)
  PCA_importance <- summary(PCA)$importance %>%
    as.data.frame() %>%
    rownames_to_column(var = 'stat') %>%
    select(stat, PC1:PC5)
  
  # pearson's r: correlations with PC1 and PC2
  PCA_table_cor <-   cor(PCA$x[,1:2], 
                         select(d_x, -site_code), 
                         method='pearson') %>% 
    t() %>%
    as.data.frame %>%
    rownames_to_column('predictor') %>%
    as_tibble() %>%
    arrange(desc(abs(PC1))) %>%
    rename(PC1_r = PC1, PC2_r=PC2)
  
  # diagnostic scree plot
  screeplot <- plot(PCA, type = "l")
  
  # Extract PCA information
  # Eigenvalues
  eig.val <- get_eigenvalue(PCA)
  # Results for Variables
  res.var <- get_pca_var(PCA)
  #res.var$coord          # Coordinates
  #res.var$contrib        # Contributions to the PCs
  #res.var$cos2           # Quality of representation 
  # Results for individuals
  res.ind <- get_pca_ind(PCA)
  #res.ind$coord          # Coordinates
  #res.ind$contrib        # Contributions to the PCs
  #res.ind$cos2           # Quality of representation   
  
  # plot-prep axis labels
  axis_labels <- PCA_importance %>%
    pivot_longer(cols=contains('PC'), names_to='PC', values_to='value') %>%
    pivot_wider(names_from=stat, values_from=value) %>%
    r_friendly_colnames() %>%
    mutate(axis_label = paste(pc, ' (', 
                              round(proportion_of_variance*100, 0),
                              "%)", sep = ''))
  
  # plot-prep variables
  pre_var <- res.var$coord %>%
    as.data.frame() %>%
    rownames_to_column(var = 'predictor') %>%
    as_tibble() %>%
    r_friendly_colnames() %>%
    select(predictor, dim_1, dim_2)
  pre_var$k=as.factor(cutree(hclust(dist(pre_var)),k=6))
  
  # plot-prep individuals
  pre_ind <- res.ind$coord %>% 
    as_tibble() %>%
    r_friendly_colnames() %>%
    select(dim_1, dim_2) %>%
    mutate(site_code=d_x$site_code)
  
  sum_ind <- pre_ind %>%
    group_by(site_code) %>%
    dplyr::summarize(dim_1=mean(dim_1),
                     dim_2=mean(dim_2))
  
  # PCA plot
  pca_site_plot <- ggplot() +
    geom_point(data = pre_ind,
               aes(x=dim_1, y=dim_2, color = site_code), size = 3,
               show.legend = F) +
    geom_polygon(data = pre_ind,
                 aes(x=dim_1, y=dim_2, color = site_code), 
                 fill=NA, alpha = .1, show.legend = F) +
    geom_label(data = pre_ind,
               aes(x=dim_1, y=dim_2, label = site_code, color = site_code),
               show.legend = F) +
    xlab(axis_labels$axis_label[1]) +
    ylab(axis_labels$axis_label[2]) +
    dark_theme_grey() +
    xlim(c(-8,4)) +
    ylim(c(-5,5)) +
    scale_color_manual(values=my_colors)
  
  pca_predictor_plot <- ggplot() +
    geom_segment(data=pre_var,
                 aes(x=0, y=0, xend=dim_1*5, yend=dim_2*5),
                 colour="grey70", linewidth=.2) +
    geom_label_repel(data=pre_var,
                     aes(label=predictor,x=dim_1*5, y=dim_2*5, color = k),
                     show.legend=F,
                     fill = 'grey10') +
    xlab(axis_labels$axis_label[1]) +
    ylab(element_blank()) +
    dark_theme_grey() +
    xlim(c(-8,4)) +
    ylim(c(-5,5)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    paletteer::scale_color_paletteer_d("ggthemes::excel_Vapor_Trail")
  
  pca_combo_plot <- pca_site_plot + pca_predictor_plot
  
  my_output <- list(pca_combo_plot, PCA_table_cor, 
                    pca_site_plot, pca_predictor_plot)
  
  return(my_output) }

remove_na <- function(df) {
  df %>% 
    rowid_to_column("ID") %>%
    pivot_longer(
    cols=-c('site_code', 'ID'), 
    names_to='predictor', 
    values_to = 'x') %>%
    filter(!is.na(x)) %>%
    pivot_wider(names_from=predictor, values_from=x) %>%
    select(-ID) %>%
    na.omit() }
#------------------------------------------------------------------------------
# PCA
#------------------------------------------------------------------------------
lte_pca <- my_pca(d_l)
ste_pca_2017 <- my_pca(d_s17)
ste_pca_2018 <- my_pca(d_s18 %>% remove_na())
ste_pca_2019 <- my_pca(d_s19 %>% remove_na())
ste_pca_2020 <- my_pca(d_s20 %>% remove_na())
ste_pca_all <- my_pca(d_s %>% remove_na())

#------------------------------------------------------------------------------
# Figures
#------------------------------------------------------------------------------
lte_pca_comboplot <- lte_pca[[3]] + ggtitle('Long-Term Averages') + lte_pca[[4]]
ste_pca_2017_comboplot <- ste_pca_2017[[3]] + ggtitle('2017')  + ste_pca_2017[[4]]
ste_pca_2018_comboplot <- ste_pca_2018[[3]] + ggtitle('2018')  + ste_pca_2018[[4]]
ste_pca_2019_comboplot <- ste_pca_2019[[3]] + ggtitle('2019')  + ste_pca_2019[[4]]
ste_pca_2020_comboplot <- ste_pca_2020[[3]] +  ggtitle('2020')  + ste_pca_2020[[4]]
ste_pca_all_comboplot <- ste_pca_all[[3]] + ggtitle('2017-2020') + ste_pca_all[[4]]

caption_pca_comboplot <- "Principal Component Analysis of short-term environmental variables for sampling locations. Variables . Sites are displayed in the left panel with hot colors indicating more arid climate. The righ panel shows scaled, environmental predictors colored according to heirarchical clusterings (k=6) informed by PC1 and PC2. Landuse categorizations were estimated using satellite imagery in 2006. Flow metrics are denoted with 'q_' and were taken from 30-years of daily discharge data. Substrate and water quality variables are average values from field collections from spring of 2017 through December of 2020."

# Export Figures
combo_plots <- list(lte_pca_comboplot,
                    ste_pca_2017_comboplot,
                    ste_pca_2018_comboplot,
                    ste_pca_2019_comboplot,
                    ste_pca_2020_comboplot,
                    ste_pca_all_comboplot )

names(combo_plots) <- c('lte_pca_comboplot',
                        'ste_pca_2017_comboplot',
                        'ste_pca_2018_comboplot',
                        'ste_pca_2019_comboplot',
                        'ste_pca_2020_comboplot',
                        'ste_pca_all_comboplot')

for (i in 1:length(combo_plots)) {
  my_place <- paste('exploration/visualization/', names(combo_plots[i]), ".png", sep='')
  my_object <- combo_plots[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 12,
         height = 6,
         units = c("in")) }

#------------------------------------------------------------------------------
# Tables
#------------------------------------------------------------------------------
lte_pca_correlation <- lte_pca[[2]]
ste_pca_2017_correlation <- ste_pca_2017[[2]]
ste_pca_2018_correlation <- ste_pca_2018[[2]]
ste_pca_2019_correlation <- ste_pca_2019[[2]]
ste_pca_2020_correlation <- ste_pca_2020[[2]]
ste_pca_all_correlation <- ste_pca_all[[2]]

caption_pca_correlation <- "Table containing correlations between the first two axes from principal component analysis and short-term environmental variables."

# Export Tables
correlation_tables <- list(lte_pca_correlation,
                           ste_pca_2017_correlation,
                           ste_pca_2018_correlation,
                           ste_pca_2019_correlation,
                           ste_pca_2020_correlation,
                           ste_pca_all_correlation )

names(correlation_tables) <- c('lte_pca_correlation',
                        'ste_pca_2017_correlation',
                        'ste_pca_2018_correlation',
                        'ste_pca_2019_correlation',
                        'ste_pca_2020_correlation',
                        'ste_pca_all_correlation')

for (i in 1:length(correlation_tables)) {
  my_place <- paste('exploration/output/', names(correlation_tables[i]), ".csv", sep='')
  my_object <- correlation_tables[[i]]
  write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# End 06_environment_PCA

