# 000_update_fish
# runs fish R scripts in order
# Sean Kinard
# last update: 2023-07-01

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('exploration/toolkit.R') # load packages and helper-functions

my_scripts <- list.files("exploration/calculation")
my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=T)]
my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
my_scripts <- my_scripts[str_detect(my_scripts,"site_map", negate=T)]
my_scripts <- sort(my_scripts, decreasing=F)

#------------------------------------------------------------------------------
# Source scripts
#------------------------------------------------------------------------------
for(i in 1:length(my_scripts)) {
  source(paste("exploration/calculation/", my_scripts[i], sep=''))
  rm(list=ls())
  my_scripts <- list.files("exploration/calculation")
  my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=T)]
  my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
  my_scripts <- my_scripts[str_detect(my_scripts,"site_map", negate=T)]
  my_scripts <- sort(my_scripts, decreasing=F)   }

#------------------------------------------------------------------------------
# End 000_update_fish