library(plyr); 
library(magrittr); 
library(tidyverse)
library(ggrepel)

prior_counts <- 0

#this is an early version of the table that will be published with the CCLE phase2 MS
M_MSI <- read_csv('~/CPDS/WRN_manuscript/data/Supplementary_Table_7_MSI_status.csv')

# #load GDSC MSI calls
# sanger_df <- read.csv('~/CPDS/data/GDSC/Cell_Lines_Details.csv') %>% 
#   dplyr::rename(sanger_msi = Microsatellite..instability.Status..MSI.) %>% 
#   mutate(CCLE_ID = CleanCellLineName(Sample.Name),
#          sanger_msi = as.vector(sanger_msi),
#          sanger_msi = ifelse(sanger_msi == '', NA, sanger_msi))

#takes all columns matching a given name pattern and creates normalized versions to correct for data-type specific intercept and slope terms
normalize_columns <- function(df, name) {
  col_set <- grep(name, colnames(df))
  x <- as.matrix(df[, col_set])
  cur_avg <- rowMeans(x, na.rm=T)
  mod <- lm(x ~ cur_avg)
  xp <- scale(x, center=mod$coefficients[1,], scale = mod$coefficients[2,])
  avg <- data.frame(x = rowMeans(xp, na.rm=T))
  colnames(avg) <- paste0(name, '_all')
  xp %<>% cbind(avg)

  colnames(xp) <- paste0(colnames(xp), '_norm')
  return(cbind(df, xp))
}

#normalize columns across data types
M_MSI %<>% normalize_columns('total_del')
M_MSI %<>% normalize_columns('msi_del')

M_MSI %<>% mutate(
  tot_dels = rowMeans(M_MSI[, c('CCLE.wes.total_del', 'GDSC.wes.total_del', 'CCLE.wgs.total_del', 'CCLE.hc.total_del')], na.rm=T),
  msi_dels = rowMeans(M_MSI[, c('CCLE.wes.msi_del', 'GDSC.wes.msi_del', 'CCLE.wgs.msi_del', 'CCLE.hc.msi_del')], na.rm=T),
  tot_dels_norm = rowMeans(M_MSI[, c('CCLE.wes.total_del_norm', 'GDSC.wes.total_del_norm', 'CCLE.wgs.total_del_norm', 'CCLE.hc.total_del_norm')], na.rm=T),
  msi_dels_norm = rowMeans(M_MSI[, c('CCLE.wes.msi_del_norm', 'GDSC.wes.msi_del_norm', 'CCLE.wgs.msi_del_norm', 'CCLE.hc.msi_del_norm')], na.rm=T)
) 
M_MSI %<>% mutate(
  msi_del_frac = pmax(0, (msi_dels + prior_counts) / (tot_dels + 2*prior_counts)),
  msi_del_frac_norm = pmax(0, (msi_dels_norm + prior_counts) / (tot_dels_norm + 2*prior_counts))
 )

# M_MSI %<>% left_join(sanger_df %>% dplyr::select(CCLE_ID, sanger_msi), by = 'CCLE_ID')
write.csv(M_MSI, '~/CPDS/WRN_manuscript/data/CCLE_MSI_info.csv', row.names = FALSE)
