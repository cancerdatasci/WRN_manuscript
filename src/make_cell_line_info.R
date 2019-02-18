library(plyr)
library(magrittr) 
library(tidyverse)
library(reshape2)
library(taigr)
library(ggrepel)
library(matrixStats)
library(cowplot)

source('~/CPDS/WRN_manuscript/src/WRN_helpers.R')

data_dir <- '~/CPDS/WRN_manuscript/data'

#PARAMETERS 
methyl_thresh <- Inf #threshold on methylation (Inf ignores methylation)
CN_loss_thresh <- -1 #threshold on log2 relative CN values for calling deletions
GE_unexpressed_thresh <- 1 #threshold on log2TPM values for calling unexpressed genes
avg_WRN_thresh <- -0.5 #dependency threshold (on the avg of CRISPR and RNAi dep scores)
RPPA_thresh <- -1 #threshold for defining protein loss from RPPA

common_MSI_lineages <- c('gastric', 'colorectal', 'endometrium', 'ovary')
MMR_genes <- c('MLH1', 'MSH2', 'MSH6', 'PMS2')

dat <- load_omics_data()
# write_rds(dat, file.path(data_dir, 'DepMap_18Q4_data.rds')) #WRITE OUT OMICS DAT FILE

sample_info <- load.from.taiga(data.name = 'master-cell-line-export-0306', data.version = '196') %>%
  dplyr::select(CCLE_ID = CCLE_name, lineage = Disease, secondary_tissue = `Disease Subtype`) %>%
  mutate(lineage = ifelse(lineage == '', NA, lineage),
         lineage = ifelse(secondary_tissue == 'uterus_endometrium', 'endometrium', lineage)) %>% 
  distinct(CCLE_ID, .keep_all=T)

#load Andrew Giacomelli's TP53 annotations (from supplement of Nat Gen TP53 paper)
#Get binary call on TP53 proficient vs non-functional
TP53_status <- read_csv('~/CPDS/data/Giacomelli/Giacomelli_TP53_status.csv') %>% 
  mutate(TP53_fnx = ifelse(grepl('Non-functional', TP53_status), 'TP53_null', 'TP53_proficient'))

#MMR loss classifications
#load MSI calls from CCLE phase 2 supplement
MSI_data <- read_csv(file.path(data_dir, 'CCLE_MSI_info.csv')) %>% 
  mutate(MSI_lab = revalue(CCLE.MSI.call, c(`inferred-MSI` = 'MSI', `inferred-MSS` = 'MSS', `undetermined` = 'indeterminate')),
         MSI_lab = factor(MSI_lab, levels = c('MSS', 'indeterminate', 'MSI')),
         MSI = ifelse(CCLE.MSI.call != 'undetermined', CCLE.MSI.call == 'inferred-MSI', NA)) %>% 
  mutate(GDSC_MSI = factor(GDSC.msi.call, levels = c('MSS/MSI-L', 'MSI-H')))

#CLASSIFY MMR LOSS
#Call damaging mutations in each MMR gene, and an OR operation (does any MMR gene have a dam mut)
MMR_dam_mut <- dat$MUT_DAM[, MMR_genes] > 0
colnames(MMR_dam_mut) %<>% paste0('_DAM_MUT')
MMR_dam_mut %<>% cbind(data.frame(CCLE_ID = rownames(dat$MUT_DAM[, MMR_genes]), 
                                  MMR_DAM_MUT = rowMaxs(dat$MUT_DAM[, MMR_genes][, MMR_genes], na.rm=T) > 0))

#Call deletions in MMR genes
MMR_DEL <- dat$CN[, MMR_genes] < CN_loss_thresh
colnames(MMR_DEL) %<>% paste0('_DEL')
MMR_DEL %<>% cbind(data.frame(CCLE_ID = rownames(dat$CN), 
                              MMR_DEL = rowMins(dat$CN[, MMR_genes], na.rm=T) < CN_loss_thresh))

#Call non-expressed MMR genes
MMR_GE <- dat$GE[, MMR_genes] < GE_unexpressed_thresh
colnames(MMR_GE) %<>% paste0('_UNEXP')
MMR_GE %<>% cbind(data.frame(CCLE_ID = rownames(dat$GE), 
                             MMR_UNEXP = rowMins(dat$GE[, MMR_genes], na.rm=T) < GE_unexpressed_thresh))

#merge MMR gene info into MSI data
MSI_data %<>% 
  left_join(MMR_dam_mut) %>% 
  left_join(MMR_DEL) %>% 
  left_join(MMR_GE) 

#MMR loss type annotation
MSI_data %<>% mutate(MMR_loss_type = 'none',
                     MMR_loss_type = ifelse(MMR_DEL, 'del', MMR_loss_type),
                     MMR_loss_type = ifelse(MMR_UNEXP, 'unexp', MMR_loss_type),
                     MMR_loss_type = ifelse(MMR_DAM_MUT, 'dam_mut', MMR_loss_type),
                     MMR_loss = MMR_DAM_MUT | MMR_DEL | MMR_UNEXP)

#per gene loss function
MSI_data %<>% mutate(
  MLH1_loss =  MLH1_DAM_MUT | MLH1_DEL | MLH1_UNEXP,
  MSH2_loss = MSH2_DAM_MUT | MSH2_DEL | MSH2_UNEXP,
  MSH6_loss = MSH6_DAM_MUT | MSH6_DEL | MSH6_UNEXP,
  PMS2_loss = PMS2_DAM_MUT | PMS2_DEL | PMS2_UNEXP,
  MMR_loss_gene = 'none',
  MMR_loss_gene = ifelse(MLH1_loss, 'MLH1', MMR_loss_gene),
  MMR_loss_gene = ifelse(MSH2_loss, 'MSH2', MMR_loss_gene),
  MMR_loss_gene = ifelse(MSH6_loss, 'MSH6', MMR_loss_gene),
  MMR_loss_gene = ifelse(PMS2_loss, 'PMS2', MMR_loss_gene),
  MMR_loss_gene = ifelse(MLH1_loss+MSH2_loss+MSH6_loss+PMS2_loss > 1, 'multiple', MMR_loss_gene)
)

MSI_data %<>% mutate(
  MMR_loss_lab = factor(ifelse(MMR_loss, 'MMR loss', 'MMR WT'), levels = c('MMR WT', 'MMR loss'))
)


#### ADD ADDITIONAL FEATURES RELEVANT FOR MSI THAT DONT GO INTO MMR LOSS CALLING######

#find cell lines with any (non-silent) MMR mutation
any_MMR_mut <- dat$MUT %>%
  filter(Hugo_Symbol %in% MMR_genes, Variant_Classification != 'Silent') %>% .[['CCLE_ID']] %>% unique()

#POLE mutants of each class
POLE_damaging_mut_set <- rownames(dat$MUT_DAM)[dat$MUT_DAM[, 'POLE'] > 0] 
POLE_other_mut_set <- rownames(dat$MUT_OTHER)[dat$MUT_OTHER[, 'POLE'] > 0] 
POLE_hotspot_mut_set <- rownames(dat$MUT_HOT)[dat$MUT_HOT[, 'POLE'] > 0] 

# Get MSH2/6 protein levels and call low-protein status
MSH2_RPPA_low <- rownames(dat$RPPA)[dat$RPPA[, 'MSH2'] < RPPA_thresh]
MSH6_RPPA_low <- rownames(dat$RPPA)[dat$RPPA[, 'MSH6_Caution'] < RPPA_thresh]
MSI_data %<>% mutate(POLE_damaging_mut = CCLE_ID %in% POLE_damaging_mut_set,
                     POLE_other_mut = CCLE_ID %in% POLE_other_mut_set,
                     POLE_hotspot_mut = CCLE_ID %in% POLE_hotspot_mut_set,
                     MSH2_RPPA_low = CCLE_ID %in% MSH2_RPPA_low,
                     MSH6_RPPA_low = CCLE_ID %in% MSH6_RPPA_low,
                     other_MMR_mut = CCLE_ID %in% any_MMR_mut
) %>% 
  mutate(POLE_damaging_mut = ifelse(CCLE_ID %in% rownames(dat$MUT_DAM), POLE_damaging_mut, NA),
         POLE_other_mut = ifelse(CCLE_ID %in% rownames(dat$MUT_OTHER), POLE_other_mut, NA),
         POLE_hotspot_mut = ifelse(CCLE_ID %in% rownames(dat$MUT_HOT), POLE_hotspot_mut, NA),
         MSH2_RPPA_low = ifelse(CCLE_ID %in% rownames(dat$RPPA), MSH2_RPPA_low, NA),
         MSH6_RPPA_low = ifelse(CCLE_ID %in% rownames(dat$RPPA), MSH6_RPPA_low, NA),
         other_MMR_mut = ifelse(CCLE_ID %in% unique(dat$MUT$CCLE_ID), other_MMR_mut, NA))

# Get MSI/MSS Differential Dependency Stats

#combine diff-dependency data with MSI data and WRN dependency
comb_data <- MSI_data %>% 
  full_join(data.frame(CCLE_ID = rownames(dat$DRIVE),
                       DRIVE_WRN_D2 = dat$DRIVE[, 'WRN']), by = 'CCLE_ID') %>% 
  full_join(data.frame(CCLE_ID = rownames(dat$CRISPR),
                       CRISPR_WRN_CERES = dat$CRISPR[, 'WRN']), by = 'CCLE_ID') %>% 
  left_join(sample_info %>% dplyr::select(lineage, CCLE_ID), by = 'CCLE_ID') %>% 
  left_join(data.frame(CCLE_ID = rownames(dat$GE), WRN_GE = dat$GE[, 'WRN']))

#compute an average WRN dependency score
comb_data$avg_WRN <- rowMeans(comb_data[, c('DRIVE_WRN_D2', 'CRISPR_WRN_CERES')], na.rm=T)
comb_data %<>% mutate(lineage = str_replace_all(lineage, '_', ' '),
                      common_MSI_lin = lineage %in% common_MSI_lineages)

comb_data %<>% left_join(TP53_status, by = 'CCLE_ID')
comb_data %<>% mutate(WRN_dep = avg_WRN < avg_WRN_thresh)


# write_rds(comb_data, '~/CPDS/WRN/data/WRN_comb_data.rds')

#Write out tables of cell line info for supplement
comb_data %>%
  dplyr::select(CCLE_ID,
                Disease = lineage,
                GDSC_MSI,
                CCLE_MSI = MSI_lab,
                DRIVE_WRN_D2,
                CRISPR_WRN_CERES,
                avg_WRN_dep = avg_WRN,
                is_WRN_dep = WRN_dep,
                TP53_status = TP53_fnx,
                common_MSI_lineage = common_MSI_lin,
                MMR_loss = MMR_loss,
                ms_deletions_normed = msi_dels_norm,
                frac_deletions_in_ms_regions = msi_del_frac_norm,
                MLH1_damaging_mut = MLH1_DAM_MUT,
                MLH1_unexpressed = MLH1_UNEXP,
                MLH1_deletion = MLH1_DEL,
                MLH1_loss,
                MMR_loss_gene,
                MSH2_damaging_mut = MSH2_DAM_MUT,
                MSH2_unexpressed = MSH2_UNEXP,
                MSH2_deletion = MSH2_DEL,
                MSH2_loss,
                MSH6_damaging_mut = MSH6_DAM_MUT,
                MSH6_unexpressed = MSH6_UNEXP,
                MSH6_deletion = MSH6_DEL,
                MSH6_loss,
                PMS2_damaging_mut = PMS2_DAM_MUT,
                PMS2_unexpressed = PMS2_UNEXP,
                PMS2_deletion = PMS2_DEL,
                PMS2_loss,
                MSH2_RPPA_low,
                MSH6_RPPA_low,
                POLE_damaging_mut,
                POLE_other_mut,
                POLE_hotspot_mut,
                any_MMR_mut
  ) %>% write_csv(file.path(data_dir, 'WRN_final_cell_line_table.csv'))
