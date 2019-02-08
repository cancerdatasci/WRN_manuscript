
# Data Loading ------------------------------------------------------------
load_omics_data <- function() {
  dat_info <- list(
    DRIVE = list(
      data.name = 'demeter2-drive-0591',
      data.file = 'gene_means_proc',
      data.version = 11,
      transpose = T
    ),
    CRISPR = list(data.name='avana-public-tentative-18q4-3406',
                  data.version=4,
                  data.file='gene_effect',
                  transpose = F),
    GE = list(
      data.name="depmap-rnaseq-expression-data-ccd0",
      data.version = 9,
      data.file = 'CCLE_depMap_18Q4_TPM_ProteinCoding',
      transpose = F),
    CN = list(data.name = 'depmap-wes-cn-data-97cc',
              data.file = 'public_18Q4_gene_cn',
              data.version = 11),
    MUT_HOT = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.file = 'hotspot_mutation',
      data.version = 5),
    MUT_OTHER = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.file = 'other_mutation',
      data.version = 5),
    MUT_DAM = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.file = 'damaging_mutation',
      data.version = 5),
    MUT = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.file = 'depmap_18Q4_mutation_calls',
      data.version = 5),
    RPPA = list(data.name='depmap-rppa-1b43', 
                data.version=2, 
                data.file='CCLE_RPPA_20180123') 
  )
  
  dat <- taigr::load.all.from.taiga(dat_info) %>% 
    cdsr::map_arxspan_to_ccle() %>% 
    cdsr::extract_hugo_symbol_colnames()
  
  #convert arxspan IDs to CCLE_ID for mut dataset
  dat$MUT %<>% mutate(CCLE_ID = celllinemapr::arxspan.to.ccle(DepMap_ID))
  return(dat)
}


# STATS -------------------------------------------------------------------

#' Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param vec: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#' @param limma_trend: Whether to fit an intensity trend with the empirical Bayes variance model
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' CRISPR = load.from.taiga(data.name='avana-2-0-1-d98f',
#' data.version=1,
#' data.file='ceres_gene_effects',
#' transpose = T)
#' is_panc <- load.from.taiga(data.name = 'ccle-lines-lineages') %>% .[, 'pancreas']
#' ulines <- intersect(rownames(CRISPR), names(is_panc))
#' lim_res <- run_lm_stats_limma(CRISPR[ulines,], is_panc[ulines])
#' @export run_lm_stats_limma
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene', limma_trend = FALSE) {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = p.adjust(p.left, method = 'BH'),
                             q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}


#' Title
#'
#' @param targ_gene Name of dependency gene
#' @param bmarker Binary vector of biomarker calls, with CCLE ID names
#' @param bmark_name Name for biomarker
#' @param use_RNAi Whether or not to use both RNAi and CRISPR data for dep binarization
#' @param dep_thresh Dependency score threshold for calling dependency (applied to the avg of CRISPR and RNAi scores if using both)
#'
#' @return
#' @export
#'
#' @examples
make_biomarker_plot <- function(targ_gene, 
                                bmarker,
                                bmark_name,
                                use_RNAi = TRUE,
                                dep_thresh = -0.5) {
  
  #merge in relevant biomarker data
  df <- data.frame(CCLE_ID = rownames(dat$DRIVE),
                   dep_R = dat$DRIVE[, targ_gene]) %>% 
    full_join(data.frame(CCLE_ID = rownames(dat$CRISPR),
                         dep_C = dat$CRISPR[, targ_gene])) %>% 
    left_join(comb_data, by = 'CCLE_ID')
  
  # if using RNAi and CRISPR compute an average dependency score
  if (use_RNAi) {
    df$avg_dep <- rowMeans(df[, c('dep_R', 'dep_C')], na.rm=T)
  } else {
    df$avg_dep <- df$dep_C
  }
  df %<>% left_join(bmarker, by = 'CCLE_ID') #merge in biomarker calls
  
  print(with(df %>% dplyr::filter(!is.na(avg_dep)), table(bmark)))
  
  #compute ppv and sensitivity
  ppv <- with(df %>% filter(bmark), mean(avg_dep < dep_thresh, na.rm=T))
  sensitivity <- with(df %>% filter(avg_dep < dep_thresh), mean(bmark, na.rm=T))
  
  #make plot
  g <- ggplot(df %>% filter(!is.na(bmark)), aes(bmark, avg_dep)) + 
    geom_violin() +
    ggbeeswarm::geom_beeswarm(alpha = 0.5, size = 1, cex = 0.75) +
    ggtitle(targ_gene) +
    xlab(bmark_name) + 
    ylab(paste0(targ_gene, ' dependency')) +
    geom_hline(yintercept = dep_thresh, linetype = 'dashed') +
    ggtitle(sprintf('PPV: %.3f, sens: %.3f', ppv, sensitivity)) +
    cdsr::theme_Publication() +
    theme(axis.text.x = element_text(angle =70, hjust = 1))
  
  res = list(targ_gene = targ_gene,
             bmark = bmark_name,
             ppv = ppv,
             sensitivity = sensitivity)
  return(list(res = res, plot = g))
}


# Plotting ----------------------------------------------------------------
make_volcano <- function(df, n_labs = 15) {
  ggplot(df, aes(EffectSize, -log10(p.value))) + 
    geom_point(alpha = 0.6, size = 1) + 
    geom_point(data = filter(df, Gene == 'WRN'), size = 3, color = 'red', alpha = 1) + 
    geom_label_repel(data = df %>%
                       filter(q.value < q_thresh) %>%
                       arrange(desc(abs(EffectSize))) %>%
                       head(n_labs),
                     aes(label = Gene), 
                     size = 2, 
                     label.size = 0.1, 
                     fontface = 'bold',
                     label.padding = 0.1) +    
    xlab('MSI-MSS mean difference') +
    ylab('-log10(q-value)') +
    theme_Publication(14) 
}


save_fig <- function(fig_dir, fname, width, height, save_as) {
  if (save_as == 'png') {
    ggsave(file.path(fig_dir, paste0(fname, '.png')), width = width, height = height)  
  } else if (save_as == 'pdf') {
    ggsave(file.path(fig_dir, paste0(fname, '.pdf')), width = width, height = height, device = grDevices::cairo_pdf) 
  }
}

# plotEnrichment_mod <- function (pathway, stats, gseaParam = 1) 
# {
#   rnk <- rank(-stats)
#   ord <- order(rnk)
#   statsAdj <- stats[ord]
#   statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
#   statsAdj <- statsAdj/max(abs(statsAdj))
#   pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
#   pathway <- sort(pathway)
#   gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
#                           returnAllExtremes = TRUE)
#   bottoms <- gseaRes$bottoms
#   tops <- gseaRes$tops
#   n <- length(statsAdj)
#   xs <- as.vector(rbind(pathway - 1, pathway))
#   ys <- as.vector(rbind(bottoms, tops))
#   toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
#   diff <- (max(tops) - min(bottoms))/16
#   x = y = NULL
#   g <- ggplot(toPlot, aes(x = x, y = y)) + 
#     geom_point(color = "#386cb0", size = 0.1) + 
#     geom_hline(yintercept = max(tops), colour = "#ef3b2c", linetype = "dashed") + 
#     geom_hline(yintercept = min(bottoms), colour = "#ef3b2c", linetype = "dashed") + 
#     geom_hline(yintercept = 0, colour = "black") + 
#     geom_line(color = "#386cb0", lwd = 1) + 
#     theme_bw() + 
#     geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
#                                                                y = -diff/2, xend = x, yend = diff/2), size = 0.2) + 
#     theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
#     labs(x = "Gene rank", y = "Enrichment score")  +
#     theme_Publication(14)
#   g
# }


theme_Publication <- function(base_size=12, base_family="Arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

col_MSI_labs <- c(MSI = "#ef3b2c", 
                  MSS = '#386cb0', 
                  indeterminate = '#fdb462',
                  `MSS/MSI-L` = '#386cb0',
                  `MSI-H` = '#ef3b2c')



print_unpaired_two_group_stats <- function(var1, var2) {
  print(sprintf(
    'Wilcox p: %.3g',
    wilcox.test(var1, var2, paired = FALSE)$p.value
  ))
  
  print(sprintf(
    'Group 1: %d, med: %.3g',
    sum(!is.na(var1)),
    median(var1, na.rm=T)
  ))
  
  print(sprintf(
    'Group 2: %d, med: %.3g',
    sum(!is.na(var2)),
    median(var2, na.rm=T)
  ))
}

print_spearman_corr <- function(var1, var2) {
  stats <- cor.test(var1, var2, method = 'spearman')
  print(sprintf(
    'Spearman rho: %.3g',
    stats$estimate
  ))
  print(sprintf(
    'P-value: %.3g',
    stats$p.value
  ))
  print(sprintf(
    'Sample Size: %d',
    sum(!is.na(var1) & !is.na(var2))
  ))

}
