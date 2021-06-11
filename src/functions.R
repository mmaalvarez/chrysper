
#### build dataframe (for one gene) for regression
buildDF = function(g, eset, conditions, refs){
  
  # extract gene from eset
  x = eset[fData(eset)$gene == g, ]
  
  # initialize data frame (counts, gRNA, gene, sample)
  df = data.frame(counts = as.vector(exprs(x)),
                  gRNA = rep(rownames(x), ncol(x)),
                  gene = g,
                  sample = rep(colnames(x),
                               rep(nrow(x),
                                   ncol(x))),
                  stringsAsFactors = F)
  
  # append conditions to the data frame
  for (condition in conditions){
    
    factor_vector = pData(x)[df$sample, condition]
    
    is_factor_missing = is.null(factor_vector)
    if (is_factor_missing == FALSE){ # factor exists in this dataset (IT SHOULD)
      df %<>%
        mutate(!!condition := factor_vector)
    }
  }
  
  # return df for NB fitting function "fitModel()"
  df
}



#### run NB regression on one gene
fitModel = function(g, eset, conditions, df){
  
  ### regression: generalized linear model, family negative binomial
  y = glm.nb(data = df,
             formula = counts ~ treatment*time + A3A_vector + HMCES + cell_line + TP53 + replicate + gRNA)
  
  ### return model summary (estimates, pvals, residual deviance) + residual deviance significance + condition means
  
  ## resdev: probability of observing by chance a residual deviance as far in value to the degrees of freedom
  # chi-square's p-value of whether the difference between the residual deviance and the degrees of freedom of the model can be due to random effects
  rd_df = pchisq(deviance(y),
                 df = df.residual(y),
                 # we don't want to reject the H0
                 lower.tail = F)
  
  ## level means at each condition
  x = eset[fData(eset)$gene == g, ]
  means = list()
  for (condition in conditions){
    condition_column = pData(x)[df$sample, condition]
    if (is.null(condition_column) == FALSE){
      factor_means = tapply(df$counts,
                            condition_column,
                            mean)
      means[[condition]] = factor_means
    }
  }
  # convert list of means into data frame
  means = melt(means,
               value.name = "count_means",
               varnames = c('level')) %>%
    rename("condition" = "L1") %>%
    select(condition, level, count_means)
  
  list("Model summary" = summary(y),
       "p-val resdev vs. d.f." = rd_df,
       "Count means by factor level" = means)
  
  ### compare real counts vs. fitted counts, we want them to be similar (stratify by experiment group "id")
  #cbind(df$counts, (df %>% mutate(sample = gsub("_.*", "", sample), sample = gsub("^[0-9]", "", sample), sample = gsub("^[0-9]", "", sample)) %>% pull(sample)), fitted(y)) %>% as.data.frame() %>% rename("real counts" = "V1", "id" = "V2", "fitted counts" = "V3") %>% ggplot(aes(x = as.numeric(`real counts`), y = as.numeric(`fitted counts`), group = id)) + geom_point() + geom_smooth(method = lm) + ggpubr::stat_cor(method = "pearson", alternative = "two.sided", cor.coef.name = "R") + facet_grid(facets = ~`id`) + xlab("real counts") + ylab("NB fitted counts") + theme_minimal()
}
