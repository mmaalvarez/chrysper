
#### build dataframe (for one gene) for regression
buildDF = function(g, eset, conditions, offset, cLine, type_treatment){
  
  # extract gene from eset
  x = eset[fData(eset)$gene == g, ]

  # initialize data frame (counts, sgRNA, gene, sample)
  df = data.frame(counts = as.vector(exprs(x)),
                  sgRNA = rep(rownames(x), ncol(x)),
                  gene = g,
                  sample = rep(colnames(x),
                               rep(nrow(x),
                                   ncol(x))),
                  stringsAsFactors = F) %>%
    ## append sample-related offset to df
    merge(offset, by = "sample")
    
  # append conditions to the data frame
  for (condition in conditions){
    
    factor_vector = pData(x)[df$sample, condition]
    
    is_factor_missing = is.null(factor_vector)
    if (is_factor_missing == FALSE){ # factor exists in this dataset (IT SHOULD)
      df %<>%
        mutate(!!condition := factor_vector)
    }
  }
  
  ## keep only cell line and type of treatment required
  
  df %<>%
    filter(cell_line == cLine &
             type == type_treatment)
  
  # return df for NB fitting function "fitModel()"
  df
}



#### run NB regression on one gene
fitModel = function(g, eset, conditions, df){
  
  ### regression: generalized linear model, family negative binomial

  formula = "counts ~ time_cat + offset(ln_sum_nontargeting)"

  ## run NB regression
  y = glm.nb(formula = formula,
             data = df,
             # see https://stackoverflow.com/questions/64002936/error-in-prettynum-internalformatx-trim-digits-nsmall-width-3l-invalid
             control = glm.control(maxit = 500)) # maximal number of IWLS iterations (default 25, not enough)
  
  ## store summary (remove possible Error)
  summary_y = try(summary(y), F)
  
  ## return list
  list("df" = df,
       "model summary" = summary_y)
}


# ## OLD stuff returned in fitModel()

# ### return list of:
# # - df with fitted values
# # - residual deviance significance (resdev)
# # - condition levels count means
# # - model (y)
# # - model summary (estimates, pvals, residual deviance...)
#
# ## merge df with fitted counts for later checking real vs. fitted counts, we want them to be similar. Also extract experiment group as "id"
# df %<>%
#   rename("real counts" = "counts") %>%
#   # extract experiment id
#   mutate(id = gsub("_.*", "", sample),
#          id = gsub("^[0-9]", "", id),
#          id = gsub("^[0-9]", "", id),
#          # add fitted counts
#          `fitted counts` = fitted(y)) %>%
#   select(-c(id))
#
# #by chance a residual deviance as far in value to the degrees of freedom
# # chi-square's p-value of whether the difference between the residual deviance and the degrees of freedom of the model can be due to random effects
# rd_df = pchisq(deviance(y),
#                df = df.residual(y),
#                # we don't want to reject the H0
#                lower.tail = F)
# 
# ## count means at each condition level
# x = eset[fData(eset)$gene == g, ]      # have to create x again, for some reason...
# means = list()
# for (condition in conditions){
#   condition_column = pData(x)[df$sample, condition]
#   if (is.null(condition_column) == FALSE){
#     factor_means = tapply(df$counts,
#                           condition_column,
#                           mean)
#     means[[condition]] = factor_means
#   }
# }
# # convert list of means into data frame
# means = melt(means,
#              value.name = "count_means",
#              varnames = c('level')) %>%
#   rename("condition" = "L1") %>%
#   select(condition, level, count_means) %>%
#   drop_na()
# 
# ## return list
# list("df" = df,
#      "p-val resdev vs. d.f." = rd_df,
#      "count means by factor level" = means,
#      "model" = y,
#      "model summary" = summary_y)
