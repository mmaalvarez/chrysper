
#### build dataframe (for one gene, one cell line, and one type of treatment) to feed into the regression
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
  
  ## keep only the type of treatment required
  df %<>%
    filter(type == type_treatment)
  
  ## and cell line, if not pooled
  if (cLine != "pooled"){
    df %<>%
      filter(cell_line == cLine)
  }

  # return df for fitting function "fitModel()"
  df
}


#### run regression on one gene, one/pooled cell line, and one type of treatment
fitModel = function(df){
  
  formula = "counts ~ time_cat + offset(ln_sum_control)"
  
  # it's necessary to wrap the regression within a tryCatch(), because otherwise the possible warning/error in glm.nb() aborts and exits the loop
               # first try negative binomial generalized linear model
  y = tryCatch(glm.nb(formula = formula, data = df),
               # if there is a failure to converge to theta, run a Poisson generalized linear model
               warning = function(w) glm(formula = formula, data = df, family = poisson),
               error = function(e) glm(formula = formula, data = df, family = poisson))
  ## return model
  y
}
