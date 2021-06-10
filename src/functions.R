
buildDF = function(g, eset, ref, cond="group"){
    x = eset[fData(eset)$gene == g,]
    df = data.frame(counts=as.vector(exprs(x)), gRNA=rep(rownames(x), ncol(x)), sample=rep(colnames(x), rep(nrow(x), ncol(x))),  stringsAsFactors=F)
    df[, cond] = pData(x)[df$sample, cond]
    df$mouse = pData(x)[df$sample,]$mouse
    df$gene = g
    df[, cond] = relevel(df[, cond], ref=ref)
    df
}


fitModel = function(df, conts, cond="group"){
    df$group = df[, cond]
    lm = lmer(counts ~ group + (1 | gRNA) + (1|mouse/sample), data=df)
    gl = glht(lm, linfct=conts)
    sm = summary(gl, test=adjusted("none"))
    pvals = as.numeric(sm$test$pvalues)
    ci = as.data.frame(confint(sm)$confint)
    data.frame(ci=ci, pvals=pvals)
}
