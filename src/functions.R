### Functions

buildDF <- function(g, eset, ref, cond="group"){
    x <- eset[fData(eset)$gene == g,]
    df <- data.frame(counts=as.vector(exprs(x)), gRNA=rep(rownames(x), ncol(x)), sample=rep(colnames(x), rep(nrow(x), ncol(x))),  stringsAsFactors=F)
    df[, cond] <- pData(x)[df$sample, cond]
    df$mouse <- pData(x)[df$sample,]$mouse
    df$gene <- g
    df[, cond] <- relevel(df[, cond], ref=ref)
    df
}

fitModel <- function(df, conts, cond="group"){
    df$group <- df[, cond]
    if(length(unique(df$gRNA)) > 1){
        lm <- lmer(counts ~ group + (1 | gRNA) + (1|mouse/sample), data=df)
        gl <- glht(lm, linfct=conts)
        sm <- summary(gl, test=adjusted("none"))
        pvals <- as.numeric(sm$test$pvalues)
        ci <- as.data.frame(confint(sm)$confint)
    }
    else {
        y <- lapply(conts$group, function(i){
            ref <- sub(".*.- ", "", sub(" = 0", "", i))
            sam <- sub(" -.*.", "", i)
            df$group <- relevel(df$group, ref =ref)
            lm <- lm(counts ~ group, data=df)
            sm <- summary(lm)
            pvals <- coef(sm)[paste0("group", sam),4]
            ci <- data.frame(Estimate=coef(sm)[paste0("group", sam), 1], confint(lm)[paste0("group", sam),, drop=F])
            list(ci=ci, pvals=pvals)
        })
        pvals <- do.call(rbind, lapply(y,"[[", "pvals")); rownames(pvals) <- conts$group
        ci <- do.call(rbind, lapply(y,"[[", "ci")); rownames(ci) <- conts$group
    }
    data.frame(ci=ci, pvals=pvals)
}
