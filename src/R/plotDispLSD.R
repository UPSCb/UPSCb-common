require("LSD") 

plotDispLSD <- function(cds, name = NULL, ymin, linecol = "#00000080", xlab = "mean of normalized counts", 
           ylab = "dispersion", log = "xy", cex = 0.45, ...) 
 {
   px = rowMeans(counts(cds, normalized = TRUE))
   sel = (px > 0)
   px = px[sel]
   py = fitInfo(cds, name = name)$perGeneDispEsts[sel]
   if (missing(ymin)) 
     ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 
                       0.1)
   heatscatter(log10(px), log10(pmax(py, ymin)), xlab = xlab, ylab = ylab, 
               pch = ifelse(py < ymin, 6, 16), cex = cex, 
               xaxt='n', yaxt='n', ...)

   # Fix logged axes labels
   atx <- axTicks(1)
   aty <- axTicks(2)
   xlabels <- sapply(atx, function (i)
                as.expression(bquote(10^ .(i)))
              )
   ylabels <- sapply(aty, function (i)
                as.expression(bquote(10^ .(i)))
              )
   axis(1, at=atx, labels=xlabels)
   axis(2, at=aty, labels=ylabels)
   xg = 10^seq(-0.5, 5, length.out = 100)
   lines(log10(xg), log10(fitInfo(cds, name = name)$dispFun(xg)), col = linecol, 
         lwd = 4, lty = 1
    )
 }

plotMALSD <- function (x, ylim, sign=0.05, col = 'forestgreen',
    linecol = "#00000080", xlab = "mean of normalized counts",
    ylab = expression(log[2] ~ fold ~ change), log = "x", cex = 0.45,
    ...)
{
    if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in%
        colnames(x))))
        stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
    x = subset(x, baseMean != 0)
    py = x$log2FoldChange
    if (missing(ylim))
        ylim = c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) *
            1.1
    heatscatter(log10(x$baseMean), pmax(ylim[1], pmin(ylim[2], py)),
        pch = ifelse(py < ylim[1], 6, ifelse(py > ylim[2], 2, 16)),
        cex = cex,
        xlab = xlab,
        ylab = ylab,
        xaxt = 'n',
        ylim = ylim, ...)

    # Fix logged x-axis
    atx <- axTicks(1)
    labels <- sapply(atx, function (i)
                as.expression(bquote(10^ .(i)))
              )
    axis(1, at=atx, labels=labels)
    abline(h = 0, lwd = 4, col = linecol, lty=1)

    # Mark the significant DEGs, quite ugly, right?
    sign <- subset(x, padj < sign)
    pointy <- sign$log2FoldChange
    points(log10(sign$baseMean),
      pmax(ylim[1], pmin(ylim[2], pointy)),
      pch = 1,
      col = col,
      cex = cex + 0.5)
}
