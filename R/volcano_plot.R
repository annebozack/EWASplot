#' Volcano plot
#'
#' This function creates a volcano plot of p-values and effect sizes using ggplot2.
#'
#' @param pvector Vector of p-values
#' @param effectvector Vector of effect sizes
#' @param FDR Should an FDR level of significance be indicated? Defaults to FALSE.
#' @param title Title text.
#' @param col Color that will be used for probes for non-significant probes. Defaults to 'gray40'.
#' @param col.sig Color of FDR or Bonferroni significant probes. Defaults to 'gray18'.
#' @param cex Size of non-significant probes. Defaults to 0.5.
#' @param cex.sig Size of significant probes. Defaults to 0.5.
#' @param alpha Alpha of non-significant probes. Defaults to 0.5.
#' @param alpha.sig alpha Alpha of non-significant probes. Defaults to 1.
#' @param Size of line indicating FDR and Bonferroni significance. Defaults to 0.5.
#' @return Volcano plot
#' @export
volcano = function(pvector, effectvector, FDR = FALSE, title = NULL, col = 'gray40', col.sig = 'gray18', cex = 0.5, cex.sig = 0.75, alpha = 0.5, alpha.sig = 1, size.line.sig = 0.5){
	# dataframe for plotting
	voldat = data.frame(effect = effectvector, P.Value = pvector, logp = -log10(pvector))
	
	# Bonferroni and FDR adjustment of p-values
	voldat$bonf = p.adjust(voldat$P.Value, method = 'bonferroni')
	if (FDR){
		voldat$FDR = p.adjust(voldat$P.Value, method = 'fdr')
	}
	
	# variables for point color, size, and alpha
	voldat$color[voldat$bonf < 0.05] = col.sig  # DMP probes
	voldat$size[voldat$bonf < 0.05] = cex.sig
	voldat$alpha[voldat$bonf < 0.05] = alpha.sig
	if (FDR){
		voldat$color[voldat$FDR < 0.05] = col.sig
		voldat$size[voldat$FDR < 0.05] = cex.sig
		voldat$alpha[voldat$FDR < 0.05] = alpha.sig
	}
	voldat$size[is.na(voldat$color)] = cex
	voldat$alpha[is.na(voldat$color)] = alpha
	voldat$color[is.na(voldat$color)] = col	
	
	volcano = ggplot(voldat, aes(x = effect, y = logp)) + 
	geom_point(size = voldat$size, alpha = voldat$alpha) + 
	geom_hline(aes(yintercept = -log10(0.05/nrow(voldat))), color = "#AB3428", size = size.line.sig) + 
	theme_minimal() + scale_y_continuous(expand = c(0, 0)) + theme(panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(size = 0.2, color = 'gray65'), panel.grid.major.x = element_line(size = 0.2, color = 'gray65')) + labs(y=expression(-log[10](italic(p))), x='Effect estimate') 
	if (!is.null(title)){
    	volcano = volcano + labs(title = title)
    }
	if (FDR == TRUE){
		volcano = volcano + geom_hline(aes(yintercept = -log10(max(P.Value[FDR < 0.05]))), color = "#AB3428", size = size.line.sig, linetype = "dashed")
	}
	return(volcano)
}