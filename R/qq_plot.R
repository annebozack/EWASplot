#' Compute genomic inflation factor lambda
#'
#' This function computes the genomic inflation factor (λ) using the median approach.
#'
#' @param pvector Vector of p-values
#' @return Value of λ
#' @export
lambda = function(p) {
	return(median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1))
	}

#' Q-Q plot
#'
#' This function creates a quantile-quantile (Q-Q) plot of p-values using ggplot2.
#'
#' @param pvector Vector of p-values
#' @return Q-Q plot
#' @export
qq_plot = function(pvector){
	l = round(lambda(pvector), 3)
	o = -log10(sort(pvector, decreasing = FALSE))
	e = -log10(ppoints(length(pvector)))
	df = data.frame(o = o, e = e)
	ggplot2::ggplot(df, ggplot2::aes(e, o)) + ggplot2::geom_point(alpha = 0.5, size = 1) + ggplot2::geom_abline(intercept = 0, slope = 1, color = '#AB3428') + ggplot2::labs(y = expression(Observed ~ ~-log[10](italic(p))), x = expression(Expected ~ ~-log[10](italic(p)))) + theme_classic() + ggplot2::annotate("text", x = 1, y = (range(o)[2]*0.9), label = paste0('lambda = ', l))
}
