#' @export
processdf1 = function(probe){
  chrdf = data.frame(chr = unique(probe$chr), chr_len = NA, tot = NA)
  chrdf = chrdf[order(chrdf$chr),]
  for (i in 1:nrow(chrdf)){
  	chrdf$chr_len[i] = as.numeric(max(probe$pos[probe$chr == chrdf$chr[i]]))
  	chrdf$tot[i] = cumsum(chrdf$chr_len[1:i])[i]-chrdf$chr_len[i]
  }
  chrdf = chrdf[,c('chr', 'tot')]
  return(chrdf)
}

#' @export
processdf2 = function(chrdf, probe){
  plotdf = merge(probe, chrdf, by = 'chr', all.x = T)
  plotdf = plotdf[order(plotdf$pos),]
  plotdf = plotdf[order(plotdf$chr),]
  plotdf$poscum = plotdf$pos + plotdf$tot
  return(plotdf)
}

#' Manhattan plot
#'
#' This function creates a Manhattan plot of EWAS results using ggplot2.
#'
#' @param probe Dataframe of EWAS results, including the columns 'cpg', 'P.Value', 'chr', and 'pos'. If 'cpg' column is not included, rownames will be used as probe names. If output from limma's topTable is provided, chromosome and position will be obtained from 450K of EPIC annotations.
#' @param region Dataframe of regional analysis results including the columns 'chr', 'start', 'end'.
#' @param array '450K' or 'EPIC'
#' @param FDR Should an FDR level of significance be indicated? Defaults to FALSE.
#' @param title Title text.
#' @param col.chr Vector of two colors that will be alternately used for probes within each chromosome. Defaults to c('gray40', 'gray55').
#' @param col.sig Color of FDR or Bonferroni significant probes. Defaults to 'gray18'.
#' @param col.dmr Color of probes within DMRs. Defaults to 'black'.
#' @param line.dmr Color of line indicating DMRs. Defaults to '#4393c3'.
#' @param size.line.dmr Size of line indicating DMRs. Defaults to 0.2.
#' @param cex Size of non-significant probes. Defaults to 0.5.
#' @param cex.sig Size of significant probes. Defaults to 0.5.
#' @param alpha Alpha of non-significant probes. Defaults to 0.5.
#' @param alpha.sig alpha Alpha of non-significant probes. Defaults to 1.
#' @param Size of line indicating FDR and Bonferroni significance. Defaults to 0.5.
#' @return Manhattan plot
#' @export
manhattan_plot = function(probe, region = NULL, array = c('450K', 'EPIC'), FDR = FALSE, title = NULL, col.chr = c('gray40', 'gray55'), col.sig = 'gray18', col.dmr = 'black', line.dmr = '#4393c3', size.line.dmr = 0.2, cex = 0.5, cex.sig = 0.75, alpha = 0.5, alpha.sig = 1, size.line.sig = 0.5){
	if (!('data.frame' %in% class(probe))){
		stop('probe must be a dataframe')
	}
	if (!('P.Value' %in% colnames(probe))){
		stop('probe must contain the column "P.Value"')
	}
	
	# Using row names as probe names if not provided
	if (!('cpg' %in% colnames(probe))){
		probe$cpg = rownames(probe)
	}
	
	# Adding chr and pos if not provided
	if (!('chr' %in% colnames(probe))){
		if (array == '450K'){
			data(Locations, package = 'IlluminaHumanMethylation450kanno.ilmn12.hg19')
			probe = cbind(probe, Locations[,c('chr')][match(probe$cpg, rownames(Locations))])
			probe = cbind(probe, Locations[,c('pos')][match(probe$cpg, rownames(Locations))])
		}
		if (array == 'EPIC'){
			data(Locations, package = 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19')
			probe = cbind(probe, Locations[,c('chr')][match(probe$cpg, rownames(Locations))])
			probe = cbind(probe, Locations[,c('pos')][match(probe$cpg, rownames(Locations))])
		}
		probe = data.frame(probe)
	}
	
	# chr as numeric
	if (!is.numeric(probe$chr)){
		probe$chr = as.numeric(gsub("chr", "", probe$chr))
	}
	
	# Bonferroni and FDR adjustment of p-values
	probe$bonf = p.adjust(probe$P.Value, method = 'bonferroni')
	if (FDR){
		probe$FDR = p.adjust(probe$P.Value, method = 'fdr')
	}
	
	# variables for point color, size, and alpha
	probe$color[probe$bonf < 0.05] = col.sig  # DMP probes
	probe$size[probe$bonf < 0.05] = cex.sig
	probe$alpha[probe$bonf < 0.05] = alpha.sig
	if (FDR){
		probe$color[probe$FDR < 0.05] = col.sig
		probe$size[probe$FDR < 0.05] = cex.sig
		probe$alpha[probe$FDR < 0.05] = alpha.sig
	}
	if (!is.null(region)){  # probes within regions
		for (i in 1:nrow(region)){
			probe$color[(probe$chr == region$chr[i]) & (probe$pos >= region$start[i]) & (probe$pos <= region$end[i])] = col.dmr
			probe$size[(probe$chr == region$chr[i]) & (probe$pos >= region$start[i]) & (probe$pos <= region$end[i])] = cex.sig
			probe$alpha[(probe$chr == region$chr[i]) & (probe$pos >= region$start[i]) & (probe$pos <= region$end[i])] = alpha.sig
		}
	}
	probe$size[is.na(probe$color)] = cex
	probe$alpha[is.na(probe$color)] = alpha
	probe$color[is.na(probe$color)] = col.chr[(probe$chr[is.na(probe$color)] %% 2) + 1]	
	
	# create dataframes for plotting probes and regions
	chrdf = processdf1(probe) 
	don = processdf2(chrdf, probe)
	if (!is.null(region)){
		don_reg = processdf2(chrdf, region)
	}
	
	axisdf = data.frame(chr = unique(probe$chr), center = NA)
	axisdf = axisdf[order(axisdf$chr),]
	for (i in 1:nrow(axisdf)){
		axisdf$center[i] = mean(c(max(don$poscum[don$chr == axisdf$chr[i]]), min(don$poscum[don$chr == axisdf$chr[i]])))
	}
			
	manhattan = ggplot2::ggplot(don, aes(x = poscum, y = -log10(P.Value))) +
	ggplot2::geom_point(color = don$color, size = don$size, alpha = don$alpha) +
	# custom axes:
	ggplot2::scale_x_continuous(expand = c(0.005, 0.005), limits = c(min(don$poscum), max(don$poscum)), label = axisdf$chr, breaks= axisdf$center) +
	ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, (max(-log10(don$P.Value)) + 0.5)), breaks = seq(from = 0, to = (max(-log10(don$P.Value)) + 0.5), by = 1)) +
	# Custom theme:
    ggplot2::theme_minimal() + ggplot2::theme( 
	legend.position="none", panel.border = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), text = element_text(size = 7.5)) + 
    ggplot2::labs(y=expression(-log[10](italic(p))), x='Chromosome') 
    if (!is.null(region)){
    	manhattan = manhattan + ggplot2::geom_vline(xintercept = don_reg$poscum, colour = line.dmr, size=size.line.dmr)
    }
    if (!is.null(title)){
    	manhattan = manhattan + ggplot2::labs(title = title)
    }
    # p-value cutoffs
    if (min(don$bonf) < 0.05){
    	manhattan = manhattan + ggplot2::geom_hline(yintercept=-log10(0.05/nrow(don)), colour = '#AB3428', size=size.line.sig)
    }	
    if (FDR == TRUE){
    	manhattan = manhattan + ggplot2::geom_hline(yintercept=-log10(max(don$P.Value[don$FDR < 0.05])), colour='#AB3428', size=size.line.sig, linetype = "dashed")
    } 
	return(manhattan)    	
}	

