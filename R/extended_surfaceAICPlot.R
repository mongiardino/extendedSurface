#' Plot AICc values of multi-OU models
#'
#' Explore model fit by visually comparing the AICc values of the results
#' obtained using \code{surfaceExtended}. The fit of OUM models obtained using
#' forward, backwards and extended phases of SURFACE, as well as that of models
#' explored with OUwie, are plotted against the number of regimes.
#'
#' @param ext_surface List of models obtained using the extended phase of
#'   SURFACE. This is always the first element in the list returned by
#'   \code{surfaceExtended}.
#' @param summary A data.frame with AICc values of OUwie models. This is always
#'   the second element in the list returned by \code{surfaceExtended}.
#' @param fwd_surface Optional. List of models obtained using the forward phase
#'   of SURFACE. If provided these will also be included in the plot.
#' @return Plot of AICc values against number of regimes.
#' @seealso \code{\link{surfaceExtended}}
#' @export
#'

extended_surfaceAICPlot <- function(ext_surface, summary, fwd_surface = NA) {

  models <- gsub('^.*_', '', as.character(summary$model))
  pch <- c(22:25, 22:23)[1:length(models)]
  col <- rep('black', length(models))
  bg <- c(rep('white', 4), rep('grey', 2))[1:length(models)]

  if(is.na(fwd_surface[1])) {
    cat('The AICc plot will only include results from backwards and extended phases of SURFACE', '\n')
    summ <- surface::surfaceSummary(ext_surface)
  } else {
    summ <- surface::surfaceSummary(fwd_surface, ext_surface)
  }

  aics <- summ$aics
  nreg <- summ$n_regimes_seq[2,]

  nreg_extended <- c()
  aics_extended <- c()
  model_type <- c()
  for(i in 1:length(models)) {
    model_results <- summary[endsWith(as.character(summary$model), models[i]),]
    model_results <- model_results[which(!is.na(model_results$AICC)),]
    if(nrow(model_results) == 0) {
      pch <- pch[-i]
      col <- col[-i]
      bg <- bg[-i]
    } else {
      nreg_extended <- c(nreg_extended, as.numeric(gsub('_.*$', '', model_results$model)))
      model_type <- c(model_type, gsub('^.*_', '', model_results$model))
      aics_extended <- c(aics_extended, model_results$AICC)
    }
  }

  model_type <- c(rep('OUM', length(aics)), model_type)
  aics <- c(aics, aics_extended)
  nreg <- c(nreg, nreg_extended)

  plot(nreg, aics, type = "n", ylim = range(aics))
  abline(h = aics[1], lty = 2)
  abline(h = min(aics[which(model_type == 'OUM')]), lty = 4)
  if(min(aics[which(model_type == 'OUM')]) > min(aics)) abline(h = min(aics), lty = 6)

  points(nreg[which(model_type == 'OUM')], aics[which(model_type == 'OUM')], type = 'l', lwd = 2)
  points(nreg[which(model_type == 'OUM')], aics[which(model_type == 'OUM')], pch = 21, bg = 'black', cex = 1.5)

  for(i in min(nreg[which(model_type != 'OUM')]):max(nreg[which(model_type != 'OUM')])) {
    to_connect_y <- c(aics[which(model_type != 'OUM' & nreg == i)],
                     min(aics[which(nreg[which(model_type == 'OUM')] == i)]))
    to_connect_x <- rep(i, 2)

    points(to_connect_x, range(to_connect_y), type = 'l', lwd = 1)

    for(j in 1:(length(to_connect_y)-1)) {
      k <- which(unique(model_type) == model_type[model_type != 'OUM' & nreg == i][j])-1
      points(i, to_connect_y[j], pch = pch[k], col = col[k], bg = bg[k], cex = 1.5)
    }
  }

  legend('topright', legend = unique(model_type), inset = 0.05, pt.bg = c('black', bg), bty = 'n',
         pch = c(21,pch), pt.cex = 1.5)
}
