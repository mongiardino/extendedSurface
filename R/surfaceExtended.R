#' Interface between SURFACE and OUwie to fit complex multi-OU models
#'
#'\code{surfaceExtended} uses the output of a SURFACE run to attempt to fit more
#'complex Ornstein-Uhlenbeck models using OUwie. These models can incorporate
#'differences between regimes in the rate of evolution (\emph{sigma^2}),
#'strength of attraction to trait optima (\emph{alpha}), as well as optimize the
#'state at the root of the tree (\emph{theta0}). Optionally, the backwards phase
#'of SURFACE is then extended by further merging regimes and attempting to use
#'these simpler models as successive inputs to OUwie. Only one morphological
#'trait can be employed. The function keeps track of successful instances of
#'model fitting, and returns the models obtained through both SURFACE and OUwie,
#'as well as an overall summary.
#'
#'Paleontological data has been shown to improve the accuracy of models
#'describing morphological evolution using OU models (Ho & Ané 2014).
#'Nonetheless, many of the methods to fit multi-OU models that do not require
#'users to specify the number and location of regime shifts work only on
#'ultrametric trees. One that does not, SURFACE (Ingram & Mahler 2013), tends to
#'favor overly complex models (Khabbazian et al. 2016), likely a consequence of
#'assuming that regimes share a common \emph{sigma^2} and \emph{alpha}
#'parameters (Mongiardino Koch & Thompson 2020). Relaxing this assumption is not
#'straightforward, as estimating these parameters for models with multiple
#'regimes is often unfeasible (Benson et al. 2017).
#'
#'\code{surfaceExtended} employs the optimal model found by the SURFACE algorithm
#'and uses OUwie (Beaulieau et al. 2012) to attempt to fit multi-OU models in
#'which rates of evolution and strengths of selection vary between regimes.
#'The function then extends the backwards phase of SURFACE to merge independent
#'regimes and find simpler multi-OU models. This is done in a stepwise fashion,
#'and every time two regimes are merged, the result is used as input for OUwie.
#'
#'The user can specify which parameters to estimate for each regime with OUwie,
#'including different rates of evolution (\code{models = 'OUMV'}), different
#'strengths of selection (\code{models = 'OUMA'}), or both (\code{models =
#''OUMVA'}). The state at the root of the tree can be further considered an
#'independent parameter by adding 'Z' at the end of the model's name (e.g.,
#'\code{models = 'OUMAZ'}), although this can destabilize parameter estimates
#'(see \code{\link[OUwie]{OUwie}} for more details). Multiple models can be
#'explored simultaneously by providing a vector with their names, or using
#'\code{models = 'all'} or \code{models = 'all_noZ'}. In the latter, only models
#'assuming the root value is distributed according to the stationary
#'distribution of the ancestral OU process are optimized. The minimum number of
#'regimes to be explored is determined by the \code{limit} parameter. By
#'default, the results will be plotted using
#'\code{\link{extended_surfaceAICPlot}}.
#'
#'Can be very time-consuming depending on the size of the phylogeny and the
#'complexity (number of regimes) of the starting model.
#'
#' @param bwd_surface List of models obtained using the backwards phase of
#'   SURFACE.
#' @param data The morphological character under investigation. Needs to be
#'   the same data.frame used to run SURFACE.
#' @param tree Phylogenetic tree in 'phylo' format. Needs to be the same used to
#'   run SURFACE.
#' @param error Optional. Measurement errors to be incorporated in the process
#'   of model fitting. Can be a data.frame or vector, but it is assumed the
#'   order matches that of \code{data}.
#' @param models Character vector specifying the models to be explored (see
#'   Details). Explores 'OUMVA' and 'OUMVAZ' models by default.
#' @param limit Minimum number of regimes to explore. The default is 2.
#' @param plot A logical indicating whether to plot the AICc of models output by
#'   SURFACE and OUwie. Default is \code{TRUE}.
#' @param fwd_surface Optional. List of models obtained using the backwards
#'   phase of SURFACE. Only used to produce a more thorough comparison of models
#'   when \code{plot=TRUE}.
#' @return A list with the following elements:
#'  \describe{
#'    \item{$ext_surface}{A list containing all the models explored by extending
#'    the backwards phase of SURFACE, identical to the one returned by
#'    \code{\link[surface]{surfaceBackward}}.}
#'    \item{$summary}{A data.frame including information on all the models
#'    explored using OUwie, including whether model fit was succesfull, and if so
#'    the AICc value.}
#'  }
#'  Additionally, if model fit was successfull, the list will also include the
#'  best option found for each of the models specified with \code{models}. If
#'  multiple models were explored, the best option for each will be returned. For
#'  example, if OUMVA and OUMVAZ models were explored (as is the default), and
#'  model fitting was successful, the returned list will also contain two more
#'  elements, 'best_OUMVA' and 'best_OUMVAZ'. Note that these might differ in
#'  the number of regimes they contain.
#' @author
#'  Nicolás Mongiardino Koch
#' @references
#'  Beaulieu J.M., Jhuwueng D.‐C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: expanding the Ornstein–Uhlenbeck model of adaptive evolution. Evolution, 66:2369–2383.
#'  Benson R.B.J., Hunt G., Carrano M.T., Campione N. (2018), Cope's rule and the adaptive landscape of dinosaur body size evolution. Palaeontology, 61:13-48.
#'  Ho L.S.T, Ané C. 2014. Intrinsic inference difficulties for trait evolution with Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 5:1133–1146.
#'  Ingram T., Mahler D.L. 2013. SURFACE: detecting convergent evolution from comparative data by fitting Ornstein‐Uhlenbeck models with stepwise Akaike Information Criterion. Methods in Ecology & Evolution, 4:416–425.
#'  Khabbazian M., Kriebel R., Rohe K., Ané, C. 2016. Fast and accurate detection of evolutionary shifts in Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 7:811–824.
#'  Mongiardino Koch N., Thompson J.R. A Total-Evidence Dated Phylogeny of Echinoids and the Evolution of Body Size across Adaptive Landscape. bioRxiv 2020.02.13.947796, doi.org/10.1101/2020.02.13.947796.
#' @seealso For details on how these models are fit visit \code{\link[surface]{surfaceBackward}} and \code{\link[OUwie]{OUwie}}.
#'   Plots of AIC values can be obtained with \code{\link{extended_surfaceAICPlot}}
#' @export
#' @examples
#'  \dontrun{
#'   data(echinoid_data)
#'   surfaceExtended(bwd_surface = echinoid_data$bwd_surface, data =
#'   echinoid_data$size, tree = echinoid_data$tree, error = echinoid_data$error,
#'   models = 'OUMVAZ', limit = 4, plot = T, fwd_surface =
#'   echinoid_data$fwd_surface)
#'  }
#'

surfaceExtended <- function(bwd_surface, data, tree, error = NA, models = c('OUMVA', 'OUMVAZ'),
                           limit = 2, plot = T, fwd_surface = NA) {

  #Initial setup
  bwd_surface2 <- bwd_surface
  number_of_regimes <- as.numeric(bwd_surface2[[length(bwd_surface2)]]$n_regimes[2])
  all_regimes <- c(number_of_regimes:limit)

  ###Models
  if(models[1] == 'all') models <- c('OUMA', 'OUMAZ', 'OUMV', 'OUMVZ', 'OUMVA', 'OUMVAZ')
  if(models[1] == 'all_noZ') models <- c('OUMA', 'OUMV', 'OUMVA')
  root.station <- !grepl('Z', models)

  ###Data
  if(class(data) != 'data.frame') data <- data.frame(data)
  if(any(!rownames(data) %in% tree$tip.label)) data <- data[-which(!rownames(data) %in% tree$tip.label),]

  ###Tree
  if(length(tree$tip.label) != dim(data)[1]) {
    tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(data)])
  }
  if(is.null(tree$node.label)) tree <- surface::nameNodes(tree)
  tree$root.time <- max(phytools::nodeHeights(tree))

  data_surface <- surface::convertTreeData(tree, data)

  ###Set up summary of results
  num_models = length(all_regimes)*length(models)
  summary <- data.frame(model = paste(rep(all_regimes, each = length(models)), models, sep = '_'),
                       worked = rep(F, num_models), AICC = rep(NA, num_models))

  while(number_of_regimes >= limit) {

    #Prepare the data to fit models to the results of surfaceBackwards
    if(number_of_regimes == as.numeric(bwd_surface[[length(bwd_surface)]]$n_regimes[2])) {
      cat('Testing', paste(models, collapse = ' and '), 'with', number_of_regimes, 'regimes','\n')
      if(is.na(error[1])) {
        data_ouwie <- data.frame(rownames(data), NA, data)
        mserr <- 'none'
      } else {
        data_ouwie <- data.frame(rownames(data), NA, data, error)
        mserr <- 'known'
      }

      #get regimes from surface object and map them on the data and the tree
      for(i in 1:nrow(data_ouwie)) {
        taxon = which(as.character(unlist(attributes(bwd_surface[[length(bwd_surface)]]$fit[[1]])[15])) ==
                        as.character(data_ouwie[i,1]))
        data_ouwie[i,2] <- as.character(unlist(attributes(bwd_surface[[length(bwd_surface)]]$fit[[1]])[6])[taxon])
      }
      tree_ouwie <- tree
      for(i in 1:length(tree_ouwie$node.label)) {
        node = which(as.character(unlist(attributes(bwd_surface[[length(bwd_surface)]]$fit[[1]])[15])) ==
                       tree$node.label[i])
        tree_ouwie$node.label[i] <- as.character(unlist(attributes(bwd_surface[[length(bwd_surface)]]$fit[[1]])[6]))[node]
      }

    } else {   #continue backwards phase of surface into suboptimal models with fewer regimes
      cat('Collapsing regimes to find simpler models', '\n')
      bwd_surface3 <- surface::surfaceBackward(data_surface[[1]], data_surface[[2]], max_steps = 2,
                                               starting_model = bwd_surface2[[length(bwd_surface2)]],
                                               verbose = T, aic_threshold = 300, only_best = T, error_skip = T)
      if(as.numeric(bwd_surface3[[length(bwd_surface3)]]$n_regimes[2]) != as.numeric(bwd_surface2[[length(bwd_surface2)]]$n_regimes[2])) {
        bwd_surface2[[(length(bwd_surface2)+1)]] <- bwd_surface3[[(length(bwd_surface3))]]
      }

      number_of_regimes <- as.numeric(bwd_surface2[[length(bwd_surface2)]]$n_regimes[2])
      cat('Testing', paste(models, collapse = ' and '), 'with', number_of_regimes, 'regimes','\n')

      #get regimes from surface object and map them on the data and the tree
      for(i in 1:nrow(data_ouwie)) {
        taxon = which(as.character(unlist(attributes(bwd_surface2[[length(bwd_surface2)]]$fit[[1]])[15])) ==
                        as.character(data_ouwie[i,1]))
        data_ouwie[i,2] <- as.character(unlist(attributes(bwd_surface2[[length(bwd_surface2)]]$fit[[1]])[6])[taxon])
      }

      for(i in 1:length(tree_ouwie$node.label)) {
        node = which(as.character(unlist(attributes(bwd_surface2[[length(bwd_surface2)]]$fit[[1]])[15])) ==
                       tree$node.label[i])
        tree_ouwie$node.label[i] <- as.character(unlist(attributes(bwd_surface2[[length(bwd_surface2)]]$fit[[1]])[6]))[node]
      }
    }

    #Fit models with OUwie and make summary of results
    for(i in 1:length(models)) {
      ouwie_result <- OUwie::OUwie(tree_ouwie, data_ouwie, model = gsub('Z', '', models[i]), simmap.tree = F,
                                   root.age = tree_ouwie$root.time, scaleHeight = F, root.station = root.station[i],
                                   clade = NULL, mserr = mserr, starting.vals = NULL, diagn = T, warn = F)

      model_name <- paste(number_of_regimes, models[i], sep = '_')
      assign(model_name, ouwie_result)

      if(all(ouwie_result$eigval > 0)) {
        summary[which(summary$model == model_name),'worked'] <- T
        summary[which(summary$model == model_name),'AICC'] <- ouwie_result$AICc
      }
    }

    number_of_regimes <- number_of_regimes - 1
  }

  #Report on best models found
  results <- rep(T, length(models))
  for(i in 1:length(models)) {
    model_results <- summary[endsWith(as.character(summary$model), models[i]),]
    model_results <- model_results[which(!is.na(model_results$AICC)),]
    if(nrow(model_results) == 0) {
      results[i] <- F
      cat('For', models[i], 'no reliable model could be fit', '\n')
    } else {
      assign(paste0('best_', models[i]),
             get(as.character(model_results[which(model_results[,'AICC'] == min(model_results[,'AICC'])),'model'])))
      cat('For', models[i], 'a', nrow(get(paste0('best_', models[i]))$theta),
          'regime solution was the best, with AICc =', get(paste0('best_', models[i]))$AICc, '\n')
      if(exists('to_output')) {
        to_output <- c(to_output, paste0('best_', models[i]))
      } else {
        to_output <- paste0('best_', models[i])
      }
    }
  }

  if(nrow(model_results) > 0 && plot) extended_surfaceAICPlot(ext_surface = bwd_surface2,
                                                              summary = summary,
                                                              fwd_surface = fwd_surface)

  #Return the extended surface object, a summary of results and the best model of each type
  if(length(to_output) == 1) {
    result <- list(bwd_surface2, summary, get(to_output[1]))
    names(result) <- c('ext_surface', 'summary', to_output[1])
    return(result)
  } else {
    if(length(to_output) > 1) {
      result <- list(bwd_surface2, summary)
      for(i in 1:length(to_output)) {
        result[[length(result)+1]] <- get(to_output[i])
      }
      names(result) <- c('ext_surface', 'summary', to_output)
      return(result)
    } else {
      cat('No reliable solution was found. The result of the backwards phase of SURFACE remains the best OUM model.', '\n')
    }
  }
}
