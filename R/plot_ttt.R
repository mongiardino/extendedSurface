#' Plot trait dynamics through time
#'
#'This function plots evolutionary dynamics of discrete traits as summarized
#'using \code{\link{transitions_through_time}}. Two types of plots can be
#'generated: the number of different states active through time (where identical
#'states that have different origins are not counted as the same one), and the
#'rates at which these states are originating (birth rate), becoming extinct
#'(death rate), or accumulating (diversification rate) thorugh time. Depending
#'on the character being investigated, these plots might or might not be
#'meaningful.
#'
#'By default, this is used by \code{\link{transitions_through_time}} to plot
#'results. However, the object returned by that function can also be used here
#'with more control on the plotting options. These include the intervals (in Ma)
#'at which the number of states are recorded, the size of the window used to
#'smooth rates, the type of plot generated, and the type of rate to plot (see
#'Arguments).
#'
#'Trends are depicted using GAM regressions (see \code{\link[mgcv]{gam}}).
#'Depending on the combination of the size of the smoothing window and the
#'number of smoothing functions used, nonsensical results can be obtained. Some
#'tuning might be necessary to correctly depict trends in the data.
#'
#'
#' @param ttt Data.frame output by \code{\link{transitions_through_time}}.
#' @param interval Numeric value that determines the temporal resolution (i.e.,
#'   the step size in Ma at which the number of active regimes is recorded).
#'   Defaults to slightly less than 100 given the time spanned by the phylogeny.
#' @param window_size Numeric value that sets the width of the window used to
#'   smoot rate estimates. Defaults to a width that includes approx. 10
#'   intervals (see above).
#' @param CI Numeric value that sets the confidence interval (expressed as
#'   percentage). Determines the amount of results that are discarded before
#'   plotting (default = 80).
#' @param trim Whether to trim a few values at the begining and end of plot that
#'   contain fewer intervals and can be noisier. Default is \code{TRUE}.
#' @param graphs Which graphs to plot. Options include \code{'active_regimes'},
#'   \code{'rate_ttt'}, and \code{'both'}.
#' @param rates Which rates to plot. Options include any combination of \code{'birth'},
#'   \code{'death'}, and \code{'diversification'}. Defaults to only the latter.
#' @param k The value of k used for gam regression. If not specified this is
#'   automatically determined (see more details in \code{\link[mgcv]{gam}}).
#' @return A plot including different visual summaries of the evolutionary
#'   dynamics of discrete traits through time.
#' @author
#'  Nicolás Mongiardino Koch
#' @references
#'  Mongiardino Koch N., Thompson J.R. A Total-Evidence Dated Phylogeny of Echinoids and the Evolution of Body Size across Adaptive Landscape. bioRxiv 2020.02.13.947796, doi.org/10.1101/2020.02.13.947796.
#' @seealso \code{\link{transitions_through_time}}
#' @export
#'

plot_ttt <- function(ttt, interval = NA, window_size = NA, CI = 80, trim = T,
                     graphs = 'both', rates = 'diversification', k = NA) {

  #total timespan of tree
  timespan <- max(ttt$Birth)

  #set interval for counting the number of active regimes
  #if not specified this defaults to an interval that guarantees slightly less than 100 points from root to tip
  if(is.na(interval)) {
    interval <- ceiling(timespan/100)
  }
  #set window size for estimating birth/death rates
  if(is.na(window_size)) {
    window_size <- ceiling(timespan/10)
  }
  #set which graphs to plot
  active_regimes <- rate_ttt <- T
  if(graphs == 'active_regimes') {
    rate_ttt <- F
  }
  if(graphs == 'ttt') {
    active_regimes <- F
  }

  #get number of replicates of stochastic character mapping that were used
  repl <- max(ttt$Repl)

  for(i in 1:repl) { #loop through the replicates

    #get one replicate at a time and sort events by their time of occurrence
    data <- ttt[which(ttt$Repl == i),]
    times <- sort(c(data[,1],data[,2]), decreasing = T)
    times <- times[-which(times == 0)]

    #place to store how the numbe of active regimes changed with each event
    n_optima <- rep(1, length(times))

    #fill it up depending on whetherthey represent regime deaths or births
    for(j in 2:length(n_optima)) {
      if(times[j] %in% data$Birth) {
        n_optima[j] <- n_optima[j-1] + 1
      } else {
        n_optima[j] <- n_optima[j-1] - 1
      }
    }

    #set up the sampling point depending on the chosen interval
    sampling_times <- rev(seq(0, timespan, interval))
    sampling_optima <- vector(length = length(sampling_times))

    for(j in 1:length(sampling_optima)) {
      if(j == length(sampling_optima)) {
        sampling_optima[j] <- n_optima[length(n_optima)]
      } else {
        sampling_optima[j] <- n_optima[which(times > sampling_times[j])[length(which(times > sampling_times[j]))]]
      }
    }

    #also express this as either gains or losses of regimes
    birth_death <- vector(length = length(n_optima))
    for(j in 1:length(n_optima)) {
      if(j == 1) {
        birth_death[j] <- 0
      } else {
        birth_death[j] <- n_optima[j] - n_optima[j-1]
      }
    }

    birth <- rep(0, length(sampling_times))
    death <- rep(0, length(sampling_times))

    times <- times[2:length(times)]
    birth_death <- birth_death[2:length(birth_death)]

    #from all of this, estimate birth and death rates over a window of chosen width
    for(j in 1:length(sampling_times)) {
      min <- sampling_times[j]+(window_size/2)
      max <- sampling_times[j]-(window_size/2)
      if(any(times < min & times > max)) {
        birth[j] <- length(which(birth_death[which(times < min & times > max)] == 1))/window_size
        death[j] <- length(which(birth_death[which(times < min & times > max)] == -1))/window_size
      } else{
        birth[j] <- 0
        death[j] <- 0
      }
    }

    #express this as a regime diversification rate as well
    diversification <- birth - death

    #put everything together into two different dataframes, one with number of active regimes through time
    #and another one with rates of change (birth, death, diversification) through time
    if(i == 1) {
      optima_to_plot <- data.frame(sampling_times, sampling_optima, repl = rep(i, length(sampling_optima)))
      diverse_to_plot <- data.frame(sampling_times = rep(sampling_times, 3), rates = c(birth, death, diversification),
                                    type = rep(c('birth', 'death', 'diversification'), each = length(sampling_times)),
                                    repl = rep(i, (length(sampling_optima)*3)))
    } else {
      optima_to_plot <- rbind(optima_to_plot,
                              data.frame(sampling_times, sampling_optima, repl = rep(i, length(sampling_optima))))
      diverse_to_plot <- rbind(diverse_to_plot,
                               data.frame(sampling_times = rep(sampling_times, 3),
                                          rates = c(birth, death, diversification),
                                          type = rep(c('birth', 'death', 'diversification'), each = length(sampling_times)),
                                          repl = rep(i, (length(sampling_optima)*3))))
    }
  }

  if(active_regimes) {
    median <- optima_to_plot %>% dplyr::group_by(sampling_times) %>%
      dplyr::summarise(median_optima = median(sampling_optima), .groups = 'drop')
    freq <- optima_to_plot %>% dplyr::group_by(sampling_times, sampling_optima) %>% dplyr::tally() %>%
      dplyr::arrange(sampling_times, n) %>% dplyr::ungroup()

    added <- rep(0, nrow(freq))
    for(i in 1:nrow(freq)) {
      if(i == 1) {
        count <- freq$n[i]
      } else {
        if(freq$sampling_times[i-1] != freq$sampling_times[i]) {
          count <- freq$n[i]
        } else {
          count <- count + freq$n[i]
        }
      }
      added[i] <- count
    }

    #remove data depending on chosen value of CI
    freq <- freq %>% dplyr::mutate(count = added) %>% dplyr::filter(count > (100-CI))

    if(is.na(k)) {
      plotA <- ggplot2::ggplot(median, aes(-sampling_times, median_optima)) +
        ggplot2::geom_point(data = freq, aes(-sampling_times, sampling_optima, color = n), shape = 16, size = 3) +
        ggplot2::scale_color_gradient(high = '#525252', low = '#f0f0f0') +
        ggplot2::stat_smooth(formula = y ~ s(x), method = "gam", se = FALSE, size = 2) +
        ggplot2::theme_bw() + ggplot2::ylab('Number of active regimes') + ggplot2::xlab('Time (Ma)') +
        ggplot2::theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
                       strip.text.x = element_text(size = 12), legend.position='none') +
        ggplot2::scale_x_continuous(breaks = pretty(-median$sampling_times),
                                    labels = abs(pretty(-median$sampling_times)))
    } else {
      plotA <- ggplot2::ggplot(median, aes(-sampling_times, median_optima)) +
        ggplot2::geom_point(data = freq, aes(-sampling_times, sampling_optima, color = n), shape = 16, size = 3) +
        ggplot2::scale_color_gradient(high = '#525252', low = '#f0f0f0') +
        ggplot2::stat_smooth(formula = y ~ s(x, k = k), method = "gam", se = FALSE, size = 2) +
        ggplot2::theme_bw() + ggplot2::ylab('Number of active regimes') + ggplot2::xlab('Time (Ma)') +
        ggplot2::theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
                       strip.text.x = element_text(size = 12), legend.position='none') +
        ggplot2::scale_x_continuous(breaks = pretty(-median$sampling_times),
                                    labels = abs(pretty(-median$sampling_times)))
    }
  }

  if(rate_ttt) {
    median_diverse <- diverse_to_plot %>% dplyr::group_by(sampling_times, type) %>%
      dplyr::summarise(median_diverse = median(rates), .groups = 'drop')
    freq_diverse <- diverse_to_plot %>% dplyr::group_by(sampling_times, rates, type) %>% dplyr::tally() %>%
      dplyr::arrange(sampling_times, type, n) %>% dplyr::ungroup()

    added <- rep(0, nrow(freq_diverse))
    for(i in 1:nrow(freq_diverse)) {
      if(i == 1) {
        count <- freq_diverse$n[i]
      } else {
        if(freq_diverse$sampling_times[i-1] != freq_diverse$sampling_times[i] ||
           freq_diverse$type[i-1] != freq_diverse$type[i]) {
          count <- freq_diverse$n[i]
        } else {
          count <- count + freq_diverse$n[i]
        }
      }
      added[i] <- count
    }

    #filter rates based on chosen CI
    freq_diverse <- freq_diverse %>% dplyr::mutate(count = added) %>% dplyr::filter(count > (100-CI))
    limits <- bind_cols(freq_diverse %>% dplyr::group_by(sampling_times, type) %>%
                          dplyr::summarise(min = min(rates), .groups = 'drop'),
                        (freq_diverse %>% dplyr::group_by(sampling_times, type) %>%
                           dplyr::summarise(max = max(rates), .groups = 'drop'))[,3])

    for(i in 1:length(rates)) {
      if(is.na(k)) {
          plot <- dplyr::bind_cols(subset(median_diverse, type == rates[i]), subset(limits, type == rates[i])[,3:4]) %>%
          ggplot2::ggplot(aes(x = -sampling_times, y = median_diverse)) +
          ggplot2::stat_smooth(formula = y ~ s(x), method = "gam", se = FALSE, size = 2) +
          ggplot2::stat_smooth(aes(y = min), formula = y ~ s(x), method = "gam", se = FALSE, size = 2) +
          ggplot2::stat_smooth(aes(y = max), formula = y ~ s(x), method = "gam", se = FALSE, size = 2)
      } else {
        plot <- dplyr::bind_cols(subset(median_diverse, type == rates[i]), subset(limits, type == rates[i])[,3:4]) %>%
          ggplot2::ggplot(aes(x = -sampling_times, y = median_diverse)) +
          ggplot2::stat_smooth(formula = y ~ s(x, k = k), method = "gam", se = FALSE, size = 2) +
          ggplot2::stat_smooth(aes(y = min), formula = y ~ s(x, k = k), method = "gam", se = FALSE, size = 2) +
          ggplot2::stat_smooth(aes(y = max), formula = y ~ s(x, k = k), method = "gam", se = FALSE, size = 2)
      }
      plot_data <- ggplot2::ggplot_build(plot)
      if(exists('plot_data_final')) {
        plot_data_final <- rbind(plot_data_final, data.frame(x = plot_data$data[[1]]$x,
                                                             y = plot_data$data[[1]]$y,
                                                             ymin = plot_data$data[[2]]$y,
                                                             ymax = plot_data$data[[3]]$y,
                                                             type = rates[i]))
      } else {
        plot_data_final <- data.frame(x = plot_data$data[[1]]$x,
                                      y = plot_data$data[[1]]$y,
                                      ymin = plot_data$data[[2]]$y,
                                      ymax = plot_data$data[[3]]$y,
                                      type = rates[i])
      }
    }

    if(trim) {
      to_trim <- (ceiling((window_size/interval)/5))*interval
      plot_data_final <- plot_data_final[-which(plot_data_final$x < min(plot_data_final$x) + to_trim),]
      plot_data_final <- plot_data_final[-which(plot_data_final$x > max(plot_data_final$x) - to_trim),]
    }

    plotB <- ggplot2::ggplot(plot_data_final, aes(x = x, y = y, color = as.factor(type))) +
      ggplot2::geom_ribbon(data = plot_data_final, aes(x = x, ymax = ymin, ymin = ymax), fill = "grey", alpha = 0.5) +
      ggplot2::geom_line(size = 2) + ggplot2::geom_hline(yintercept = 0, linetype="dashed") + ggplot2::theme_bw() +
      ggplot2::ylab('Rate per Ma') + ggplot2::xlab('Time (Ma)') +
      ggplot2::theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
                     strip.text.x = element_text(size = 12), legend.position="bottom") +
      ggplot2::labs(color = "Type of rate") + ggplot2::scale_color_manual(values = c('#F95335','#50A3A4','#674A40'))

    if(exists('plotA')) {
      plotB <- plotB + ggplot2::scale_x_continuous(breaks = pretty(plot_data_final$x),
                                                   labels = abs(pretty(plot_data_final$x)),
                                                   limits = layer_scales(plotA)$x$range$range)
    } else {
      plotB <- plotB + ggplot2::scale_x_continuous(breaks = pretty(plot_data_final$x),
                                                   labels = abs(pretty(plot_data_final$x)))
    }
  }

  if(active_regimes & rate_ttt) {
    print(plot_grid(plotA, plotB, ncol = 1, align = "v"))
  } else {
    if(exists('plotA')) {
      plot(plotA)
    } else {
      plot(plotB)
    }
  }
}
