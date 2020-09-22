#' Summarize trait dynamics through time
#'
#'This function provides a way of summarizing the temporal dynamics of traits by
#'registering the times when there are transitions between character states.
#'First, a number of replicates of stochastic character mappings are created
#'(see \code{\link[phytools]{make.simmap}} for more details), and the phylogeny
#'is traversed in order to register the moments in time when there are
#'transitions between states. For every time a transition occurs, the function
#'will follow all descendants and also register the point in time where that
#'instance of a given state goes extinct (either by all descendants
#'transitioning to a different state or by all tips going extinct).
#'
#'The funtion only requires a tree and a character, and allows the user to
#'modify the number of replicates, the model used for character mapping, as well
#'as whether to return a plot summarizing the results. More control on what and
#'how to plot these results is attained by providing the output to the function
#'\code{\link{plot_ttt}}.
#'
#' @param tree Phylogenetic tree in 'phylo' format.
#' @param char Named vector including character states for all tips. Names need
#'   to correspond to tips in the phylogeny.
#' @param repl Number of replicates of stochastic character mapping (default = 100)
#' @param model A character specifying the model of evolution used to
#'   reconstruct the evolutionary history of the character. Options include
#'   \code{'ER'} (default), \code{'SYM'} and \code{'ARD'}, see
#'   \code{\link[ape]{ace}}.
#' @param plot A logical indicating whether to plot the results. Default is
#'   \code{TRUE}.
#' @return A data.frame including the times of birth and death of each instance
#'   of a character across replicates. This object can be passed to
#'   \code{\link{plot_ttt}} to obtain different visualizations.
#' @author
#'  Nicolás Mongiardino Koch
#' @references
#'  Mongiardino Koch N., Thompson J.R. A Total-Evidence Dated Phylogeny of Echinoids and the Evolution of Body Size across Adaptive Landscape. bioRxiv 2020.02.13.947796, doi.org/10.1101/2020.02.13.947796.
#' @seealso For details on stochastic character mapping visit
#'   \code{\link[phytools]{make.simmap}}. Plots can be explored with
#'   \code{\link{plot_ttt}}
#' @export
#'

transitions_through_time <- function(tree, char, repl = 100, model = 'ER', plot = T) {

  #generate repl number of stochastic character mappings using the specified model
  simmap <- phytools::make.simmap(tree, char, model, repl)

  for(i in 1:repl) { #loop through the mappings
    #branches_transitions stores the number of branches that contain events (births/deaths)
    num_transitions <- unlist(lapply(simmap[[i]]$maps, length)) - 1
    branches_transitions <- rep.int(1:length(num_transitions), times = num_transitions)

    #before_ and after_nodes store nodes at the beginning and end of branches with events
    before_node <- as.numeric(sapply(strsplit(rownames(simmap[[i]]$mapped.edge[branches_transitions,]), ','), '[', 1))
    after_node <- as.numeric(sapply(strsplit(rownames(simmap[[i]]$mapped.edge[branches_transitions,]), ','), '[', 2))

    #add ancestral regime
    after_node <- c(after_node, (length(simmap[[i]]$tip.label)+1))

    #get ages for all internal nodes and all tips
    node_ages <- max(phytools::nodeHeights(simmap[[i]])) -
      ape::dist.nodes(simmap[[i]])[(length(simmap[[i]]$tip.label)+1),(length(simmap[[i]]$tip.label)+1):ncol(ape::dist.nodes(simmap[[i]]))]
    tip_ages <- max(phytools::nodeHeights(simmap[[i]])) -
      ape::dist.nodes(simmap[[i]])[(length(simmap[[i]]$tip.label)+1),(1:length(simmap[[i]]$tip.label))]

    #which transitions occur in terminal branches?
    which_terminal <- c(after_node <= length(simmap[[i]]$tip.label))

    #data.frame to store the results
    data <- data.frame(Birth = rep(NA, length(after_node)), Death = rep(NA, length(after_node)), Repl = i)

    #store time spent in each regime for each branch
    times_in_regimes <- simmap[[i]]$maps

    for(j in 1:nrow(data)) { #loop through the events
      transient <- F
      if(j < (nrow(data) - 1)) {
        if(before_node[j+1] == before_node[j]) {  #transient event, changes again in same branch
          transient <- T
        }
      }

      if(transient) {
        data$Birth[j] <- node_ages[which(names(node_ages) == before_node[j])] -
          unlist(times_in_regimes[branches_transitions[j]])[which(which(before_node == before_node[j]) == j)]
        data$Death[j] <- data$Birth[j] -
          unlist(times_in_regimes[branches_transitions[j]])[which(which(before_node == before_node[j]) == j) + 1]
      } else {
        if(which_terminal[j]) { #it is a terminal branch and the regimes ends with the death of the terminal
          data$Death[j] <- tip_ages[which(names(tip_ages) == after_node[j])]
          data$Birth[j] <- data$Death[j] +
            unlist(times_in_regimes[branches_transitions[j]])[length(unlist(times_in_regimes[branches_transitions[j]]))]
        } else { #events that happen in internal branches and are not reversed within the same branch
          if(j != nrow(data)) {
            data$Birth[j] <- node_ages[which(names(node_ages) == before_node[j])] -
              sum(unlist(times_in_regimes[branches_transitions[j]])[1:(length(unlist(times_in_regimes[branches_transitions[j]])) - 1)])
          } else {
            data$Birth[j] <- max(phytools::nodeHeights(simmap[[i]]))
          }

          #store all tip descendants, this will help traverse the tree and see up until which point in time
          #the regime is still alive
          desc_tips <- unlist(phangorn::Descendants(simmap[[i]], after_node[j], type = 'tips'))

          if(j != nrow(data)) {
            time_spent <- unname(rep(unlist(times_in_regimes[branches_transitions[j]])[length(unlist(times_in_regimes[branches_transitions[j]]))], length(desc_tips)))
          } else {
            time_spent <- rep(0, length(desc_tips))
          }

          #get the two nodes directly descending from the branch with the event
          direct_desc <- simmap[[i]]$edge[which(simmap[[i]]$edge[,1] == after_node[j]), 2]

          k <- 1
          while(k <= length(direct_desc)) { #crawl across its descendant branches
            if(direct_desc[k] <= length(simmap[[i]]$tip.label)) { #if it is a terminal
              if(direct_desc[k] %in% after_node) { #check if there is an event in the branch
                #if there is, only add the time spent in ancestral regime
                time_spent[which(desc_tips == direct_desc[k])] <- time_spent[which(desc_tips == direct_desc[k])] +
                  unlist(times_in_regimes[branches_transitions[which(after_node == direct_desc[k])]])[1]
              } else { #if there are not transitions in this tip
                #add the time until the terminal goes extinct
                time_spent[which(desc_tips == direct_desc[k])] <- time_spent[which(desc_tips == direct_desc[k])] +
                  simmap[[i]]$edge.length[which(simmap[[i]]$edge[,2] == direct_desc[k])]
              }

              #advance
              k <- k + 1

            } else { #if it is an internal branch
              if(direct_desc[k] %in% after_node) { #if there is an event in the branch
                #find all descendants who are no longer in the regime and add the time spent
                descendants_this_branch <- unlist(phangorn::Descendants(simmap[[i]], direct_desc[k], type = 'tips'))
                time_spent[which(desc_tips %in% descendants_this_branch)] <- time_spent[which(desc_tips %in% descendants_this_branch)] +
                  unlist(times_in_regimes[branches_transitions[which(after_node == direct_desc[k])]])[1]

                #advance
                k <- k + 1
              } else { #if there are no events, find all direct descendants of this branch
                descendants_this_branch <- simmap[[i]]$edge[which(simmap[[i]]$edge[,1] == direct_desc[k]), 2]
                time_spent[which(desc_tips %in% unlist(phangorn::Descendants(simmap[[i]], descendants_this_branch, type = 'tips')))] <-
                  time_spent[which(desc_tips %in% unlist(phangorn::Descendants(simmap[[i]], descendants_this_branch, type = 'tips')))] +
                  simmap[[i]]$edge.length[which(simmap[[i]]$edge[,2] == direct_desc[k])]

                #add new descendant nodes
                direct_desc <- c(direct_desc, descendants_this_branch)

                k <- k + 1
              }
            }
          }

          #complete dataset with time of death of regime
          data$Death[j] <- data$Birth[j] - max(time_spent)
        }
      }
    }

    #combine data into final dataframe
    if(i == 1) {
      all_data <- round(data, 2)
    } else {
      all_data <- rbind(all_data, round(data, 2))
    }
  }

  #plot and return data
  if(plot) plot_ttt(ttt = all_data, interval = NA, window_size = NA, CI = 80, trim = T,
                    graphs = 'both', rates = 'diversification', k = NA)
  return(all_data)
}
