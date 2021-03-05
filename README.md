# extendedSurface
extendedSurface is an approach to fitting complex multi Ornstein-Uhlenbeck models to comparative datasets. It combines methods implemented in [SURFACE](https://www.otago.ac.nz/ecoevotago/code/surface.html) and [OUwie](https://www.jeremybeaulieu.org/r.html). Additionally, several methods to visualize the temporal dynamics of discrete traits in pylogenies are also implemented.

## Installation
```R
install.packages("devtools")
devtools::install_github("mongiardino/extendedSurface")
```

## Description
### Macroevolutionary modeling
Multi-peak Ornstein-Uhlenbeck (OUM) processes are used to model the evolution of traits across a macroevolutionary adaptive landscape. Many approaches exist to fit these types of models, but they generally suffer from one of several shortcomings:
 1) They only work on ultrametric trees, even when paleontological data has been shown to improve accuracy of OU models (Ho & Ané 2014).
 2) They require *a priori* hypotheses on the number and location of regimes shifts (Beaulieau et al. 2012).
 3) They favor overly complex models (i.e., too many adaptive peaks; Khabbazian et al. 2016).
 
extendedSurface was developed in an attempt to solve these different issues. The method starts by employing the output of surface (Ingram & Mahler 2013), an approach that can employ trees with fossils and does not require any *a priori* input besides the trait data and tree topology. However, the method has been found to favor models that are too complex, a result that might stem from the fact that it assumes all regimes share a common rate of evolution (&sigma;<sup>2</sup>) and strength of attraction (&alpha;; Mongiardino Koch & Thompson 2020). This assumption is not only problematic from a model-fitting perspective, but is also likely biologically unrealistic. However, relaxing this assumption is not always feasible, as model fit with other approaches sush as OUwie (Beaulieu et al. 2012) might fail when attempting to optimize all of these parameters independently (Benson et al. 2018).

extendedSurface provides an easy interface between surface and OUwie, so that parameters of OUM models found with the former can be estimated with the latter. This allows to test the relative fit of models in which the regimes have different rates of evolution (OUMV), different strengths of attraction (OUMA), or both (OUMVA). Furthermore, it also allows optimizing the state at the root of the tree (Z<sub>0</sub>), which is otherwise assumed to be the optimal value of the ancestral regime (models OUMVZ, OUMAZ, OUMVAZ). All of this can be done both with the model preferred by the backwards phase of surface, as well as with simpler models obtained after independent regimes are merged (the extended phase of surface). Even though these simpler OUM models are suboptimal by surface standards, they have less number of parameters and model fitting with OUwie is more likely to be successful. Furthermore, simpler OUMV/OUMA/OUMVA(Z) models can fit the data drastically better that more complex OUM models. For example, Mongiardino Koch & Thompson (2020) found a 4-peak OUMVAZ model to strongly outperform the 7-peak OUM model favored by surface (see Fig. 1 below). These functions therefore provide a way to explore OUMV/OUMA/OUMVA(Z) models in non-ultrametric trees, without specifying *a prior* regime-shift hypotheses, and reducing problems of over-fitting.


![model_comparison](https://github.com/mongiardino/extendedSurface/blob/master/images/model_comparison.jpg)
**Fig. 1:** Comparison of the relative fit of multi Ornstein-Uhlenbeck models obtained for a comparative dataset of echinoid body size using Surface and extendedSurface (highlighted in yellow). extendedSurface explores suboptimal OUM models and uses them as starting points for the optimization of different &alpha; and &sigma;<sup>2</sup> parameters for each regime. In this case, exploration of a suboptimal 4-peak OUM model (&Delta;AICc = 5.76 with respect to the optimal 7-peak OUM model returned by Surface) allows the method to find an much better 4-peak OUMVAZ (&Delta;AICc = -27.81). Modified from Mongiardino Koch & Thompson (2020).


Although extendedSurface was developed to explore adaptive landscapes in non-ultrametric trees, as these attain higher levels of model accuracy than topologies that sample only extant lineages (Ho & Ané 2014), the SURFACE algorithm is still routinely employed to study macroevolutionary dynamics of living clades. Preliminary analyses show that extendedSurface is also able to find better fitting models when using trees containing only extant lineages (Fig. 2). 

![other_clades](https://github.com/mongiardino/extendedSurface/blob/master/images/other_clades.jpg)
**Fig. 2:** extendedSurface is able to improve upon the fit of multi-OU models supported by SURFACE, even when no fossils tips are sampled. The five clades explored  represent different radiations of vertebrate clades within Australia. Black dots show the progression of the forward and backward phase of SURFACE, while white shapes show models explored using extendedSurface. Even in the absense of fossil terminals, extendedSurface finds better fitting multi-OU models for all five clades, often times including as few as half of the number of regimes supported by SURFACE. Datasets are from Brennan & Keogh (2018).


### Visualizing discrete trait evolution
This package has been extended to provide novel ways of summarizing and visualizing the evolution of discrete traits in phylogenies. Although these approaches were designed to capture certain aspects of the evolution of selective regimes, they can be useful to capture the temporal dynamics of discrete traits in general (although this depends on the type of trait and the goals of the researcher). For similar types of plots see [this](http://blog.phytools.org/2017/11/visualizing-rate-of-change-in-discrete.html) entry of the phytools blog by Liam Revell.

These approaches rely on stochastic character mappings (SCMs; Bolback 2006, Revell 2012) to provide a sample of plausible evolutionary histories conditioned on observed states at the tips, and can be used with either ultrametric or non-ultrametric trees. Both the number of SCMs and the model of trait evolution can be specified. These replicates are then summarized using the function transitions_through_time, which records the time (in Ma) of all transitions between states. Given an event of transition between states at a branch, the tree topology is then traversed towards the tips and the last point in time when the derived state is still active is recorded. This can represent the moment when all descendants have transitioned to a different state, when they have all gone extinct (if fossils are sampled), or otherwise the present day. Thus, for every time any character originates, it's time of birth and death is recorded.

Summarizing and visualizing this data is handled by function plot_ttt (a plot is also returned by transitions_through_time, but more control on its graphical aspects can be gained by saving the function's output and passing it to plot_ttt). Two general types of graphics can be generated (see Fig. 3): A. The number of active regimes through time ('active_regimes'), and B. The rate at which regimes are originating, becoming extinct or accumulating ('rate_ttt'). The first of these is probably most relevant in the context of multi-peak OU models, as the number of active regimes through times can be linked to expected levels of morphological disparity or the breadth of occupied ecospace. The second plot however can be used for any type of discrete trait, as changes in the rate at which states originate and diversify can reveal unexpected patterns of trait dynamics through time that can be difficult to capture otherwise. By default both plots are produced, but the user can choose to output only one. Most attributes of these plots can be tuned, including the temporal resolution, the width of the sliding window used to smooth rates, the desired confidence interval, etc. For 'rate_ttt', the default is to plot the rate at which traits are accumulating (i.e., birth - death rates, analogous to the net diversification rate of species), but all 3 rates, or any one in particular, can be selected. For both types of plots, trends are smoothed using general additive models (GAMs).


![plot_examples](https://github.com/mongiardino/extendedSurface/blob/master/images/plot_examples.jpg)
**Fig. 3:** Example of visualization approach. **A.** The optimal OUMVAZ model for the dataset included in this package is coded as a discrete trait and its evolution explored using replicates of stochastic character mapping. **B.** In the top plot, the number of active regimes through time in 100 replicates are summarized. Grey dots show possible outcomes and their color scales with their relative frequency (darker = more likely). In the bottom plot, the rate at which regimes are accumulating through time (i.e., the difference between the birth rate and the death rate) is plotted. Blue (top) and red (bottom) lines are median values smoothed using GAM, confidence intervals span 80% of replicates.


## Usage
The package comes with exemplary data that includes a time-calibrated phylogeny of echinoids, body size and measurement errors for all tips, as well as the outputs of running the forward and backwards phases of Surface on this data.
```R
data(echinoid_data)
OUmodels <- surfaceExtended(bwd_surface = echinoid_data$bwd_surface,
                            data = echinoid_data$size,
                            tree = echinoid_data$tree,
                            error = echinoid_data$error,
                            models = 'OUMVAZ', limit = 4, plot = T,
                            fwd_surface = echinoid_data$fwd_surface)

regimes <- OUmodels$best_OUMVAZ$data[,1]
names(regimes) <- rownames(OUmodels$best_OUMVAZ$data)
ttt <- transitions_through_time(echinoid_data$tree, regimes, repl = 100, model = 'ER', plot = T)
```

## Author
Nicolás Mongiardino Koch. Department of Earth & Planetary Sciences, Yale University.

**_Citation:_** Mongiardino Koch, N & Thompson, J. R. 2020. A Total-Evidence Dated Phylogeny of Echinoids and the Evolution of Body Size across Adaptive Landscape. bioRxiv, https://doi.org/10.1101/2020.02.13.947796.

## References
Beaulieu J.M., Jhuwueng D.‐C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: expanding the Ornstein–Uhlenbeck model of adaptive evolution. Evolution, 66:2369–2383.

Benson R.B.J., Hunt G., Carrano M.T., Campione N. (2018), Cope's rule and the adaptive landscape of dinosaur body size evolution. Palaeontology, 61:13-48.

Bollback J.P. 2006. Stochastic character mapping of discrete traits on phylogenies. BMC Bioinformatics, 7:88.

Brennan I.G., Keogh J.S. 2018. Miocene biome turnover drove conservative body size evolution across Australian vertebrates. Proceedings of the Royal Society B 285:20181474.

Ho L.S.T, Ané C. 2014. Intrinsic inference difficulties for trait evolution with Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 5:1133–1146.

Ingram T., Mahler D.L. 2013. SURFACE: detecting convergent evolution from comparative data by fitting Ornstein‐Uhlenbeck models with stepwise Akaike Information Criterion. Methods in Ecology & Evolution, 4:416–425.

Khabbazian M., Kriebel R., Rohe K., Ané, C. 2016. Fast and accurate detection of evolutionary shifts in Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 7:811–824.

Revell, L.J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology & Evolution, 3:217-223
