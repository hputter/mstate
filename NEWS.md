# Version 0.3.3 (2024-07-03)
- Fixed issue with msprep() that occurred with very long vectors/mismatch between numeric id columns and integer id levels - thanks to @fumi-github for fixing #26 
- Fixed issue #16 inside msfit (dimension of hlp) - thanks to @MetzgerSK
- Fixed issue #27 with msfit() and vartype=="aalen", similar to an issue in previous version - thanks to @MetzgerSK
- Other fixes: #31 (dimensions preserved when transitions made by only one individual), #30 (two-state model with vartype=="aalen" now works), #13 (proper output returned for mssample() when output='data'),a minor issue with vis.mirror.pt(), Cuminc() update (coerces data to data.frame, and removes defunct failcode argument)

# Version 0.3.2 (2021-11-08)
- Multi-state models can now be run accounting for population mortality,
  by adding relsurv functionality, see Manevski et al. (2021)
- plot.msfit() allows confidence intervals around the cumulative hazards
- Various plotting fixed
- msprep() now accepts tibbles and data.tables
- Fixed issues #5 and #9 from github

# Version 0.3.1 (2020-12-17)
- Added test for Markov assumption from Titman & Putter (2020) in
  Biostatistics
- Added ggplot # Versions of plot.msfit and plot.probtrans, as
  well as visualisation of multiple probtrans objects and mirrored
  plots, see the vignette on visualisation
- Improved summary functions for msfit and probtrans objects
- Added help functions of transition matrices
- Fixed bug in msprep
- Oxygenized all R functions

# Version 0.2.12 (2019-12-10)
- fixed bug in Cuminc caused by new survival (>=3) package

# Version 0.2.11 (2018-04-06)
- fixed bug in msprep

# Version 0.2.10 (2016-12-03)
- Changed call to survfit in Cuminc function
- Changed to.trans2 function in misc.R

# Version 0.2.9 (2016-02-28)
- Repaired small bug in msprep introduced in # Version 0.2.8, in case
  of one covariate.
- Added LMAJ function to estimate non-parametrically transition
  probabilities for possibly non-Markov multi-state models. LMAJ
  stands for landmark Aalen-Johansen.

# Version 0.2.8 (2015-12-05)
- crprep function. Object orientation is introduced. Function is
  renamed to crprep.default. Code that calculates the weights has
  been rewritten. It contains an extra argument "strata" that
  computes value-specific weights. Fixed some errors that affect the
  way in which the id and keep columns are presented. Default value
  of prec.factor changed to 1000 (with 100, values may not be
  correct in case of ties). Many changes in help file. Creation of
  the counting process data set in function create.wData.omega.
- Added aidssi2 data set, which	contains data for multi-state
  analyses. Many changes in help file for aidssi. 
- expand.covs function. Object orientation is introduced. Renamed
  to expand.covs.mstate. expand.covs.default is used for competing
  risks data sets. 
- Repaired bug in events function. Rows corresponding to absorbing
  states were empty. The numbers entering absorbing states are no
  longer 0.
- Added Tentry and time in resulting data frame of cutLMms function.
- Warning added to ELOS, in case of negative predicted proabilities.
- Dimension check for beta.state in mssample repaired.
- In msprep function, id and keep are now correctly handled.
- Time points for fixedhor prediction corrected in probtrans. 
- In plot.probtrans, correct handling now of lty, lwd, col etc.
- Cuminc function is effectively made obsolete. For backward compatibility,
  it has been "redefined" by calling survfit with type="mstate" and
  translating it to the output it used to have. Print, plot and summary
  methods are the same as for survfit. Subsetting like survfit objects
  is not possible. Some of the functionality of defining events and
  failcodes has gone.

# Version 0.2.7 (2014-05-29)
- Added function cutLMms (mirrored after cutLM from dynpred) that cuts a multi-
  state data frame (object of type "msdata") at a landmark time point LM;
  administrative censoring at one common time point is also supported
- Added function xsect that takes a cross-section at a particular time point
- Added function ELOS that calculates expected length of stay (ELOS) from
  a probtrans object
- Added disclaimers to the EBMT data sets
- Repaired a bug in the example of events()
- Made a much user-friendlier # Version of the covariance matrix (list element
  K+1 with K states) of probtrans. NOTE: the dimension has changed from
  nt x K^2 x K^2 (with nt the distinct transition time points) to
  K^2 x K^2 x (nt+1); the extra time point (time 0) is also present in the column
  "time" of each of the data frames in [[1]] up to [[K]], containing the
  transition probabilities and SE's.
- plot.msfit now uses names rather than numbers for transitions (set legend=1:K
  if you don't want that  (K is no of transitions))
- Changed default colors for plot.probtrans
- Fixed a bug in plot.probtrans that misplaced text with filled or stacked
  probability plots when xlim has been set
- Fixed lwd bugs in plot.msfit and plot.probtrans
- Citation to Geskus, Biometrics 2010, updated

# Version 0.2.6 (2010-12-04)
- Citation updated
- Function crprep added to implement the method of Geskus, Biometrics 2010

# Version 0.2.5 (2010-10-12)
- Function events had been omitted from # Version 0.2.4; this has been restored
- Error fixed in handling of xlim in plot.probtrans
- Reference to the mstate paper in Computer Methods and Programs in Biomedicine
  updated

# Version 0.2.4 (2010-05-12)
- The ebmt4 dataset has been adapted; recovery and adverse event are abstracted.
  Simultaneous occurrences of recovery and adverse events are removed; these
  have been adapted so that recoveries occur half a day earlier. Also simultaneous
  relapse and death have been set to relapse, the death status has been set to 0.
  A new variable recae has been created, defined as the time at which both recovery
  and AE have been obtained. The corresponding status variable recae.s
  is only 1 in case the patient has experienced both recovery and AE.
- In probtrans, the argument direction now can take the values of "forward" and
  "fixedhorizon". The option "fixedhorizon" was previously called "backward".
- In the summary method for probtrans, a bug has beeen fixed. This summary used
  to show head and again head, this has been changed to head and tail.
- In msfit, the newdata argument is now first sorted according to trans
- The function redrank now asks for two formulas, rather than two character strings
  containing column names in the data for reduced rank and full covariates
- xlim has been added as explicit argument to plot.probtrans
- Error fixed in handling last time point in msfit

# Version 0.2.3 (2009-11-06)

- The results of calls to msfit and probtrans are now 'msfit' and 'probtrans'
  S3 objects
- Plot and summary methods for 'msfit' and 'probtrans' objects have been added
- The result of call to msprep is now an 'msdata' S3 object
- The first argument of expand.covs and the result of a call to expand.covs
  is now an 'msdata' object
- A print method for 'msdata' objects has been added
- A convenient function, transMat, to define transition matrices, kindly
  provided by Steven McKinney, has been added

# Version 0.2.2 (2009-09-15)

- Changed a small error in the vignette on page 9, where the written text was
  incompatible with the R output of c3
- mstate had a call to "model.newframe" in msfit from the survival package,
  which was removed in R # Version 2.9.2. In msfit,
  survival:::model.newframe(Terms, newdata, response=FALSE) has been replaced by
  model.frame(delete.response(Terms), newdata, xlev=object$xlevels)
- Replaced the example in redrank by a much shorter one, the original is now
  only mentioned in \dontrun{}


# Version 0.2.1 (2009-07-10)

- Added vignette showing how to do the analyses of Putter H, Fiocco M,
  Geskus RB (2007). Tutorial in biostatistics: Competing risks and multi-state
  models. Statistics in Medicine 26, 2389-2430. Usage: vignette("Tutorial")
- Changed name of cuminc to Cuminc to avoid confusion with the cuminc function
  of cmprsk
- Repaired bug in Cuminc that caused error when no events of certain cause
  were present
- Repaired bug in Cuminc that caused error when only single cause of failure
  was present
- Added option longnames to expand.covs to have the possibility of shorter
  covariate names
- Implemented "aalen" option in msfit for models without covariates
- Implemented "efron" method for ties in msfit; the coxph object that is used
  as input for msfit may have used either "breslow" or "efron" as method for
  ties. Previously only "breslow" was possible.
