#' Data from the Amsterdam Cohort Studies on HIV infection and AIDS
#' 
#' These data sets give the times (in years) from HIV infection to AIDS, SI
#' switch and death in 329 men who have sex with men (MSM). Data are from the
#' period until combination anti-retroviral therapy became available (1996).
#' For more background information on the cohort, ccr5 and SI, see Geskus
#' \emph{et al.} (2000, 2003)
#' 
#' \code{aidssi} contains follow-up data until the first of AIDS and SI switch.
#' It was used as example for the competing risks analyses in Putter, Fiocco,
#' Geskus (2007) and in Geskus (2016).
#' 
#' \code{aidssi2} extends the \code{aidssi} data set in three ways. First, it
#' considers events after the initial one. Second, it includes the entry times
#' of the individuals that entered the study after HIV infection. Third, age at
#' HIV infection has been added as extra covariable. Numbers differ slightly
#' from the \code{aidssi} data set. Some individuals were diagnosed with AIDS
#' only when they died and others had their last follow-up at AIDS diagnosis.
#' In order to prevent two transitions to happen at the same time, their time
#' to AIDS was shortened by 0.25 years. This data set was used as example for
#' the multi-state analyses in Geskus (2016).
#' 
#' 
#' @name aidssi
#' @aliases aidssi aidssi2
#' @docType data
#' @format aidssi \tabular{ll}{ patnr:\tab Patient identification number\cr
#' time:\tab Time from HIV infection to first of SI appearance and AIDS, or
#' last follow-up\cr status:\tab Event indicator; 0 = censored, 1 = AIDS, 2 =
#' SI appearance\cr cause:\tab Failure cause; factor with levels "event-free",
#' "AIDS", "SI"\cr ccr5:\tab CCR5 genotype; factor with levels "WW" (wild type
#' allele on both chromosomes),\cr \tab "WM" (mutant allele on one
#' chromosome)\cr } aidssi2 \tabular{ll}{ patnr:\tab Patient identification
#' number\cr entry.time:\tab Time from HIV infection to cohort entry. Value is
#' zero if HIV infection occurred while in follow-up.\cr aids.time:\tab Time
#' from HIV infection to AIDS, or last follow-up if AIDS was not observed\cr
#' aids.stat:\tab Event indicator with respect to AIDS, with values 0
#' (censored) and 1 (AIDS)\cr si.time:\tab Time from HIV infection to SI
#' switch, or last follow-up if SI switch was not observed\cr si.stat:\tab
#' Event indicator with respect to SI switch, with values 0 (no switch) and 1
#' (switch)\cr death.time:\tab Time from HIV infection to death, or last
#' follow-up if death was not observed\cr death.stat:\tab Event indicator with
#' respect to death; 0 = alive, 1 = dead\cr age.inf:\tab Age at HIV
#' infection\cr ccr5:\tab CCR5 genotype; factor with levels "WW" (wild type
#' allele on both chromosomes),\cr \tab "WM" (mutant allele on one
#' chromosome)\cr }
#' @references Geskus, Ronald B. (2016). \emph{Data Analysis with Competing
#' Risks and Intermediate States.} CRC Press, Boca Raton.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' @source Geskus RB (2000). On the inclusion of prevalent cases in HIV/AIDS
#' natural history studies through a marker-based estimate of time since
#' seroconversion.  \emph{Statistics in Medicine} \bold{19}, 1753--1769.
#' 
#' Geskus RB, Miedema FA, Goudsmit J, Reiss P, Schuitemaker H, Coutinho RA
#' (2003).  Prediction of residual time to AIDS and death based on markers and
#' cofactors.  \emph{Journal of AIDS} \bold{32}, 514--521.
#' @keywords datasets
NULL





#' BMT data from Klein and Moeschberger
#' 
#' A data frame of 137 rows (patients) and 22 columns. The included variables
#' are \describe{ \item{group}{ Disease group; 1 = ALL, 2 = AML Low Risk, 3 =
#' AML High Risk } \item{t1}{ Time in days to death or last follow-up }
#' \item{t2}{ Disease-free survival time in days (time to relapse, death or
#' last follow-up) } \item{d1}{ Death indicator; 1 = dead, 0 = alive }
#' \item{d2}{ Relapse indicator; 1 = relapsed, 0 = disease-free } \item{d3}{
#' Disease-free survival indicator; 1 = dead or relapsed, 0 = alive and
#' disease-free) } \item{ta}{ Time in days to Acute Graft-versus-Host Disease
#' (AGVHD) } \item{da}{ Acute GVHD indicator; 1 = Acute GVHD, 0 = No Acute GVHD
#' } \item{tc}{ Time (days) to Chronic Graft-vrsus-Host Disease (CGVHD) }
#' \item{dc}{ Chronic GVHD indicator; 1 = Chronic GVHD, 0 = No Chronic GVHD }
#' \item{tp}{ Time (days) to platelet recovery } \item{dp}{ Platelet recovery
#' indicator; 1 = platelets returned to normal, 0 = platelets never returned to
#' normal } \item{z1}{ Patient age in years } \item{z2}{ Donor age in years }
#' \item{z3}{ Patient sex; 1 = male, 0 = female } \item{z4}{ Donor sex; 1 =
#' male, 0 = female } \item{z5}{ Patient CMV status; 1 = CMV positive, 0 = CMV
#' negative } \item{z6}{ Donor CMV status; 1 = CMV positive, 0 = CMV negative }
#' \item{z7}{ Waiting time to transplant in days } \item{z8}{ FAB; 1 = FAB
#' grade 4 or 5 and AML, 0 = Otherwise } \item{z9}{ Hospital; 1 = The Ohio
#' State University, 2 = Alferd , 3 = St. Vincent, 4 = Hahnemann } \item{z10}{
#' MTX used as a Graft-versus-Host prophylactic; 1 = yes, 0 = no } }
#' 
#' @name bmt
#' @aliases bmt
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Klein and Moeschberger (1997). \emph{Survival Analysis
#' Techniques for Censored and Truncated Data}, Springer, New York.
#' @keywords datasets
NULL





#' Data from the European Society for Blood and Marrow Transplantation (EBMT)
#' 
#' A data frame of 1977 patients transplanted for CML. The included variables
#' are \describe{ \item{patid}{Patient identification number} \item{srv}{Time
#' in days from transplantation to death or last follow-up}
#' \item{srvstat}{Survival status; 1 = death; 0 = censored} \item{rel}{Time in
#' days from transplantation to relapse or last follow-up}
#' \item{relstat}{Relapse status; 1 = relapsed; 0 = censored}
#' \item{yrel}{Calendar year of relapse; factor with levels "1993-1996","
#' 1997-1999", "2000-"} \item{age}{Patient age at transplant (years)}
#' \item{score}{Gratwohl score; factor with levels "Low risk", "Medium risk",
#' "High risk"} }
#' 
#' 
#' @name EBMT year of relapse data
#' @aliases ebmt1
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @source We acknowledge the European Society for Blood and Marrow
#' Transplantation (EBMT) for making available these data. Disclaimer: these
#' data were simplified for the purpose of illustration of the analysis of
#' competing risks and multi-state models and do not reflect any real life
#' situation. No clinical conclusions should be drawn from these data.
#' @keywords datasets
NULL





#' Data from the European Society for Blood and Marrow Transplantation (EBMT)
#' 
#' A data frame of 8966 patients transplanted at the EBMT. The included
#' variables are \describe{ \item{id}{Patient identification number}
#' \item{time}{Time in months from transplantation to death or last follow-up}
#' \item{status}{Survival status; 0 = censored; 1,...,6 = death due to the
#' following causes: Relapse (1), GvHD (2), Bacterial infections (3), Viral
#' infections (4), Fungal infections (5), Other causes (6)} \item{cod}{Cause of
#' death as factor with levels "Alive", "Relapse", "GvHD", "Bacterial",
#' "Viral", "Fungal", "Other"} \item{dissub}{Disease subclassification; factor
#' with levels "AML", "ALL", "CML"} \item{match}{Donor-recipient gender match;
#' factor with levels "No gender mismatch", "Gender mismatch"}
#' \item{tcd}{T-cell depletion; factor with levels "No TCD", "TCD", "Unknown"}
#' \item{year}{Year of transplantation; factor with levels "1985-1989",
#' "1990-1994", "1995-1998"} \item{age}{Patient age at transplant; factor with
#' levels "<=20", "20-40", ">40"} }
#' 
#' 
#' @name EBMT cause of death data
#' @aliases ebmt2
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Fiocco M, Putter H, van Houwelingen JC (2005). Reduced rank
#' proportional hazards model for competing risks. \emph{Biostatistics}
#' \bold{6}, 465--478.
#' @source We acknowledge the European Society for Blood and Marrow
#' Transplantation (EBMT) for making available these data. Disclaimer: these
#' data were simplified for the purpose of illustration of the analysis of
#' competing risks and multi-state models and do not reflect any real life
#' situation. No clinical conclusions should be drawn from these data.
#' @keywords datasets
NULL





#' Data from the European Society for Blood and Marrow Transplantation (EBMT)
#' 
#' A data frame of 2204 patients transplanted at the EBMT between 1995 and
#' 1998. These data were used in Section 4 of the tutorial on competing risks
#' and multi-state models (Putter, Fiocco & Geskus, 2007). The included
#' variables are \describe{ \item{id}{Patient identification number}
#' \item{prtime}{Time in days from transplantation to platelet recovery or last
#' follow-up} \item{prstat}{Platelet recovery status; 1 = platelet recovery, 0
#' = censored} \item{rfstime}{Time in days from transplantation to relapse or
#' death or last follow-up (relapse-free survival time)}
#' \item{rfsstat}{Relapse-free survival status; 1 = relapsed or dead, 0 =
#' censored} \item{dissub}{Disease subclassification; factor with levels "AML",
#' "ALL", "CML"} \item{age}{Patient age at transplant; factor with levels
#' "<=20", "20-40", ">40"} \item{drmatch}{Donor-recipient gender match; factor
#' with levels "No gender mismatch", "Gender mismatch"} \item{tcd}{T-cell
#' depletion; factor with levels "No TCD", "TCD"} }
#' 
#' 
#' @name EBMT platelet recovery data
#' @aliases ebmt3
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
#' Competing risks and multi-state models. \emph{Statistics in Medicine}
#' \bold{26}, 2389--2430.
#' @source We acknowledge the European Society for Blood and Marrow
#' Transplantation (EBMT) for making available these data. Disclaimer: these
#' data were simplified for the purpose of illustration of the analysis of
#' competing risks and multi-state models and do not reflect any real life
#' situation. No clinical conclusions should be drawn from these data.
#' @keywords datasets
NULL





#' Data from the European Society for Blood and Marrow Transplantation (EBMT)
#' 
#' A data frame of 2279 patients transplanted at the EBMT between 1985 and
#' 1998. These data were used in Fiocco, Putter & van Houwelingen (2008), van
#' Houwelingen & Putter (2008, 2012) and de Wreede, Fiocco & Putter (2011). The
#' included variables are \describe{ \item{id}{Patient identification number}
#' \item{rec}{Time in days from transplantation to recovery or last follow-up}
#' \item{rec.s}{Recovery status; 1 = recovery, 0 = censored} \item{ae}{Time in
#' days from transplantation to adverse event (AE) or last follow-up}
#' \item{ae.s}{Adverse event status; 1 = adverse event, 0 = censored}
#' \item{recae}{Time in days from transplantation to both recovery and AE or
#' last follow-up} \item{recae.s}{Recovery and AE status; 1 = both recovery and
#' AE, 0 = no recovery or no AE or censored} \item{rel}{Time in days from
#' transplantation to relapse or last follow-up} \item{rel.s}{Relapse status; 1
#' = relapse, 0 = censored} \item{srv}{Time in days from transplantation to
#' death or last follow-up} \item{srv.s}{Relapse status; 1 = dead, 0 =
#' censored} \item{year}{Year of transplantation; factor with levels
#' "1985-1989", "1990-1994", "1995-1998"} \item{agecl}{Patient age at
#' transplant; factor with levels "<=20", "20-40", ">40"}
#' \item{proph}{Prophylaxis; factor with levels "no", "yes"}
#' \item{match}{Donor-recipient gender match; factor with levels "no gender
#' mismatch", "gender mismatch"} }
#' 
#' 
#' @name EBMT data
#' @aliases ebmt4
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
#' proportional hazards regression and simulation-based prediction for
#' multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' 
#' van Houwelingen HC, Putter H (2008). Dynamic predicting by landmarking as an
#' alternative for multi-state modeling: an application to acute lymphoid
#' leukemia data. \emph{Lifetime Data Anal} \bold{14}, 447--463.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analaysis. Chapman & Hall/CRC Press, Boca Raton.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2011). mstate: An R Package for the
#' Analysis of Competing Risks and Multi-State Models. \emph{Journal of
#' Statistical Software}, Volume 38, Issue 7.
#' @source We acknowledge the European Society for Blood and Marrow
#' Transplantation (EBMT) for making available these data. Disclaimer: these
#' data were simplified for the purpose of illustration of the analysis of
#' competing risks and multi-state models and do not reflect any real life
#' situation. No clinical conclusions should be drawn from these data.
#' @keywords datasets
NULL





#' Data preparation, estimation and prediction in multi-state models
#' 
#' Functions for data preparation, descriptives, (hazard) estimation and
#' prediction (Aalen-Johansen) in competing risks and multi-state models.
#' 
#' \tabular{ll}{ Package: \tab mstate\cr Type: \tab Package\cr Version: \tab
#' 0.2.10\cr Date: \tab 2016-12-03\cr License: \tab GPL 2.0\cr }
#' 
#' @name mstate-package
#' @aliases mstate-package mstate
#' @docType package
#' @author Liesbeth de Wreede, Marta Fiocco, Hein Putter. Maintainer: Hein
#' Putter <H.Putter@@lumc.nl>
#' @references Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
#' Competing risks and multi-state models. \emph{Statistics in Medicine}
#' \bold{26}, 2389--2430.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2010). The mstate package for
#' estimation and prediction in non- and semi-parametric multi-state and
#' competing risks models. \emph{Computer Methods and Programs in Biomedicine}
#' \bold{99}, 261--274.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2011). mstate: An R Package for the
#' Analysis of Competing Risks and Multi-State Models. \emph{Journal of
#' Statistical Software}, Volume 38, Issue 7.
#' @keywords package
NULL





#' Abnormal prothrombin levels in liver cirrhosis
#' 
#' A data frame of 488 liver cirrhosis patients from a randomized clinical
#' trial concerning prednisone treatment in these patients. The dataset is in
#' long format. The included variables are \describe{ \item{id}{Patient
#' identification number} \item{from}{Starting state} \item{to}{Receiving
#' state} \item{trans}{Transition number} \item{Tstart}{Starting time}
#' \item{Tstop}{Transition time} \item{status}{Status variable; 1=transition,
#' 0=censored} \item{treat}{Treatment; factor with levels "Placebo",
#' "Prednisone"} }
#' 
#' This data was kindly provided by Per Kragh Andersen. It was introduced in
#' Andersen, Borgan, Gill & Keiding (1993), Example 1.3.12, and used as
#' illustration for computation of transition probabilities in multi-state
#' models, see Sections IV.4 (Example IV.4.4) and VII.2 (Example VII.2.10).
#' 
#' @name Liver cirrhosis data
#' @aliases prothr
#' @docType data
#' @format A data frame, see \code{\link{data.frame}}.
#' @references Andersen PK, Borgan O, Gill RD, Keiding N (1993).
#' \emph{Statistical Models Based on Counting Processes}. Springer, New York.
#' @keywords datasets
NULL





#' Help functions for transition matrix
#' 
#' Help functions to get insight into the structure of a transition matrix.
#' 
#' Function \code{to.trans2} simply lists the transitions in \code{trans} in a
#' data frame; function \code{trans2Q} converts \code{trans} to a \code{Q}
#' matrix, the (j,k)th element of which contains the (shortest) number of
#' transitions needed to travel from the jth to the kth state; function
#' \code{absorbing} returns a vector (named if \code{trans} contains row or
#' columnc names) with the state numbers that are absorbing; function
#' \code{is.circular} returns (a Boolean) whether the transition matrix
#' specified in \code{trans} is circular or not.
#' 
#' @name transhelp
#' @aliases to.trans2 trans2Q absorbing is.circular
#' @param trans Transition matrix, for instance produced by \code{transMat}),
#' \code{trans.comprisk}, or \code{trans.illdeath}
#' @return See details.
#' @author Hein Putter <H.Putter@@lumc.nl>
#' @keywords univar
#' @examples
#' 
#' # Irreversible illness-death model
#' tmat <- trans.illdeath(c("Healthy", "Illness", "Death"))
#' tmat
#' to.trans2(tmat)
#' trans2Q(tmat)
#' absorbing(tmat)
#' is.circular(tmat)
#' # Reversible illness-death model
#' tmat <- transMat(x = list( c(2, 3), c(1, 3), c() ),
#'                  names = c("Healthy", "Illness", "Death"))
#' tmat
#' to.trans2(tmat)
#' trans2Q(tmat)
#' absorbing(tmat)
#' is.circular(tmat)
#' 
NULL



