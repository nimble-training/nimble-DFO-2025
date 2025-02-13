# nimble-DFO-2025

This repository contains materials for the NIMBLE workshop for the 2025 Statistical Computing in NIMBLE Fisheries and Oceans workshop.

This is scheduled as a 3-day workshop from 9am-1230pm. The schedule will provide plenty of time for breaks.

To prepare for the workshop:

 - Install NIMBLE (see below)
 - Install additional packages (see below)
 - Download the materials provided in this repository

All materials for the workshop will be in this GitHub repository. If you're familiar with Git/GitHub, you already know how to get all the materials on your computer. If you're not, go [here](https://github.com/nimble-training/nimble-DFO-2025), click the (green) "Code" button, and choose the "Download ZIP" option.

## Background for the workshop

This workshop will focus on the `nimble` R package, not on statistical methodology per se.  The material assumes attendees have basic knowledge of fisheries and ecological statistical models (e.g., basic stock recruit models and generalize linear mixed models). You will still be able to follow the workshop without this background, but the workshop is geared towards participants already familiar with standard statistical topics.

## Tentative Schedule (Under construction)

Tuesday February 18th:

1. (9:00 am - 9:30 am) Introductions and quick review of module 0, basic introduction to NIMBLE.
2. (9:30 am - 10:20 am) Introduction to Bayesian Statistics: Basic theory 
3. (10:30 am - 11:20 am) Introduction to NIMBLE
4. (11:30 am - 12:30 pm) Continued Introduction to NIMBLE

Wednesday February 19th:

7. (9:00 am - 10:20 am) User defined distributions.
8. (10:30 am - 11:50 am) Improving MCMC sampling
9. (12:00 pm - 12:30 pm) Posterior Predictions

Thursday February 20th:

10. (9:00 am - 9:50 am) Special Topics - Bayesian nonparametrics
11. (10:00 am - 10:50 am) Special Topics - Bayesian variable selection
12. (11:00 am - 11:50 am) Special Topics - Numerical methods
13. (12:00 pm - 12:30 pm) If time, writing your own MCMC sampler

## Help with NIMBLE

The NIMBLE web site is [here](https://r-nimble.org).

The NIMBLE user manual is [here](https://r-nimble.org/html_manual/cha-welcome-nimble.html).

A NIMBLE "cheatsheet" is available [here](https://r-nimble.org/documentation).

## Installing NIMBLE

NIMBLE is an R package available on CRAN, so in general it will be straightforward to install as with any R package. However, NIMBLE does require a compiler and related tools installed on your system.

The steps to install NIMBLE are:

1. Install compiler tools on your system. [https://r-nimble.org/download](https://r-nimble.org/download) will point you to more details on how to install *Rtools* on Windows and how to install the command line tools of *Xcode* on a Mac. Note that if you have packages requiring a compiler (e.g., *Rcpp*) on your computer, you should already have the compiler tools installed.

2. Install the *nimble* package from CRAN in the usual fashion of installing an R package (e.g. `install.packages("nimble")`). More details (including troubleshooting tips) can also be found in Section 4 of the [NIMBLE manual](https://r-nimble.org/html_manual/cha-installing-nimble.html).

3) Test that the installation is working, by running the following code in R:

```
library(nimble)
code <- nimbleCode({ x ~ dnorm(0, 1) })
model <- nimbleModel(code)
cModel <- compileNimble(model)
```

If the above code runs without error, you're all set. If not, please see the troubleshooting tips.  The most common problems have to do with proper installation of the compiler tools.  On Windows, the `PATH` system variable must be set (see link to Rtools installation details from our download linked above).  On Mac OSX, command line tools must be installed as part of Xcode.  If you still have problems, please email the [nimble-users group](https://r-nimble.org/more/issues-and-groups) for help.

In general, we encourage you to update to the most recent version of NIMBLE (version 1.2.1).

## Installing additional packages

Prior to the workshop, you should also install the following R packages (beyond those automatically installed with `nimble`), which can be installed as follows:

```
install.packages(c("nimbleHMC", "mcmcplots", "coda", "nimbleEcology", "compareMCMCs"))
```

We will use some other packages to set up various examples. To be able to run everything in the workshop material, you will also want:

```
install.packages(c("CARBayesdata","sp","spdep","classInt", "glmmTMB"))
```
