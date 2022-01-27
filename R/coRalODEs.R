################################################################################
#' Solve coral-\emph{Symbiodinium} ODEs

#' Function solveCoral uses the deSolve package to generate values similar to
#' run_coral from the coRal package.
#'
#' @param times The time values at which output is desired. Default c(0,500).
#' @param pars The paramaters of the model, given as a list. See function
#' 	\code{\link{defPars}}. Default defPars().
#' @param lambda The sensitivity of the runs. High values are more sensitive
#' 	and small values are less sensitive. Default 5.
#' @param method The character method argument of ode() from deSolve desired.
#' 	Default "vode".
#' @param func The function to find the rates of change of the state variables.
#'  Default coralODEs.
#' @param ... Any other arguments to be passed to ode().
#' @return Matrix of values for fluxes, biomass, and host growth rate at
#' 	explicitly desired time values.
#' @examples
#' solveCoral()
#' solveCoral(times = seq(0,365,0.1), pars = defPars(), lambda = 10,atol = 0.01,
#' 	rtol = 0.01)
#' @seealso \code{\link{defPars}}, \code{\link{initState}},
#'  \code{\link{coralODEs}}, \code{\link{coralODEs.noROS}}
#' @import deSolve
#' @export
solveCoral <- function(times = c(0,500), pars = defPars(), lambda = 5,
		       method = "vode", func = coralODEs,
		       init = initState, ...) {
  # Solve system and return
  return(ode(y = init(pars), times = times, func = func,
	     parms = append(pars, list(lambda = lambda)), method = method, ...))
}  # End function solveCoral

#' Wrapper function to coRal
#'
#' Has the same output structure as run_coral() from the coRal package
#' @param time The time values at which output is desired.
#' @param env list of functions of time describing the environment, can
#'  be obtained with initEnv
#' @param pars list of parameters in the structure of def_pars from coRal
#' @examples
#' time <- seq(0,365,0.1)
#' env <- initEnv(time, c(30,30,0), c(1e-6,1e-6,0), c(1e-7,1e-7,0))
#' pars <- coRal::def_pars()
#' run_coral_ode(time, env, pars)
#' @seealso \code{\link{initEnv}}
#' @import deSolve
#' @export
run_coral_ode <- function(time, env, pars) {
	# time is equivalent to times
  	# can only be used for constant env values at this point
  	newPars <- list(vals = c(pars,jHG0 = 0.25, jeC0 = 10, cROS0 = 1),
			  env = env)
  	run <- solveCoral(times = time, pars = newPars)
  	with(as.data.frame(run),{
    		data.frame(time,L,N,X,jX,jN,rNH,rhoN,jeC,jCO2,jHG,rCH,
                      dH.Hdt=dH.Hdt,H,rNS,jL,jCP,jeL,jNPQ,jSG,rhoC,jST,rCS,
                      cROS,dS.Sdt=dS.Sdt,S)
  	})
}

# Helper methods
# ==============

#' Determine the rates of change of the state variables
#'
#' Default func argument for solveCoral.
#'
#' @export
coralODEs <- function(t, y, pars) {
	# I think the way we're looking at aux functions right now requires 1
	# symb

	# Target Fluxes: jCP and rhoC

	# auxOrder <- c('jL',jeL','jNPQ','cROS','jHG','jeC','rCH','rhoN',
	# 'jCO2','jSG','rCS')

	# Test necessity of max(0,stuff)

	names(y) <- c('jCP','rhoC','H','S') # Need to establish the names and
	#				    # order of the state variables
	othervars <- c('L', 'N', 'X', 'jX', 'jN', 'rNH', 'rhoN', 'jeC', 'jCO2',
		       'jHG', 'rCH', 'dH.Hdt', 'rNS', 'jL', 'jeL', 'jNPQ',
		       'jSG', 'jST', 'rCS', 'cROS', 'dS.Sdt')

	aims <- with(as.list(c(y, pars$vals,L=pars$env$L(t), N=pars$env$N(t),
			       X=pars$env$X(t))), {
		# Perform path through flux graph
    		jX <- mmk(X, KX, jXm)
		jNm <- jNm
      		jN <- mmk(N, KN, jNm)
      		jHT <- jHT0
      		rNH <- jHT0 * nNH * sigmaNH
      		rNS <- jST0 * nNS * sigmaNS
		jL <- (1.26 + 1.39 * exp(-6.48 * S / H)) * L * astar
		jeL <- max(0, jL - jCP / yCL)
		jNPQ <-((kNPQ)^-1 + (jeL)^-1)^-1
		cROS <- 1 + (max(0, jeL - jNPQ) / kROS)^k
		jHG <- synth(jHGm, yC * (rhoC * S / H + jX),
			     (jN + nNX * jX + rNH) / nNH)
		jeC <- max(0, jX + rhoC * S / H - jHG / yC)
		rCH <- sigmaCH * (jHT + (1 - yC) * jHG / yC)
		rhoN <- max(0, jN + nNX * jX + rNH - nNH * jHG)
		jCO2 <- kCO2 * jeC
		jSG <- synth(jSGm, yC * jCP, (rhoN * H / S + rNS) / nNS)
		rCS <- sigmaCS * (jST0 + (1 - yC) * jSG / yC)
		jCP <- synth(jCPm, yCL * jL, (jCO2 + rCH) * H / S + rCS) / cROS
		rhoC <- max(0, jCP - jSG / yC)
		jST <- jST0 * (1 + b * (cROS - 1))
		dS <- (jSG - jST) * S
		dH <- (jHG - jHT) * H
		# final product is a vector containing the aims of all fluxes
		# and the deltas for S and H
		c(jCP = jCP, rhoC = rhoC, dS = dS, dH = dH,
		  L=L, N=N, X=X, jX=jX, jN=jN, rNH=rNH, rhoN=rhoN, jeC=jeC,
		  jCO2=jCO2, jHG=jHG, rCH=rCH, dH.Hdt=dH / H, rNS=rNS, jL=jL,
		  jeL=jeL, jNPQ=jNPQ, jSG=jSG, jST=jST, rCS=rCS, cROS=cROS,
		  dS.Sdt = dS / S)
	})

	auxes <- c('jCP','rhoC') # names of target fluxes, sees repeated use

	# Return a list for which the first element is a vector containing the
	# deltas for each state variable jCP, rhoC, H, and S. The second element
	# of the list is a vector containing all other fluxes in the model.
	# eq(...) produces the deltas for jCP and rhoC
	# the deltas for H and S are already stored in aims vector
	return(list(c(eq(y[auxes], aims[auxes], pars$lambda),
		      aims[c('dH','dS')]), aims[othervars]))
}  # end function coralODEs

#' Determine the rates of change of the state variables in a system without ROS
#'
#' Can be used as func argument in solveCoral
#'
#' @export
coralODEs.noROS <- function(t, y, pars) {

  # Target Fluxes: jCP and rhoC

  # auxOrder <- c('jL',jeL','jNPQ','cROS','jHG','jeC','rCH','rhoN',
  # 'jCO2','jSG','rCS')

  names(y) <- c('jCP','rhoC','H','S') # Need to establish the names and
  #				    # order of the state variables
  othervars <- c('L', 'N', 'X', 'jX', 'jN', 'rNH', 'rhoN', 'jeC', 'jCO2',
                 'jHG', 'rCH', 'dH.Hdt', 'rNS', 'jL', 'jeL', 'jNPQ',
                 'jSG', 'jST', 'rCS', 'dS.Sdt')

  aims <- with(as.list(c(y, pars$vals,L=pars$env$L(t), N=pars$env$N(t),
                         X=pars$env$X(t))), {
                           # Perform path through flux graph
                           jX <- mmk(X, KX, jXm)
                           jNm <- jNm
                           jN <- mmk(N, KN, jNm)
                           jHT <- jHT0
                           rNH <- jHT0 * nNH * sigmaNH
                           rNS <- jST0 * nNS * sigmaNS
                           jL <- (1.26 + 1.39 * exp(-6.48 * S / H)) * L * astar
                           jeL <- max(0, jL - jCP / yCL)
                           jNPQ <-((kNPQ)^-1 + (jeL)^-1)^-1
                           cROS <- 1 + (max(0, jeL - jNPQ) / kROS)^k
                           jHG <- synth(jHGm, yC * (rhoC * S / H + jX),
                                        (jN + nNX * jX + rNH) / nNH)
                           jeC <- max(0, jX + rhoC * S / H - jHG / yC)
                           rCH <- sigmaCH * (jHT + (1 - yC) * jHG / yC)
                           rhoN <- max(0, jN + nNX * jX + rNH - nNH * jHG)
                           jCO2 <- kCO2 * jeC
                           jSG <- synth(jSGm, yC * jCP, (rhoN * H / S + rNS) / nNS)
                           rCS <- sigmaCS * (jST0 + (1 - yC) * jSG / yC)
                           jCP <- synth(jCPm, yCL * jL, (jCO2 + rCH) * H / S + rCS)
                           rhoC <- max(0, jCP - jSG / yC)
                           jST <- jST0
                           dS <- (jSG - jST) * S
                           dH <- (jHG - jHT) * H
                           # final product is a vector containing the aims of all fluxes
                           # and the deltas for S and H
                           c(jCP = jCP, rhoC = rhoC, dS = dS, dH = dH,
                             L=L, N=N, X=X, jX=jX, jN=jN, rNH=rNH, rhoN=rhoN, jeC=jeC,
                             jCO2=jCO2, jHG=jHG, rCH=rCH, dH.Hdt=dH / H, rNS=rNS, jL=jL,
                             jeL=jeL, jNPQ=jNPQ, jSG=jSG, jST=jST, rCS=rCS,
                             dS.Sdt = dS / S)
                         })

  auxes <- c('jCP','rhoC') # names of target fluxes, sees repeated use

  # Return a list for which the first element is a vector containing the
  # deltas for each state variable jCP, rhoC, H, and S. The second element
  # of the list is a vector containing all other fluxes in the model.
  # eq(...) produces the deltas for jCP and rhoC
  # the deltas for H and S are already stored in aims vector
  return(list(c(eq(y[auxes], aims[auxes], pars$lambda),
                aims[c('dH','dS')]), aims[othervars]))
}  # end function coralODEs.noROS

# Synthesis rate of a product given maximum rate m and substrates A and B.
synth <- function(m,A,B) {A*B * (A + B)*m /
	(A^2 * B + A * B^2 + A^2 * m + A*B * m + B^2 * m)}

# Uptake of a substrate A given half-saturation constant half and
# maximum rate max.
mmk <- function(A, half, max) {max * A / (A + half)}

# Form of the ODEs
eq <- function(x, aim, lambda) {lambda * (aim - x)} # dy.dt = lambda * (aim - y)

#' Default parameters for solveCoral
#'
#' Helper function, based off of def_pars() from package coRal.
#' @examples defPars()
#' @return Named list of model parameters. Includes environment.
#' @seealso \code{\link{solveCoral}}
#' @export
defPars <- function() {
	# comments with units and stuff from the coRal package
  return(list(vals = c(
    jHT0=0.03,  # Host specific biomass turnover rate (d^-1)
    nNH=0.18,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCH=0.1,  # Proportion of host carbon turnover recycled (-)
    jXm=0.13,  # Maximum specific host feeding rate (molX/CmolH/d)
    jNm=0.035,  # Maximum specific host DIN uptake rate (molN/CmolH/d)
    jHGm=1,  # Maximum specific host growth rate (CmolH/CmolH/d)
    jHG0 = 0.25,  # initial host growth rate
    jeC0 = 10,  # initial excess carbon flux
    kCO2=10,  # Rate of host CCM's (molCO2/molC/d)
    KN=1.5e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    KX=1e-6,  # Half-saturation constant for host feeding (CmolX/L)
    initH=1,  # Initial host biomass (CmolH)
    yC=0.8,
    jST0=0.03,  # Symbiont specific biomass turnover rate (d^-1)
    nNS=0.13,  # N:C ratio in symbiont biomass (-)
    yCL=0.1,  # L:C ratio in fixed carbon (=quantum yield) (molC/mol ph)
    kNPQ=112,  # capacity of non-photochemical quenching (mol ph/CmolS/d)
    # calculated as 4x max. photochemical quenching (Gorbunov et al. 2001)
    kROS=80,  # amount of excess light beyond NPQ capacity (e.g., jeL-jNPQ)
    	      # that doubles ROS production relative to
    	      # baseline (mol ph/CmolS/d)
    cROS0 = 1,  # capacity for NPQ
    k=1,  # exponent on ROS production (-)
    astar=1.34,  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=0.9,  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=0.9,  # Proportion of symbiont carbon turnover recycled (-)
    jCPm=2.8,  # Maximum specific photosynthate production rate (Cmol/CmolS/d)
    jSGm=0.25,  # Maximum specific symbiont growth rate (CmolS/CmolS/d)
    initS=1,  # Initial symbiont biomass (CmolS)
    b=5  # Scaling parameter for bleaching response
  ),
   env = c(
    L=Vectorize(function(x) return(30)),  # Environmental light (mol photons)
    N=Vectorize(function(x) return(1e-6)),  # Environmental DIN (mol N)
    X=Vectorize(function(x) return(1e-7))  # Environmental prey (mol X)
	   )))
}  # End function defpars

#' Create initial state for solveCoral with target fluxes jCP and rhoC
#'
#' Helper function, based off of run_coral() from the coRal package.
#' @param pars A list of model parameters in the same structure as defPars().
#' @return A named numeric vector containing the initial flux and biomass
#' 	values.
#' @seealso \code{\link{solveCoral}}, \code{\link{defPars}}
#' @examples initState(defPars())
#' @export
initState <- function(pars) {
  with(as.list(c(pars$vals,
	       L=pars$env$L(0), N=pars$env$N(0), X=pars$env$X(0))), {
    # Initial Host fluxes
    rhoN <- mmk(N, KN, jNm)
    jeC <- jeC0
    jCO2 <- kCO2 * jeC
    jHG <- jHG0
    rCH <- jHT0 * sigmaCH
    dH.Hdt <- jHGm
    H <- initH
    # Initial Symbiont fluxes
    jL <- L * astar
    jCP <- max(0, synth(jL * yCL, jCO2 * H / initS, jCPm))
    jeL <- max(0, jL - jCP / yCL)
    jNPQ <- kNPQ
    jSG <- jSGm/10
    rhoC <- jCP
    jST <- jST0
    rCS <- jST0 * sigmaCS
    cROS <- cROS0
    dS.Sdt <- jSGm
    S <- initS
    # Return named numeric vector
    return(c(# Initial flux Values
      jCP = jCP,
      rhoC = rhoC,
      # Host and Symbiont biomass
      H = H,
      S = S))
  })
}  # End function initState


#' Create functions of the time for the environment, for use with run_coral_ode
#'
#' @param times same as times for run_coral_ode
#' @param L same as L for init_env from coRal
#' @param N same as N for init_env from coRal
#' @param X same as X for init_env from coRal
#' @return list of functions of t
#' @examples init_env(seq(0,365,0.1), c(30,30,0), c(1e-6,e-6,0), c(1e-7,1e-7,0))
#' @export
# Want to use this function for the wrapper, return functions for L, N, i
# and X at an arbitrary time t
initEnv <- function(times, L, N, X) {

  # Light
  Lf <- function(L) {
	  if (L[3]==0) {
	  # linear from 1 to 2
	  function(t) (L[2]-L[1])/(times[length(times)] - times[1]) * t +
		  L[1] - (L[2]-L[1])/(times[length(times)] - times[1])*times[1]
  } else if (L[3]==1) {
	  # linear from 2 to 1
  	  function(t) (L[1]-L[2])/(times[length(times)] - times[1]) * t +
		  L[2] - (L[2]-L[1])/(times[length(times)] - times[1])*times[1]
  } else {
    	  function(t) 0.5*(L[2]-L[1]) * sin(0.0172*t) + L[1] + 0.5 * (L[2]-L[1])
  }
	}
  # DIN
  Nf <- function(N) {
	  if (N[3]==0) {
	  # linear from 1 to 2
	  function(t) (N[2]-N[1])/(times[length(times)] - times[1]) * t +
		  N[1] - (N[2]-N[1])/(times[length(times)] - times[1])*times[1]

  } else if (N[3]==1) {
	  # linear from 2 to 1
	  function(t) (N[1]-N[2])/(times[length(times)] - times[1]) * t +
		  N[2] - (N[1]-N[2])/(times[length(times)] - times[1])*times[1]
  } else {
    	  function(t) 0.5*(N[2]-N[1]) * sin(0.0172*t) + N[1] + 0.5 * (N[2]-N[1])
  }
	}
  # Prey
  Xf <- function(X) {
	  if (X[3]==0) {
	  # linear from 1 to 2
	  function(t) (X[2]-X[1])/(times[length(times)] - times[1]) * t +
		  X[1] - (X[2]-X[1])/(times[length(times)] - times[1])*times[1]
  } else if (X[3]==1) {
	  # linear from 2 to 1
	  function(t) (X[1]-X[2])/(times[length(times)] - times[1]) * t +
		  X[2] - (X[1]-X[2])/(times[length(times)] - times[1])*times[1]
  } else {
    	  function(t) 0.5*(X[2]-X[1]) * sin(0.0172*t) + X[1] + 0.5 * (X[2]-X[1])
  }
  }
  # Set environment specifications
  return(list(L=Lf(L), N=Nf(N), X=Xf(X)))
}
