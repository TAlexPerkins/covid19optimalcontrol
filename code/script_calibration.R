# clear the workspace
rm(list=ls())

# load libraries
library(deSolve)

# read data
# https://github.com/nytimes/covid-19-data
d = read.csv('../data/us-states_20200408.csv')
d.all = read.csv('../data/us-states.csv')

# logistic function
Lfun = function(x,h1,h2){
  (1+exp(-h1*(x-h2)))^-1
}

# functions defining RHS of state variable ODEs
stateEqns = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dS = mu - (delta + beta*(1-u(t))*(alpha*A+I+H) + iota + ifelse(t>tau_nu,nu,0)) * S
    dE = beta*(1-u(t))*(alpha*A+I+H) * (S + (1-epsilon)*V) + iota * S - (delta+rho) * E
    dA = (1-sigma)*rho * E - (delta+gamma) * A
    dI = sigma*rho * E - (delta+gamma) * I
    dH = gamma*kappa * I - (delta+eta) * H
    dV = ifelse(t>tau_nu,nu,0) * S - (delta + beta*(1-u(t))*(alpha*A+I+H)*(1-epsilon)) * V
    list(c(dS,dE,dA,dI,dH,dV))
  })
}

# function for predicting new deaths over time
Dfun = function(H,t,eta,Delta.min,Delta.max,h){
  ifelse(
    H <= Hmax(t),
    eta * H * Delta.min,
    eta * H *
      (Delta.min + (Delta.max-Delta.min) *
         (1 - exp(h * (H - Hmax(t))))))
}

# function for maximum hospital capacity
Hmax = function(t){
  parameters['Hmax'] * (1 + parameters['Hmax.slope'] * t / 365)
}



# different calibration settings
calib = expand.grid(
  R0 = c(3, 3.5, 4),
  umax = c(0.5, 0.7, 0.9))

# allocate storage for death results
mod = list()

# loop across parameter scenarios
for(jj in 1:nrow(calib)){

  # parameter values
  N = 331e6
  parameters = c(
    alpha = 0.602,
    beta = NA,
    gamma = 0.31,
    delta = 0.0116/365,
    Delta.min = 0.104,
    Delta.max = 0.28,
    epsilon = 0.8,
    eta = 0.075,
    iota = NA,
    kappa = 0.26,
    mu = 0.0116/365,
    nu = 0.00197,
    rho = 0.2,
    sigma = 0.82,
    tau_nu = 456,
    h = 2*log(0.5)/9.879154e-4,
    Hmax = 9.879154e-4,
    Hmax.slope = 0,
    umax = calib$umax[jj],
    c = NA,
    N = N,
    R0 = calib$R0[jj])
  
  alpha = parameters['alpha']
  gamma = parameters['gamma']
  delta = parameters['delta']
  kappa = parameters['kappa']
  rho = parameters['rho']
  eta = parameters['eta']
  sigma = parameters['sigma']
  R0 = parameters['R0']
  beta = R0 /
    ((1-sigma)*alpha*rho*(delta+eta)+rho*sigma*(delta+eta)+kappa*gamma*rho*sigma) *
    ((delta+gamma)*(delta+rho)*(delta+eta))
  parameters['beta'] = beta
  
  # initial conditions of state variables
  stateState = c(
    S = 1,
    E = 0,
    A = 0,
    I = 0,
    H = 0,
    V = 0)
  
  # other parameters
  days = 100
  timestep = 1
  
  # times over which to solve the model
  times = seq(32,days,by=timestep)
  
  # allow control to go down to zero at any point
  umin = rep(0,length(times))
  
  # maximum control prior to declaration of national emergency
  umax = parameters['umax'] * Lfun(times,0.3,95)
  
  # initial control function
  u = approxfun(times, umax)
  
  # simulate deaths as a function of importation detection probability
  iota = 10 ^ seq(-12,-4,length.out=300)
  deaths = numeric(length=length(iota))
  for(ii in 1:length(iota)){
    
    # set importation rate value
    parameters['iota'] = iota[ii]
    
    # solve for state variables forward in time
    outState = ode(y = stateState, times = times, func = stateEqns, parms = parameters)
    total.deaths = 
      round(sum(parameters['N'] * timestep *
                Dfun(outState[,'H'],times,parameters['eta'],
                     parameters['Delta.min'],parameters['Delta.max'],
                     parameters['h']),na.rm=T))
    # print(c(iota[ii],total.deaths))
    
    # store deaths
    deaths[ii] = total.deaths
  }
  
  # get value of importation multiplier that matches deaths
  mod[[jj]] = approxfun(log(iota,10),log(deaths,10))
}



# expand calibration results for different scenarios about death reporting
omega = rep(c(0.5,0.8),each=nrow(calib))
calib = cbind(rbind(calib,calib),omega=omega)
calib$iota = NA
for(cc in 1:nrow(calib)){
  calib$iota[cc] =
    10^optimize(f=function(par){
      (mod[[(cc-1)%%6+1]](par)-
       log(sum(d$deaths)/calib$omega[(cc-1)%%6+1],10))^2},
      interval=c(-12,-5))$minimum
}
save(parameters,calib,file='../data/calibrated_params.RData')

# R0 umax omega         iota
# 1  3.0  0.5   0.5 9.350472e-07
# 2  3.5  0.5   0.5 3.159165e-07
# 3  4.0  0.5   0.5 1.070143e-07
# 4  3.0  0.7   0.5 9.692553e-07
# 5  3.5  0.7   0.5 3.322363e-07
# 6  4.0  0.7   0.5 1.137230e-07
# 7  3.0  0.9   0.5 9.350472e-07
# 8  3.5  0.9   0.5 3.159165e-07
# 9  4.0  0.9   0.5 1.070143e-07
# 10 3.0  0.5   0.8 9.692553e-07
# 11 3.5  0.5   0.8 3.322363e-07
# 12 4.0  0.5   0.8 1.137230e-07
# 13 3.0  0.7   0.8 9.350472e-07
# 14 3.5  0.7   0.8 3.159165e-07
# 15 4.0  0.7   0.8 1.070143e-07
# 16 3.0  0.9   0.8 9.692553e-07
# 17 3.5  0.9   0.8 3.322363e-07
# 18 4.0  0.9   0.8 1.137230e-07



# plot deaths over time based on data and calibrated models
pdf('../figures/calibration.pdf',width=3.25,height=2)
  
  par(mar=c(1.75,4.1,0.5,0.5))

  plot(23:100,diff(aggregate(d.all$deaths,by=list(d.all$date),sum)[,2]),
       xlab='',ylab='Daily deaths',xlim=c(23,100),ylim=c(0,3000),las=1,
       type='l',col=2,lwd=2,xaxs='i',yaxs='i',xaxt='n')
  axis(1,at=365/12*(0:24),labels=rep('',25))
  mtext(c('February','March'),1,
        at=-365/12/2+365/12*(2:3),line=0.5)
  legend('topleft',legend=c('Reported deaths','Simulated deaths'),col=c(2,1),lwd=c(2,1),bty='n')
  
  # loop across parameter scenarios
  deaths.total = numeric()
  for(jj in 1:nrow(calib)){
    
    # parameter values
    N = 331e6
    parameters = c(
      alpha = 0.602,
      beta = NA,
      gamma = 0.31,
      delta = 0.0116/365,
      Delta.min = 0.104,
      Delta.max = 0.28,
      epsilon = 0.8,
      eta = 0.075,
      iota = calib$iota[jj],
      kappa = 0.26,
      mu = 0.0116/365,
      nu = 0.00197,
      rho = 0.2,
      sigma = 0.82,
      tau_nu = 456,
      h = 2*log(0.5)/9.879154e-4,
      Hmax = 9.879154e-4,
      Hmax.slope = 0,
      umax = calib$umax[jj],
      c = NA,
      N = N,
      R0 = calib$R0[jj])
    
    alpha = parameters['alpha']
    gamma = parameters['gamma']
    delta = parameters['delta']
    kappa = parameters['kappa']
    rho = parameters['rho']
    eta = parameters['eta']
    sigma = parameters['sigma']
    R0 = parameters['R0']
    beta = R0 /
      ((1-sigma)*alpha*rho*(delta+eta)+rho*sigma*(delta+eta)+kappa*gamma*rho*sigma) *
      ((delta+gamma)*(delta+rho)*(delta+eta))
    parameters['beta'] = beta
    
    # solve for state variables forward in time
    outState = ode(y = stateState, times = times, func = stateEqns, parms = parameters)
  
    # calculate deaths over time
    deaths = 
      parameters['N'] * timestep *
      Dfun(outState[,'H'],times,parameters['eta'],
           parameters['Delta.min'],parameters['Delta.max'],
           parameters['h'])
    lines(times, deaths * calib[jj,'omega'])
  }
  
  lines(23:100,diff(aggregate(d.all$deaths,by=list(d.all$date),sum)[,2]),
        col=2,lwd=3)

dev.off()
