# clear the workspace
rm(list=ls())

# load library
library(deSolve)
library(doParallel)
library(fields)
library(foreach)

# read calibrated parameter values
load('../data/calibrated_params.RData')



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

# functions defining RHS of adjoint variable ODEs
adjointEqns = function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dlambdaS =
      lambdaS * (delta + beta*(1-u(t))*(alpha*Afun(t)+Ifun(t)+Hfun(t)) + iota + nu) -
      lambdaE * (beta*(1-u(t))*(alpha*Afun(t)+Ifun(t)+Hfun(t)) + iota) -
      lambdaV * (nu)
    dlambdaE =
      lambdaE * (delta + rho) -
      lambdaA * ((1-sigma)*rho) -
      lambdaI * (sigma*rho)
    dlambdaA = 
      lambdaS * (beta*(1-u(t))*alpha*Sfun(t)) -
      lambdaE * (beta*(1-u(t))*alpha*(Sfun(t)+(1-epsilon)*Vfun(t))) +
      lambdaA * (delta + gamma) +
      lambdaV * (beta*(1-u(t))*alpha*(1-epsilon)*Vfun(t))
    dlambdaI = 
      lambdaS * (beta*(1-u(t))*Sfun(t)) -
      lambdaE * (beta*(1-u(t))*(Sfun(t)+(1-epsilon)*Vfun(t))) +
      lambdaI * (delta + gamma) -
      lambdaH * (gamma*kappa) +
      lambdaV * (beta*(1-u(t))*(1-epsilon)*Vfun(t))
    dlambdaH =
      -dD2dH(t) +
      lambdaS * (beta*(1-u(t))*Sfun(t)) -
      lambdaE * (beta*(1-u(t))*(Sfun(t)+(1-epsilon)*Vfun(t))) +
      lambdaH * (delta + eta) +
      lambdaV * (beta*(1-u(t))*(1-epsilon)*Vfun(t))
    dlambdaV =
      -lambdaE * (beta*(1-u(t))*(alpha*Afun(t)+Ifun(t)+Hfun(t))*(1-epsilon)) +
      lambdaV * (delta + beta*(1-u(t))*(alpha*Afun(t)+Ifun(t)+Hfun(t))*(1-epsilon))
    list(c(dlambdaS,dlambdaE,dlambdaA,dlambdaI,dlambdaH,dlambdaV))
  })
}

# function defining control that maximizes Hamiltonian
control = function(t, parameters){
  with(as.list(c(parameters)),{
    u = -beta*(alpha*Afun(t)+Ifun(t)+Hfun(t)) *
      (lambdaSfun(t)*Sfun(t) + lambdaVfun(t)*(1-epsilon)*Vfun(t) - lambdaEfun(t)*(Sfun(t)+(1-epsilon)*Vfun(t))) / (2 * c)
    list(c(u))
  })
}



# specify the value of c
parameters['c'] = 1e-12

# specify jobs that need to be run
ustart = seq(-21,21,7)

# specify jobs that need to be run
jobs = expand.grid(pp = 1:nrow(calib), cc = 1:length(ustart))

# define place to store the results
toSave = list()

# run jobs
cl = makeCluster(48)
registerDoParallel(cl)
toSave = foreach(ii = 1:nrow(jobs), .packages='deSolve') %dopar% {
  pp = jobs$pp[ii]
  cc = jobs$cc[ii]

  # initial control function
  u = function(t){
    rep(0,length(t))
  }
  
  # parameter values
  parameters['iota'] = calib$iota[pp]
  parameters['umax'] = calib$umax[pp]
  parameters['R0'] = calib$R0[pp]
  parameters['omega'] = calib$omega[pp]
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
  
  # final conditions of adjoint variables
  stateAdjoint = c(
    lambdaS = 0,
    lambdaE = 0,
    lambdaA = 0,
    lambdaI = 0,
    lambdaH = 0,
    lambdaV = 0)
  
  # other parameters
  years = 2
  timestep = 1
  smoothingWindow = 20
  maxit = 2000
  weights = 1 / smoothingWindow
  
  # times over which to solve the model
  times = seq(32,years*365,by=timestep)

  # maximum control prior to declaration of national emergency
  umax = parameters['umax'] * Lfun(times,0.3,95+ustart[cc])
  
  # allow control to go down to zero at any point after April 30
  umin = c(umax[times<=(121+ustart[cc])],rep(0,sum(times>(121+ustart[cc]))))
  
  # track the objective function value across iterations
  J = rep(1e5,maxit)
  
  # iterate through the forward-backward sweep for a certain number of iterations
  umat = matrix(NA,maxit,length(times))
  iteration = 1
  repeat{
    # solve for state variables forward in time
    outState = ode(
      y = stateState, times = times, func = stateEqns,
      parms = parameters, method = 'rk4')
    
    # get functions approximating state variable solutions
    Sfun = approxfun(times, outState[,'S'], rule = 2)
    Efun = approxfun(times, outState[,'E'], rule = 2)
    Afun = approxfun(times, outState[,'A'], rule = 2)
    Ifun = approxfun(times, outState[,'I'], rule = 2)
    Hfun = approxfun(times, outState[,'H'], rule = 2)
    Vfun = approxfun(times, outState[,'V'], rule = 2)
    
    # calculate objective function
    D2 = integrate(
      function(t){(Dfun(
        Hfun(t),t,parameters['eta'],parameters['Delta.min'],
        parameters['Delta.max'],parameters['h']))^2},
      min(outState[,'time']),max(outState[,'time']),stop.on.error=F)$value
    u2 = integrate(
      function(t){u(t)^2},min(outState[,'time']),max(outState[,'time']),stop.on.error=F)$value
    J[iteration] = D2 + parameters['c'] * u2
    
    # function for second partial derivative of D with respect to H
    eta = parameters['eta']
    Delta.max = parameters['Delta.max']
    Delta.min = parameters['Delta.min']
    h = parameters['h']
    dD2dH = approxfun(
      times,
      ifelse(
        outState[,'H'] < Hmax(times),
        2 * (eta * Delta.min) ^ 2 * outState[,'H'],
        (2 / outState[,'H'] * Dfun(outState[,'H'],times,eta,Delta.min,Delta.max,h) ^ 2 -
           2 * eta * h * outState[,'H'] * Dfun(outState[,'H'],times,eta,Delta.min,Delta.max,h) *
           (Delta.max - Delta.min) * exp(h * (outState[,'H'] - Hmax(times))))),
      rule=2)
    
    # solve for adjoint variables backward in time
    outAdjoint = ode(
      y = stateAdjoint, times = rev(times), func = adjointEqns,
      parms = parameters, method = 'rk4')
  
    # get functions approximating adjoint variable solutions
    lambdaSfun = approxfun(rev(times), outAdjoint[,'lambdaS'], rule = 2)
    lambdaEfun = approxfun(rev(times), outAdjoint[,'lambdaE'], rule = 2)
    lambdaAfun = approxfun(rev(times), outAdjoint[,'lambdaA'], rule = 2)
    lambdaIfun = approxfun(rev(times), outAdjoint[,'lambdaI'], rule = 2)
    lambdaHfun = approxfun(rev(times), outAdjoint[,'lambdaH'], rule = 2)
    lambdaVfun = approxfun(rev(times), outAdjoint[,'lambdaV'], rule = 2)
    
    # get function approximating control that maximizes Halmitonian
    u = approxfun(times, pmax(umin,pmin(umax,unlist(control(times, parameters)))),rule=2)

    if(iteration < (10+2*smoothingWindow)){
      # record the control
      umat[iteration,] = u(times)
    } else {
      # record a convex combination of controls
      umat[iteration,] =
        colSums(weights * rbind(
          umat[iteration-(1:(smoothingWindow-1)),], u(times)))
      # update the control function to this convex combination
      u = approxfun(times, umat[iteration,],rule=2)
    }
    
    # exit once the control stops changing sufficiently
    if(iteration == maxit){
      break
    }
    
    # update iteration
    iteration = iteration + 1
  }

  # save outputs from this parameter combination
  list(
    umat = umat,
    J = J,
    parameters = parameters)
}

# save results
save(toSave,file='../data/sweep_objectives_ustart.RData')
stopCluster(cl)



# figure out which control best minimizes the objective and is well converged
best = foreach(ii = 1:length(toSave),.combine='rbind') %do% {
  J = toSave[[ii]]$J
  best = 1900+rev(sort(order(tail(J,100),decreasing=F)[1:10],decreasing=T))
  umat.sub = toSave[[ii]]$umat[best,]
  convergence = sapply(1:nrow(umat.sub),function(jj)
    sum(abs(umat.sub[jj,]-umat.sub[jj-1,]))/sum(umat.sub[jj-1,]))
  c(iteration=tail(best,1),objective=J[tail(best,1)],convergence=tail(convergence,1))
}



# add additional results to the saved outputs
for(ii in 1:length(toSave)){
  # retrieve the best control for this set of parameters
  jj = best[ii,'iteration']

  # times over which to solve the model
  years = 2; timestep = 1
  times = seq(32,years*365,by=timestep)

  # parameter values
  parameters = toSave[[ii]]$parameters

  # control function
  u = approxfun(times,toSave[[ii]]$umat[jj,],rule=2)

  # initial conditions of state variables
  stateState = c(
    S = 1,
    E = 0,
    A = 0,
    I = 0,
    H = 0,
    V = 0)

  # solve for state variables forward in time
  outState = ode(
    y = stateState, times = times, func = stateEqns,
    parms = parameters, method = 'rk4')

  # get functions approximating state variable solutions
  Hfun = approxfun(times, outState[,'H'], rule = 2)

  # calculate objective function
  D = integrate(
    function(t){Dfun(
      Hfun(t),t,parameters['eta'],parameters['Delta.min'],
      parameters['Delta.max'],parameters['h'])},
    min(outState[,'time']),max(outState[,'time']),stop.on.error=F)$value
  u = integrate(
    function(t){u(t)},min(outState[,'time']),max(outState[,'time']),stop.on.error=F)$value

  # add these outputs
  toSave[[ii]]$outState = outState
  toSave[[ii]]$D = D
  toSave[[ii]]$u = u
}



# plot trade-off between control and deaths as a function of c
pdf('../figures/objective_sensitivity_ustart.pdf',width=6.5,height=5.5)

  layout(matrix(1:9,3,3))

  par(oma=c(2,2,0,2),mar=c(2.5,3,2.5,3))

  R0.vec = c(3,3.5,4)
  umax.vec = c(0.5,0.7,0.9)

  for(ll in 1:3){
    for(kk in 1:3){
      which.plot = which(
        sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['R0']})==R0.vec[ll] &
        sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['umax']})==umax.vec[kk] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['omega']})==0.5)
      plot(ustart[jobs[which.plot,'cc']],
           sapply(which.plot,function(ii)toSave[[ii]]$u),type='b',pch=4,col=4,
           xlab='',ylab='',las=1,yaxt='n')
      axis(2,las=1,col=4,col.ticks=4,col.axis=4)
      which.plot = which(
        sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['R0']})==R0.vec[ll] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['umax']})==umax.vec[kk] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['omega']})==0.8)
      lines(ustart[jobs[which.plot,'cc']],
            sapply(which.plot,function(ii)toSave[[ii]]$u),type='b',pch=3,lty=3,col=4)
      par(new=T)
      which.plot = which(
        sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['R0']})==R0.vec[ll] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['umax']})==umax.vec[kk] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['omega']})==0.5)
      plot(ustart[jobs[which.plot,'cc']],
           sapply(which.plot,function(ii)toSave[[ii]]$D),type='b',pch=4,col=2,
           xaxt='n',yaxt='n',xlab='',ylab='')
      axis(4,las=1,col=2,col.ticks=2,col.axis=2)
      which.plot = which(
        sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['R0']})==R0.vec[ll] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['umax']})==umax.vec[kk] &
          sapply(1:length(toSave),function(jj){toSave[[jj]]$parameters['omega']})==0.8)
      lines(ustart[jobs[which.plot,'cc']],
            sapply(which.plot,function(ii)toSave[[ii]]$D),type='b',pch=3,lty=3,col=2)
      eval(parse(text=paste(
        'mtext(expression(R[0]*" = ',R0.vec[ll],', "*u["max"]*" = ',umax.vec[kk],'"),3,cex=0.7,line=0.25)',sep='')))
      if(kk==3 & ll==2){
        mtext('Timing of initial increase in control relative to default timing (days)',1,cex=0.7,line=2.75)
      }
      if(kk==2 & ll==1){
        mtext('Days spent fully under control, u',2,cex=0.7,line=3.25,col=4)
      }
      if(kk==2 & ll==3){
        mtext('Cumulative deaths, D',4,cex=0.7,line=3.25,col=2)
      }
    }
  }

dev.off()



# make a plot for each parameter scenario
for(ii in 1:length(toSave)){

  # retrieve info for each plot
  label = paste(paste(c('R0','umax','omega'),round(10*toSave[[ii]]$parameters[c('R0','umax','omega')]),sep='_',collapse='_'),'ustart',ustart[jobs[ii,'cc']],sep='_')
  outState = toSave[[ii]]$outState
  parameters = toSave[[ii]]$parameters
  u = approxfun(times,toSave[[ii]]$umat[best[ii,'iteration'],],rule=2)

  # plot outputs
  png(file=paste('../figures/output_',label,'_ustart.png',sep=''),width=3.25,height=3.25,units='in',res=200)

    layout(1:3)
    par(mar=c(2,5,1,0.5),oma=c(0.5,0,0.25,0))

    ## Plot susceptibles
    plot(outState[,'time'], outState[,'S'], type = "l", lwd = 1, ylim = c(0,1.0),
         xaxs = "i", yaxs = "i", las = 1,
         col = "#7570b3", ylab = "", xlab = "", xaxt='n')
    eval(parse(text=paste("mtext(expression(
      R[0]*' = '*",toSave[[ii]]$parameters['R0'],"*', '*
        u['max']*' = '*",toSave[[ii]]$parameters['umax'],"*', '*
        omega*' = '*",toSave[[ii]]$parameters['omega'],"*', '*
        u['start']*' = '*",ustart[jobs[ii,'cc']],"),
      3,cex=0.7)",sep='')))
    mtext('Susceptibles, S(t)',2,cex=0.7,line=4)
    polygon(c(0,121,121,0),c(0,0,1,1),
            border=NA,col="#c0c0c095")
    polygon(x = c(outState[,'time'],rev(outState[,'time'])),
            y = c(rep(0,nrow(outState)), rev(outState[,'S'])),
            border=NA,col="#7570b380")
    polygon(x = c(outState[,'time'],rev(outState[,'time'])),
            y = c(outState[,'S']+outState[,'V'], rev(outState[,'S'])),
            border=NA,col=rgb(217/255,95/255,2/255,0.8))
    lines(outState[,'time'], outState[,'S']+outState[,'V'], type = "l", lwd = 1,
          col = rgb(217/255,95/255,2/255,0.8))
    lines(outState[,'time'], outState[,'S'], type = "l", lwd = 1,
          col = "#7570b3")
    abline(v=365)
    arrows(parameters['tau_nu'],1/5*1.0,parameters['tau_nu'],0,length=0.05,lwd=1)
    axis(1,at=365/12*(0:24),labels=rep('',25))
    mtext(rep(c('J','F','M','A','M','J','J','A','S','O','N','D'),2),1,
          at=-365/12/2+365/12*(1:24),line=0.5,cex=0.7)
    mtext(c(2020,2021),1,at=365/2+c(0,1)*365,line=1.5,cex=0.7)

    ## Plot hospitalizations
    hosp = outState[,'H']
    hospmax = Hmax(times)
    plot(outState[,'time'], hosp, type = "l", lwd = 1, las=1,xaxt='n',
         ylab = "", xlab = "",xaxs="i",yaxs="i", col = rgb(231/255,41/255,138/255,0.3), ylim = c(0,1.02*max(c(hosp,hospmax))))
    mtext('Hospitalizations, H(t)',2,cex=0.7,line=4)
    polygon(c(0,121,121,0),c(0,0,1,1),
            border=NA,col="#c0c0c095")
    polygon(c(outState[,'time'],rev(outState[,'time'])),
            c(rep(0,nrow(outState)),
              rev(hosp)),
            border=NA,col=rgb(231/255,41/255,138/255,0.3))
    lines(outState[,'time'], hosp, lwd = 1,col = rgb(231/255,41/255,138/255,1))
    abline(v=365)
    lines(times,hospmax,lwd=1,lty=2)
    arrows(parameters['tau_nu'],1/5*1.02*max(c(hosp,hospmax)),parameters['tau_nu'],0,length=0.05,lwd=1)
    axis(1,at=365/12*(0:24),labels=rep('',25))
    mtext(rep(c('J','F','M','A','M','J','J','A','S','O','N','D'),2),1,
          at=-365/12/2+365/12*(1:24),line=0.5,cex=0.7)
    mtext(c(2020,2021),1,at=365/2+c(0,1)*365,line=1.5,cex=0.7)

    ## Plot control
    plot(-100,-100,xlim=range(times),
         ylim=c(0,1),type='l',xlab="",las=1,
         ylab='',xaxs='i',yaxs='i', col = "#1b9e77", lwd = 1,xaxt='n')
    mtext('Control, u(t)',2,cex=0.7,line=4)
    polygon(c(0,121,121,0),c(0,0,1,1),
            border=NA,col="#c0c0c095")
    polygon(c(times,rev(times)),
            c(rep(0,length(times)),
              rev(u(times))),
            border=NA,col="#1b9e7750")
    lines(times, u(times), lwd = 1, col = "#1b9e77")
    abline(v=365)
    arrows(parameters['tau_nu'],1/5*1,parameters['tau_nu'],0,length=0.05,lwd=1)
    axis(1,at=365/12*(0:24),labels=rep('',25))
    mtext(rep(c('J','F','M','A','M','J','J','A','S','O','N','D'),2),1,
          at=-365/12/2+365/12*(1:24),line=0.5,cex=0.7)
    mtext(c(2020,2021),1,at=365/2+c(0,1)*365,line=1.5,cex=0.7)

  dev.off()
}



# convert -delay 30 output_R0_40_umax_9_omega_8_ustart_-21_ustart.png output_R0_40_umax_9_omega_8_ustart_-14_ustart.png output_R0_40_umax_9_omega_8_ustart_-7_ustart.png output_R0_40_umax_9_omega_8_ustart_0_ustart.png output_R0_40_umax_9_omega_8_ustart_7_ustart.png output_R0_40_umax_9_omega_8_ustart_14_ustart.png output_R0_40_umax_9_omega_8_ustart_21_ustart.png output_R0_40_umax_9_omega_8_ustart.gif
