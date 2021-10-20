# Helper functions for rShiny App for exploring the impact of isolation on an epidemic using an SEIR model
# By Hannah Meredith
# Last updated: Sept 21, 2021

## Create an SIR function
seiqr <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    S[S < 0] <- 0                       ## Ensures population numbers stay positive
    E[E < 0] <- 0
    I[I < 0] <- 0
    Q[Q < 0] <- 0
    R[R < 0] <- 0
    V[V < 0] <- 0
    
    N <- S + E + I + Q + R +V                 ## Total population 
    
    dS <- mu * N - (beta * I / N) * S  + alpha * E                                  - mu * S + omega * (R + V) - eta * S      ## Change in Susceptibles         
    dE <-          (beta * I / N) * S  - alpha * E - lambda * E                     - mu * E                            + (1 - rho) * (beta * I / N) * V ## Change in exposed
    dI <-                              (1 - sigma) * lambda * E  - gamma * I        - mu * I                            ## Change in Infected
    dQ <-                                  (sigma) * lambda * E  - gamma * Q        - mu * Q                            ## Change in Isolated
    dR <-                                                          gamma * (I + Q)  - mu * R - omega * R                ## Change in Recovered
    dV <-                                                                           - mu * V - omega * V      + eta * S - (1 - rho) * (beta * I / N) * V 
    dCases <- lambda * E                                                                                                    ## Number of new cases each time step - used to generate epi curve, not part of SIR model
    dVaccinations <- eta * S                                                                                                      ## Number of new vaccinations each time step 
    return(list(c(dS, dE, dI, dQ, dR, dV, dCases, dVaccinations)))
  })
}

generate.data.frames <- function(alpha, beta, gamma, mu, lambda, sigma, omega, rho, eta, time.dur){
  ### Set parameters
  init       <- c(S = 9999, E = 0, I = 1, Q = 0, R = 0, V = 0, Cases = 1, Vaccinations = 0)              ## Counts in each compartment
  parameters <- c(alpha = alpha, beta = beta, gamma = gamma, mu = mu, lambda = lambda, sigma = sigma, omega = omega, rho = rho, eta = eta, time.dur = time.dur)   ## beta: infection parameter; gamma: recovery parameter
  times      <- seq(0, time.dur, by = 0.5)                                  ## Define time frame
  
  ### Solve using ode (General Solver for Ordinary Differential Equations)
  out <- ode(y = init, times = times, func = seiqr, parms = parameters)
  out <- as.data.frame(out)
  out <- round(out,2)
  
  ## Calculate new cases that arise daily for Epi-Curve
  out$day <- floor(out$time)
  out <- out %>% mutate(new.cases = c(init[2], diff(Cases, lag = 1)))
  
  daily.counts <- out %>% group_by(day)%>% summarise(daily.cases = sum(new.cases))
  
  ## Calculate R_effective
  R_0 = beta / gamma
  out <- out %>% group_by(time) %>% mutate(N = sum(S, E, I, Q, R, V),
                                           R_effective = round(S / N * R_0,2))
  ## Wrangle dataframe for Time course plot
  time_course <- out[ , !colnames(out) %in% c("C", "N")] %>% tidyr::pivot_longer(
    cols = c("S": "V"),
    names_to = "Population",
    values_to = "Counts"
  )
  
  time_course$Population <- factor(time_course$Population,
                                   levels = c("S", "E", "R", "I", "Q", "V"))
  
  return(list(as.data.frame(time_course), as.data.frame(out), as.data.frame(daily.counts)))
}

SEIQR.table <- function(alpha, beta, gamma, mu, lambda, sigma, omega, rho, eta, time.dur){
  ## Import datasets for making table
  seiqr.data <- generate.data.frames(alpha, beta, gamma, mu, lambda, sigma, omega, rho, eta, time.dur)
  out <- seiqr.data[[2]]
  out <- subset(out, out$time%%1==0)
  colnames(out) <- c("time", "S", "E", "I", "Q", "R", "V", "C", "Cummulative Vaccinations", "Time (days)", "New cases", "N", "R_effective")
  table.out <- out[ , c("Time (days)", "S", "E", "I", "Q", "R", "V", "New cases", "Cummulative Vaccinations", "R_effective")]
  return(table.out)
}

SEIQR.plots <- function(alpha, beta, gamma, mu, lambda, sigma, omega, rho, eta, time.dur){
  ## import datasets for plotting
  
  seiqr.data <- generate.data.frames(alpha, beta, gamma, mu, lambda, sigma, omega, rho, eta, time.dur)
  time_course <- seiqr.data[[1]]
  out <- seiqr.data[[2]]
  daily.counts <- seiqr.data[[3]]
  
  ## Plots
  cols = c("S"= "blue", "E" = "gold", "I" = "red", "Q" =  "pink","R" = "dark grey", "V" = "black")
  
  ## Plot 1: Temporal dynamics of S, E, I, Q, and R
  pop.plot <- ggplot(time_course, aes(time, Counts, color = Population)) + 
    geom_line(size = 2) +
    scale_color_manual(values= cols)+
    labs(x = "Time (days)", y = "Counts", 
         title = "Figure 1. Temporal dynamics of stages.")+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18), 
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "bottom",
      plot.title = element_text(size = 20)
    )
  
  
  ## Plot 2: Epidemic curve 
  epi.curve <- ggplot(daily.counts, aes(day, daily.cases))+
    geom_bar(stat = "identity")+
    labs(x = "Time (days)", y = "New cases",
         title = "Figure 2. Epidemic curve.")+
    ylim(min = 0, max = 1.2*max(daily.counts$daily.cases))+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size = 20)
    )
  
  ## Plot 3: Vaccination curve 
  vax.curve <- ggplot(out, aes(time, Vaccinations))+
    geom_bar(stat = "identity")+
    labs(x = "Time (days)", y = "Vaccinated",
         title = "Figure 3. Cumulative Vaccination curve.")+
    ylim(min = 0, max = 1.2*max(out$Vaccinations))+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size = 20)
    )
  
  # Plot 4: temporal dynamics of effective reproductive number 
  
  R_eff <- ggplot(out, aes(time, R_effective))+
    geom_line(size = 2)+
    labs(x = "Time (days)", y = expression(R[effective]),
         title = expression(paste("Figure 4. Temporal dynamics of ", R[effective], sep = "")))+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18), 
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(size = 20)
    )+
    ylim(c(0, NA))+
    # xlim(c(0, t_end))+
    geom_hline(yintercept = 1, linetype='dotted', size = 2)
  
  ggarrange(pop.plot, epi.curve, vax.curve, R_eff, 
            ncol = 1,
            align = c("v"))
  
}


### If you want to run this without running the app, uncomment the lines below. 
### Make sure the libraries called at the top of the app.R file are loaded.
### Instead of using the sliders in the app to set the parameter values, you can define the parameters directly.

# test.table <- SEIQR.table(alpha = 1/100,
#                         beta = 0.1429, #transmission rate
#                         gamma = 0.0714, #recovery rate
#                         mu = (1/76)/365, # birth/death rate
#                         lambda = 1/5, # maturation rate of pathogen
#                         sigma = 0.3, # proportion who isolate
#                         omega = (1/10)/365, # waning immunity rate
#                         rho = 0.95, # vaccine efficacy
#                         eta = 3.5E-4, # vaccination rate
#                         time.dur = 1500)
# 
# test.plot <- SEIQR.plots(alpha = 1/100,
#                          beta = 0.1429, #transmission rate
#                          gamma = 0.0714, #recovery rate
#                          mu = (1/76)/365, # birth/death rate
#                          lambda = 1/5, # maturation rate of pathogen
#                          sigma = 0, # proportion who isolate
#                          omega = (1/10)/365, # waning immunity rate
#                          rho = 0.95, # vaccine efficacy
#                          eta = 3.5E-3, #3.5E-3, # vaccination rate
#                          time.dur = 1500)
