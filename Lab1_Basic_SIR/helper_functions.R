# R Shiny App for exploring parameters of a basic SIR Model
# By Hannah Meredith
# Last updated: August 9, 2021

## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    S[S < 0] <- 0                        ## Ensures population numbers stay positive
    I[I < 0] <- 0
    R[R < 0] <- 0
    
    N <- S + I + R                       ## Total population
    
    dS <- -beta * S * I / N              ## Change in Susceptibles  
    dI <-  beta * S * I / N  - gamma * I ## Change in Infected
    dR <-                      gamma * I ## Change in Recovered  
    dC <-  beta * S * I / N              ## Number of new cases each time step - used to generate epi curve, not part of SIR model
    
    return(list(c(dS, dI, dR, dC)))
  })
}

generate.data.frames <- function(beta, gamma, time.dur){
  ### Set parameters
  init       <- c(S = 9999, I = 1, R = 0, C = 1)       ## Set initial conditions for each compartment
  parameters <- c(beta = beta, gamma = gamma)          ## beta: infection parameter; gamma: recovery parameter
  times      <- seq(0, time.dur, by = 0.5)             ## Define time frame  

  ### Solve using ode (General Solver for Ordinary Differential Equations)
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  out <- as.data.frame(out)
  out <- round(out,2)                                   ## Round to make nicer for tables in App

  ## Calculate new cases that arise daily for Epi- curve
  out$day <- floor(out$time)
  out <- out %>% mutate(new.cases = c(init[2], diff(C, lag = 1)))
  daily.cases <- out %>% group_by(day)%>% summarise(daily.cases = sum(new.cases))


  ## Calculate R_effective
  R_0 = beta / gamma
  out <- out %>% group_by(time) %>% mutate(N = sum(S,I,R),
                                           R_effective = round(S / N * R_0,2))
  
  ## Wrangle dataframe for Time Course Plot
  time_course <- out[ , !colnames(out) %in% c("C", "N")] %>% tidyr::pivot_longer(
    cols = c("S": "R"),
    names_to = "Population",
    values_to = "Counts"
  )

  time_course$Population <- factor(time_course$Population,
                                   levels = c("S", "R", "I"))

  return(list(as.data.frame(time_course), as.data.frame(out), as.data.frame(daily.cases)))
}


SIR.table <- function(beta, gamma, time.dur){
  ## Import datasets for making table
  sir.data <- generate.data.frames(beta, gamma, time.dur)
  out <- sir.data[[2]]
  out <- subset(out, out$time%%1==0)
  colnames(out) <- c("time", "S", "I", "R", "C", "Time (days)", "New cases", "N", "R_effective")
  table.out <- out[ , c("Time (days)", "S", "I", "R", "New cases", "R_effective")]
  return(table.out)
}

SIR.plots <- function(beta, gamma, time.dur){
  ## import datasets for plotting
  sir.data <- generate.data.frames(beta, gamma, time.dur)
  time_course <- sir.data[[1]]
  out <- sir.data[[2]]
  daily.cases <- sir.data[[3]]
  
  ## Plots
  cols = c("S"= "blue", "I" = "red", "R" = "dark grey")
  
  ## Plot 1: Temporal dynamics of S, I, and R
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
  epi.curve <- ggplot(daily.cases, aes(day, daily.cases))+
    geom_bar(stat = "identity")+
    labs(x = "Time (days)", y = "New cases",
         title = "Figure 2. Epidemic curve.")+
    ylim(min = 0, max = 1.2*max(daily.cases$daily.cases))+
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
  
  # Plot 3: temporal dynamics of effective reproductive number 
  R_eff <- ggplot(out, aes(time, R_effective))+
    geom_line(size = 2)+
    labs(x = "Time (days)", y = expression(R[effective]),
         title = expression(paste("Figure 3. Temporal dynamics of ", R[effective], sep = "")))+
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
    geom_hline(yintercept = 1, linetype='dotted', size = 2)
  
  ggarrange(pop.plot, epi.curve, R_eff, 
            ncol = 1,
            align = c("v"))
  
}

### If you want to run this without running the app, uncomment the lines below. 
### Make sure the libraries called at the top of the app.R file are loaded.
### Instead of using the sliders in the app to set the parameter values, you can define the parameters directly.

# test.table <- SIR.table(beta = 0.25,
#                         gamma = 1/14,
#                         time.dur = 500)
# 
# test.plots <- SIR.plots(beta = 0.21,
#                         gamma = 1/14,
#                         time.dur = 500)
