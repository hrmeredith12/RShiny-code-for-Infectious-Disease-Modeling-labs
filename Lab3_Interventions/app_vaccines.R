# R Shiny App for exploring the impact of isolation on an epidemic using an SEIR model
# By Hannah Meredith
# Last updated: September 21, 2021

# Load packages ----
library(shiny)
library(maps)
library(mapproj)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggrepel)
library(mathjaxr)
library(ggpubr)
library(deSolve)
library(ggplot2)

# Source helper functions -----
source("helper_functions_vaccines.R")

# User interface ----
ui <- fluidPage(theme=("css/style.css"),
                shinyjs::useShinyjs(),
                htmlOutput("masthead"),
                
                navbarPage("SEIR epidemic model with Isolation",id="tabs",
                           tabPanel("Model output",value=1,id=1,
                                    sidebarLayout(
                                      sidebarPanel(
                                        sliderInput("sigma_slider",
                                                    label = HTML("Proportion of infected who isolate, <i>&sigma;</i>"),  ## average global birth rate in 2016
                                                    min = 0, max = 1, value = 0.3, step = 0.05),
                                        sliderInput("rho_slider",
                                                    label = HTML("Vaccine efficacy, <i>&rho;</i>"),  ## vaccine efficacy
                                                    min = 0, max = 1, value = 1, step = 0.05),
                                        sliderInput("eta_slider",
                                                    label = HTML("Vaccination rate, <i>&eta;</i>, (1/day)"),  ## vaccination rate
                                                    min = 0, max = 3.5E-3, value = 3.5E-4, step = 5E-4),
                                        sliderInput("beta_slider",
                                                    label = HTML("Effective contact rate, <i>&beta;</i> (1/day)"),
                                                    min = 0, max = 5, value = 0.25, step = 0.05),
                                        sliderInput("omega_slider",
                                                    label = HTML("Rate of waning immunity, <i>&omega;</i> (1/year)"),
                                                    min = 0, max = 1, value = 1/2, step = 0.005),
                                        sliderInput("time_slider",
                                                    label = HTML("Time frame (days)"),
                                                    min = 0, max = 5000, value = 1500, step = 50)
                                      ),
                                      mainPanel(h3("Explore the impact of interventions on the output of an SEIR model:"),
                                                plotOutput("plots"))
                                    )
                           ),
                           tabPanel("Equations",value=2,id=2,
                                    p("Here, we use an SEIR model to explore the impact of self-isolation when infectious and a vaccination campaign. 
                                      (S)usceptibles become (E)xposed depending on the transmission rate (beta) and the number of (I)nfected 
                                      individuals they could come into contact with. (S)usceptibles are (V)accinated at rate (eta). (E)xposed individuals 
                                      become infectious depending on the pathogen's incubation rate (lambda). A small proportion (alpha) eliminate the 
                                     infection and return to (S)usceptible. When an individual becomes infectious, a proportion of them (sigma) self-
                                     isolate (Q) where they have no contact with the rest of the population. The remaining proportion of the infectious 
                                     form the (I)nfectious population capable of infecting others. We assume that people isolate immediately upon becoming 
                                      infectious - this would be representative of an infection that develops symptoms on the same day that the 
                                     host becomes infectious. When might this assumption not be a good approximation? Infected individuals (both Q 
                                     and I) (R)ecover at the same rate (gamma) to a phase where they have immunity against reinfection. The (R)ecovered 
                                     lose their immunity and return back to (S)usceptible at a rate (omega). (V)accinated individuals are less likely to be 
                                      infected upon contact with an (I)nfected persion; however, the lower the vaccine's efficacy (rho), the more likely 
                                     a (V)accinated person is to become infected. (V)accinated individuals lose their immunity at the same rate as someone who 
                                     became infected an recovered (omega). Individuals, regardless of their category, have a natural birth/death rate (mu). "),
                                    p("Equations"),
                                    withMathJax(
                                      helpText("Susceptible $$\\frac{dS}{dt} = \\mu N + \\alpha E + \\omega (R + V) - \\frac{\\beta I}{N} S - (\\eta + \\mu) S$$"),
                                      helpText("Exposed $$\\frac{dE}{dt} = \\frac{\\beta I}{N} (S + (1 - \\rho) V) - (\\lambda + \\alpha + \\mu) E$$"),
                                      helpText("Infectious $$\\frac{dI}{dt} = (1 - \\sigma)\\lambda E - (\\gamma + \\mu) I$$"),
                                      helpText("Isolated $$\\frac{dQ}{dt} =         \\sigma\\lambda E - (\\gamma + \\mu) Q$$"),
                                      helpText("Recovered $$\\frac{dR}{dt} = \\gamma (I + Q) - (\\omega + \\mu) R$$"),
                                      helpText("Vaccinated $$\\frac{dV}{dt} = \\eta S - \\frac{\\beta I}{N} (1 - \\rho) V - (\\omega + \\mu) V$$"),
                                      helpText("Basic reproduction number $$R_0 =  \\frac{\\beta}{\\gamma}$$"),
                                      helpText("Effective reproduction number $$R_e =  R_0 \\frac{S}{N}$$")
                                    ),
                                    p("Initial Conditions"),
                                    withMathJax(
                                      helpText("$$S(0) = 9999$$"),
                                      helpText("$$E(0) = 0$$"),
                                      helpText("$$I(0) = 1$$"),
                                      helpText("$$Q(0) = 0$$"),
                                      helpText("$$R(0) = 0$$"),
                                      helpText("$$V(0) = 0$$"),
                                      helpText("$$S + E + I + Q + R + V = N$$")
                                    ),
                                    p("Parameters values. To focus on the impact of interventions, these parameters have been set:"),
                                    withMathJax(
                                      helpText("Recovery rate [1/day] $$\\lambda = 1/14$$"),
                                      helpText("Pathogen maturation rate [1/day] $$\\lambda = 1/5.5$$"),
                                      helpText("Proportion of E who recover to S $$\\alpha = 1/100$$"),
                                      helpText("Birth/death rate [1/year] $$\\mu  = 1/76$$")
                                    )
                           ),
                           tabPanel("Time course data", value = 3, id = 3,
                                    p("You can download the table seen below by clicking the Download button. 
                                      When you download, if you save it with .csv at the end (e.g. MyFileName.csv), 
                                      you can open it in Excel for further analysis needed to complete the assignment.
                                      Tip - make sure to record the parameter values (beta, gamma, etc.) you used to 
                                      create this dataset."),
                                    downloadButton('download',"Download the data"),
                                    tableOutput("table"))
                )
)


# Server logic ----
server <- function(input, output) {
  
  output$plots <- renderPlot({ 
    SEIQR.plots(alpha = 1/100,
                beta = input$beta_slider[1], 
                gamma = 1/14, 
                mu = (1/76)/365, # birth/death rate
                lambda = 1/5, # maturation rate of pathogen
                sigma = input$sigma_slider[1],
                omega = input$omega_slider[1]/365, # waning immunity rate
                rho = input$rho_slider[1], # vaccine inefficacy
                eta = input$eta_slider[1], #input$eta_slider[1], #3.5E-3, # vaccination rate
                time.dur = input$time_slider[1])
  },
  width = "auto",
  height = 800)
  
  output$table <- renderTable({
    SEIQR.table(alpha = 1/100,
                beta = input$beta_slider[1], 
                gamma = 1/14, 
                mu = (1/76)/365, # birth/death rate
                lambda = 1/5, # maturation rate of pathogen
                sigma = input$sigma_slider[1],
                omega = input$omega_slider[1]/365, # waning immunity rate
                rho = input$rho_slider[1], # vaccine inefficacy
                eta = input$eta_slider[1], #input$eta_slider[1], #3.5E-3, # vaccination rate
                time.dur = input$time_slider[1])
  })
  
  data <- reactive({
    out <- SEIQR.table(alpha = 1/100,
                       beta = input$beta_slider[1], 
                       gamma = 1/14, 
                       mu = (1/76)/365, # birth/death rate
                       lambda = 1/5, # maturation rate of pathogen
                       sigma = input$sigma_slider[1],
                       omega = input$omega_slider[1]/365, # waning immunity rate
                       rho = input$rho_slider[1], # vaccine inefficacy
                       eta = input$eta_slider[1], #input$eta_slider[1], #3.5E-3, # vaccination rate
                       time.dur = input$time_slider[1])
  })
  
  output$download <- downloadHandler(
    filename = "Vaccinate_Isolate_SEIR_output.csv",
    content = function(file) {
      write.csv(data(), file)
    }
  )
  
}

# Run app ----
shinyApp(ui, server)
