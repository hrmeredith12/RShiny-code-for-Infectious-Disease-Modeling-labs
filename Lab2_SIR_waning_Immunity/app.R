# R Shiny App for exploring parameters of a SIR Model with births, deaths, and waning immunity
# By Hannah Meredith
# Last updated: August 11, 2021

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
source("helper_functions.R")

# User interface ----
ui <- fluidPage(theme=("css/style.css"),
                shinyjs::useShinyjs(),
                htmlOutput("masthead"),
                
                navbarPage("SIR epidemic model with births, deaths, and waning immunity ",id="tabs",
                           tabPanel("Model output",value=1,id=1,
                                    sidebarLayout(
                                      sidebarPanel(
                                        sliderInput("beta_slider",
                                                    label = HTML("Effective contact rate, <i>&beta;</i> (1/day)"),
                                                    min = 0, max = 5, value = 0.25, step = 0.05),
                                        sliderInput("gamma_slider",
                                                    label = HTML("Recovery rate, <i>&gamma;</i> (1/day)"),
                                                    min = 1/20, max = 1/2, value = 1/7, step = 0.05),
                                        sliderInput("mu_slider",
                                                    label = HTML("Death/birth rate, <i>&mu;</i> (1/year)"),  ## average global birth rate in 2016
                                                    min = 0, max = 1/20, value = 1/76, step = 0.001),
                                       sliderInput("omega_slider",
                                                    label = HTML("Rate of waning immunity, <i>&omega;</i> (1/year)"),
                                                    min = 0, max = 1, value = 1/2, step = 0.005),
                                       sliderInput("time_slider",
                                                   label = HTML("Time frame (days)"),
                                                   min = 0, max = 5000, value = 1500, step = 50)
                                      ),
                                      mainPanel(h3("Explore the impact of life expectancy and waning immunity the output of an SIR model:"),
                                                plotOutput("plots"))
                                    )
                           ),
                           tabPanel("Equations",value=2,id=2,
                                    withMathJax(
                                      helpText("Susceptible $$\\frac{dS}{dt} = \\mu N + \\omega R - (\\frac{\\beta I}{N} + \\mu) S$$"),
                                      helpText("Infectious $$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - (\\gamma + \\mu) I$$"),
                                      helpText("Recovered $$\\frac{dR}{dt} = \\gamma I - (\\omega + \\mu) R$$"),
                                      helpText("Recovery rate $$\\gamma = \\frac{1}{\\text{infectious period}}$$"),
                                      helpText("Basic reproduction number $$R_0 =  \\frac{\\beta}{\\gamma}$$"),
                                      helpText("Effective reproduction number $$R_e =  R_0 \\frac{S}{N}$$"),
                                      helpText("Waning immunity rate $$\\omega = \\frac{1}{\\text{immunity duration}}$$")
                                    ),
                                    p("Initial parameters"),
                                    withMathJax(
                                      helpText("$$S(0) = 9999$$"),
                                      helpText("$$I(0) = 1$$"),
                                      helpText("$$R(0) = 0$$"),
                                      helpText("$$S + I + R = N$$")
                                    )),
                           tabPanel("Time course data", value = 3, id = 3,
                                    tableOutput("table"))
                )
)


# Server logic ----
server <- function(input, output) {
  
  output$plots <- renderPlot({ 
    SIR.plots(beta = input$beta_slider[1],
              gamma = input$gamma_slider[1],     # take inverse 
              mu = (input$mu_slider[1])/365,     # take inverse and convert to days
              omega = (input$omega_slider[1])/365,  # take inverse and convert to days
              time.dur = input$time_slider[1])
  },
  width = "auto",
  height = 800)
  
  output$table <- renderTable({
    SIR.table(beta = input$beta_slider[1],
              gamma = input$gamma_slider[1],     # take inverse and convert to days
              mu = (input$mu_slider[1])/365,     # take inverse and convert to days
              omega = (input$omega_slider[1])/365,  # take inverse and convert to days  
              time.dur = input$time_slider[1])
  })
}

# Run app ----
shinyApp(ui, server)
