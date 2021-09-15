# R Shiny App for exploring parameters of a basic SIR Model
# By Hannah Meredith
# Last updated: August 9, 2021

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

# Source helper functions -----
source("helper_functions.R")

# User interface ----
ui <- fluidPage(theme=("css/style.css"),
                shinyjs::useShinyjs(),
                htmlOutput("masthead"),
                
                navbarPage("Basic SIR epidemic model",id="tabs",
                           tabPanel("The SIR model output",value=1,id=1,
                                    sidebarLayout(
                                      sidebarPanel(
                                        sliderInput("beta_slider",
                                                    label = HTML("Effective contact rate, <i>&beta;</i> (1/day)"),
                                                    min = 0, max = 5, value = 0.25, step = 0.05),
                                        sliderInput("gamma_slider",
                                                    label = HTML("Recovery rate, <i>&gamma;</i> (1/day)"),
                                                    min = 1/20, max = 1/2, value = 1/7, step = 0.05),
                                        sliderInput("time_slider",
                                                    label = HTML("Time frame (days)"),
                                                    min = 0, max = 500, value = 100, step = 1)
                                      ),
                                      mainPanel(h3("Explore the impact of parameters on epidemic size, duration, etc.:"),
                                                plotOutput("plots"))
                                    )
                           ),
                           tabPanel("Equations",value=2,id=2,
                                    withMathJax(
                                      helpText("Susceptible $$\\frac{dS}{dt} = - \\frac{\\beta I S}{N}$$"),
                                      helpText("Infectious $$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - \\gamma I$$"),
                                      helpText("Recovered $$\\frac{dR}{dt} = \\gamma I$$"),
                                      helpText("Recovery rate $$\\gamma = \\frac{1}{\\text{infectious period}}$$"),
                                      helpText("Basic reproduction number $$R_0 =  \\frac{\\beta}{\\gamma}$$"),
                                      helpText("Effective reproduction number $$R_e =  R_0 \\frac{S}{N}$$")
                                    ),
                                    p("Initial parameters"),
                                    withMathJax(
                                      helpText("$$S(0) = 9,999$$"),
                                      helpText("$$I(0) = 1$$"),
                                      helpText("$$R(0) = 0$$"),
                                      helpText("$$S + I + R = N = 10,000$$")
                                    )),
                           tabPanel("Time course data", value = 3, id = 3,
                                    tableOutput("table"))
                )
)


# Server logic ----
server <- function(input, output) {
  
  output$plots <- renderPlot({ 
    SIR.plots(beta = input$beta_slider[1],
              gamma = input$gamma_slider[1],
              time.dur = input$time_slider[1])
  },
  width = "auto",
  height = 800)
  
  output$table <- renderTable({
    SIR.table(beta = input$beta_slider[1], 
              gamma = input$gamma_slider[1],
              time.dur = input$time_slider[1])
  })
  
  
}

# Run app ----
shinyApp(ui, server)
