# RShiny-code-for-Infectious-Disease-Modeling-labs

Here, you can find the code used to generate the Shiny apps for the labs featured in Global Health 778 - Global Health Programming, Policy, and Response: Approaches to and Use of Infectious Disease Models. Each lab has two files: 
- app.R : this file creates the user interface, communicates user input to the functions, and updates the Shiny app after running the functions with the new user input. 
- helper_functions.R : this file has the code for running different functions. In the case of our course, it codes for solving the SIR equations over time, generating plots, and compiling tables.

If you want to learn more about Shiny apps, check out this tutorial: https://shiny.rstudio.com/tutorial/

To run the apps, a couple of things are needed: 
1. Install R: https://www.r-project.org/
2. Install RStudio: https://www.rstudio.com/
3. Install and load libraries called in app.R 
-- The ones currently called are listed at the top of the app.R file. If you don't have them already installed, use the command "install.library("libraryname")"
4. Download the files from github. Make sure you keep the relevant app.R and helper_functions.R files together. 
5. Open the app.R file in R. Make sure the working directory is set correctly. (Determine current working directory using command "getwd()". Set working directory to location with your files with the command "setwd(directory_location_here)")
6. Then click "Run App" at the top right of the Rstudy window, and the Shiny interface should automatically pop up on your screen. 
