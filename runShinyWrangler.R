# If you don't have the following packages installed, the app will do so for you:
# shiny, shinydashboard, tidyverse, DT, minpack.lm, broom, modelr


# Just run these two lines
# -> You need to install the shiny package first with install.packages("shiny")

library(shiny)
runGitHub("ShinyLipidDataWrangler", "jannikbuhr", launch.browser = T)
