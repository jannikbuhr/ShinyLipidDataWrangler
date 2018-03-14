
# Handling packaged installing and dependencies -------------------------------------------------------------------

#list of packages required
list.of.packages <- c("shiny","shinydashboard","tidyverse", "readxl")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)


# Load packages --------------------------------------------------------------------------------------------------------

library(shiny)
library(shinydashboard)
library(tidyverse)
library(readxl)


# UI --------------------------------------------------------------------------------------------------------------

ui <- dashboardPage(skin = "green",
                    dashboardHeader(title = "Shiny Data Wrangler for ShinyLipids"),


                    #** Sidebar ---------------------------------------------------------------------------------------------------------

                    dashboardSidebar(
                        h5("Jannik Buhr",
                           a(href= "mailto:jannik.buhr@stud.uni-heidelberg.de",
                             icon("envelope", lib = "font-awesome")
                           )),
                        hr(),
                        fileInput("files", label = "Files of raw data", multiple = T, accept = c(".txt")),
                        fileInput("meta", label = "Metadata with smaple information", multiple = T, accept = c(".txt")),
                        actionButton("go","Go!")
                    ),


                    # ** Body ------------------------------------------------------------------------------------------------------------

                    dashboardBody(
                        tableOutput("debug")
                    )


)

# SERVER ----------------------------------------------------------------------------------------------------------

server <- function(input, output){
    raw <- eventReactive(input$go,{
        input$files

        # Extract the class from filenames to use as .id in the complete dataframe
        filenames <- str_split(files, "_", simplify = T)
        names(files) <- str_replace_all(filenames[,ncol(filenames)],".txt", "")

        data <- map_df(files, read_tsv, skip = 2, .id = "class") %>%
            gather(matches("\\d+"), key = "sample_num", value = "intensity", convert = T)
    })

    meta <- eventReactive(input$go, {
        input$meta
    })





    output$debug <- renderTable({
        raw()
    })


}
