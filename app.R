
# Handling packaged installing and dependencies -------------------------------------------------------------------

#list of packages required
list.of.packages <- c("shiny","shinydashboard","tidyverse", "readxl","DT")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)


# Load packages --------------------------------------------------------------------------------------------------------

library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(readxl)

# Make console output less cluttered
options(readr.num_columns = 0)


# Functions -------------------------------------------------------------------------------------------------------

# For DAGs
DAG_fun <- function(data, x) {
    data %>%
        slice(x) %>%
        mutate(
            C1 = first(C),
            C2 = last(C),
            D1 = first(D),
            D2 = last(D)
        ) %>%
        select(-C, -D,-index) %>%
        mutate_at(.vars = vars(starts_with("Sample"), starts_with("Kon")),
                  .funs = sum) %>% distinct()
}

# For TAGs
TAG_fun <- function(data, x) {
    data %>%
        slice(x) %>%
        mutate(
            C1 = first(C),
            C2 = nth(C,2),
            C3 = last(C),
            D1 = first(D),
            D2 = nth(D,2),
            D3 = last(D)
        ) %>%
        select(-C, -D,-index) %>%
        mutate_at(.vars = vars(starts_with("Sample"), starts_with("Kon")),
                  .funs = sum) %>% distinct()
}


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
                        fileInput("meta", label = "Metadata with sample information", multiple = F, accept = c(".xlsx")),
                        actionButton("go","Load Data!"),
                        radioButtons("radio", "Calculate DAGs and TAGs?",
                                     choices = list("Only simple sums" = 1, "All allowed permutations" = 2),
                                     selected = 2),
                        checkboxInput("drop0","Drop lipids with 0 intensity for all samples in wide form?",
                                      value = T),
                        textInput("dataset_name", "Name your dataset for export")
                    ),


                    # ** Body ------------------------------------------------------------------------------------------------------------

                    dashboardBody(
                        p("TODO: "),
                        tabsetPanel(
                            tabPanel("Raw Data",
                                     # Table Output with all measurements, raw
                                     DTOutput("rawtbl"),
                                     textOutput("debugtxt")
                            ),
                            tabPanel("Processed Data",
                                     DTOutput("datatbl"),
                                     downloadButton("downloadexp", "Download (wide form)"),
                                     downloadButton("downloadfinal", "Download (long form)")
                            )
                        )

                    )
)

# SERVER ----------------------------------------------------------------------------------------------------------

server <- function(input, output){


    # * Data Input ---------------------------------------------------------------------------------------------------

    # ** Secondary Names ----------------------------------------------------------------------------------------------

    secondary_names <- reactive({
        meta <- read_xlsx(input$meta$datapath, skip = 2, n_max = 1) %>%
            select(-c(1:3)) %>% gather(key = "secondary_name", value = "sample") %>%
            mutate(
                sample = case_when(
                    str_detect(sample, "^Sample") ~ sample,
                    str_detect(sample, "^Kon") ~ sample,
                    TRUE ~ paste("Sample", sample)
                )
            )
        cat(file=stderr(), "Reached Secondary Names \n")
        return(meta)
    })


    # ** Raw data ------------------------------------------------------------------------------------------------------

    raw <- eventReactive(input$go,{

        files <- input$files$datapath

        # Extract the class from filenames to use as .id in the complete dataframe
        filenames <- str_split(input$files$name, "_", simplify = T)
        names(files) <- str_replace(filenames[,ncol(filenames)],".txt", "")

        data <- map_df(files, read_tsv, skip = 2, .id = "class") %>%
            gather(matches("\\d+"), key = "sample_num", value = "intensity", convert = T)

        # augment the data with sample names from the raw data
        sample_names <- map_df(files, read_tsv, skip = 1, n_max = 1, .id = "class")[-c(2:7)] %>%
            gather(-class, key = "sample", value = "sample_num") %>%
            mutate(
                sample = case_when(
                    str_detect(sample, "^Sample") ~ sample,
                    str_detect(sample, "^Kon") ~ sample,
                    TRUE ~ paste("Sample", sample)
                )
            )

        cat(file=stderr(), "Reached Sample Names \n")
        data <- left_join(data, sample_names)

        # Delete unnecessary columns
        data <- data %>% select(-Polarity, -`View Type`) %>%
            rename(ID = `Sample ID`, PIS = `PIS m/z`, scan = `(ScanName)`)

        # remove rows with non-existent samples / intensities
        data <- data %>% filter(!is.na(sample))

        # Now the template with meta data
        meta_path <- input$meta$datapath

        # The rest of the meta data, first the amount of standard used per sample and class (in pmol)
        meta <- read_xlsx(meta_path, skip = 3) %>%
            select(-3) %>% filter(!is.na(Rf)) %>%
            gather(-1, -2, key = "sample", value = "standard_input")

        cat(file=stderr(), "Reached Metadata \n")

        # Then the volumes of the samples
        sample_volumes <- read_xlsx(meta_path, skip = 3, n_max = 1) %>%
            select(-1, -2, -3) %>%
            gather(key = "sample", value = "sample_volume")

        cat(file=stderr(), "Reached Sample Volumes \n")

        data <- left_join(data, meta)

        data <- left_join(data, sample_volumes) %>%
            select(-sample_num, -PIS)

        data <- data %>% mutate(LipidName = str_replace(LipidName, "\\+NH4", "")) %>%
            rename(lipid = LipidName)

        cat(file=stderr(), "Raw data is fine. \n")
        return(data)
    })


    # Data Processing -------------------------------------------------------------------------------------------------

    final <- reactive({

        if (input$radio == 2){

            withProgress(message = "Calculating DAG and TAG combinations",
                         {

                             # ** DAGs TAGs ----------------------------------------------------------------------------------------------------

                             # Extracting information about the side chains from the scan names and filtering out internal standards as `AGs_IS`.
                             # TAGs and DAGs
                             AGs <- raw() %>% filter(class %in% c("DAG", "TAG"))

                             # Their respective internal standards
                             AGs_IS <- AGs %>% filter(str_detect(lipid, "^IS")) %>% select(-scan, -ID)

                             # Information about sidechain lenghts and desaturation is in the scan name
                             AGs <- AGs %>%
                                 filter(!str_detect(lipid, "^IS")) %>%
                                 group_by(class) %>%
                                 separate(lipid, into = c("first","rest"), sep = " ", remove = F) %>%
                                 separate(rest, into = c("C_total","D_total"), sep = ":") %>% select(-first) %>%
                                 mutate(C_total = parse_number(C_total),
                                        D_total = parse_number(D_total)) %>%
                                 separate(scan, into = c("C","D"), sep = ":", remove = F) %>%
                                 mutate(C = parse_number(str_replace(C,"-","")),
                                        D = parse_number(D)) %>% select(-scan) %>%
                                 ungroup()

                             # Data and meta information for DAGs and TAGs
                             DAGs <- filter(AGs, class == "DAG")
                             meta_DAGs <- raw() %>% filter(class == "DAG") %>% select(class, sample, Rf, standard_input, sample_volume) %>% distinct()

                             TAGs <- filter(AGs, class == "TAG")
                             meta_TAGs <- raw() %>% filter(class == "TAG") %>% select(class, sample, Rf, standard_input, sample_volume) %>% distinct()

                             cat(file=stderr(), "TAGs and DAGs separated \n")
                             incProgress(0.05)

                             # *** DAG TAG calculations ----------------------------------------------------------------------------------------

                             # DAGs
                             DAGcombos <- DAGs %>% select(-Rf, - ID, -standard_input, -sample_volume) %>%
                                 spread(key = "sample", value = "intensity") %>%
                                 split(.$lipid) %>%
                                 map(function(x){
                                     x <- x %>%
                                         mutate(index = row_number())
                                     x %>% bind_rows(x) %>%
                                         combn(x = .$index, m = 2, simplify = F, FUN = DAG_fun, data = .) %>%
                                         bind_rows() %>%
                                         distinct() %>%
                                         select(-lipid)
                                 }) %>%
                                 bind_rows(.id = "lipid") %>%
                                 filter(C_total == C1 + C2 & D_total == D1 + D2) %>%
                                 mutate(
                                     key1 = paste0(C1, ":", D1),
                                     key2 = paste0(C2,":", D2),
                                     sort1 = paste0(C1, D1),
                                     sort2 = paste0(C2, D2)
                                 )

                             incProgress(0.2)

                             DAGsums <- DAGcombos %>%
                                 filter(sort1 >= sort2) %>%  # larger value first in the lipid identifyier, rest is duplicates
                                 gather(starts_with("Sample"), starts_with("Kon"), key = "sample", value = "intensity") %>%
                                 mutate(
                                     intensity = if_else(C1 == C2 & D1 == D2, intensity/2, intensity),
                                     key = paste0(key1, "_", key2)
                                 ) %>%
                                 select(-C1, -C2, -D1, -D2) %>%
                                 distinct() %>%
                                 mutate(lipid = paste0("DAG ", key)) %>%
                                 select(-key, -C_total, -D_total, -key1, -key2, -sort1, -sort2)

                             newDAGs <- left_join(DAGsums, meta_DAGs) %>%
                                 distinct() %>%
                                 bind_rows(filter(AGs_IS, class == "DAG"))


                             incProgress(0.05)

                             # TAGs
                             TAGcombos <- TAGs %>% select(-Rf, - ID, -standard_input, -sample_volume) %>%
                                 spread(key = "sample", value = "intensity") %>%
                                 split(.$lipid) %>%
                                 map(function(x){
                                     x <- x %>% mutate(index = row_number())
                                     x %>% bind_rows(x) %>%
                                         combn(x = .$index, m = 2, simplify = F, FUN = TAG_fun, data = .) %>%
                                         bind_rows() %>%
                                         distinct() %>%
                                         select(-lipid)
                                 }) %>%
                                 bind_rows(.id = "lipid") %>%
                                 filter(C_total == C1 + C2 + C3 & D_total == D1 + D2 + D3) %>%
                                 mutate(
                                     key1 = paste0(C1, ":", D1),
                                     key2 = paste0(C2,":", D2),
                                     key3 = paste0(C3,":", D3),
                                     sort1 = paste0(C1, D1),
                                     sort2 = paste0(C2, D2),
                                     sort3 = paste0(C3, D3)
                                 )

                             incProgress(0.4)

                             TAGsums <- TAGcombos %>%
                                 filter(sort1 >= sort2 & sort2 >= sort3) %>%  # larger or same value first in the lipid identifier, rest is duplicates
                                 gather(starts_with("Sample"), starts_with("Kon"), key = "sample", value = "intensity") %>%
                                 mutate(
                                     intensity = case_when(
                                         C1 == C2 & C1 == C3 & D1 == D2 & C1 == D3 ~ intensity/3,
                                         C1 == C2 & D1 == D2 ~ intensity/2,
                                         C1 == C3 & D1 == D3 ~ intensity/2,
                                         C2 == C3 & D2 == D3 ~ intensity/2,
                                         TRUE ~ intensity
                                     ),
                                     key = paste0(key1, "_", key2, "_", key3)
                                 ) %>%
                                 select(-C1, -C2, -D1, -D2, -C3,-D3) %>% distinct() %>%
                                 mutate(lipid = paste0("TAG ", key)) %>%
                                 select(-key, -C_total, -D_total, -key1, -key2, -key3, -sort1, -sort2, -sort3)

                             incProgress(0.05)

                             newTAGs <- left_join(TAGsums, meta_TAGs) %>%
                                 distinct() %>%
                                 bind_rows(filter(AGs_IS, class == "TAG"))

                             cat(file=stderr(), "TAGs and DAGs calculated \n")

                             rest <- raw() %>%
                                 filter(!class %in% c("DAG", "TAG")) %>%
                                 select(-scan, -ID) # we won't need the scan names after augmenting the data
                             data <- bind_rows(rest, newDAGs, newTAGs)

                             incProgress(0.05)

                         })  # End of withProgress
        } # end of if


        # *** simple sums -------------------------------------------------------------------------------------------------
        else{
            data <- raw() %>% select(-scan, -ID) %>% group_by_at(.vars = vars(-intensity)) %>%
                summarise(intensity = sum(intensity)) %>% ungroup()
            cat(file=stderr(), "Data was only summed, no permutations used. \n")
        }

        # ** Numbercrunching ----------------------------------------------------------------------------------------------

        # Internal Standards exist per lipid class and start with IS.
        data <- data %>%
            group_by(class, sample) %>%
            mutate(standard_mean = mean(intensity[str_detect(lipid, "^IS ")]),
                   standard_normalised = standard_input/standard_mean) %>%
            ungroup()

        # calculate pmol based on intensity and standards
        data <- data %>%
            mutate(pmol = intensity * standard_normalised) %>% filter(!str_detect(lipid, "^IS "))

        # blanc substraction
        data <- data %>%
            mutate(pmol = pmol - mean(pmol[str_detect(sample, "^Kon")])) %>%
            mutate(pmol = if_else(pmol < 0, 0, pmol)) %>% filter(!str_detect(sample, "^Kon"))

        # Response Factor
        data <- data %>%
            mutate(pmol = pmol * Rf)

        # molar µmolar
        data <- data %>%
            mutate(molar = pmol/sample_volume)

        # final data, cleaned of unneccessary columns used in earlier calculations
        final <- data %>% ungroup() %>%
            select(class, lipid, sample, molar) %>%
            left_join(secondary_names())

        return(final)

    })

    # ** Export Ready Data --------------------------------------------------------------------------------------------

    export <- reactive({
        df <- final() %>%
            select(-class, -sample) %>%
            spread(key = "secondary_name", value = "molar") %>%
            rename(Lipid = lipid)

        if (input$drop0){
            df <- df[rowSums(df[,-1] == 0) != ncol(df[,-1]),] # Delete rows with all 0s.
        }

        return(df)
    })

    # ** Output --------------------------------------------------------------------------------------------------------

    # Raw Datatable
    output$rawtbl <- DT::renderDT({
        raw()
    }, filter="top", options = list(sDom  = '<"top">lrt<"bottom">ip'))

    # Processed Data Datatable
    output$datatbl <- DT::renderDT({
        final()
    }, filter="top", options = list(sDom  = '<"top">lrt<"bottom">ip'))


    output$downloadexp <- downloadHandler(
        filename = function() {
            paste0(input$dataset_name, "_wide" ,".csv")
        },
        content = function(con){
            write_csv(export(), con)
        }
    )

    output$downloadfinal <- downloadHandler(
        filename = function() {
            paste0(input$dataset_name, "_long" ,".csv")
        },
        content = function(con){
            write_csv(final(), con)
        }
    )


    # Debug
    # output$debugtxt <- renderPrint({
    #     raw()
    # })



}



# Run App ---------------------------------------------------------------------------------------------------------

shinyApp(ui, server)
