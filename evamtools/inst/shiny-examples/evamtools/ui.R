library(DT)
# library(shinydashboard)
library(markdown)
library(shinyjs)

# source("ui/user_input_csd.R")
# source("ui/see_results_simple.R")

cpm_info <- function(){
  tags$div(id = "background",
    tags$head(
      tags$style(HTML("
        html, body{
          margin: 0;
        }
        #background{
          background-color: #f1f5f8;
          margin: 0;
          padding: 0;
          width: 100%;
          height: 100vh;
          top: 0;
          left: 0;
          position: absolute;
        }

        #cpm_info{
          background: white;
          padding: 10px 50px;
          width: 75%;
          # height: 100vh;
          margin: auto;
          margin-top: 100px;
          border: 5px solid white;
          border-radius: 5px;
          box-shadow: 0 0 10px 5px rgba(160,160,160, 0.5);
        }

        .container-fluid{
          padding:0;
        }
        "
    ))),
    tags$div(id = "cpm_info",
      includeMarkdown("test.md")
    )
  )
}

results_simple <- function(){
    fluidPage(
        tags$head(
            tags$style(
                HTML("
                    .row{
                      width: 98%;
                      margin: auto;
                        # background-color:red;
                        # height: 45vh;
                        # overflow: auto;
                    }

                    .max_height{
                      max-height: 70vh;
                      overflow: auto;
                    }
                    #sims > .col-sm-3{
                      padding: 0;
                    }

                    #noprogress input.form-control{
                      visibility: hidden;
                    }

                    .form-control{
        height: 42px !important;
      }

      #freq2label-label{
        width: 0%;
      }

      #freq2label-wrap .form-group{
        display: block !important;
        margin-bottom: 0;
      }

       #freq2label-wrap .irs-single, .irs-min, .irs-max{
        visibility: visible !important;
      }

       #freq2label-wrap .irs-grid-text{
        visibility: hidden !important;
      }

                   
                #select_cpm div.radio{
                  background-color: rgba(200,200,200, 0.5);
                  text-align: center;
                  border: 2px solid gray;
                  border-radius: 3px;
                  margin-top: 5px;
                  # margin-left: -20px;
                  white-space: nowrap;
                  overflow: hidden;
                  text-overflow: ellipsis;
                }

                @media only screen and (min-width: 1400px) {
                  #select_cpm div.radio{
                    max-width: 150px;
                  }
                  #select_cpm div.radio{
                    max-width: 150px;
                  }
                }

                @media only screen and (min-width: 1900px) {
                  #select_cpm div.radio{
                    max-width: 200px;
                  }
                  #select_cpm div.radio{
                    max-width: 200px;
                  }

                }

                #select_cpm div.radio:hover {
                  background-color: rgba(160,160,160, 0.5);
                  box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
                  cursor: pointer;
                }

                #select_cpm div.radio{
                  background-color: rgba(200,200,200, 0.5);
                  text-align: center;
                  border: 2px solid gray;
                  border-radius: 3px;
                  margin-top: 5px;
                  margin-left: -20px;
                  white-space: nowrap;
                  overflow: hidden;
                  text-overflow: ellipsis;

                    }
                    "
                )
            )
        ),
        tags$div(class = "row",
          column(1,
            tags$h3("Outputs"),
            uiOutput("cpm_list")),
          column(11,
            column(2,
              tags$div(class = "frame",
                tags$h3("1. Load  & download"),
                tags$div(id = "noprogress",
                fileInput("output_cpms", "Load your results"
                  , multiple = FALSE,
                  accept = c(".Rdata", ".rds", ".RDS"))
                ),
                tags$div(class = "download_button",
                  downloadButton("download_cpm", "Download!")
                )
              ),
              uiOutput("customize")
              
              
            ),
            column(10,
              column(12, uiOutput("sims")),
              column(12, uiOutput("sims2"))
            )
            ,
            column(4,
              uiOutput("original_data")
              
            ),
            column(6,
              uiOutput("tabular_data")
            )

          )
        )
    )

}

user_input <- function(){
  fluidPage(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
      body{
        font-size: 15px;
      }
      
      .irs-grid-pol.small{
        height: 0px !important; 
      }

      #define_genotype>*>*>label{
        margin-right: 10px;
        text-align: right;
        margin-bottom: 0;
      }
      input[type='checkbox']{
        zoom: 1.5;
        margin-right: 4px;
      }

      .checkbox{
        vertical-align: middle;
        margin-top: 10px !important; 
      }

      #genotype  input[type=checkbox] {
        zoom: 0;
      }

      #genotype>.shiny-options-group{
        display: flex;
        margin-right: 15px;
      }

      label{
        margin-right: 10px;
      }
      .frame .form-group{
        display: flex;
      }

      .form-control{
        height: 42px !important;
      }

      .form-group{
        display: flex;
        margin-bottom: 0;
      }
      .upload_file .form-group{
        width: 500px !important; 
      }
      input.formc-control{
        height: 40px !important;
      }

      #genotype>label{
        margin-right: 20px;
      }
      #add_genotype{
        width: 35%;
        height: 40px;
      }
      #fr{
        display: flex;
        justify-content: space-between;
      }

      #genotype_freq-label{
        margin-right: 20px;
      }

      div.inlin label { 
        width: 15%; 
      }

      div.inlin2 label { 
        width: 100%; 
      }
      .inlin2>div{ 
        width: 100% !important; 
      }

      input{
        z-index: 100;
      }
      .inlin label{ 
        display: table-cell; 
        text-align: left; 
        vertical-align: middle; 
      }

      .inlin>.form-group>.irs--shiny.irs-with-grid{
        margin-left: 10%;
        width: 75%;
      }

      .btn, input.form-control{
        font-size: 15px;
      }
      .inlin .form-group { 
        display: table-row;
        font-size: 15px;
      }

      .irs > span{
        font-size: 15px;
      }

      span [class*='irs'] { 
        font-size: 15px !important; 
      }

      .irs-single, .irs-min, .irs-max{
        visibility: hidden !important;
      }

      .download_button{
        heigth: auto;
        display: block;
        display:flex;
        # height: 80px;
        margin-top: 15px;
        justify-content: center;
      }

      .submit_button>button{
        font-size: 30px !important;
        padding: 15px 30px;
        background: #FAA0A0;
      }

      .submit_button>button:hover{
        font-size: 30px !important;
        padding: 15px 30px;
        background: #F20D0D;
      }

      .shiny-download-link{
        margin: auto;
        # font-size: 18px;
      }
      .not_show{
        color: rgba(0,0,0,0);
      }

      .frame{
        margin-bottom:15px;
        padding: 15px;
        border: 3px solid rgba(100, 100, 100, 0.5);
        border-radius: 5px;
      }
      .upload_file{
        width: 100%;
        display: table-cell;
        # justify_content: center
      }
      h3{
        margin-top:0;
      }

      #csd_table{
        height: 40vh;
        overflow: auto;
      }
      .row{
        padding-top: 20px;
      }

      .flex>*{
        display: inline;
      }
      .flex>h3{
        width: 50%;
        margin-right: 0;
      }

      #move div.radio:hover {
        background-color: rgba(160,160,160, 0.5);
        box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
        cursor: pointer;
      }

      #select_csd div.radio{
        background-color: rgba(200,200,200, 0.5);
        text-align: center;
        border: 2px solid gray;
        border-radius: 3px;
        margin-left: -20px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }

      @media only screen and (min-width: 1200px) {
        #input2build div.radio{
          width: 140px;
        }
        #select_csd div.radio{
          width: 140px;
        }
      }

      @media only screen and (min-width: 1900px) {
        #input2build div.radio{
          width: 150px;
        }
        #select_csd div.radio{
          width: 150px;
        }
      }

      #input2build div.radio:hover {
        background-color: rgba(160,160,160, 0.5);
        box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
        cursor: pointer;
      }

      #input2build div.radio{
        background-color: rgba(200,200,200, 0.5);
        text-align: center;
        border: 2px solid gray;
        border-radius: 3px;
        margin-left: -20px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }

      #input2build div label{
        padding-left: 20px;
      }
        ") # end HTML
      ) # end tags$style
    ),
    column(width=12,
      
      sidebarLayout(
        column(width = 1,
          tags$div(
            tags$h3("Input to build"),
            tagList(
              radioButtons(inputId = "input2build", label = "", 
                choiceNames = c("Cross sectional data", "DAG builder", "Matrix Builder"),
                choiceValues = c("csd", "dag", "matrix"),
                selected = "csd"
              )
            ),
            tags$h3("Some examples"),
            uiOutput("csd_list")
          )
          ## do it with a render UI
        ),

        column(width = 11,
          titlePanel("Cross sectional data input"),
          column(width = 6,
              column(width = 12,
            tags$div(class = "frame",
                tags$div(class = "flex",
                    tags$h3("1. Set the number of genes"),
                    actionButton("change_gene_names", "Change gene names")
                  ),
                uiOutput("genes_number")),
            tags$div(class = "frame",
              uiOutput("define_genotype"),
            ),
            uiOutput("change_freqs"),

            tags$div(class = "frame",
              tags$h3("Save & Download data"),
              uiOutput("dataset_name"),
              tags$div(
                tags$div(class = "download_button",
                  disabled(actionButton("save_csd_data", "Save Data")),
                ),
                tags$div(class = "download_button",
                  downloadButton("download_csd", "Download your data")
                )
              )
            ),
            tags$div(class = "frame",
              tags$h3("... or load your own data"),
              tags$div(class = "upload_file",
                  fileInput("csd", "Load Data",
                    multiple = FALSE,
                    accept = c(
                      ".rds", ".RDS",
                      "text/csv",
                      ".csv")) 
            ))
            )
          ),
          column(width = 6,
          
            tags$div(class = "download_button submit_button",
              actionButton("analysis", "Run evamtools!")
            ),
            tags$div(class = "download_button",
              actionButton("advanced_options", "Advanced Options")
            ),
            plotOutput("plot")
            ,
            plotOutput("dag_plot")

          )
        )
      )
      )
  )
}
ui <- 
  navbarPage( id = "navbar",
    title = "evamtools",
    tabPanel("User input", 
      value = "csd_builder",
      user_input()
    ),
    tabPanel("Results",
      value = "result_viewer",
      results_simple()
    ),
    tabPanel("Example Library",
      "WIP"
    ),
    tabPanel("About CPMs",
      cpm_info()
    ),
    tags$head(
      tags$style(
        HTML(
          "
        h2{
          margin: 0;
          padding: 0;
        }
        .navbar{
          padding-left: 20px;
          margin-bottom: 0px;
        }
        .dropdown-menu > li > a{
          font-size: 15px;
          padding: 10px 10px;
        }
        "
      )))
)

