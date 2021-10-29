library(DT)
library(shinydashboard)
library(markdown)

source("ui/user_input_csd.R")
source("ui/see_results_simple.R")

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
                    "
                )
            )
        ),
        tags$div(class = "row",
          column(1,
            "list here all CPM outputs"),
          column(11,
            column(3,
              # fileInput("output_cpms", "Load your results"
              #   , multiple = FALSE,
              #   accept = c(".Rdata", ".rds", ".RDS")),
              uiOutput("cpm_names"),
              numericInput("top_paths", "Paths to show",
                min = 0, step = 1, value = NULL)
            ),
            column(5,
                plotOutput("sims")
            ),
            column(4,
                plotOutput("csd"),
                tags$div(class = "download_button",
                    actionButton("modify_data", "Modify data")
                )
            )
          )
        )
    )

}

ui <- 
  navbarPage( id = "navbar",
    title = "guloMAM",
    tabPanel("Results",
      results_simple()
    ),
    navbarMenu("User Input",
      tabPanel("Cross sectional Data Builder", value = "csd_builder",
        user_input()
      ),
      tabPanel("DAG Builder (CBN)",
      ),
      tabPanel("Matrix Builder (MHN)",
      )
    ),
    tabPanel("About CPMS",
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
          font-size: 20px;
          padding: 10px 10px;
        }
        "
      )))
)

