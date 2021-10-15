library(DT)

user_input <- function(){
  fluidPage(
    tags$head(
      tags$style(HTML("
       #define_genotype {
          border: 3px solid red;
          border-radius: 3px;
          padding: 5px 15px;
          width: auto;
      }
      #define_genotype>*{
        padding:0;
      }
      #define_genotype>*>*{
        width: auto;
        display: inline-flex;
        margin-right: 15px;
        margin-bottom: 0;
      }

      #define_genotype>*>*>label{
        margin-right: 10px;
        text-align: right;
        margin-bottom: 0;
      }

      #genotype>.shiny-options-group>*{
        display: inline-block;
        margin-right: 15px;
      }
      #genotype_freq{padding:0}

      #genes_number{
        border: 3px solid red;
        border-radius: 3px;
        padding: 5px 15px;
        width: auto;
      }
      #genes_number>div{
        display: inline-flex;
        padding:0;
      }
      #genes_number>div>*{
        margin-right: 20px;
        display: inline-flex;
      }
      #gene_number-label{
        width: 100%;
      }
        ") # end HTML
      ) # end tags$style
    ),
    column(width=12,
      column(width = 6
        ,fileInput("csd", "Load a CSV"),
        "Double click in a cell to edit",
        "Set genotype to zero to remove it",
        "shift + enter to save changes when editing column",
        column(width = 12,
          uiOutput("genes_number"),
          uiOutput("define_genotype")
        ),
        DTOutput("csd_freqs"),
        ),
      column(width = 6,
        plotOutput("plot"),
        plotOutput("plot2")
      ),
    ),
    actionButton("run_cpms", "Run analysis!")
  )

}

define_CPM_usage <- function(){
  "Define CPM usage"
}

display_output <- function(){
  fluidPage(
    textOutput("out_cpms")
  )
}

ui <- fluidPage(
  tabsetPanel(id = "inTabSet",
    tabPanel("Input", fluid = TRUE, user_input()),
    tabPanel("Run CPMs", fluid = TRUE, define_CPM_usage()),
    tabPanel("Loading", fluid = TRUE, "Loading"),
    tabPanel("Output", fluid = TRUE, display_output())
  )
)