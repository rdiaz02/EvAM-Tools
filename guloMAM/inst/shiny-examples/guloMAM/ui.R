library(DT)

user_input <- function(){
  fluidPage(
    column(width=12,
      column(width = 6
        ,fileInput("csd", "Load a CSV"),
        "Double click in a cell to edit",
        "Set genotype to zero to remove it",
        "shift + enter to save changes when editing column",
        column(width = 12,
          uiOutput("gene_number"),
          uiOutput("define_genotype"),
          numericInput("genot_freq", "Frequency", NA, min = 0, width = 100),
          actionButton("add_genotype", "Add")
        ),
        DTOutput("csd_freqs"),
       
        ),
      column(width = 6,
        plotOutput("plot")
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