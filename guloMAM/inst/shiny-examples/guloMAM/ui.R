library(DT)

user_input <- function(){
  fluidPage(

  fileInput("csd", "Load a CSV"),
  DT::dataTableOutput("csd"),
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