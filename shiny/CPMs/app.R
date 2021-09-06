library(shiny)
library(DT)



  csd <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0), 300) #A
      , rep(c(0, 0, 1, 0, 0), 300) #C
      , rep(c(1, 1, 0, 0, 0), 200) #AB
      , rep(c(0, 0, 1, 1, 0), 200) #CD
      , rep(c(1, 1, 1, 0, 0), 100) #ABC
      , rep(c(1, 0, 1, 1, 0), 100) #ACD
      , rep(c(1, 1, 0, 0, 1), 100) #ABE
      , rep(c(0, 0, 1, 1, 1), 100) #CDE
      , rep(c(1, 1, 1, 0, 1), 100) #ABCE
      , rep(c(1, 0, 1, 1, 1), 100) #ACDE
      , rep(c(1, 1, 1, 1, 0), 50) # ABCD
      , rep(c(0, 0, 0, 0, 0), 10) # WT
    ), ncol = 5, byrow = TRUE
  )
  colnames(csd) <- LETTERS[1:5]
ui <- basicPage(
  h2("Cross sectional data"),
  DT::dataTableOutput("mytable")
)

server <- function(input, output) {
  output$mytable = DT::renderDataTable({
    csd
  })
}

shinyApp(ui, server)
