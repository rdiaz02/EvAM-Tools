library(DT)
library(guloMAM)

server <- function(input, output, session) {
    dB_c1 <- matrix(
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
        colnames(dB_c1) <- LETTERS[1:5]

    
    output$csd <- DT::renderDataTable({data.frame(dB_c1)})

    output_cpms <- reactiveValues(data = NULL)

    # output_cpms2 <- readRDS("/home/pablo/CPM-SSWM-Sampling/guloMAM/inst/shiny-examples/cpm_out_with_simulations.rds")
    output$out_cpms <- renderText({paste(list("a" = "no hay nada"))})
    
    observeEvent(input$run_cpms, {
        updateTabsetPanel(session, "inTabSet",selected = "Loading")
    })

    observeEvent(input$run_cpms, {
        output_cpms$data <- all_methods_2_trans_mat(dB_c1)
        output$out_cpms <- renderText({paste(output_cpms$data)})
        updateTabsetPanel(session, "inTabSet",selected = "Output")
    })
}
