library(DT)
library(guloMAM)
library(OncoSimulR)

plot_genotypes_freqs <- function(data){
    par(las = 2)
    barplot(data[, 2]
        , names = data$Genotype
        , ylab="Counts", xlab="Genotype"
        , horiz = FALSE
        , panel.first=grid())
    grid(nx = NA, ny = NULL, col='gray', lwd = 2)
}

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

    csd_freqs <- data.frame(sampledGenotypes(dB_c1))
    rownames(csd_freqs) <- csd_freqs$Genotype
    n_genes <- ncol(dB_c1)
    gene_names <- LETTERS[1: n_genes]

    ## Define number of genes
    output$gene_number <- renderUI({
             numericInput("gene_number", "Number of genes", n_genes, min = 2, width = 100)
    })

    observeEvent(input$gene_number, {
        ## Change number of genes to show
        ## Change plot
        ## Filter csd
        ## Recalculate freqs
    })

    ## Define new genotype
    output$define_genotype <- renderUI({
         checkboxGroupInput(inputId = "genotype", label = "Mutations", choices =  lapply(1:n_genes, function(i)gene_names[i]))
    })

    observeEvent(input$genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- csd_freqs[, 2][csd_freqs[, 1] == genotype]
        updateNumericInput(session, "genot_freq", value = genot_freq)
    })

    observeEvent(input$add_genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- input$genot_freq
        if(!is.na(genot_freq)){
            csd_freqs[genotype, ] <<- c(genotype, genot_freq)
            csd_freqs[, 2] <<- as.numeric(csd_freqs[, 2])
            replaceData(proxy_csd, csd_freqs, resetPaging = FALSE, rownames = FALSE)
            output$plot <- renderPlot({
                plot_genotypes_freqs(csd_freqs)
            })
        }
    })
    ## Genotypes table
    output$csd_freqs <-  DT::renderDT(csd_freqs, selection =    'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
        , rownames = FALSE,
        options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
    )
    output_cpms <- reactiveValues(data = NULL)

    proxy_csd <- dataTableProxy("csd_freqs")

    observeEvent(input$csd_freqs_cell_edit, {
        info = input$csd_freqs_cell_edit
        # str(info)
        info[ , "col"] <- 2
        csd_freqs <<- editData(csd_freqs, info, "csd_freqs") 
        ## Filtering out non-positive counts
        csd_freqs <<- csd_freqs[csd_freqs[,2] > 0,]
        # info$value <- as.numeric(info$value)
        replaceData(proxy_csd, csd_freqs, resetPaging = FALSE, rownames = FALSE)
        output$plot <- renderPlot({
            plot_genotypes_freqs(csd_freqs)
        })
    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        # Add a little noise to the cars data
        plot_genotypes_freqs(csd_freqs)
    })

    ## Run CPMS
    observeEvent(input$run_cpms, {
        updateTabsetPanel(session, "inTabSet",selected = "Loading")
    })

    observeEvent(input$run_cpms, {
        output_cpms$data <- all_methods_2_trans_mat(dB_c1)
        output$out_cpms <- renderText({paste(output_cpms$data)})
        updateTabsetPanel(session, "inTabSet",selected = "Output")
    })

    # output_cpms2 <- readRDS("/home/pablo/CPM-SSWM-Sampling/guloMAM/inst/shiny-examples/cpm_out_with_simulations.rds")
    output$out_cpms <- renderText({paste(list("a" = "no hay nada"))})
}
