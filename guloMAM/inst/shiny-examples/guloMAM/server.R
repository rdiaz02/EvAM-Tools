library(DT)
library(guloMAM)
library(OncoSimulR)

plot_genotypes_freqs <- function(data){
    par(las = 2, cex.main=1.6, cex.lab=1.5, cex.axis=1.2)
    barplot(data[, 2]
        , names = data$Genotype
        , ylab="Counts", main="Genotype Frequencies"
        , horiz = FALSE
        , panel.first=grid())
    grid(nx = NA, ny = NULL, col='gray', lwd = 2)
}

freqs2csd <- function(freqs, gene_names){
    csd <- apply(
        freqs, 1
        # mut = freqs$Genotype, freq = freqs$Freq
        , function(x){
            mut <- x[[1]]
            freq <- as.numeric(x[[2]])
            genot <- rep(0, length(gene_names))
            if(mut != "WT"){
                mut <-  strsplit(mut, ", ")[[1]]
                genot[which( gene_names %in% mut)] <- 1
            }
            csd <- matrix(rep(genot, freq), ncol= length(gene_names), byrow = TRUE)
            return(csd)
        })
    csd <- do.call(rbind, csd)
    colnames(csd) <- gene_names
    return(csd)
}

get_display_freqs <- function(freqs, n_genes, gene_names){
    valid_gene_names <- c("WT", gene_names[1:n_genes])

    selected_rows <- sapply(freqs$Genotype, function(x){
        genes <- strsplit(x, ", ")[[1]]
        return(all(genes %in% valid_gene_names))
    })

    return(freqs[selected_rows, ])
}

get_csd <- function(complete_csd){
    csd <- data.frame(sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    return(csd)
}

server <- function(input, output, session) {
    complete_csd <- matrix(
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
        colnames(complete_csd) <- LETTERS[1:5]

    min_genes <- 2
    max_genes <- 10

    data <- reactiveValues(csd_freqs =  get_csd(complete_csd)
        , complete_csd =  complete_csd
        , n_genes = ncol(complete_csd))
    
    gene_names <- reactive({LETTERS[1: input$gene_number]})
    display_freqs <- reactive({get_display_freqs(data$csd_freqs, input$gene_number, gene_names())})

    ## Upload csv
    observeEvent(input$csd, {
        # TODO hanlde corrupt files
        data$complete_csd <- read.csv(input$csd$datapath)
        data$csd_freqs <- sampledGenotypes(data$complete_csd)
        data$n_genes <- ncol(data$complete_csd)
    })

    ## Define number of genes
    output$genes_number <- renderUI({
        tags$div(id="inlin",
            sliderInput("gene_number", "Number of genes",
                value = data$n_genes, max = max_genes, min = min_genes, step = 1)
        )
    })

    observeEvent(input$display_help, {
      showModal(modalDialog(
        title = tags$h3("How does it work?"),
         tags$div(
                tags$p("1. Double click in a Frequency cell to edit it"),
                tags$p("2. Press Tab to move to the next row"),
                tags$p("3. Use Shift + Enter to save changes"),
                tags$p("4. Set a frequency to 0 to remove a genotype"),
                tags$p("5. Type in the Search bar to filter genotypes")
                )
              )
      )
    })
    ## Define new genotype
    output$define_genotype <- renderUI({
        tags$div(
            tags$div(class = "inline",
                checkboxGroupInput(inputId = "genotype", 
                    label = "Mutations", 
                    choices =  lapply(1:input$gene_number, function(i) gene_names()[i]))
            ),
            tags$div(id="fr",
                numericInput(label="Frequency", value = NA, min = 0, inputId = "genotype_freq",width = NA),
                actionButton("add_genotype", "Add Genotype")
                )
        )
    })

    observeEvent(input$genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- data$csd_freqs[, 2][data$csd_freqs[, 1] == genotype]
        updateNumericInput(session, "genotype_freq", value = genot_freq)
    })
    
    observeEvent(input$add_genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- input$genotype_freq
        if(!is.na(genot_freq)){
            data$csd_freqs[genotype, ] <- c(genotype, genot_freq)
            data$csd_freqs[, 2] <- as.numeric(data$csd_freqs[, 2])
            rownames(data$csd_freqs) <- data$csd_freqs$Genotype
        }
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices =  lapply(1:input$gene_number, function(i)gene_names()[i]), selected = NULL)
    })
    
    ## Genotypes table
    output$csd_freqs <-  DT::renderDT(display_freqs(), selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
        , rownames = FALSE,
        options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
    )
    # output_cpms <- reactiveValues(data = NULL)

    # proxy_csd <- dataTableProxy("csd_freqs")

    observeEvent(input$csd_freqs_cell_edit, {
        info = input$csd_freqs_cell_edit
        info[ , "col"] <- 2
        data$csd_freqs <- editData(data$csd_freqs, info, "csd_freqs") 
        ## Filtering out non-positive counts
        data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]

    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        plot_genotypes_freqs(display_freqs())
    })

    ## Download button
    observeEvent(input$download, {
    })
    output$download <- downloadHandler(
        filename = "cross_section_data.csv",
        content = function(file) {
            print(file)
            data2download <- freqs2csd(display_freqs(), gene_names())
            write.csv(data2download, file, row.names = FALSE)
        }
    )

    ## Run CPMS
    # observeEvent(input$run_cpms, {
    #     updateTabsetPanel(session, "inTabSet",selected = "Loading")
    # })

    # observeEvent(input$run_cpms, {
    #     output_cpms$data <- all_methods_2_trans_mat(dB_c1)
    #     output$out_cpms <- renderText({paste(output_cpms$data)})
    #     updateTabsetPanel(session, "inTabSet",selected = "Output")
    # })

    # # output_cpms2 <- readRDS("/home/pablo/CPM-SSWM-Sampling/guloMAM/inst/shiny-examples/cpm_out_with_simulations.rds")
    # output$out_cpms <- renderText({paste(list("a" = "no hay nada"))})
}
