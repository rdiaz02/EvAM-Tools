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

    n_genes <- ncol(complete_csd)
    min_genes <- 2
    max_genes <- 10

    csd <- data.frame(sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    data <- reactiveValues(csd_freqs = csd)
    gene_names <- reactive({LETTERS[1: input$gene_number]})
    display_freqs <- reactive({get_display_freqs(data$csd_freqs, input$gene_number, gene_names())})
    

    ## Define number of genes
    output$genes_number <- renderUI({
             numericInput("gene_number", "Number of genes",
                value = ncol(complete_csd), max = max_genes, min = min_genes)
    })

    ## Define new genotype
    output$define_genotype <- renderUI({
        fluidPage(
            checkboxGroupInput(inputId = "genotype", label = "Mutations", choices =  lapply(1:input$gene_number, function(i) gene_names()[i])),
            numericInput(label="Frequency", value = NA, min = 0, inputId = "genotype_freq",width = NA),
            actionButton("add_genotype", "Add")
        )
    })

    observeEvent(input$add_genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- input$genotype_freq
        if(length(genot_freq) > 0){
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
