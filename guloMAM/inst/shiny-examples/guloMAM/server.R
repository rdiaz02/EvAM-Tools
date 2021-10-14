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
    display_csd <- complete_csd
    min_genes <- 2
    max_genes <- 10

    csd_freqs <- data.frame(sampledGenotypes(display_csd))
    rownames(csd_freqs) <- csd_freqs$Genotype
    gene_names <- LETTERS[1: n_genes]
    display_freqs <- get_display_freqs(csd_freqs, n_genes, gene_names)
    # complete_csd <- freqs2csd(csd_freqs, gene_names)

    ## Define number of genes
    output$gene_number <- renderUI({
             numericInput("gene_number", "Number of genes",
                n_genes, max = max_genes, min = min_genes, 
                width = 100)
    })

    observeEvent(input$gene_number, {
        ## Change number of genes to show
        old_gene_number <- n_genes
        if (input$gene_number >= 2) n_genes <<- input$gene_number

        ## Update Labels 
        gene_names <<- LETTERS[1:n_genes]
        output$define_genotype <- renderUI({
            checkboxGroupInput(inputId = "genotype", label = "Mutations", choices =  lapply(1:n_genes, function(i)gene_names[i]))
        })

        ## Recalculate freqs
        display_freqs <<- get_display_freqs(csd_freqs, n_genes, gene_names)
        # rownames(csd_freqs) <<- csd_freqs$Genotype
        ##Update table 
        replaceData(proxy_csd, display_freqs, resetPaging = FALSE, rownames = FALSE)
        
        ## Change plot 
        output$plot <- renderPlot({
            plot_genotypes_freqs(display_freqs)
        })
    })

    ## Define new genotype
    output$define_genotype <- renderUI({
         checkboxGroupInput(inputId = "genotype", inline = TRUE, label = "Mutations", choices =  lapply(1:n_genes, function(i)gene_names[i]))
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
            rownames(csd_freqs) <- csd_freqs$Genotype
            replaceData(proxy_csd, csd_freqs, resetPaging = FALSE, rownames = FALSE)
            display_freqs <- get_display_freqs(csd_freqs, n_genes, gene_names)
            output$plot <- renderPlot({
                plot_genotypes_freqs(display_freqs)
            })
        }
         updateNumericInput(session, "genot_freq", value = NA)
         updateCheckboxGroupInput(session, "genotype",label = "Mutations", choices =  lapply(1:n_genes, function(i)gene_names[i]), selected = NULL)
        #  complete_csd <<- freqs2csd(csd_freqs, gene_names)
    })
    ## Genotypes table
    output$csd_freqs <-  DT::renderDT(display_freqs, selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
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
        # complete_csd <<- freqs2csd(csd_freqs, gene_names)
        # display_csd <<- get_display_csd(complete_csd, n_genes)
        # info$value <- as.numeric(info$value)
        replaceData(proxy_csd, csd_freqs, resetPaging = FALSE, rownames = FALSE)
        display_freqs <- get_display_freqs(csd_freqs, n_genes, gene_names)
        output$plot <- renderPlot({
            plot_genotypes_freqs(display_freqs)
        })
    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        # Add a little noise to the cars data
        plot_genotypes_freqs(display_freqs)
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
