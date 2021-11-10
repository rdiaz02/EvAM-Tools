library(DT)
library(guloMAM)
library(OncoSimulR)
library(shinyjs)

source("../../../data/toy_datasets.R")

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
    # browser()
    csd <- data.frame(sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    return(csd)
}


available_cpms <- function(data){
    data$csd_data <- NULL
    cpm_names <- unique(sapply(names(data), function(x) str_split(x, "_")[[1]][[1]]))
    return(cpm_names)
}


plot_model <- function(cpm_output, mod){
    ## DAG relationships colors 
    standard_relationship <- "gray73"
    colors_relationships <- c(standard_relationship, "coral2", "cornflowerblue", "darkolivegreen3")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    
    
    model_data2plot <- process_data(cpm_output, mod)
    ## Plotting data
    if(!is.null(model_data2plot$dag_tree)) {
        if(!is.null(model_data2plot$parent_set)){
            for(i in names(model_data2plot$parent_set)){
                E(model_data2plot$dag_tree)[.to(i)]$color <- colors_relationships[model_data2plot$parent_set[[i]]]
            }
        } else E(model_data2plot$dag_tree)$color <- standard_relationship
        plot(model_data2plot$dag_tree
            , layout = layout.reingold.tilford
            , vertex.size = 55 
            , vertex.label.color = "black"
            , vertex.label.family = "Helvetica"
            , vertex.label.cex = 1.5
            , font.best = 2
            , vertex.frame.width = 0.5
            , vertex.color = "white"
            , vertex.frame.color = "black" 
            , vertex.label.cex = 1
            , edge.arrow.size = 0
            , edge.width = 5
            )
        # par(mar=c(0, 0, 0, 0))
        if(!is.null(model_data2plot$parent_set)){
            legend("topleft", legend = names(colors_relationships),
                col = colors_relationships, lty = 1, lwd = 5, bty = "n")
        }
    }else if(!is.null(model_data2plot$theta)) {
        op <- par(mar=c(3, 3, 5, 3), las = 1)
        plot(model_data2plot$theta, cex = 1.4, digits = 2, key = NULL
            , axis.col = list(side = 3)
            , xlab = "Effect of this (effector)"
            , ylab = " on this (affected)"
            , main = ""
            , mgp = c(2, 1, 0))
        par(op)
    }
    title(mod, cex.main = 1.8)
}

compare_cpm_freqs <- function(data){
    all_counts <- data.frame(Genotype = data[["MHN_genotype_freqs"]]$Genotype)
    for(name in names(data)){
        if(grepl("_genotype_freqs", name)){
            method_name <- strsplit(name, "_")[[1]][[1]]
            all_counts[[method_name]] <- data[[name]]$Counts
        }
    }

    order_by_counts <- sort(rowSums(all_counts[-1]), 
        decreasing = TRUE, index.return = TRUE)$ix
    return(all_counts[order_by_counts, ])
}

server <- function(input, output, session) {
    ## CSD input
    ## Load many examples

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
            , rep(c(0, 0, 0, 0, 0), 100) # WT
        ), ncol = 5, byrow = TRUE
        )
        colnames(complete_csd) <- LETTERS[1:5]
    
    all_csd_data <- c(
        list(user = list(data = complete_csd, name = "User Data"))
        , all_examples_csd_2)

    

    min_genes <- 2
    max_genes <- 10

    data <- reactiveValues(
        csd_freqs =  get_csd(complete_csd)
        , all_csd = all_csd_data
        , complete_csd = complete_csd
        , n_genes = complete_csd
        )
    
    gene_names <- reactive({LETTERS[1: input$gene_number]})
    display_freqs <- reactive({get_display_freqs(data$csd_freqs, input$gene_number, gene_names())})

    ## Display List of availabe CSD 

    output$csd_list <- renderUI({
        all_names <- unname(sapply(data$all_csd, function(dataset) dataset$name))
        
        tagList(
            radioButtons(
                inputId = "select_csd",
                label = "",
                selected = "se",
                choiceNames = all_names,
                choiceValues = names(data$all_csd)
            )
      )
    })

    observeEvent(input$select_csd, {
        data$complete_csd <- data$all_csd[[input$select_csd]]$data
        data$csd_freqs <- sampledGenotypes(data$complete_csd)
        data$n_genes <- ncol(data$complete_csd)
    })


    ## Upload csv
    observeEvent(input$csd, {
        # TODO hanlde corrupt files
        data$complete_csd <- read.csv(input$csd$datapath)
        data$csd_freqs <- sampledGenotypes(data$complete_csd)
        data$n_genes <- ncol(data$complete_csd)
        updateRadioButtons(session, "select_csd", selected = "user")
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
        if(input$input2build == "CSD"){
            tags$div(
                tags$h3("2. Add new genotypes"),
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
        } else if (input$input2build == "DAG"){
            tags$div(
                tags$h3("2. Define edges"),
                tags$h3("WIP")
            )
        }
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
            rownames(data$csd_freqs) <- data$csd_freqs$Genotype
            data$csd_freqs[genotype, ] <- c(genotype, genot_freq)
            data$csd_freqs[, 2] <- as.numeric(data$csd_freqs[, 2])
            ## Filtering out non-positive counts
            data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
            data$all_csd[[input$select_csd]]$data <- freqs2csd(data$csd_freqs, gene_names())
        }
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices =  lapply(1:input$gene_number, function(i)gene_names()[i]), selected = NULL)
    })
    
    ## Genotypes table
    output$csd_freqs <- DT::renderDT(display_freqs(), selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
        , rownames = FALSE,
        options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
    )

    observeEvent(input$csd_freqs_cell_edit, {
        info = input$csd_freqs_cell_edit
        info[ , "col"] <- 2
        data$csd_freqs <- editData(data$csd_freqs, info, "csd_freqs") 
        ## Filtering out non-positive counts
        data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
        # browser()
        data$all_csd[[input$select_csd]]$data <- freqs2csd(data$csd_freqs, gene_names())
    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        plot_genotypes_freqs(display_freqs())
    })

    ## Download csd button
    output$download_csd <- downloadHandler(
        filename = "cross_section_data.csv",
        content = function(file) {
            data2download <- freqs2csd(display_freqs(), gene_names())
            write.csv(data2download, file, row.names = FALSE)
        }
    )

    ## Run CPMS
    # observeEvent(input$run_cpms, {
    #     updateTabsetPanel(session, "inTabSet",selected = "Loading")
    # })

    cpm_out <- readRDS("/home/pablo/CPM-SSWM-Sampling/guloMAM/inst/shiny-examples/sims.RDS")
    cpm_out$MHN_f_graph <- cpm_out$MHN_transitionRateMatrix
    
    all_cpm_out <- reactiveValues(output = list(default = cpm_out))

    observeEvent(input$analysis, {
        shinyjs::disable("analysis")
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Running GuloMAM", value = 0)

        data2run <- freqs2csd(display_freqs(), gene_names())
        progress$inc(1/4, detail = "Setting up data")
        Sys.sleep(0.5)
        progress$inc(2/4, detail = "Running CPMs")
        cpm_output <- all_methods_2_trans_mat(data2run)
        n_samples <- 10000 
        progress$inc(3/4, detail = paste("Running ", n_samples, " samples"))
        new_data <- sample_all_CPMs(cpm_output, n_samples, input$gene_number)
        progress$inc(4/4, detail = "Post processing data")
        Sys.sleep(0.5)
        new_data$MHN_f_graph <- new_data$MHN_transitionRateMatrix
        all_cpm_out$output[["user_input"]] <- new_data
        updateTabsetPanel(session, "navbar", selected = "result_viewer")
        shinyjs::enable("analysis")
    })
    
    genotype_freq_df <- reactive(compare_cpm_freqs(all_cpm_out$output[[
        length(all_cpm_out$output)
        ]]))

    output$cpm_freqs <- DT::renderDT(genotype_freq_df(),  
        selection = 'none', server = TRUE
        , rownames = FALSE
        , options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE)
    )

    output$sims <- renderUI({
        column_models2show <- floor(12 / length(input$cpm2show)) 

        lapply(input$cpm2show, function(mod){
            output[[sprintf("plot_sims_%s", mod)]] <- renderPlot({
                pl <- plot_model(all_cpm_out$output[["user_input"]], mod)
            })
            return(
                column(3,
                    plotOutput(sprintf("plot_sims_%s", mod)))
            )
        }
        )
    })

    output$sims2 <- renderUI({
        column_models2show <- floor(12 / length(input$cpm2show)) 

        attribute_name <- c("f_graph", "genotype_transitions", "trans_mat")

        names(attribute_name) <- c("Probabilities", "Transitions", "Transition Matrix")

        selected_plot_type <- attribute_name[input$data2plot]

        top_paths <- input$top_paths
        if(top_paths == 0) top_paths <- NULL
        
        lapply(input$cpm2show, function(mod){
            data2plot <- all_cpm_out$output[["user_input"]][[sprintf("%s_%s", mod, 
                # "trans_mat"
                selected_plot_type
                )]]
            output[[sprintf("plot_sims2_%s", mod)]] <- renderPlot({
                pl <- plot_genot_fg(data2plot, 
                    all_cpm_out$output[["user_input"]]$csd_data, 
                    all_cpm_out$output[["user_input"]][[sprintf("%s_genotype_%s", mod, "freqs")]], 
                    top_paths = top_paths)
            })
            return(
                column(3,
                    plotOutput(sprintf("plot_sims2_%s", mod)))
            )
        })
    })

    output$csd <- renderPlot({
        plot_genotypes_freqs(get_csd(all_cpm_out$output[["user_input"]]$csd_data))
    })

    observeEvent(input$modify_data, {
        data$csd_freqs <- sampledGenotypes(all_cpm_out$output[["user_input"]]$csd_data)
        rownames(data$csd_freqs) <- data$csd_freqs$Genotype
        data$n_genes <- ncol(all_cpm_out$output[["user_input"]]$csd_data)
        updateTabsetPanel(session, "navbar",
            selected = "csd_builder"
        )
        updateRadioButtons(session, "select_csd", selected = "user")
    })

    ## Download button
    output$download_cpm <- downloadHandler(
        filename = "cpm.RDS",
        content = function(file) {
            saveRDS(all_cpm_out$output[["user_input"]], file)
        }
    )
}
