library(DT)
library(guloMAM)
library(OncoSimulR)
library(shinyjs)
library(igraph)

source("../../../data/toy_datasets.R")

check_if_csd <- function(data){
    tmp_names <- c("data", "dag", "dag_parent_set", "name", "type", "trm", "thetas")
    types <- c("csd", "dag", "matrix")
    if(all(names(data) %in% tmp_names)){
        if((data$type %in% types 
        & all(unique(c(data$data)) %in% c(0, 1)))){
            return(TRUE)
        } else {return(FALSE)}
    } else { ## We assume a data.frame
        tmp_data <- unlist(data)
        names(tmp_data) <- NULL
        return(all(unique(tmp_data) %in% c(0, 1)))
    }
}

plot_genotypes_freqs <- function(data){
    if(is.null(data)) return()
    par(las = 2, cex.main=1.6, cex.lab=1.5, cex.axis=1.2)
    barplot(data[, 2]
        , names = data$Genotype
        , ylab="Counts", main="Genotype Frequencies"
        , horiz = FALSE
        , panel.first=grid())
    grid(nx = NA, ny = NULL, col='gray', lwd = 2)
    ## TODO short genotypes
}

plot_dag <- function(dag, parent_set){
    if (is.null(dag)) return()
    standard_relationship <- "gray73"
    colors_relationships <- c(standard_relationship, "coral2", "cornflowerblue", "darkolivegreen3")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")

    dag <- graph_from_adjacency_matrix(dag, mode = "directed")
    dag <- igraph::decompose(dag)[[1]]
    ## Plotting data
    if(!is.null(parent_set)){
        for(i in setdiff(names(V(dag)), "WT")){
            E(dag)[.to(i)]$color <- colors_relationships[parent_set[[i]]]
        }
    } else E(dag)$color <- standard_relationship
        
    plot(dag
        , layout = layout.reingold.tilford
        , vertex.size = 30 
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
    if(!is.null(parent_set)){
        legend("topleft", legend = names(colors_relationships),
            col = colors_relationships, lty = 1, lwd = 5, bty = "n")
    }
    title("Direct acyclic graph", cex.main = 1.8)
}

freqs2csd <- function(freqs, gene_names){
    csd <- apply(
        freqs, 1
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
    if(is.null(freqs)) return(NULL)
    valid_gene_names <- c("WT", gene_names[1:n_genes])

    selected_rows <- sapply(freqs$Genotype, function(x){
        genes <- strsplit(x, ", ")[[1]]
        return(all(genes %in% valid_gene_names))
    })

    return(freqs[selected_rows, ])
}

get_csd <- function(complete_csd){
    if(is.null(complete_csd)) return(NULL)
    csd <- data.frame(sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    return(csd)
}

available_cpms <- function(data){
    data$csd_data <- NULL
    cpm_names <- unique(sapply(names(data), function(x) str_split(x, "_")[[1]][[1]]))
    return(cpm_names)
}

get_mhn_data <- function(n_genes, n_samples, gene_names, thetas = NULL){
    if(is.null(thetas)) thetas <- Random.Theta(n=n_genes)
    rownames(thetas) <- colnames(thetas) <- gene_names
    samples <- floor(Finite.Sample(Generate.pTh(thetas), n_samples)*n_samples)
    trm <- theta_to_trans_rate_3_SM(thetas,
                                    inner_transition = inner_transitionRate_3_1)
    state_names <- vapply(1:(ncol(trm)), function(x){
        x <- x - 1
        if(x == 0) state_name <- "WT"
        else state_name <- paste(gene_names[which(int2binary(x, n_genes) == 1)], collapse = ", ")
        return(state_name)
    }, character(1))
    rownames(trm) <- colnames(trm) <- state_names
    samples <- data.frame("Genotype" = state_names, "Freq" = samples)
    rownames(samples) <- samples$Genotype
    return(list(thetas = thetas, trm = trm, samples = samples))
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

dataModal <- function(error_message) {
      modalDialog(
        easyClose = TRUE,
        title = tags$h3("There was an error"),
        tags$div(
            error_message
        )
      )
    }

server <- function(input, output, session) {
    all_csd_data <- all_examples_csd_2
    min_genes <- 2
    max_genes <- 10
    default_mhn_samples <- 5000
    keep_dataset_name <- FALSE
    all_gene_names <- LETTERS[1: max_genes]
    error_message <- NULL
    template_dag <- matrix(0, ncol= max_genes + 1, nrow = max_genes + 1)
    rownames(template_dag) <- colnames(template_dag) <- c("WT", all_gene_names)
    template_parent_set <- rep("Single", max_genes)
    names(template_parent_set) <- all_gene_names
    template_lambdas <- rep(1, max_genes)
    names(template_lambdas) <- all_gene_names

    data <- reactiveValues(
        csd_freqs =  NULL
        , all_csd = all_csd_data
        , complete_csd = NULL
        , dag = template_dag
        , dag_parent_set = template_parent_set
        , lambdas = template_lambdas
        , thetas = NULL
        , trm = NULL
        , n_genes = 3
        , gene_names = LETTERS[1: max_genes]
    )

    last_visited_pages <- list(csd = "user", dag = "user", matrix = "user")
    
    display_freqs <- reactive({
        get_display_freqs(data$csd_freqs, input$gene_number, data$gene_names
    )})

    ## Upload data
    observeEvent(input$csd, {
        # TODO hanlde corrupt files
        if(grepl(".csv", input$csd$datapath)){
            dataset_name <- strsplit(strsplit(input$csd$name, ".csv")[[1]], "_")[[1]][[1]]
            tmp_data <- read.csv(input$csd$datapath)
            if(check_if_csd(tmp_data)){
                data$all_csd[[dataset_name]]$data <- tmp_data
                data$all_csd[[dataset_name]]$name <- dataset_name
                data$all_csd[[dataset_name]]$type <- "csd"
                # data$complete_csd <- tmp_data
                # data$csd_freqs <- sampledGenotypes(data$complete_csd)
                # data$gene_names <- colnames(data$complete_csd)
                keep_dataset_name <<- dataset_name 
                updateRadioButtons(session, "input2build", selected = "csd")
                updateRadioButtons(session, "select_csd", selected = dataset_name)
                error_message <<- "Your csv data can not be loaded. Make sure it only contains 0 and 1."
                
            }else {
                showModal(dataModal(error_message))
            }
        }else if(grepl(".rds", input$csd$datapath, ignore.case = TRUE)){
            tmp_data <- readRDS(input$csd$datapath)
            if(check_if_csd(tmp_data)){
                data$all_csd[[tmp_data$name]] <- tmp_data
                updateRadioButtons(session, "input2build", selected = tmp_data$type)
                updateRadioButtons(session, "select_csd", selected = tmp_data$name)
            } else {
                error_message <<- "There was a problem when checking your .rds file. Make sure it containis $type (either 'csd', 'dag', or 'matrix'), $data only with 0 and 1"
                showModal(dataModal(error_message))
            }
        }
    })

    observeEvent(input$input2build, {
        if(keep_dataset_name %in% names(data$all_csd)){
            updateRadioButtons(session, "select_csd", selected = keep_dataset_name)
            keep_dataset_name <- FALSE
        }else{
            updateRadioButtons(session, "select_csd", selected = last_visited_pages[[input$input2build]])
        }
    })

    ## Download csd button
    output$download_csd <- downloadHandler(
        filename = function() sprintf("%s_csd.RDS", input$select_csd),
        content = function(file) {
            saveRDS(data$all_csd[[input$select_csd]], file)
        }
    )

    ## Display List of availabe CSD 
    output$csd_list <- renderUI({
        all_names <- c()
        all_choice_names <- c()
        for (i in names(data$all_csd[[input$input2build]])){
            tmp_data <- data$all_csd[[input$input2build]][[i]]
            # if(tmp_data$name == "User Data" || tmp_data$type == input$input2build){
            all_names <- c(all_names, i)
            all_choice_names <- c(all_choice_names, tmp_data$name)
            # }
        }

        tagList(
            radioButtons(
                inputId = "select_csd",
                label = "",
                selected = last_visited_pages[[input$input2build]],
                choiceNames = all_choice_names,
                choiceValues = all_names
            )
        )
    })

    observeEvent(input$select_csd, {
        ## Cleaning stuf
        # browser()
        data$csd_freqs <-  NULL
        data$complete_csd <- NULL
        data$dag <- template_dag
        data$dag_parent_set <- template_parent_set
        data$lambdas <- template_lambdas
        data$thetas <- NULL
        data$trm <- NULL
        data$gene_names <- LETTERS[1: max_genes]
        # browser()

        last_visited_pages[[input$input2build]] <<- input$select_csd

        tmp_complete_csd <- data$all_csd[[input$input2build]][[input$select_csd]]$data
        tmp_thetas <- data$all_csd[[input$input2build]][[input$select_csd]]$thetas
        tmp_dag <- data$all_csd[[input$input2build]][[input$select_csd]]$dag
        tmp_parent_set <- data$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set
        
        if(input$input2build == "dag" & !is.null(tmp_dag)){
            ## Filtering
            to_keep <- which(colSums(tmp_dag)>0 | rowSums(tmp_dag)>0)
            tmp_dag <- tmp_dag[to_keep, to_keep]
            n_genes <- ncol(tmp_dag)
            new_dag <- template_dag 
            new_dag[1:n_genes, 1:n_genes] <- tmp_dag
            rownames(new_dag) <- colnames(new_dag) <- c("WT", data$gene_names)
            data$dag <- new_dag
            data$dag_parent_set[names(tmp_parent_set)] <- tmp_parent_set
            data$lambdas[names(tmp_parent_set)] <- data$all_csd[[input$input2build]][[input$select_csd]]$lambdas
            n_genes <- n_genes - 1
            tmp_gene_names <- colnames(new_dag)[-1]
        }else if(input$input2build == "thetas" & !is.null(tmp_thetas)){
            data$thetas <- tmp_thetas 
            data$trm <- data$all_csd[[input$input2build]][[input$select_csd]]$trm
            n_genes <- ncol(tmp_thetas)
            tmp_gene_names <- colnames(tmp_thetas)
        }else if(input$input2build == "csd" & !is.null(tmp_complete_csd)){
            data$complete_csd <- tmp_complete_csd
            data$csd_freqs <- sampledGenotypes(tmp_complete_csd) 
            n_genes <- ncol(data$complete_csd)
            tmp_gene_names <- colnames(data$complete_csd)
        }else { 
            n_genes <- 3 
            tmp_gene_names <- data$gene_names    
        }

        data$gene_names <- c(tmp_gene_names, 
                                LETTERS[(length(tmp_gene_names) + 1): max_genes])
        
        updateNumericInput(session, "gene_number", value = n_genes)
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices =  lapply(1:input$gene_number, function(i)data$gene_names[i]), selected = NULL)
    })

    ## Define number of genes
    output$genes_number <- renderUI({
        val <- ifelse(is.null(data$n_genes), 3, data$n_genes)
        tags$div(class="inlin",
            sliderInput("gene_number", "Number of genes",
                value = val, max = max_genes, min = min_genes, step = 1)
        )
    })

    ## Define dataset name
    output$dataset_name <- renderUI({
        tags$div(class="inlin2",
            textInput(inputId="dataset_name", "Give your dataset a name", value = input$select_csd)
        )
    })

    observeEvent(input$change_gene_names, {
        showModal(modalDialog(
            easyClose = TRUE,

            title = tags$h3("Modify gene names"),
            tags$div(class = "inlin2",
                textInput(inputId = "new_gene_names", "Genes names", 
                    value = paste(data$gene_names[1:input$gene_number], 
                    collapse = ", ")
                    # value = "goku, vegeta, krilin", 
                ),
            tags$h3("Please, separate you gene names with a ','"),    
            tags$div(class = "download_button",
                actionButton("action_gene_names", "Change genes names"),
                )
            )
        ))
    })

    ## Updating gene names
    observeEvent(input$action_gene_names,{
        new_gene_names <- unique(strsplit(gsub(" ", "", input$new_gene_names), ",")[[1]])
        all_gene_names <- data$gene_names
        all_gene_names[1:length(new_gene_names)] <- new_gene_names
        
        data$gene_names <- all_gene_names
        colnames(data$complete_csd) <- data$gene_names[1:ncol(data$complete_csd)]
        data$csd_freqs <- sampledGenotypes(data$complete_csd)
        if(!is.null(data$dag)){
            names(data$dag_parent_set) <- all_gene_names
            rownames(data$dag) <- colnames(data$dag) <- c("WT", all_gene_names)
        } else if (!is.null(data$thetas)){
            rownames(data$thetas) <- colnames(data$thetas) <- data$gene_names[1:input$gene_number]
            state_names <- vapply(1:(ncol(data$trm)), function(x){
                x <- x - 1
                if(x == 0) state_name <- "WT"
                else state_name <- paste(data$gene_names[which(int2binary(x, input$gene_number) == 1)], collapse = ", ")
                return(state_name)
            }, character(1))
            rownames(data$trm) <- colnames(data$trm) <- state_names
        }
        if(!is.null(input$select_csd) & input$select_csd != "user"){
            data$all_csd[[input$input2build]][[input$select_csd]]$data <- data$complete_csd
            data$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
            data$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set <- data$dag_parent_set
            data$all_csd[[input$input2build]][[input$select_csd]]$lambdas <- data$lambdas
            data$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
            data$all_csd[[input$input2build]][[input$select_csd]]$trm <- data$trm
        }
    })

    ## Saving dataset
    observeEvent(input$save_csd_data,{
        ## 1 save dataset to list after user data

        if(!(input$dataset_name %in% names(data$all_csd))){

            data$all_csd[[input$input2build]][[input$dataset_name]]$data <- freqs2csd(data$csd_freqs, data$gene_names)
            data$all_csd[[input$input2build]][[input$dataset_name]]$name <- input$dataset_name
            data$all_csd[[input$input2build]][[input$dataset_name]]$type <- input$input2build
            data$all_csd[[input$input2build]][[input$dataset_name]]$dag <- data$dag
            data$all_csd[[input$input2build]][[input$dataset_name]]$dag_parent_set <- data$dag_parent_set
            data$all_csd[[input$input2build]][[input$dataset_name]]$thetas <- data$thetas
            data$all_csd[[input$input2build]][[input$dataset_name]]$trm <- data$trm

            data$all_csd[[input$input2build]] <- 
                c(data$all_csd[[input$input2build]]["user"],
                data$all_csd[[input$input2build]][input$dataset_name],
                data$all_csd[[input$input2build]][which(!(names(data$all_csd) %in% c("user", input$dataset_name, names(all_examples_csd_2))))],
                data$all_csd[[input$input2build]][which(names(data$all_csd) %in% names(all_examples_csd_2))]
            )
            ## 2 restore default values
            data$all_csd[[input$input2build]][[input$select_csd]]$data <- all_examples_csd_2[[input$select_csd]]$data

            ## 3 update selected entry
            updateRadioButtons(session, "select_csd", selected = input$dataset_name)
        }
    })

    observeEvent(input$display_help, {
      showModal(modalDialog(
        easyClose = TRUE,
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

    observeEvent(input$advanced_options, {
        showModal(modalDialog(
            size = "l",
            easyClose = TRUE,
            title = tags$h3("Advanced options?"),
            tags$div(
                numericInput("num_steps", "Sampling steps", 10000, min=0, max=1000000,step=1000, width="100%"),
                checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("HyperTRAPS", "MCCBN"), choiceValues = c("hypertraps", "mccbn")),
                tags$h4("DISCLAIMER: Both HyperTraps and MCCBN may take hours to run")
            )
            )
        )
    })
    ## Define new genotype
    output$define_genotype <- renderUI({
        options <- data$gene_names[1:input$gene_number]
        if(input$input2build == "csd"){
            tags$div(
                tags$h3("2. Add new genotypes"),
                tags$div(class = "inline",
                    checkboxGroupInput(inputId = "genotype", 
                        label = "Mutations", 
                        choices =  options)
                ),
                tags$div(id="fr",
                    numericInput(label="Frequency", value = NA, min = 0, inputId = "genotype_freq",width = NA),
                    actionButton("add_genotype", "Add Genotype")
                )
            )
        } else if (input$input2build == "dag"){
            tags$div(
                tags$div(class = "flex",
                    tags$h3("2. Define a Direct Acyclic Graph (DAG)"),
                    actionButton("how2build_dag", "Help")
                ),
                if(!is.null(data$lambdas)){
                    tags$div(
                        tags$h3("New Edge"),
                        tags$div(class = "inline",
                            radioButtons(inputId = "dag_from", 
                            label = "From", 
                            inline = TRUE,
                            choices =  c("Root", options))
                        ),
                        tags$div(class = "inline",
                            radioButtons(inputId = "dag_to", 
                            label = " To ", 
                            inline = TRUE,
                            choices =  options)
                        ),
                        actionButton("add_edge", "Add edge"),
                        tags$h3("DAG table"),
                        DTOutput("dag_table"),
                        numericInput("dag_samples", "Total genotypes to sample", value = default_mhn_samples, min= 100, max= 100000, width = "50%"),
                        actionButton("resample_dag", "Sample from DAG")
                    )
                }
            )
        } else if (input$input2build == "matrix"){
            tags$div(
                tags$div(class = "flex",
                    tags$h3("2. Define input with a Matrix"),
                    actionButton("how2build_matrix", "Help")
                ),
                if(!is.null(data$thetas)){
                    tags$div(
                        tags$h3("Thetas table"),
                        DTOutput("thetas_table"),
                        numericInput("mhn_samples", "Total genotypes to sample", value = default_mhn_samples, min= 100, max= 100000, width = "50%"),
                        actionButton("resample_mhn", "Sample from MHN")
                    )
                }
            )
        }
    })
    output$change_freqs <- renderUI({
        if(input$input2build == "csd"){
            tags$div(class = "frame",
                tags$div(class = "flex",
                    tags$h3("3. Change frequencies"),
                    actionButton("display_help", "Help"),
                ),
                tags$div(id = "csd_table", 
                    DTOutput("csd_freqs")
                )
            )
        }
    })

    ## Controling dag builder
    dag_data <- reactive({
        all_gene_names <- c("Root", data$gene_names)
        edges <- which(data$dag == 1, arr.ind = TRUE)
        tmp_dag_parent_set <- data$dag_parent_set
        x <- length(tmp_dag_parent_set)
        ## I have to this weird thing because using data$gene_names does not work for some unkown reason
        names(tmp_dag_parent_set) <- all_gene_names[seq(2, x + 1)]
        dag_data <- data.frame(From = all_gene_names[edges[, "row"]]
            , To = all_gene_names[edges[, "col"]]
            , Relationship = tmp_dag_parent_set[edges[, "col"] - 1]
            , Lambdas = data$lambdas[edges[, "col"] - 1])
        dag_data
    })

    output$dag_table = DT::renderDT(
        dag_data(), escape = FALSE, selection = 'none', server = FALSE,
        rownames = FALSE,
        editable = list(target = "all", disable = list(columns = c(0, 1))),
        options = list(dom = 't', paging = FALSE, ordering = FALSE, 
            columnDefs = list(list(className = 'dt-center', 
            targets = "_all")))
    )
    
    ## Adding new edge
    observeEvent(input$add_edge, {
        from_gene <- input$dag_from
        to_gene <- input$dag_to
        if(is.null(input$dag_from) | is.null(input$dag_to)){
            error_message <<- "Both From and To options has to be defined and must be different"
            showModal(dataModal(error_message))
        } else if (input$dag_from == input$dag_to){
            error_message <<- "Both From and To options must be different"
            showModal(dataModal(error_message))
        } else{
            from_node <- ifelse(input$dag_from == "Root", "WT", input$dag_from)
            if(data$dag[from_node, input$dag_to] == 1){
                error_message <<- "That edge is already present"
                showModal(dataModal(error_message))
            } else if(data$dag[input$dag_to, from_node] == 1){
                error_message <<- "Relathionships cannot be bidirectional"
                showModal(dataModal(error_message))
            } else{
                tmp_dag <- data$dag
                tmp_dag [from_node, input$dag_to] = 1
                g <- graph_from_adjacency_matrix(tmp_dag, mode = "directed")
                if(is_dag(g)){
                    data$dag <- tmp_dag
                } else {
                    error_message <<- "This relathionship breaks the DAG. Revise it."
                showModal(dataModal(error_message))
                }
            }
        }
    })

    observeEvent(input$dag_table_cell_edit, {
        names(data$dag_parent_set) <- data$gene_names[1:length(data$dag_parent_set)]
        names(data$lambdas) <- data$gene_names[1:length(data$dag_parent_set)]
        info <- input$dag_table_cell_edit
        all_genes <- dag_data()$To

        ## Different lambdas
        new_lambdas <- as.integer(info[info["col"] == 3,"value"])
        old_lambdas <- dag_data()$Lambdas
        changed_genes <- all_genes[new_lambdas != old_lambdas
            & new_lambdas > 0]
        changed_lambdas <- new_lambdas[new_lambdas != old_lambdas 
            & new_lambdas > 0]

        tmp_lambdas <- data$lambdas
        tmp_lambdas[changed_genes] <- changed_lambdas 
        data$lambdas <- tmp_lambdas

        ## Remove relationships 
        new_relationships <- info[info["col"] == 3,]
        rels2remove <- new_relationships[new_relationships$value == 0, ]
        orig_data <- dag_data()
        if(nrow(rels2remove)){
            for(i in 1:nrow(rels2remove)){
                from2remove <- orig_data$From[rels2remove[i,]$row]
                to2remove <- orig_data$To[rels2remove[i,]$row]
                data$dag[from2remove, to2remove] <- 0
            }
        }
        
        ## Restructure the DAG
        number_of_parents2 <- colSums(data$dag)
        number_of_children2 <- rowSums(data$dag)
        for(i in colnames(data$dag)[-1]){
            if(number_of_children2[i] > 0 & number_of_parents2[i] == 0){
                ## We add link to WT
                data$dag["WT", i] <- 1
            }
        }

        ## Different relationships
        new_relationships <- info[info["col"] == 2,"value"]
        old_relationships <- orig_data$Relationship
        changed_genes <- all_genes[new_relationships != old_relationships]
        changed_relationships <- new_relationships[new_relationships != old_relationships]

        tmp_parent_set <- data$dag_parent_set
        tmp_parent_set[changed_genes] <- toupper(changed_relationships)
        ## Replace relationships that are not AND, OR, XOR
        ## Genes with only one parent with single relationship
        ## Default relationship for several parents is OR

        number_of_parents <- colSums(data$dag)
        for(i in unique(all_genes)){
            if(number_of_parents[[i]] <= 1){
                tmp_parent_set[[i]] <- "Single"
            } else if (!(tmp_parent_set[i] %in% c("AND", "OR", "XOR"))){
                tmp_parent_set[i] <- "OR"
            }
        }
        data$dag_parent_set <- tmp_parent_set

        # Resampling
        shinyjs::click("resample_dag")
    })

    observeEvent(input$how2build_matrix, {
      showModal(modalDialog(
        easyClose = TRUE,
        title = tags$h3("How to build a matrix"),
        tags$div(
            tags$p("A positive theta indicates that the presence of one the row gene makes column gene more likely. A negative means the opposite."),
            tags$p("Diagonal thetas indicates how likely is that event to be the first one (positive values likely, negative unlikely."),
            tags$p("Once the thetas are defined hit the 'Sample from MHN' to generate a sample."),
            tags$p("To make a sample we take into account multiplicative effects of all thetas"),
            tags$h3("How to modify the table"),
            tags$p("1. Double click in a Frequency cell to edit it"),
            tags$p("2. Press Tab to move to the next row"),
            tags$p("3. Use Shift + Enter to save changes"),
            tags$p("4. Set a frequency to 0 to remove a genotype"),
            tags$p("5. Type in the Search bar to filter genotypes")
            )
        )
      )
    })
    ## Building trm from dag

    observeEvent(input$resample_dag, {
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Running GuloMAM", value = 0)

        progress$inc(1/2, detail = "Doing sampling")
        shinyjs::disable("resample_dag")
        tmp_data <- list(edges = dag_data(), parent_set = data$dag_parent_set)
        trm <- cpm_access_genots_paths_w_simplified_relationships(tmp_data)$weighted_fgraph
        samples <- population_sample_from_trm(trm, input$dag_samples)
        process_data <- process_samples(samples, input$gene_number)
        tmp_csd <- process_data$frequencies
        tmp_csd <- tmp_csd[tmp_csd$Counts > 0, ]
        rownames(tmp_csd) <- tmp_csd$Genotype
        data$csd_freqs <- tmp_csd
        data$complete_csd <- freqs2csd(tmp_csd, data$gene_names)
        shinyjs::enable("resample_dag")
        progress$inc(1/2, detail = "Sampling Finished")

        if(input$select_csd != "user"){
            data$all_csd[[input$input2build]][[input$select_csd]]$data <- data$complete_csd
            data$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
            data$all_csd[[input$input2build]][[input$select_csd]]$trm <- trm
            data$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set <- data$dag_parent_set
            data$all_csd[[input$input2build]][[input$select_csd]]$type <- input$input2build
        }
    })
    ## Controling thetas

    output$thetas_table <- DT::renderDT(data$thetas, selection = 'none', server = TRUE, editable = list(target = "column")
        , rownames = TRUE,
        options = list(
            searching = FALSE, columnDefs = list(list(className = 'dt-center', 
            targets = "_all")), info = FALSE, paginate= FALSE),
    )

    observeEvent(input$thetas_table_cell_edit, {
        info <-input$thetas_table_cell_edit
        data$thetas <- editData(data$thetas, info, "thetas") 
        
        ## Resample based on changes
        shinyjs::click("resample_mhn")
        # mhn_data <- get_mhn_data(input$gene_number, input$mhn_samples, data$gene_names[1:input$gene_number], thetas = data$thetas)
        # data$csd_freqs <- mhn_data$samples
        # data$complete_csd <- freqs2csd(data$csd_freqs, data$gene_names)
        # data$trm <- mhn_data$trm


        # if(input$select_csd != "user"){
        #     data$all_csd[[input$select_csd]]$data <- data$complete_csd
        #     data$all_csd[[input$select_csd]]$thetas <- data$thetas
        #     data$all_csd[[input$select_csd]]$trm <- data$trm
        #     data$all_csd[[input$select_csd]]$type <- input$input2build
        # }
    })

    observeEvent(input$resample_mhn, {
        mhn_samples <- input$mhn_samples
        if(is.null(input$mhn_samples)) mhn_samples <- default_mhn_samples
        mhn_data <- get_mhn_data(input$gene_number, mhn_samples, data$gene_names[1:input$gene_number])
        data$csd_freqs <- mhn_data$samples
        data$complete_csd <- freqs2csd(data$csd_freqs, data$gene_names)
        data$thetas <- mhn_data$thetas
        data$trm <- mhn_data$trm

        if(input$select_csd != "user"){
            data$all_csd[[input$input2build]][[input$select_csd]]$type <- input$input2build
            data$all_csd[[input$input2build]][[input$select_csd]]$data <- data$complete_csd
            data$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
            data$all_csd[[input$input2build]][[input$select_csd]]$trm <- data$trm
        }
    })

    observe({
        ncol_trm <- ifelse(is.null(data$thetas), 0, ncol(data$thetas))
        tmp_ngene <- ifelse(is.null(input$gene_number), 0, input$gene_number)
        tmp_dataset_name <- ifelse(is.null(input$dataset_name), 0, input$dataset_name)
        tmp_select_csd <- ifelse(is.null(input$select_csd), 1, input$select_csd)
        if(tmp_dataset_name == tmp_select_csd & input$input2build == "matrix" 
            & tmp_ngene > ncol_trm 
            # & !is.null(data$thetas)
            ){
            mhn_samples <- input$mhn_samples
            if(is.null(input$mhn_samples)) mhn_samples <- default_mhn_samples
            mhn_data <- get_mhn_data(input$gene_number, mhn_samples, data$gene_names[1:input$gene_number])
            data$csd_freqs <- mhn_data$samples
            data$complete_csd <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
            data$thetas <- mhn_data$thetas
            data$trm <- mhn_data$trm
        } 
    })

    observeEvent(input$how2build_matrix, {
      showModal(modalDialog(
        easyClose = TRUE,
        title = tags$h3("How to build a matrix"),
        tags$div(
            tags$p("A positive theta indicates that the presence of one the row gene makes column gene more likely. A negative means the opposite."),
            tags$p("Diagonal thetas indicates how likely is that event to be the first one (positive values likely, negative unlikely."),
            tags$p("Once the thetas are defined hit the 'Sample from MHN' to generate a sample."),
            tags$p("To make a sample we take into account multiplicative effects of all thetas"),
            tags$h3("How to modify the table"),
            tags$p("1. Double click in a Frequency cell to edit it"),
            tags$p("2. Press Tab to move to the next row"),
            tags$p("3. Use Shift + Enter to save changes"),
            tags$p("4. Set a frequency to 0 to remove a genotype"),
            tags$p("5. Type in the Search bar to filter genotypes")
            )
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
            rownames(data$csd_freqs) <- data$csd_freqs$Genotype
            data$csd_freqs[genotype, ] <- c(genotype, genot_freq)
            data$csd_freqs[, 2] <- as.numeric(data$csd_freqs[, 2])
            ## Filtering out non-positive counts
            data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
            data$all_csd[[input$select_csd]]$data <- freqs2csd(data$csd_freqs, data$gene_names)
        }
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices =  lapply(1:input$gene_number, function(i)data$gene_names[i]), selected = NULL)
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
        data$all_csd[[input$select_csd]]$data <- freqs2csd(data$csd_freqs, data$gene_names)
    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        plot_genotypes_freqs(display_freqs())
    })

    ## Plot dag of dataset
    output$dag_plot <- renderPlot({
        if(input$input2build %in% c("csd", "dag")){
            tmp_dag <- NULL
            if(!is.null(data$dag) & !is.null(input$gene_number)) {
                num_genes <- ncol(data$dag) - 1
                if(num_genes >= input$gene_number){
                    genes2show <- input$gene_number
                } else {
                    genes2show <- num_genes
                }
                print(genes2show)
                tmp_dag <- data$dag[c("WT", data$gene_names[1:genes2show]), 
                         c("WT", data$gene_names[1:genes2show])]
            }
            plot_dag(tmp_dag
                , data$dag_parent_set[1:genes2show])
        }else if(input$input2build %in% c("matrix")){
            if(!is.null(data$thetas)){
                op <- par(mar=c(3, 3, 5, 3), las = 1)
                plot(data$thetas[1:input$gene_number, 1:input$gene_number], cex = 1.8, digits = 2, key = NULL
                    , axis.col = list(side = 3)
                    , xlab = "Effect of this (effector)"
                    , ylab = " on this (affected)"
                    , main = ""
                    , mgp = c(2, 1, 0))
                par(op)
                title("Theta Matrix", cex.main = 1.8)
            }
        }
    })

    ## Run CPMS
    observeEvent(input$analysis, {
        shinyjs::disable("analysis")
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Running GuloMAM", value = 0)

        data2run <- freqs2csd(display_freqs(), data$gene_names)
        progress$inc(1/4, detail = "Setting up data")
        Sys.sleep(0.5)
        progress$inc(2/4, detail = "Running CPMs")
        cpm_output <- all_methods_2_trans_mat(data2run)
        n_samples <- 100
        progress$inc(3/4, detail = paste("Running ", n_samples, " samples"))
        new_data <- sample_all_CPMs(cpm_output, n_samples, input$gene_number)
        progress$inc(4/4, detail = "Post processing data")
        Sys.sleep(0.5)
        new_data$MHN_f_graph <- new_data$MHN_transitionRateMatrix
        new_data$name <- input$select_csd
        all_cpm_out$output[[input$select_csd]] <- new_data
        shinyjs::enable("analysis")
        updateTabsetPanel(session, "navbar", selected = "result_viewer")
        updateRadioButtons(session, "select_cpm", selected = input$select_csd)
    })
    
    output$cpm_freqs <- DT::renderDT(genotype_freq_df(),  
        selection = 'none', server = TRUE
        , rownames = FALSE
        , options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE)
    )

    cpm_out <- readRDS("/home/pablo/CPM-SSWM-Sampling/guloMAM/inst/shiny-examples/sims.RDS")
    cpm_out$MHN_f_graph <- cpm_out$MHN_transitionRateMatrix
    
    all_cpm_out <- reactiveValues(output = list(user = cpm_out))

    output$sims <- renderUI({
        column_models2show <- floor(12 / length(input$cpm2show)) 

        lapply(input$cpm2show, function(mod){
            output[[sprintf("plot_sims_%s", mod)]] <- renderPlot({
                pl <- plot_model(all_cpm_out$output[[input$select_cpm]], mod)
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

        names(attribute_name) <- c("Transition Rate matrix", "Transitions", "Transition Probability Matrix")

        selected_plot_type <- attribute_name[input$data2plot]

        top_paths <- input$top_paths
        if(top_paths == 0) top_paths <- NULL
        
        lapply(input$cpm2show, function(mod){
            data2plot <- all_cpm_out$output[[input$select_cpm]][[sprintf("%s_%s", mod, 
                # "trans_mat"
                selected_plot_type
                )]]
            output[[sprintf("plot_sims2_%s", mod)]] <- renderPlot({
                pl <- plot_genot_fg(data2plot, 
                    all_cpm_out$output[[input$select_cpm]]$csd_data, 
                    all_cpm_out$output[[input$select_cpm]][[sprintf("%s_genotype_%s", mod, "freqs")]], 
                    top_paths = top_paths)
            })
            return(
                column(3,
                    plotOutput(sprintf("plot_sims2_%s", mod)))
            )
        })
    })

    output$csd <- renderPlot({
        plot_genotypes_freqs(get_csd(all_cpm_out$output[[input$select_cpm]]$csd_data))
    })

    observeEvent(input$modify_data, {
        data$csd_freqs <- sampledGenotypes(all_cpm_out$output[[input$select_cpm]]$csd_data)
        rownames(data$csd_freqs) <- data$csd_freqs$Genotype
        data$n_genes <- sum(colSums(all_cpm_out$output[[input$select_cpm]]$csd_data) > 0)
        updateTabsetPanel(session, "navbar",
            selected = "csd_builder"
        )
        selected <- ifelse(input$select_cpm %in% names(data$all_csd), input$select_cpm, "user")
        updateRadioButtons(session, "select_csd", selected = selected)
    })

    output$cpm_list <- renderUI({
        all_names <- unname(sapply(all_cpm_out$output, function(dataset) dataset$name))
        selected <- ifelse(is.null(input$select_csd), "user", input$select_csd)
        selected <- ifelse(input$select_csd %in% names(all_cpm_out$output),input$select_csd, "user")
        tagList(
            radioButtons(
                inputId = "select_cpm",
                label = "",
                selected = selected,
                choiceNames = names(all_cpm_out$output),
                choiceValues = names(all_cpm_out$output)
            )
      )
    })

    genotype_freq_df <- reactive({
        compare_cpm_freqs(all_cpm_out$output[[
        input$select_cpm
        ]])})

    ## Download button
    output$download_cpm <- downloadHandler(
        filename = "cpm.RDS",
        content = function(file) {
            saveRDS(all_cpm_out$output[[input$select_cpm]], file)
        }
    )

    ## Upload button
    observeEvent(input$output_cpms, {
        cpm_out <- readRDS(input$output_cpms$datapath)
        cpm_out$MHN_f_graph <- cpm_out$MHN_transitionRateMatrix
        if(is.null(cpm_out$name)) cpm_out$name <- "User Data"
        all_cpm_out$output[[cpm_out$name]] <- cpm_out
        updateRadioButtons(session, "select_cpm", selected = cpm_out$name)
    })
}
