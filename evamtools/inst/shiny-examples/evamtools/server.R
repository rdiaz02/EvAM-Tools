#'  I followed this link to structure the shiny app whithin the package
#'  https://deanattali.com/2015/04/21/r-package-shiny-app/

default_genes <- 3
max_genes <- 10
all_gene_names <- LETTERS[1: max_genes]
template_dag <- matrix(0, ncol= max_genes + 1, nrow = max_genes + 1)
rownames(template_dag) <- colnames(template_dag) <- c("WT", all_gene_names)
template_parent_set <- rep("Single", max_genes)
names(template_parent_set) <- all_gene_names
template_lambdas <- rep(1, max_genes)
names(template_lambdas) <- all_gene_names
template_thetas <- matrix(0, ncol = max_genes, nrow = max_genes)
rownames(template_thetas) <- colnames(template_thetas) <- all_gene_names
template_csd_freqs <- data.frame(Genotype = character(), Freq = integer())
template_csd_data <- matrix(0, ncol=3, nrow=0)

check_if_csd <- function(data){
    tmp_names <- c("data", "dag", "dag_parent_set", "gene_names", "name", "type", "thetas")
    types <- c("csd", "dag", "matrix")
    if(is.null(data)) return(FALSE)
    if(all(tmp_names %in% names(data))){
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
    if(nrow(data) == 0) return()
    par(las = 2, cex.main=1.6, cex.lab=1.5, cex.axis=1.2)
    barplot(data[, 2]
        , names = data$Genotype
        , ylab="Counts", main="Genotype Frequencies"
        , horiz = FALSE
        , panel.first=grid())
    grid(nx = NA, ny = NULL, col='gray', lwd = 2)
    ## TODO sort genotypes
}

plot_dag <- function(dag, parent_set){
    if (is.null(dag)) return()
    standard_relationship <- "gray73"
    colors_relationships <- c(standard_relationship, standard_relationship, "cornflowerblue", "coral2")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    dag <- igraph::graph_from_adjacency_matrix(dag, mode = "directed")
    dag <- igraph::decompose(dag)[[1]]
    ## Plotting data
    if(!is.null(parent_set)){
        for(i in igraph::E(dag)){
            igraph::E(dag)[i]$color <- colors_relationships[parent_set[[igraph::head_of(dag, igraph::E(dag)[i])$name]]]
        }
    } else igraph::E(dag)$color <- standard_relationship
        
    plot(dag
        , layout = igraph::layout.reingold.tilford
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
    csd <- NULL
    if(nrow(freqs) > 0){
        csd2 <- apply(
            freqs, 1
            , function(x){
                mut <- x[[1]]
                freq <- as.numeric(x[[2]])
                genot <- rep(0, length(gene_names))
                if(mut != "WT"){
                    mut <-  strsplit(mut, ", ")[[1]]
                    genot[which(gene_names %in% mut)] <- 1
                }
                csd <- matrix(rep(genot, freq), ncol= length(gene_names), byrow = TRUE)
                return(list(csd))
            })
        
        csd <- csd2[[1]][[1]]
        if(nrow(freqs) > 1){ 
            for (i in 2:nrow(freqs)){
                csd <- rbind(csd, csd2[[i]][[1]])
            }
        } 
        colnames(csd) <- gene_names
    }

    return(csd)
}

get_display_freqs <- function(freqs, n_genes, gene_names){
    if(nrow(freqs) == 0) return(template_csd_freqs)
    valid_gene_names <- c("WT", gene_names[1:n_genes])

    selected_rows <- sapply(freqs$Genotype, function(x){
        genes <- strsplit(x, ", ")[[1]]
        return(all(genes %in% valid_gene_names))
    })

    return(freqs[selected_rows, ])
}

get_csd <- function(complete_csd){
    if(is.null(complete_csd)) return(NULL)
    csd <- data.frame(OncoSimulR::sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    return(csd)
}

available_cpms <- function(data){
    data$csd_data <- NULL
    cpm_names <- unique(sapply(names(data), function(x) str_split(x, "_")[[1]][[1]]))
    return(cpm_names)
}

get_mhn_data <- function(n_genes, n_samples, gene_names, thetas = NULL){
    if(is.null(thetas)) thetas <- evamtools:::Random.Theta(n=n_genes)
    rownames(thetas) <- colnames(thetas) <- gene_names
    samples <- floor(evamtools:::Finite.Sample(evamtools:::Generate.pTh(thetas), n_samples)*n_samples)
    trm <- evamtools:::theta_to_trans_rate_3_SM(thetas,
                                    inner_transition = evamtools:::inner_transitionRate_3_1)
    state_names <- vapply(1:(ncol(trm)), function(x){
        x <- x - 1
        if(x == 0) state_name <- "WT"
        else state_name <- paste(gene_names[which(evamtools:::int2binary(x, n_genes) == 1)], collapse = ", ")
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
    colors_relationships <- c(standard_relationship, standard_relationship, "cornflowerblue", "coral2")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    
    model_data2plot <- evamtools:::process_data(cpm_output, mod)
    ## Plotting data
    if(!is.null(model_data2plot$dag_tree)) {

        if(!is.null(model_data2plot$parent_set)){
            for(i in names(model_data2plot$parent_set)){
                igraph::E(model_data2plot$dag_tree)[.to(i)]$color <- colors_relationships[model_data2plot$parent_set[[i]]]
            }
        } else igraph::E(model_data2plot$dag_tree)$color <- standard_relationship
        plot(model_data2plot$dag_tree
            , layout = igraph::layout.reingold.tilford
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

create_tabular_data <- function(data, type){
    available_methods <- c("Source", "OT", "CBN", "MHN", "HESBCN")
    # , "DBN", "MCCBN")
    if(type %in% c("freqs")){
        all_counts <- data.frame(Genotype = data[["MHN_genotype_freqs"]]$Genotype)
        for(name in names(data)){
            if(grepl("_genotype_freqs", name)
            #  & !grepl("^OT", name) ##Now we have genotypes frequencies for OT
             ){
                method_name <- strsplit(name, "_")[[1]][[1]]
                all_counts[[method_name]] <- data[[name]]$Counts
            } 
        }
        
        order_by_counts <- sort(rowSums(all_counts[-1]), 
        decreasing = TRUE, index.return = TRUE)$ix
    
        all_counts[order_by_counts, ]
        return(all_counts[order_by_counts, ])

    } else if(type %in% c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")){
        var2var <- c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")
        names(var2var) <- c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")
        
        var2use <- var2var[type]
        ## 1 Methods to compute
        methods2compute <- vapply(available_methods, function(x){
            if(!is.null(data[[sprintf("%s_%s", x, var2use)]])){
                return(!any(is.na(data[[sprintf("%s_%s", x, var2use)]])))
            }
            return(FALSE)
        }, logical(1))

        methods2compute <- names(methods2compute)[methods2compute]
        if(type == "genotype_transitions"){
            methods2compute <- setdiff(methods2compute, "OT")
        }
            
        ## 2 Fill the data frame
        all_genotypes <- matrix(0, nrow = 0, ncol = 2)
        all_methods <- matrix(0, nrow = 0, ncol = length(methods2compute))
        n_methods2compute <- length(methods2compute)
        base_vector <- rep(0, n_methods2compute)
        names(base_vector) <- methods2compute
        ## I use MHN because it has all the genotypes
        for(i in rownames(data[[sprintf("MHN_%s", var2use)]])){
            for(j in rownames(data[[sprintf("MHN_%s", var2use)]])){
                row_data <- sapply(methods2compute, function(x){
                    tryCatch({
                        return(data[[sprintf("%s_%s", x, var2use)]][i, j])
                    }, error = function(e){
                        return(0)
                    })
                })

                all_genotypes <- rbind(all_genotypes, c(i, j))
                all_methods <- rbind(all_methods, unname(row_data))
            }
        }
        
        selected_rows <- rowSums(abs(all_methods))>0
        all_methods <- round(all_methods[selected_rows, ], 2)
        all_genotypes <- all_genotypes[selected_rows, ]

        all_the_data <- data.frame(From = all_genotypes[, 1]
            , To = all_genotypes[, 2])
        for(i in 1:n_methods2compute){
            all_the_data[[methods2compute[i]]] <- all_methods[, i]
        }

        colnames(all_the_data) <- c("From", "To", methods2compute)

        order_by_counts <- sort(rowSums(all_methods), 
        decreasing = TRUE, index.return = TRUE)$ix
    
        return(all_the_data[order_by_counts, ])

    } else if(type %in% c("lambdas")){
        lambda_field <- c("Lambdas", "OT_edgeWeight", "rerun_lambda", "Lambdas", "lambda", "Thetas")
        names(lambda_field) <- c("Source", "OT", "CBN", "HESBCN", "MCCBN")

        gene_names <- sort(unique(data$OT_model$To))
        all_counts <- data.frame(Gene = gene_names)
        for(name in names(data)){
            if(grepl("_model", name)){
                method_name <- strsplit(name, "_")[[1]][[1]]
                if(!is.null(data[[name]]) & !is.na(data[[name]])){
                    tmp_data <- data[[name]][[lambda_field[method_name]]]
                    names(tmp_data) <- data[[name]]$To
                    all_counts[[method_name]] <- round(tmp_data[all_counts$Gene], 2)
                }
            }
        }

        return(all_counts)
    }

    order_by_counts <- sort(rowSums(all_counts[-1]), 
        decreasing = TRUE, index.return = TRUE)$ix
    
    all_counts[order_by_counts, ]
    return(to_return)
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

standarize_dataset<-function(data){
    new_data <- list()

    if(!is.null(colnames(data$data))){
        new_data$gene_names <- c(colnames(data$data)
            , LETTERS[(ncol(data$data) + 1) : max_genes ])
    } else {
        new_data$gene_names <- all_gene_names
    }
    new_data$name <- data$name  

    if(is.null(data$lambdas)) {
        new_data$lambdas <- template_lambdas
    }else{
        new_lambdas <- template_lambdas
        new_lambdas[1:length(data$lambdas)] <- data$lambdas
        names(new_lambdas) <- new_data$gene_names
        new_data$lambdas <- new_lambdas
    }

    if(is.null(data$dag_parent_set)) {
        new_data$dag_parent_set <- template_parent_set
    }else {
        new_parent_set <- template_parent_set
        new_parent_set[1:length(data$dag_parent_set)] <- data$dag_parent_set
        names(new_parent_set) <- new_data$gene_names
        new_data$dag_parent_set <- new_parent_set
    }

    if(is.null(data$dag)) {
        new_data$dag <- template_dag
    } else {
        to_keep <- which(colSums(data$dag)>0 | rowSums(data$dag)>0)
        tmp_dag <- data$dag[to_keep, to_keep]
        n_genes <- ncol(data$dag)
        new_dag <- template_dag 
        new_dag[to_keep, to_keep] <- tmp_dag
        new_data$dag <- new_dag
    }

    if(is.null(data$thetas)) {
        new_data$thetas <- template_thetas
    } else {
        new_data$thetas <- data$thetas
        rownames(new_data$thetas) <- colnames(new_data$thetas) <-
            new_data$gene_names[1:ncol(new_data$thetas)]
    }

    if(is.null(data$data)) {
        new_data$data <- template_csd_data
    } else{
        new_data$data <- data$data
        colnames(new_data$data) <- new_data$gene_names[1:ncol(new_data$data)]
    }
    return(new_data)
}

standarize_all_datasets <- function(datasets){
    all_new_data <- list()
    for(i in c("csd", "dag", "matrix")){
        tmp_data <- datasets[[i]] 
        for(j in names(tmp_data)) all_new_data[[i]][[j]] <-     
            standarize_dataset(tmp_data[[j]])
    }
    return(all_new_data)
}

server <- function(input, output, session) {
    examples_csd$csd <- examples_csd$csd[1:5]
    all_csd_data <- standarize_all_datasets(examples_csd)
    min_genes <- 2
    max_genes <- 10
    default_mhn_samples <- 5000
    keep_dataset_name <- FALSE

    error_message <- NULL
    last_visited_pages <- list(csd = "User", dag = "User", matrix = "User")
    last_visited_cpm <- "user"

    datasets <- reactiveValues(
        all_csd = all_csd_data
    )
    data <- reactiveValues(
        csd_freqs =  template_csd_freqs
        , data = NULL
        , dag = template_dag
        , dag_parent_set = template_parent_set
        , lambdas = template_lambdas
        , thetas = template_thetas
        , trm = NULL
        , n_genes = 3
        , gene_names = LETTERS[1: max_genes]
    )
    
    display_freqs <- reactive({
        get_display_freqs(data$csd_freqs, input$gene_number, data$gene_names)
    })

    ## Upload data
    observeEvent(input$csd, {
        if(grepl(".csv", input$csd$datapath)){
            dataset_name <- strsplit(strsplit(input$csd$name, ".csv")[[1]], "_")[[1]][[1]]
            tmp_data <- list()
            tmp_data$data <- read.csv(input$csd$datapath)
            if(check_if_csd(tmp_data$data)){
                datasets$all_csd[["csd"]][[dataset_name]] <- standarize_dataset(tmp_data)
                datasets$all_csd[["csd"]][[dataset_name]]$name <- dataset_name
                datasets$all_csd[["csd"]][[dataset_name]]$gene_names <- c(colnames(tmp_data$data)
                , LETTERS[(ncol(tmp_data$data) + 1): max_genes])
                # keep_dataset_name <<- dataset_name 
                last_visited_pages["csd"] <<- dataset_name
                updateRadioButtons(session, "input2build", selected = "csd")
                updateRadioButtons(session, "select_csd", selected = dataset_name)
            } else {
                error_message <<- "Your csv data can not be loaded. Make sure it only contains 0 and 1."
                showModal(dataModal(error_message))
            }
        } else if(grepl(".rds", input$csd$datapath, ignore.case = TRUE)){
            tmp_data <- readRDS(input$csd$datapath)
            if(check_if_csd(tmp_data$data)){
                last_visited_pages[tmp_data$type] <<- tmp_data$name
                datasets$all_csd[[tmp_data$type]][[tmp_data$name]] <- tmp_data
                updateRadioButtons(session, "input2build", selected = tmp_data$type)
                updateRadioButtons(session, "select_csd", selected = tmp_data$name)
            } else {
                error_message <<- "There was a problem when checking your .rds file. Make sure it containis $type (either 'csd', 'dag', or 'matrix'), $data only with 0 and 1"
                showModal(dataModal(error_message))
            }
        }
    })

    observeEvent(input$input2build, {
        if(keep_dataset_name %in% names(datasets$all_csd)){
            updateRadioButtons(session, "select_csd", selected = keep_dataset_name)
            keep_dataset_name <- FALSE
        }else{
            updateRadioButtons(session, "select_csd", selected = last_visited_pages[[input$input2build]])
        }
    })

    ## Define dataset name
    output$dataset_name <- renderUI({
        tags$div(class="inlin2",
            textInput(inputId="dataset_name", "Give your dataset a name", value = input$select_csd)
        )
    })

    listen2dataset_name <- reactive({
        list(input$dataset_name, input$save_csd_data)
    })

    ## Saving dataset
    observeEvent(input$save_csd_data,{
        ## 1 save dataset to list after user data
        if(!(input$dataset_name %in% names(datasets$all_csd[[input$input2build]]))){
            datasets$all_csd[[input$input2build]][[input$dataset_name]]$name <- input$dataset_name

            if(nrow(data$csd_freqs) > 0){
                datasets$all_csd[[input$input2build]][[input$dataset_name]]$data <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
            }

            tmp_data <- list(
                data = data$data
                , dag = data$dag
                , gene_names = data$gene_names
                , dag_parent_set = data$dag_parent_set
                , lambdas = data$lambdas
                , thetas = data$thetas
                , trm = data$trm
                , name = input$dataset_name) 
            datasets$all_csd[[input$input2build]][[input$dataset_name]] <- tmp_data

            tmp_data_2 <- datasets$all_csd[[input$input2build]]
            datasets$all_csd[[input$input2build]] <- 
                c(tmp_data_2["User"]
                , tmp_data_2[input$dataset_name]
                , tmp_data_2[which(!(names(tmp_data_2) %in% c("User",           
                    input$dataset_name, names(all_csd_data))))]
                , tmp_data_2[which(names(datasets$all_csd[[input$input2build]]) %in% names(all_csd_data))]
            )

            ## 2 restore default values
            try({
                datasets$all_csd[[input$input2build]][[input$select_csd]] <- all_csd_data[[input$input2build]][[input$select_csd]]
            })

            ## 3 update selected entry
            updateRadioButtons(session, "select_csd", selected = input$dataset_name)

            shinyjs::disable("save_csd_data")
        }
    })
    
    # observeEvent(listen2dataset_name(), {
    observeEvent(input$dataset_name, {
        dataset_name <- ifelse(is.null(input$dataset_name), "", input$dataset_name)
        if(dataset_name != "" &
            !(dataset_name %in% names(datasets$all_csd[[input$input2build]]))) {
            shinyjs::enable("save_csd_data")
        }else if(dataset_name %in% names(datasets$all_csd[[input$input2build]])){
            shinyjs::disable("save_csd_data")
        }
    })
    
    ## Download csd button
    output$download_csd <- downloadHandler(
        filename = function() sprintf("%s_csd.RDS", input$select_csd),
        content = function(file) {
            tmp_data <- datasets$all_csd[[input$input2build]][[input$select_csd]] 
            tmp_data$type <- input$input2build
            saveRDS(tmp_data, file)
        }
    )

    observeEvent(input$select_csd, {
        last_visited_pages[[input$input2build]] <<- input$select_csd
    })

    ## Display List of availabe CSD 
    output$csd_list <- renderUI({
        all_names <- c()
        all_choice_names <- c()
        for (i in names(datasets$all_csd[[input$input2build]])){
            tmp_data <- datasets$all_csd[[input$input2build]][[i]]
            all_names <- c(all_names, i)
            all_choice_names <- c(all_choice_names, tmp_data$name)
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

    toListen <- reactive({
        list(input$select_csd, input$input2build)
    })

    observeEvent(toListen(), {
        ## Cleaning stuf
        selected <- last_visited_pages[[input$input2build]]
        tmp_data <- datasets$all_csd[[input$input2build]][[selected]]
        data$gene_names <- tmp_data$gene_names
        data$data <- tmp_data$data

        shinyjs::disable("analysis")

        if(nrow(data$data) > 0){
            data$csd_freqs <- get_csd(data$data) 
            shinyjs::enable("analysis")
        } else{
            data$csd_freqs <- template_csd_freqs
        }

        data$dag <- tmp_data$dag
        data$dag_parent_set <- tmp_data$dag_parent_set
        data$lambdas <- tmp_data$lambdas
        data$thetas <- tmp_data$thetas
        data$trm <- tmp_data$trm
        data$name <- tmp_data$name

        if(input$input2build == "dag"){
            to_keep <- length(which(colSums(data$dag)>0 | rowSums(data$dag)>0)) - 1
            n_genes <- ifelse(to_keep < 1 , default_genes, to_keep)
        } else if(input$input2build == "matrix"){
            n_genes <- length(which(colSums(abs(data$thetas))>0 
                | rowSums(abs(data$thetas))>0))
            n_genes <- ifelse(n_genes <= 0, 3, n_genes)
        } else if (input$input2build == "csd"){
            n_genes <- ncol(data$data)
        }

        updateNumericInput(session, "gene_number", value = n_genes)
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices = lapply(1:n_genes, function(i)data$gene_names[i]), selected = NULL)
    })

    observeEvent(input$change_gene_names, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("Modify gene names"),
            tags$div(class = "inlin2",
                textInput(inputId = "new_gene_names", "Genes names", 
                    value = paste(data$gene_names[1:input$gene_number], 
                    collapse = ", ")
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
        data$gene_names <- c(
            new_gene_names
            , LETTERS[(length(new_gene_names) + 1):max_genes]
        )

        ## Rename stuff
        colnames(data$data) <- data$gene_names[1:ncol(data$data)]
        if(nrow(data$data)>0) {
            data$csd_freqs <- get_csd(data$data)
        }
        names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
        rownames(data$dag) <- colnames(data$dag) <- c("WT", data$gene_names)
        rownames(data$thetas) <- colnames(data$thetas) <- data$gene_names[1:ncol(data$thetas)]

        tmp_data <- list(data = data$data, dag = data$dag
            , name = data$name
            , dag_parent_set = data$dag_parent_set, lambdas = data$lambdas
            , thetas = data$thetas, trm = data$trm)

        datasets$all_csd[[input$input2build]][[input$select_csd]] <- tmp_data
        
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

    ## Advanced option for running evamtools
    observeEvent(input$advanced_options, {
        showModal(modalDialog(
            size = "l",
            easyClose = TRUE,
            title = tags$h3("Advanced options?"),
            tags$div(
                numericInput("num_steps", "Sampling steps", 10000
                    , min = 0, max = 1000000, step = 1000, width="100%"),
                # checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("HyperTRAPS", "MCCBN"), choiceValues = c("hypertraps", "mccbn")),
                checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("MCCBN"), choiceValues = c("mccbn")),
                tags$h4("DISCLAIMER: MCCBN may take hours to run")
                # tags$h4("DISCLAIMER: Both HyperTraps and MCCBN may take hours to run")
            )
            )
        )
    ## TODO this has no effect so far 
    })
    
    ## Define number of genes
    output$genes_number <- renderUI({
        val <- ifelse(is.null(data$n_genes), 3, data$n_genes)
        tags$div(class="inlin",
            sliderInput("gene_number", "Number of genes",
                value = val, max = max_genes, min = min_genes, step = 1)
        )
    })

    ## Define new genotype
    output$define_genotype <- renderUI({
        n_genes <- ifelse(is.null(input$gene_number), 3, input$gene_number)
        options <- data$gene_names[1:n_genes]
        if(input$input2build == "csd"){
            tags$div(
                tags$h3("2. Add new genotypes"),
                    tags$div(class = "inline",
                        checkboxGroupInput(inputId = "genotype", 
                            label = "Mutations", 
                            choices = options)
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
                        DT::DTOutput("dag_table"),
                        numericInput("dag_samples", "Total genotypes to sample", value = default_mhn_samples, min= 100, max = 10000, step = 100, width = "50%"),
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
                        DT::DTOutput("thetas_table"),
                        numericInput("mhn_samples", "Total genotypes to sample", value = default_mhn_samples, min= 100, max= 10000, step = 100, width = "50%"),
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
                    DT::DTOutput("csd_freqs")
                )
            )
        }
    })

    ## ÐAG builder
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
        } else if(data$dag["WT", to_gene] == 1){
            error_message <<- "A direct children of Root cannot have multiple parents"
            showModal(dataModal(error_message))
        } else{
            from_node <- ifelse(input$dag_from == "Root", "WT", input$dag_from)
            if(data$dag[from_node, input$dag_to] == 1){
                error_message <<- "That edge is already present"
                showModal(dataModal(error_message))
            } else if(data$dag[input$dag_to, from_node] == 1){
                error_message <<- "Relationships cannot be bidirectional"
                showModal(dataModal(error_message))
            } else{
                tmp_dag <- data$dag
                tmp_dag [from_node, input$dag_to] = 1
                g <- igraph::graph_from_adjacency_matrix(tmp_dag, mode = "directed")
                if(igraph::is_dag(g)){
                    data$dag <- tmp_dag
                } else {
                    error_message <<- "This relationship breaks the DAG. Revise it."
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

    })

    listen2dag_change <- reactive({
        list(input$dag_table_cell_edit, input$add_edge)
    })

    observeEvent(listen2dag_change(), {
        info <- input$dag_table_cell_edit
        orig_data <- dag_data()
        all_genes <- dag_data()$To

        ## Restructure the DAG
        number_of_parents <- colSums(data$dag)
        number_of_children <- rowSums(data$dag)
        for(i in colnames(data$dag)[-1]){
            if(number_of_children[i] > 0 & number_of_parents[i] == 0){
                ## We add link to WT
                data$dag["WT", i] <- 1
            }
        }

        ## Different relationships
        ## Genes with only one parent: "Single" relationship
        ## Default relationship for several parents is OR
        number_of_parents <- colSums(data$dag)[-1]
        tmp_parent_set <- number_of_parents
        new_relationships <- info[info["col"] == 2,"value"]
        names(new_relationships) <- info[info["col"] == 1,"value"]
        for(idx in c(1:length(new_relationships))){
            tmp_parent_set[names(new_relationships[idx])] <- new_relationships[idx]
        }

        tmp_parent_set[!(tmp_parent_set %in% c("Single", "AND", "OR", "XOR"))] <- "OR"
        tmp_parent_set[number_of_parents <= 1] <- "Single"

        data$dag_parent_set <- tmp_parent_set
    })

    ## Building trm from dag
    observeEvent(input$resample_dag, {
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Running evamtools", value = 0)

        progress$inc(1/2, detail = "Doing sampling")
        shinyjs::disable("resample_dag")
        tmp_data <- list(edges = dag_data(), parent_set = data$dag_parent_set)
        trm <- evamtools:::cpm2tm(tmp_data)$weighted_fgraph
        samples <- evamtools:::population_sample_from_trm(trm, input$dag_samples)
        process_data <- evamtools:::process_samples(samples, input$gene_number, data$gene_names[1:input$gene_number])
        tmp_csd <- process_data$frequencies
        tmp_csd <- tmp_csd[tmp_csd$Counts > 0, ]
        rownames(tmp_csd) <- tmp_csd$Genotype
        data$csd_freqs <- tmp_csd
        data$data <- freqs2csd(tmp_csd,data$gene_names[1:input$gene_number])
        shinyjs::enable("resample_dag")
        progress$inc(1/2, detail = "Sampling Finished")

        datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
        datasets$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
        datasets$all_csd[[input$input2build]][[input$select_csd]]$trm <- trm
        datasets$all_csd[[input$input2build]][[input$select_csd]]$lambdas <- data$lambdas
        datasets$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set <- data$dag_parent_set
    
        shinyjs::enable("analysis")
    })
    
    ## Help for DAG building
    observeEvent(input$how2build_dag, {
      showModal(modalDialog(
        easyClose = TRUE,
        title = tags$h3("How to build a DAG"),
        tags$div(
            tags$p("Select a parent node and child node and hit 'Add edge'."),
            tags$p("But, TAKE CARE"),
            tags$p("Some edge wont be allowed if: "),
            tags$li("That edge is already present"),
            tags$li("It introduces cycles"),
            tags$p("TO REMOVE EDGES: set a lambda of the relationship to 0"),
            tags$p("But, TAKE CARE"),
            tags$p("Breaking edges might reestructure the DAG:"),
            tags$li("If a node has no parent, it will be assigned as descendent of WT")
            )
        )
      )
    })

    ## MHN
    ## Controling thetas
    output$thetas_table <- DT::renderDT(data$thetas[1:input$gene_number, 1:input$gene_number]
        , selection = 'none', server = TRUE
        , editable = list(target = "all", disable = list(columns = c(0)))
        , rownames = TRUE,
        options = list(
            searching = FALSE, columnDefs = list(list(className = 'dt-center', 
            targets = "_all")), info = FALSE, paginate= FALSE),
    )

    observeEvent(input$thetas_table_cell_edit, {
        info <-input$thetas_table_cell_edit
        data$thetas[1:input$gene_number, 1:input$gene_number] <- 
            DT::editData(data$thetas[1:input$gene_number, 1:input$gene_number], info, "thetas") 
        
        datasets$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
        ## Resample based on changes
        shinyjs::click("resample_mhn")
    })

    observeEvent(input$resample_mhn, {
        mhn_data <-get_mhn_data(input$gene_number, input$mhn_samples, 
            data$gene_names[1:input$gene_number], thetas = data$thetas[1:input$gene_number, 1:input$gene_number])
        data$csd_freqs <- mhn_data$samples
        data$data <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
        # data$thetas <- mhn_data$thetas
        # data$trm <- mhn_data$trm

        datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
        # datasets$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
        datasets$all_csd[[input$input2build]][[input$select_csd]]$trm <- data$trm
        shinyjs::enable("analysis")
    })

    observeEvent(input$how2build_matrix, {
      showModal(modalDialog(
        easyClose = TRUE,
        title = tags$h3("How to build a matrix"),
        tags$div(
            tags$p("Positive theta: gene i in row makes gene j in column more likely. A negative means the opposite."),
            tags$p("Diagonal theta: likelihood of that event i to be the first one (positive values likely, negative unlikely."),
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

    ## Working with raw CSD
    observeEvent(input$genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genot_freq <- data$csd_freqs[, 2][data$csd_freqs[, 1] == genotype]
        updateNumericInput(session, "genotype_freq", value = genot_freq)
    }, ignoreNULL = FALSE)

    observeEvent(input$add_genotype, {
        genotype <- paste(input$genotype, collapse = ", ")
        genotype <- ifelse(genotype == "", "WT", genotype)
        genot_freq <- ifelse(is.na(input$genotype_freq), -1, input$genotype_freq)

        if(genot_freq >= 0){
            data$csd_freqs[genotype, ] <- c(genotype, genot_freq)
            rownames(data$csd_freqs) <- data$csd_freqs$Genotype
            data$csd_freqs[, 2] <- as.numeric(data$csd_freqs[, 2])
            ## Filtering out non-positive counts
            data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
            data$data <- datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
            shinyjs::enable("analysis")
        }
        updateNumericInput(session, "genotype_freq", value = NA)
        updateCheckboxGroupInput(session, "genotype", label = "Mutations", 
            choices = lapply(1:input$gene_number, function(i)data$gene_names[i]), selected = NULL)
    })
    
    ## Genotypes table
    output$csd_freqs <- DT::renderDT(display_freqs(), selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
        , rownames = FALSE,
        options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
    )

    observeEvent(input$csd_freqs_cell_edit, {
        info <- input$csd_freqs_cell_edit
        info[ , "col"] <- 2
        data$csd_freqs <- DT::editData(data$csd_freqs, info, "csd_freqs") 
        ## Filtering out non-positive counts
        data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
        data$data <- datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
    })

    ## Plot histogram of genotypes
    output$plot <- renderPlot({
        plot_genotypes_freqs(display_freqs())
    })

    ## Plot dag of dataset
    output$dag_plot <- renderPlot({
        if(input$input2build %in% c("csd", "dag")){
            tmp_dag <- NULL
            if(sum(data$dag)>0 & !is.null(input$gene_number)) {
                num_genes <- ncol(data$dag) - 1
                if(num_genes >= input$gene_number){
                    genes2show <- input$gene_number
                } else {
                    genes2show <- num_genes
                }
                tmp_dag <- data$dag[c("WT", data$gene_names[1:genes2show]), 
                         c("WT", data$gene_names[1:genes2show])]
            }
            plot_dag(tmp_dag
                , data$dag_parent_set[1:genes2show])
        }else if(input$input2build %in% c("matrix")){
            if(!is.null(data$thetas) 
                && length(data$thetas[1:input$gene_number, 1:input$gene_number])>0){
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

    ## Run CPMs
    observeEvent(input$analysis, {
        ## Calculate TRM for DAG and for matrices
        
        source_trm <- NULL
        
        if(input$input2build == "dag"){
            tmp_data <- list(edges = dag_data(), parent_set = data$dag_parent_set)
            source_trm <- evamtools:::cpm2tm(tmp_data)$weighted_fgraph
        }else if(input$input2build == "matrix"){
            source_trm <- evamtools:::theta_to_trans_rate_3_SM(data$thetas[1:input$gene_number, 1:input$gene_number],
                                    inner_transition = evamtools:::inner_transitionRate_3_1)
        }

        ## TODO: Put everything with a try and display a modal if I have an error
        shinyjs::disable("analysis")
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Running evamtools", value = 0)

        data2run <- freqs2csd(display_freqs(), data$gene_names[1:input$gene_number])
        progress$inc(1/5, detail = "Setting up data")
        Sys.sleep(0.5)
        progress$inc(2/5, detail = "Running CPMs")
        do_MCCBN <- "mccbn" %in% input$more_cpms
        methods <- c("CBN", "OT", "HESBCN", "MHN")
        # , "DBN"
        # if(do_MCCBN) methods <- c(methods, "MCCBN")
        cpm_output <- evam(data2run, methods = methods)
        ## To see Source data in the results section
        if(input$input2build != "csd"){
            ## FIXME: rename f_graph to trans_rate_mat: ???
            cpm_output$Source_trans_rate_mat <- source_trm
            cpm_output$Source_trans_mat <- evamtools:::rowScaleMatrix(source_trm)
        }
        if(input$input2build == "dag"){
            cpm_output$Source_model <- dag_data()
            cpm_output$Source_parent_set <- data$dag_parent_set[1:input$gene_number]
        } else if(input$input2build == "matrix"){
            cpm_output$Source_theta <- data$thetas[1:input$gene_number
                , 1:input$gene_number]
        }
        # n_samples <- 100 ## To be replaced by input$something
        n_samples <- 100 
        progress$inc(3/5, detail = paste("Running ", n_samples, " samples"))
        new_data <- sample_all_CPMs(cpm_output, n_samples
            , input$gene_number, data$gene_names[1:input$gene_number])
        progress$inc(4/5, detail = "Post processing data")
        Sys.sleep(0.5)
        ## FIXME: rename f_graph to trans_rate_mat?
        # new_data$MHN_f_graph <- new_data$MHN_trans_rate_mat
        new_data$OT_f_graph <- NULL
        orig_data <- list(data = data2run, name = data$name
            , type = input$input2build, gene_names = data$gene_names
            , thetas = data$thetas, lambdas = data$lambdas
            , dag = data$dag, dag_parent_set = data$dag_parent_set)
        new_data$name <- input$select_csd
        new_data$source_data <- orig_data

        ##CPM output name
        result_index <- length(grep(sprintf("^%s", input$select_csd), names(all_cpm_out$output)))
        result_name <- ifelse(result_index == 0
            , input$select_csd
            , sprintf("%s__%s", input$select_csd, result_index))

        all_cpm_out$output[[result_name]] <- new_data
        last_visited_cpm <<- result_name
        updateRadioButtons(session, "select_cpm", selected = result_name)
        progress$inc(5/5, detail = "You can see your result by going to the Results tab")
        Sys.sleep(1)
        shinyjs::enable("analysis")
        
        ## Maybe I do no want this
        updateTabsetPanel(session, "navbar", selected = "result_viewer")
        updateRadioButtons(session, "select_cpm", selected = result_name)

        # selected <- ifelse(is.null(input$select_cpm), last_visited_cpm, input$select_cpm)
        # if(is.null(all_cpm_out$output[[selected]]$Source_genotype_transitions)){
        #     updateCheckboxGroupInput(session, "cpm2show", 
        #         selected = setdiff(c(input$cpm2show), "Source"))
        #     shinyjs::disable(selector = "#cpm2show input[value='Source']")
        # }else{
        #     shinyjs::enable(selector = "#cpm2show input[value='Source']")
        # }
    })
    
    output$cpm_freqs <- DT::renderDT(genotype_freq_df(),  
        selection = 'none', server = TRUE
        , rownames = FALSE
        , options = list(
            columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE)
    )

    ## FIXME get sample data from namespace
    cpm_out <- and_cpm_with_simulations
    ## FIXME: rename to trans_rate_mat??
    # cpm_out$MHN_f_graph <- cpm_out$MHN_trans_rate_mat
    
    all_cpm_out <- reactiveValues(output = list(user = cpm_out))

    # all_cpm_out <- reactiveValues(output = list())

    output$sims <- renderUI({
        if(length(all_cpm_out$output) > 0){
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
        }
    })

    output$sims2 <- renderUI({
        if(length(all_cpm_out$output) > 0){
            ## Enabling donwload button
            shinyjs::enable(selector = "#download_cpm")

            ## Main display
            column_models2show <- floor(12 / length(input$cpm2show)) 
            selected_plot_type <- input$data2plot

            ## To make all plots of the same type comparable
            max_edge <- 0
            min_edge <- 0
            if(!(is.null(selected_plot_type))){
                if (selected_plot_type %in% c("trans_mat", "td_trans_mat")){
                    min_edge <- 0
                    max_edge <- 1
                } else if (selected_plot_type == "genotype_transitions") {
                    for(i in input$cpm2show){
                        if (i != "OT"){
                            tmp_data <- all_cpm_out$output[[input$select_cpm]][[sprintf("%s_%s", i, selected_plot_type)]]
                            
                            tmp_max_edge <- max(tmp_data)
                            max_edge <- ifelse(tmp_max_edge > max_edge
                                , tmp_max_edge
                                , max_edge)
                            
                            tmp_min_edge <- min(tmp_data)
                            mom_edge <- ifelse(tmp_min_edge > min_edge
                                , tmp_min_edge
                                , min_edge)
                        }
                    }
                } else if (selected_plot_type == "trans_rate_mat") {
                    max_edge <- min_edge <- NULL
                }
            } else {
                max_edge <- min_edge <- NULL
            }

            lapply(input$cpm2show, function(mod){

                data2plot <- all_cpm_out$output[[input$select_cpm]][[
                    sprintf("%s_%s", mod, selected_plot_type)]]
                
                if(selected_plot_type == "td_trans_mat"){
                    diag(data2plot) <- 0
                    data2plot <- drop0(data2plot, 0)
                }

                if(is.null(data2plot)){
                    output[[sprintf("plot_sims2_%s", mod)]] <- renderPlot({})
                } else{
                    output[[sprintf("plot_sims2_%s", mod)]] <- renderPlot({
                        pl <- evamtools:::plot_genot_fg(data2plot 
                            , all_cpm_out$output[[input$select_cpm]]$csd_data 
                            , all_cpm_out$output[[input$select_cpm]][[sprintf("%s_genotype_%s", mod, "freqs")]] 
                            , freq2label = input$freq2label
                            , top_paths = 5
                            , max_edge = max_edge
                            , min_edge = min_edge)
                    })
                }
                return(
                    column(3,
                        plotOutput(sprintf("plot_sims2_%s", mod)))
                )
            })
        } else {
            ## Disabling donwload button
            shinyjs::disable(selector = "#download_cpm")

            return(tags$h3("There are not results to show yet. Go to the input tab, select a dataset and hit the 'Run evamtools!' button"))
        }
    })

    ## Go back to input to work again with the data
    observeEvent(input$modify_data, {
        if(length(all_cpm_out$output) > 0){
            tmp_data <- all_cpm_out$output[[input$select_cpm]]$source_data
            dataset_name <- strsplit(input$select_cpm, "__")[[1]][[1]]
            dataset_type <- tmp_data$type
            last_visited_pages[[tmp_data$type]] <<- dataset_name
            tmp_data <- datasets$all_csd[[tmp_data$type]][[dataset_name]] <- standarize_dataset(tmp_data)
            
            data <- tmp_data 
            data$csd_freqs <- get_csd(tmp_data$data)
            data$n_genes <- ncol(data$data)
            
            updateNumericInput(session, "gene_number", value = data$n_genes)
            updateTabsetPanel(session, "navbar",
                selected = "csd_builder")
            updateRadioButtons(session, "input2build", selected = dataset_type)
            updateRadioButtons(session, "select_csd", selected = dataset_name)
        }
    })

    output$customize <- renderUI({
            tags$div(class = "frame",
                tags$h3("2. Customize the visualization"),
                tags$div(class = "inline",
                  checkboxGroupInput(inputId = "cpm2show", 
                      label = "CPMs to show", 
                      choices = c("Source", "OT", "CBN", "MHN", "HESBCN"),
                    #   , "DBN")
                    # , "MCCBN"),
                      selected = c("CBN", "MHN", "HESBCN")),
                      
                tags$div(class = "inline",
                  radioButtons(inputId = "data2plot", 
                      label = "Data to show", 
                      choiceNames =  c("Transition rate matrix", "Transitions",
                                       "Transition matrix", "Time-discretized transition matrix"),
                      choiceValues = c("trans_rate_mat", "genotype_transitions",
                                       "trans_mat", "td_trans_mat"),
                      selected = "genotype_transitions"
                      )
                    ),
                ),
              tags$p("Label genotypes with frequency bigger than:"),
              tags$div(id="freq2label-wrap",
                sliderInput("freq2label", "", width = "500px",
                  value = 0.05, max = 1, min = 0, step = 0.05)
              )
            )
        # }
    })

    output$cpm_list <- renderUI({
        all_names <- unname(sapply(all_cpm_out$output, function(dataset) dataset$name))

        if(length(all_names) > 0){
            selected <- ifelse(is.null(input$select_csd), "user", input$select_csd)
            selected <- ifelse(input$select_csd %in% names(all_cpm_out$output),input$select_csd, "user")
        
            tagList(
                radioButtons(
                    inputId = "select_cpm",
                    label = "",
                    selected = last_visited_cpm,
                    choiceNames = names(all_cpm_out$output),
                    choiceValues = names(all_cpm_out$output)
                )
            )
        }
    })

    output$csd <- renderPlot({
        plot_genotypes_freqs(get_csd(all_cpm_out$output[[input$select_cpm]]$csd_data))
    })

    output$original_data <- renderUI({
        ## To see if I disable original data
        if(length(all_cpm_out$output) > 0){
            tags$div(class="frame max_height",
                tags$h3("3. The original data"),
                plotOutput("csd"),
                tags$div(class = "download_button",
                    actionButton("modify_data", "Modify data")
                )
            )
        }
    })

    ## Tabular data
    genotype_freq_df <- reactive({
        create_tabular_data(all_cpm_out$output[[
            input$select_cpm]]
            ,  input$data2table )
        }
    )

    output$tabular_data <- renderUI({
        if(length(all_cpm_out$output) > 0){
            tags$div(class="frame max_height",
                tags$h3("4. Tabular data"),
                radioButtons(inputId = "data2table", 
                        label = "", 
                        inline = TRUE,
                        choiceNames =  c("Transition rates",
                                         "Genotype transitions counts",
                                         "Genotype frequencies",
                                         "Transition probabilities",
                                         "Lambdas/probabilities",
                                         "Time-discretized transition matrix"),
                        choiceValues =  c("trans_rate_mat",
                                          "genotype_transitions",
                                          "freqs",
                                          "trans_mat",
                                          "lambdas",
                                          "td_trans_mat"),
                        selected =  "freqs"
                        ),
                tags$div( 
                    DT::DTOutput("cpm_freqs")
                )
            )
        }
    })

    ## Download button
    output$download_cpm <- downloadHandler(
        filename = function() sprintf("%s_cpm.RDS", input$select_cpm),
        content = function(file) {
            saveRDS(all_cpm_out$output[[input$select_cpm]], file)
        }
    )

    ## We only want the "Source" option enable if we have the data to show it
    observeEvent(input$select_cpm, {
        selected <- ifelse(is.null(input$select_cpm), last_visited_cpm, input$select_cpm)
        if(is.null(all_cpm_out$output[[selected]]$Source_genotype_transitions)){
            updateCheckboxGroupInput(session, "cpm2show", 
                selected = setdiff(c(input$cpm2show), "Source"))
            shinyjs::disable(selector = "#cpm2show input[value='Source']")
        }else{
            shinyjs::enable(selector = "#cpm2show input[value='Source']")
        }
        if(is.na(all_cpm_out$output[[selected]]$MCCBN_model)){
            updateCheckboxGroupInput(session, "cpm2show", 
                selected = setdiff(c(input$cpm2show), "Source"))
            shinyjs::disable(selector = "#cpm2show input[value='MCCBN']")
        }else{
            shinyjs::enable(selector = "#cpm2show input[value='MCCBN']")
        }
    })

    ## Upload button
    observeEvent(input$output_cpms, {
        cpm_out <- readRDS(input$output_cpms$datapath)
        ## FIXME: f_graph to trans_rate_mat?
        # cpm_out$MHN_f_graph <- cpm_out$MHN_trans_rate_mat
        if(is.null(cpm_out$name)) cpm_out$name <- "User_Data"
        all_cpm_out$output[[cpm_out$name]] <- cpm_out
        updateRadioButtons(session, "select_cpm", selected = cpm_out$name)
        ## To see if I disable original data
        shinyjs::disable(selector = "#variable input[value='cyl']")
    })
}
