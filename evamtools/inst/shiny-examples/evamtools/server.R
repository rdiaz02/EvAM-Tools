#'  I followed this link to structure the shiny app whithin the package
#'  https://deanattali.com/2015/04/21/r-package-shiny-app/

dataModal <- function(error_message, type="Error: ") {
  modalDialog(
    easyClose = TRUE,
    title = tags$h3(type),
    tags$div(
      error_message
    )
  )
  }



server <- function(input, output, session) {
  examples_csd$csd <- examples_csd$csd[1:5]
  all_csd_data <- evamtools:::standarize_all_datasets(examples_csd)
  min_genes <- SHINY_DEFAULTS$min_genes
  max_genes <- SHINY_DEFAULTS$max_genes
  default_csd_samples <- SHINY_DEFAULTS$csd_samples
  default_cpm_samples <- SHINY_DEFAULTS$cpm_samples
  default_dag_model <- SHINY_DEFAULTS$dag_model
  keep_dataset_name <- FALSE

  last_visited_pages <- list(csd = "User", dag = "User", matrix = "User")

  ## FIXME: No default "user"
  ## last_visited_cpm <- "user" ## And why "user" and not "User"?
  last_visited_cpm <- ""
  

  ## random_name <- paste(sample(c(LETTERS, letters), size = 8, replace = TRUE),
  ##                    collapse = "")
  ## last_visited_pages <- list(csd = random_name,
  ##                            dag = random_name,
  ##                            matrix = random_name)
  ## last_visited_cpm <- random_name

  datasets <- reactiveValues(
    all_csd = all_csd_data
  )

   data <- reactiveValues(
    csd_counts = SHINY_DEFAULTS$template_data$csd_counts
    , data = SHINY_DEFAULTS$template_data$data
    , dag = SHINY_DEFAULTS$template_data$dag
    , dag_parent_set = SHINY_DEFAULTS$template_data$dag_parent_set
    , lambdas = SHINY_DEFAULTS$template_data$lambdas
    , thetas = SHINY_DEFAULTS$template_data$thetas
    , n_genes = SHINY_DEFAULTS$ngenes
    , gene_names = LETTERS[1: max_genes]
  )

  display_freqs <- reactive({
    evamtools:::get_display_freqs(data$csd_counts, input$gene_number, data$gene_names)
  })

  ## Upload data
  observeEvent(input$csd, {
    
    if(grepl(".csv", input$csd$datapath)){
      dataset_name <- strsplit(strsplit(input$csd$name, ".csv")[[1]], "_")[[1]][[1]]
      tmp_data <- list()
      tmp_data$data <- read.csv(input$csd$datapath)

      tryCatch({
        datasets$all_csd[["csd"]][[dataset_name]] <- evamtools:::standarize_dataset(tmp_data)
        datasets$all_csd[["csd"]][[dataset_name]]$name <- dataset_name
        last_visited_pages["csd"] <<- dataset_name
        updateRadioButtons(session, "input2build", selected = "csd")
        updateRadioButtons(session, "select_csd", selected = dataset_name)
      }, error = function(e){
        showModal(dataModal(e[[1]]))
      })
    } else if(grepl(".rds", input$csd$datapath, ignore.case = TRUE)){
        tmp_data <- readRDS(input$csd$datapath)
        tryCatch({
          new_data <- evamtools:::standarize_dataset(tmp_data)
          datasets$all_csd[[tmp_data$type]][[new_data$name]] <- new_data
          last_visited_pages[tmp_data$type] <<- tmp_data$name
          # datasets$all_csd[[tmp_data$type]][[tmp_data$name]] <- tmp_data
          updateRadioButtons(session, "input2build", selected = tmp_data$type)
          updateRadioButtons(session, "select_csd", selected = tmp_data$name)
        }, error = function(e){
          showModal(dataModal(e[[1]]))
        })
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

  # ## Define dataset name
  output$dataset_name <- renderUI({
    tags$div(class = "inlin2",
             textInput(inputId = "dataset_name",
                       "Give your dataset a name",
                       value = input$select_csd)
    )
  })

  # listen2dataset_name <- reactive({
  #   list(input$dataset_name, input$save_csd_data)
  # })

  ## Saving dataset
  observeEvent(input$save_csd_data, {
    tryCatch({
      ## 1 save dataset to list after user data
      # if (!(input$dataset_name %in% names(datasets$all_csd[[input$input2build]]))){
          datasets$all_csd[[input$input2build]][[input$dataset_name]]$name <-
              input$dataset_name

          if (nrow(data$csd_counts) > 0) {
              ## FIXME: comment_out_freqs2csd
              ## datasets$all_csd[[input$input2build]][[input$dataset_name]]$data <-
              ##     evamtools:::freqs2csd(data$csd_counts,
              ##                           data$gene_names[1:input$gene_number])
            datasets$all_csd[[input$input2build]][[input$dataset_name]]$data <-
                evamtools:::genotypeCounts_to_data(data$csd_counts,
                                                   e = 0)
        }

        tmp_data <- list(
          data = data$data
          , dag = data$dag
          , gene_names = data$gene_names
          , dag_parent_set = data$dag_parent_set
          , lambdas = data$lambdas
          , thetas = data$thetas
          , trm = data$trm
          , n_genes = input$gene_number
          , name = input$dataset_name)
        datasets$all_csd[[input$input2build]][[input$dataset_name]] <- tmp_data

        # tmp_data_2 <- datasets$all_csd[[input$input2build]]
        # datasets$all_csd[[input$input2build]] <-
        #   c(tmp_data_2["User"]
        #   , tmp_data_2[input$dataset_name]
        #   , tmp_data_2[which(!(names(tmp_data_2) %in% c("User",
        #     input$dataset_name, names(all_csd_data))))]
        #   , tmp_data_2[which(names(datasets$all_csd[[input$input2build]]) %in% names(all_csd_data))]
        # )

        # ## 2 restore default values
        # try({
        #     datasets$all_csd[[input$input2build]][[input$select_csd]] <-
        #         all_csd_data[[input$input2build]][[input$select_csd]]
        # })

        ## 3 update selected entry
        updateRadioButtons(session, "select_csd", selected = input$dataset_name)

        # shinyjs::disable("save_csd_data")
      # }
      }, error=function(e){
        showModal(dataModal(e[[1]]))
      })
  })

  # observeEvent(listen2dataset_name(), {
  # observeEvent(input$dataset_name, {
    
  #   tryCatch({
  #     dataset_name <- ifelse(is.null(input$dataset_name), "", input$dataset_name)
  #     if(dataset_name != "" &
  #       !(dataset_name %in% names(datasets$all_csd[[input$input2build]]))) {
  #       shinyjs::enable("save_csd_data")
  #     }else if(dataset_name %in% names(datasets$all_csd[[input$input2build]])){
  #       shinyjs::disable("save_csd_data")
  #     }
  #   }, error = function(e){
  #     showModal(dataModal(e[[1]]))
  #   })
  # })

  ## Download csd button
  output$download_csd <- downloadHandler(
    filename = function() sprintf("%s_csd.csv", input$select_csd),
    content = function(file) {
      tmp_data <- datasets$all_csd[[input$input2build]][[input$select_csd]]$data
      write.csv(tmp_data, file, row.names = FALSE)
      # tmp_data$type <- input$input2build
      # saveRDS(tmp_data, file)
    }
  )

  observeEvent(input$select_csd, {
    tryCatch({
      last_visited_pages[[input$input2build]] <<- input$select_csd
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  # ## Display List of availabe CSD
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
  # observe({
    tryCatch({
      input$select_csd
      input$input2build
      ## Cleaning stuf
      selected <- last_visited_pages[[input$input2build]]
      tmp_data <- datasets$all_csd[[input$input2build]][[selected]]
      data$gene_names <- tmp_data$gene_names
      data$data <- tmp_data$data

      shinyjs::disable("analysis")
      shinyjs::hide("all_advanced_options")
      if (!is.null(data$data)) {
        data$csd_counts <- evamtools:::get_csd(data$data)
        shinyjs::enable("analysis")
      } else {
        data$csd_counts <- SHINY_DEFAULTS$template_data$csd_counts
      }

      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$dag_parent_set
      data$lambdas <- tmp_data$lambdas
      data$thetas <- tmp_data$thetas
      data$name <- tmp_data$name
      data$n_genes <- tmp_data$n_genes
      
      if (input$input2build == "dag") {
          # to_keep <- length(which(colSums(data$dag) > 0 |
          #                         rowSums(data$dag) > 0)) - 1
          ## to_keep <- sum(colSums(data$dag) > 0)
          ## if (to_keep < 2) {
          ##     showModal(dataModal("You have only defined one relationship"))
          ##     updateRadioButtons(session, "dag_model", selected = "HESBCN")
          ## }
          # ngenedats <- to_keep
          number_of_parents <- colSums(data$dag)
          to_keep <- sum(number_of_parents > 0)
          n_genes <- ifelse(to_keep < 1, SHINY_DEFAULTS$ngenes, to_keep)
          # if (input$dag_model == "OT" && any(number_of_parents > 1)) {
          #     showModal(dataModal("This DAG cannot be transformed into a tree"))
          #     # updateRadioButtons(session, "dag_model", selected = "HESBCN")
          # }
          updateRadioButtons(session, "dag_model", selected = "HESBCN")
      } else if (input$input2build == "matrix") {
          # n_genes <- length(which(colSums(abs(data$thetas)) > 0
          #                         | rowSums(abs(data$thetas)) > 0))
          # n_genes <- ifelse(n_genes <= 0, 3, n_genes)
          n_genes <- data$n_genes
          if(is.null(n_genes)){
            n_genes <- SHINY_DEFAULTS$ngenes
          }
      } else if (input$input2build == "csd" && !is.null(data$data)) {
          n_genes <- ncol(data$data)
      } else if (input$input2build == "csd" && is.null(data$data)) {
          n_genes <- SHINY_DEFAULTS$ngenes
      }

      updateNumericInput(session, "gene_number", value = n_genes)
      updateNumericInput(session, "genotype_freq", value = NA)
      updateCheckboxGroupInput(session, "genotype", label = "Mutations",
                               choices = lapply(1:n_genes, function(i)data$gene_names[i]),
                               selected = NULL)
    }, error = function(e) {
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$gene_number, {
    datasets$all_csd[[input$input2build]][[input$select_csd]]$n_genes <- input$gene_number
  })

  observeEvent(input$change_gene_names, {
    showModal(modalDialog(
      easyClose = TRUE,
      title = tags$h3("Change gene names"),
      tags$div(class = "inlin2",
        textInput(inputId = "new_gene_names", "Gene names",
          value = paste(data$gene_names[1:input$gene_number],
          collapse = ", ")
          ),
        tags$h3(HTML("<br/>")),
        tags$h4("Separate you gene names with a ','. ",
                "Use only alphanumeric characters ",
                "(of course, do not use comma as part of a gene name), ",
                "and do not start ",
                "a gene name with a number; ",
                "keep gene names short (for figures)." 
                ),
      tags$div(class = "download_button",
        actionButton("action_gene_names", "Change genes names"),
        )
      )
    ))
  })

  ## Updating gene names
  observeEvent(input$action_gene_names,{
    tryCatch({
      new_gene_names <- unique(strsplit(gsub(" ", "", input$new_gene_names), ",")[[1]])
      data$gene_names <- c(
        new_gene_names
        , LETTERS[(length(new_gene_names) + 1):max_genes]
      )
      ## Rename stuff
      new_data <- evamtools:::standarize_dataset(data)
      data$data <- new_data$data
      data$dag <- new_data$dag
      data$dag_parent_set <- new_data$dag_parent_set
      data$thetas <- new_data$thetas
      data$lambdas <- new_data$lambdas
      data$csd_counts <- new_data$csd_counts

      datasets$all_csd[[input$input2build]][[input$select_csd]] <- new_data
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$display_help, {
    showModal(modalDialog(
      easyClose = TRUE,
      title = tags$h3("Changing genotype's counts"),
      tags$div(
        tags$p("1. Double click in a Counts cell to edit it"),
        tags$p("2. Press Tab to move to the next row"),
        tags$p("3. Use Ctrl + Enter to save changes"),
        tags$p("4. Set a frequency to 0 to remove a genotype"),
        tags$p("5. Type in the Search bar to filter genotypes")
        )
        )
      )
  })

  # ## Advanced option for running evamtools
  observeEvent(input$advanced_options, {
    shinyjs::toggle("all_advanced_options")
  })

  # ## Define number of genes
  output$gene_number <- renderUI({
      val <- ifelse(is.null(data$n_genes), 3, data$n_genes)
      tags$div(class="inlin",
               tags$h3(HTML("<br/>")),
               sliderInput("gene_number", "Number of genes",
                           value = val, max = max_genes, min = min_genes,
                           step = 1)
    )
  })

  ## Define new genotype
  observeEvent(input$dag_model, {
    number_of_parents <- colSums(data$dag)
    if(input$dag_model == "OT" && any(number_of_parents > 1)){
      showModal(dataModal("This DAG cannot be transformed into a tree"))
      updateRadioButtons(session, "dag_model", selected = "HESBCN")
    }else{
      default_dag_model <<- input$dag_model
    }
  })

  output$define_genotype <- renderUI({
    n_genes <- ifelse(is.null(input$gene_number), 3, input$gene_number)
    options <- data$gene_names[1:n_genes]
    if(input$input2build == "csd"){
      tags$div(
               tags$h3("2. Add genotypes"),
               tags$h5("WT is added by not clicking on any mutations; ",
                       "but the WT genotype should not be the first one added ",
                       "(or you'll get an innocuous error message)."),
               tags$h5("Any gene without mutations is excluded from the data, ",
                       "regardless of the setting for number of genes. "),
               tags$h5("If any gene is always observed mutated ",
                       "(i.e., has a constant value of 1 for all observations), ",
                       "one observation with no genes mutated is added ",
                       "to the sample."),
               tags$h5("Genes that have identical patterns",
                       "(i.e., that lead to identical columns in the data matrix), ",
                       "are fused into",
                       "a single gene."),
               tags$h5("Genotypes are always shown with gene names ",
                       "sorted alphanumerically. "),
               tags$h5(" "),
            tags$div(class = "inline",
              checkboxGroupInput(inputId = "genotype",
                label = "Mutations",
                choices = options)
            ),
            tags$div(id="fr",
                     numericInput(label = "Counts", value = NA, min = 0,
                                  inputId = "genotype_freq", width = NA),
                     actionButton("add_genotype", "Add genotype")
                     ),
            tags$h5(HTML("<br>")),
            actionButton("clear_genotype",
                         "Clear all genotype data")
        )
    } else if (input$input2build == "dag"){
      tags$div(
        tags$div(class = "flex",
          tags$h3("2. Define a Directed Acyclic Graph (DAG)"),
          actionButton("how2build_dag", "Help")
        ),
        if(!is.null(data$lambdas)){
          tags$div(
            tags$h4("Type of model"),
            tags$div(class = "inline",
              radioButtons(inputId = "dag_model",
              label = "Model: ",
              inline = TRUE,
              ## choices =  c("OT", "OncoBN", "HESBCN"),
              choiceNames = list("OT", "OncoBN", "CBN/HESBCN"),
              choiceValues = list("OT", "OncoBN", "HESBCN"),
              selected = default_dag_model)
            ),
            tags$h4("New Edge"),
                   tags$h5(HTML("<p></p>")),
            tags$div(class = "inline",
              radioButtons(inputId = "dag_from",
              label = "From (parent node)",
              inline = TRUE,
              choices =  c("Root", options))
            ),
            tags$div(class = "inline",
              radioButtons(inputId = "dag_to",
              label = " To (child node)",
              inline = TRUE,
              choices =  options)
              ),
            tags$h5(HTML("<p></p>")),
            actionButton("add_edge", "Add edge"),
            actionButton("remove_edge", "Remove edge"),
            actionButton("clear_dag", "Clear DAG"),
            tags$h3(HTML("<br/>DAG table")),
            DT::DTOutput("dag_table"),
            tags$h3(HTML("<br/>")),
            numericInput("dag_samples", HTML("Total genotypes to sample"),
                         value = default_csd_samples, min = 100, max = 10000,
                         step = 100, width = "50%"),
            numericInput("dag_noise", HTML("Noise"),
                         value = 0.01, min = 0, max = 1,
                         step = 0.05, width = "50%"),
            numericInput("dag_epos", HTML("Epos/Epsilon"),
                         value = 0.01, min = 0, max = 1,
                         step = 0.05, width = "50%"),
            tags$h3(HTML("<br/>")),
            actionButton("resample_dag", "Sample from DAG")
            )
          }
        )
      } else if (input$input2build == "matrix"){
        tags$div(
          tags$div(class = "flex",
                   ## tags$h3("2. Define input with a Matrix"),
                   tags$h3("2. Define MHN's log-Theta",
                           HTML("matrix (log-&Theta;):")),
            actionButton("how2build_matrix", "Help")
          ),
          if(!is.null(data$thetas)){
            tags$div(
                     tags$h3("Entries are ",
                             "lower case thetas, ",
                             HTML("&theta;s, range &plusmn; &infin;"),),
              DT::DTOutput("thetas_table"),
              tags$h3(HTML("<br/>")),
              numericInput("mhn_samples", "Total genotypes to sample",
                           value = default_csd_samples, min = 100, max = 10000,
                           step = 100, width = "50%"),
              numericInput("mhn_noise", HTML("Noise"),
                         value = 0.01, min = 0, max = 1,
                         step = 0.05, width = "50%"),
              tags$h3(HTML("<br/>")),
              actionButton("resample_mhn", "Sample from MHN"),
              actionButton("clear_mhn", "Clear matrix data")
            )
          }
        )
      }
  })

  output$change_counts <- renderUI({
    if(input$input2build == "csd"){
      tags$div(class = "frame",
        tags$div(class = "flex",
                 tags$h3("3. Change genotype's counts"),
                 actionButton("display_help", "Help"),
                 tags$h3(HTML("<br/>")),
        ),
        tags$div(id = "csd_table",
          DT::DTOutput("csd_counts")
        )
      )
    }
  })


  # ## DAG builder
  # ## Controling dag builder
  dag_data <- reactive({
    input$dag_model
    input$dag_table_cell_edit
    all_gene_names <- c("Root", data$gene_names)
    edges <- which(data$dag == 1, arr.ind = TRUE)
    tmp_dag_parent_set <- data$dag_parent_set
    x <- length(tmp_dag_parent_set)
    ## I have to this weird thing because using data$gene_names does not work for some unkown reason
    names(tmp_dag_parent_set) <- all_gene_names[seq(2, x + 1)]

    # lambdas <- data$lambdas
    # if(default_dag_model %in% c("OT")) {
    #   lambdas <- rep(0.9, 10)
    #   names(lambdas) <- data$gene_names
    # }
    dag_data <- data.frame(From = all_gene_names[edges[, "row"]]
      , To = all_gene_names[edges[, "col"]]
      , Relation = tmp_dag_parent_set[edges[, "col"] - 1]
      , Lambdas = data$lambdas[edges[, "col"] - 1])
    
    if((default_dag_model %in% c("OT", "OncoBN"))
      & (any(dag_data$Lambdas < 0) | any(dag_data$Lambdas > 1))){
        showModal(dataModal("Thetas/probabilities should be between 0 and 1"))
        updateRadioButtons(session, "dag_model", selected = "HESBCN")
    }

    if(default_dag_model %in% c("OT")) {
      colnames(dag_data) <- c("From", "To", "Relation", "Weight")
      dag_data$Relation <- NULL
      # data$lambdas <- data$lambdas - 0.1
      # dag_data$Lambdas <- data$lambdas[edges[, "col"] - 1]
     
    } else if(default_dag_model %in% c("OncoBN")){
      if(length(unique(dag_data$Relation))>2){
        showModal(dataModal("OncoBN model must only include one type of relationship"))
        updateRadioButtons(session, "dag_model", selected = "HESBCN")
      }
      colnames(dag_data) <- c("From", "To", "Relation", "Thetas")
     
    }

    return(dag_data)
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
    tryCatch({

      tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene,
                                          operation = "add",
                                          parent_set = data$dag_parent_set,
                                          dag_model = default_dag_model)
      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$parent_set
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
      if(default_dag_model != "OT"){
        shinyjs::click("resample_dag")
      }
    },error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Remove edge
  observeEvent(input$remove_edge, {
    from_gene <- input$dag_from
    to_gene <- input$dag_to
    tryCatch({
        tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene,
                                           operation = "remove",
                                           parent_set = data$dag_parent_set,
                                           dag_model = default_dag_model)
      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$parent_set
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
      if(sum(data$dag) == 0) {
        data$csd_counts <- datasets$all_csd[[input$input2build]][[input$select_csd]]$csd_counts
      } else if (default_dag_model != "OT"){
        shinyjs::click("resample_dag")
      }
    },error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Clear DAG
  observeEvent(input$clear_dag, {
    tryCatch({
      tmp_data <- evamtools:::modify_dag(data$dag, NULL, NULL, operation = "clear")
      tmp_dag <- tmp_data$dag
      colnames(tmp_dag) <- rownames(tmp_dag) <- c("WT", data$gene_names)
      tmp_dag["WT", data$gene_names[1]] <- 1
      data$dag <- tmp_dag
      data$csd_counts <- SHINY_DEFAULTS$template_data$csd_counts
      data$data <- SHINY_DEFAULTS$template_data$data
      data$dag_parent_set <- tmp_data$dag_parent_set
      data$lambdas <- SHINY_DEFAULTS$template_data$lambdas
      names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
      shinyjs::disable("analysis")
    }, error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$dag_table_cell_edit, {
    tryCatch({
      names(data$dag_parent_set) <- data$gene_names[1:length(data$dag_parent_set)]
      names(data$lambdas) <- data$gene_names[1:length(data$dag_parent_set)]
      info <- input$dag_table_cell_edit
      tmp_data <- evamtools:::modify_lambdas_and_parent_set_from_table(dag_data(), info, data$lambdas
        , data$dag, data$dag_parent_set
        , dag_model= default_dag_model)
      data$lambdas <- tmp_data$lambdas
      data$dag_parent_set <- tmp_data$parent_set
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
      shinyjs::click("resample_dag")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Building trm from dag
  observeEvent(input$resample_dag, {
    tryCatch({
      gene_names <- setdiff(unique(c(dag_data()$From, dag_data()$To)), "Root")
      number_of_genes <- length(gene_names)
      tmp_dag_data <- evamtools:::get_dag_data(dag_data()
        , data$dag_parent_set[gene_names]
        , noise = input$dag_noise
        , N = input$dag_samples
        , dag_model = default_dag_model
        , epos = input$dag_epos)
      data$csd_counts <- tmp_dag_data$csd_counts
      data$data <- tmp_dag_data$data

      datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
      datasets$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
      datasets$all_csd[[input$input2build]][[input$select_csd]]$lambdas <- data$lambdas
      datasets$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set <- data$dag_parent_set

      shinyjs::enable("analysis")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Help for DAG building
  observeEvent(input$how2build_dag, {
    showModal(modalDialog(
      easyClose = TRUE,
      title = tags$h3("How to build a DAG and generate a sample"),
      tags$div(
               tags$p("1. Select 'From' (parent node) and 'To' (child node) ",
                      HTML("and hit 'Add edge' or 'Remove edge'.<ul>")),
               tags$li(HTML("An edge wont be allowed if: <ul>")),
               tags$li("it is already present;"),
               tags$li("it introduces cycles."),
               tags$p(HTML("</ul>")),
               tags$li(HTML("To remove edge you can also "),
                       "set the lambda of the relationship to 0."),
               tags$li(HTML("Removing edges might restructure the DAG."),
                       "If a node has no parent, " ,
                       "it will be assigned as descendant of WT."),
               tags$p(HTML("</ul>")),
               tags$p(HTML("2. To <strong>change the value of a lambda</strong> "),
                       "click on the cell, ",
                      "edit the cell's content and press Ctrl+Enter."),
               tags$p(HTML("3. Set the value of <strong>Relation</strong> "),
                       "to one of 'Single' (single parent), ",
                      "AND, OR, XOR.",
                      "Edit the cell's content and press Ctrl+Enter. ",
                      "All incoming edges to a node must have the same ",
                      "Relation (the program will force this)."),
               tags$p(HTML("4. Modify, if you want, the <strong>size of the sample</strong> "),
                      "('Total genotypes to sample') and ",
                      "click on 'Sample from DAG' to generate a sample. "),
               tags$p("5. Possible random noise when sampling is controlled under 'Advanced options'.")
           )
    )
    )
  })

  ## MHN
  ## Controling thetas
  output$thetas_table <- DT::renderDT(data$thetas[1:input$gene_number,
                                                  1:input$gene_number]
    , selection = 'none', server = TRUE
    , editable = list(target = "all", disable = list(columns = c(0)))
    , rownames = TRUE,
    options = list(
      searching = FALSE, columnDefs = list(list(className = 'dt-center',
      targets = "_all")), info = FALSE, paginate= FALSE),
  )

  observeEvent(input$thetas_table_cell_edit, {
    tryCatch({
      info <-input$thetas_table_cell_edit
      data$thetas[1:input$gene_number, 1:input$gene_number] <-
        DT::editData(data$thetas[1:input$gene_number, 1:input$gene_number], info, "thetas")

      datasets$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
      datasets$all_csd[[input$input2build]][[input$select_csd]]$n_genes <- input$gene_number
      ## Resample based on changes
      shinyjs::click("resample_mhn")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$resample_mhn, {
    tryCatch({
      mhn_data <- evamtools:::get_mhn_data(data$thetas[1:input$gene_number
                                                  , 1:input$gene_number]
                                          , noise = input$mhn_noise 
                                          , N = input$mhn_samples)
      data$csd_counts <- mhn_data$csd_counts
      data$data <- mhn_data$data
      datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- mhn_data$data
      shinyjs::enable("analysis")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$clear_mhn, {
    tryCatch({
      tmp_thetas <- SHINY_DEFAULTS$template_data$thetas
      colnames(tmp_thetas) <- rownames(tmp_thetas) <- data$gene_names
      data$thetas <- tmp_thetas
      data$dag <- NULL
      data$csd_counts <- SHINY_DEFAULTS$template_data$csd_counts
      data$data <- SHINY_DEFAULTS$template_data$data
      data$dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
      data$lambdas <- SHINY_DEFAULTS$template_data$lambdas
      names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
      shinyjs::disable("analysis")
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
    }, error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$how2build_matrix, {
    showModal(modalDialog(
      easyClose = TRUE,
      title = tags$h3(HTML("How to input &theta;s and generate a sample")),
      tags$div(
        ## tags$p("Positive theta: gene i in row makes gene j in column more likely. A negative means the opposite."),
        ## tags$p("Diagonal theta: likelihood of that event i to be the first one (positive values likely, negative unlikely."),
               ## tags$p("Once the thetas are defined hit the",
               ##        " 'Sample from MHN' to generate a sample."),
               tags$p("0. Select the number of genes with the slider, above."
                     ## , " (Even if you set the diagonal and all off-diagonal",
                     ## " entries as 0, that gene is still part of the data set",
                     ## HTML(", with a value &Theta; = 1 for all its contributing terms).")
                     ),
               ## tags$p("To make a sample we take into account multiplicative effects of all thetas"),
               ## tags$h3("How to modify the table"),
               tags$p(HTML("1. &Theta;<sub>i,j</sub> ",
                           "(i.e., <em>e<sup>&theta;<sub>i,j</sub></sup></em>) ",
                           "is the multiplicative ",
                           "effect of gene in column <em>j</em> on ",
                           "gene in row <em>i</em>. ",
                           "&Theta;<sub>i,i</sub> is the baseline hazard rate ",
                           "of event <em>i</em>. "## ,
                           ## "See the ",
                           ## "figure on the right."
                           )),
               tags$p("2. Double click in a cell to edit it."),
               tags$p("3. Press Tab to move to the next row."),
               tags$p("4. Use Ctrl + Enter to save changes. ",
                      HTML("You <strong>must</strong> save the changes.")),
               tags$p("5. Modify, if you want, the size of the sample ",
                      "('Total genotypes to sample') and ",
                      "click on 'Sample from MHN' to generate a sample. ",
                      "The sample is also updated as soon as you save an entry."),
               tags$p("6. Possible random noise when sampling is controlled under 'Advanced options'."),
               tags$p(HTML("Make sure <b>the &theta;s have been updated</b> "),
                      "by checking the figure of the matrix on the right.")
        ## tags$p("4. Set a frequency to 0 to remove a genotype"),
        ## tags$p("5. Type in the Search bar to filter genotypes")
        )
      )
    )
  })

  ## Working with raw CSD
  observeEvent(input$genotype, {
    tryCatch({
      genotype <- paste(input$genotype, collapse = ", ")
      genot_count <- data$csd_counts[, 2][data$csd_counts[, 1] == genotype]
      updateNumericInput(session, "genotype_freq", value = genot_count)

    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  }, ignoreNULL = FALSE)

  observeEvent(input$clear_genotype, {
    tryCatch({
      tmp_dag <- SHINY_DEFAULTS$template_data$dag
      colnames(tmp_dag) <- rownames(tmp_dag) <- c("WT", data$gene_names)
      tmp_dag["WT", data$gene_names[1]] <- 1
      data$dag <- NULL
      data$csd_counts <- SHINY_DEFAULTS$template_data$csd_counts
      data$data <- SHINY_DEFAULTS$template_data$data
      data$dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
      data$lambdas <- SHINY_DEFAULTS$template_data$lambdas
      names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
      shinyjs::disable("analysis")
      datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::standarize_dataset(data)
    }, error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$add_genotype, {
      tryCatch({
      genotype <- paste(sort(input$genotype), collapse = ", ")
      genotype <- ifelse(genotype == "", "WT", genotype)
      genot_count <- ifelse(is.na(input$genotype_freq), -1,
                            input$genotype_freq)

      if (genot_count >= 0) {
        data$csd_counts[genotype, ] <- c(genotype, genot_count)
        rownames(data$csd_counts) <- data$csd_counts$Genotype
        data$csd_counts[, 2] <- as.numeric(data$csd_counts[, 2])
        ## Filtering out non-positive counts
        data$csd_counts <- data$csd_counts[data$csd_counts[, 2] > 0,]
        ## FIXME: comment_out_freqs2csd
        ## data$data <- datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
        ##              evamtools:::freqs2csd(data$csd_counts,
        ##                                    data$gene_names[1:input$gene_number])
        data$data <-
            datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
                evamtools:::genotypeCounts_to_data(data$csd_counts, e = 0)
        
        
        shinyjs::enable("analysis")
      }
      updateNumericInput(session, "genotype_freq", value = NA)
      updateCheckboxGroupInput(session, "genotype", label = "Mutations",
                               choices = lapply(1:input$gene_number,
                                                function(i)data$gene_names[i]),
                               selected = NULL)
    }, error = function(e) {
      showModal(dataModal(e[[1]]))
    })
  })

  ## Genotypes table
  output$csd_counts <- DT::renderDT(display_freqs(), selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
    , rownames = FALSE,
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
  )

  observeEvent(input$csd_counts_cell_edit, {
    tryCatch({
      info <- input$csd_counts_cell_edit
      info[ , "col"] <- 2
      data$csd_counts <- DT::editData(data$csd_counts, info, "csd_counts")
      ## Filtering out non-positive counts
      data$csd_counts <- data$csd_counts[data$csd_counts[,2] > 0,]

      ## FIXME: comment_out_freqs2csd
      ## data$data <-
      ##     datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
      ##         evamtools:::freqs2csd(data$csd_counts,
      ##                               data$gene_names[1:input$gene_number])
      data$data <-
          datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
              evamtools:::genotypeCounts_to_data(data$csd_counts, e = 0)

    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Plot histogram of genotypes
  output$plot <- renderPlot({
    tryCatch({
      evamtools:::plot_genotype_counts(display_freqs())
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Plot dag of dataset
  output$dag_plot <- renderPlot({
    data2plot <- NULL
    edges <- NULL
    if(input$input2build %in% c("csd", "dag")
      && sum(data$dag)>0
      && !is.null(input$gene_number)
      ){
      data2plot <- igraph::graph_from_adjacency_matrix(data$dag)
      data2plot <- igraph::decompose(data2plot)[[1]]
      edges <- igraph::as_data_frame(data2plot)
      colnames(edges) <- c("From", "To")
      if(!is.null(data$dag_parent_set)) edges$Relation <- data$dag_parent_set[edges$To]
    }else if(input$input2build %in% c("matrix") 
      && !is.null(data$thetas)
      && length(data$thetas[1:input$gene_number, 1:input$gene_number])>0
      ){
        data2plot <- data$thetas[1:input$gene_number, 1:input$gene_number]
    }
    evamtools:::plot_method(data2plot, data$dag_parent_set, edges)
  })

  # ## Run CPMs
  observeEvent(input$analysis, {
    ## Calculate TRM for DAG and for matrices

    # source_trm <- NULL
    ## FIXME the following has to be moved to R/
    # if(input$input2build == "dag"){
    #   tmp_data <- list(edges = dag_data(), parent_set = data$dag_parent_set)
    #   source_trm <- evamtools:::cpm2tm(tmp_data)$weighted_fgraph
    # }else if(input$input2build == "matrix"){
    #   source_trm <- evamtools:::theta_to_trans_rate_3_SM(data$thetas[1:input$gene_number, 1:input$gene_number],
    #                       inner_transition = evamtools:::inner_transitionRate_3_1)
    # }

    tryCatch({
      if(input$gene_number >= 7){
        showModal(dataModal("Beware! You are running a dataset with 7 genes or more. This can take longer than usual and plots may be crowded. We recommend using top_paths options in the Results' tab.", type = "Warning: "))
      }
      shinyjs::disable("analysis")
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Running evamtools", value = 0)

      mhn_opts <- list()
      if(!is.na(input$MHN_lambda)) mhn_opts$lambda <- input$MHN_lambda
      
      ot_opts <- list()
      if(input$OT_with_error == "TRUE"){
        ot_opts$with_errors_dist_ot <- TRUE
      } else ot_opts$with_errors_dist_ot <- FALSE

      cbn_opts <- list(init_poset = input$CBN_init_poset)
      hesbcn_opts <- list(
        steps = input$HESBCN_steps,
        reg = input$HESBCN_reg
      ) 
      if(!is.na(input$HESBCN_seed)) hesbcn_opts$seed <- input$HESBCN_seed

      oncobn_opts <- list(
        model = input$OncoBN_model,
        algorithm = input$OncoBN_algorithm,
        k = input$OncoBN_k
      ) 
      if(!is.na(input$OncoBN_epsilon)) oncobn_opts$epsilon <- input$OncoBN_epsilon

      mccbn_opts <- list(
        model = input$MCCBN_model,
        L = input$MCCBN_L,
        sampling = input$MCCBN_sampling,
        max.iter = input$MCCBN_max_iter,
        update.step.size = input$MCCBN_update_step_size,
        tol = input$MCCBN_tol,
        max.lambda.val = input$MCCBN_max_lambda_val,
        T0 = input$MCCBN_T0,
        adap.rate = input$MCCBN_adapt_rate,
        max.iter.asa = input$MCCBN_max_iter_asa,
        neighborhood.dist = input$MCCBN_neighborhood_dist
      )
      if(!is.na(input$MCCBN_seed)) mccbn_opts$seed <- input$MCCBN_seed
      if(!is.na(input$MCCBN_acceptance_rate)) mccbn_opts$acceptance.rate <-  input$MCCBN_acceptance_rate
      if(!is.na(input$MCCBN_step_size)) mccbn_opts$step.size <-  input$MCCBN_acceptance_rate
      if(input$MCCBN_adaptive == "TRUE"){
        mccbn_opts$adaptive <- TRUE
      } else mccbn_opts$adaptive <- FALSE

      ## FIXME: comment_out_freqs2csd
      ## data2run <- evamtools:::freqs2csd(display_freqs(),
      ##                                   data$gene_names[1:input$gene_number])
      data2run <- evamtools:::genotypeCounts_to_data(display_freqs(),
                                                     e = 0)
        
      progress$inc(1/5, detail = "Setting up data")
      Sys.sleep(0.5)
      progress$inc(2/5, detail = "Running CPMs")

      methods <- SHINY_DEFAULTS$cpms2run
      if(!is.null(input$more_cpms)){
        methods <- unique(c(methods, input$more_cpms[input$more_cpms %in% SHINY_DEFAULTS$all_cpms]))
      }
      cpm_output <- evam(data2run, methods = methods
        , mhn_opts = mhn_opts
        , ot_opts = ot_opts
        , cbn_opts = cbn_opts
        , hesbcn_opts = hesbcn_opts
        , oncobn_opts = oncobn_opts
        , mccbn_opts = mccbn_opts)

      ## To see Source data in the results section
      # if(input$input2build != "csd"){
      #   ## FIXME: rename f_graph to trans_rate_mat: ???
      #   cpm_output$Source_trans_rate_mat <- source_trm
      #   cpm_output$Source_trans_mat <- evamtools:::rowScaleMatrix(source_trm)
      # }
      # if(input$input2build == "dag"){
      #   cpm_output$Source_model <- dag_data()
      #   cpm_output$Source_parent_set <- data$dag_parent_set[1:input$gene_number]
      # } else if(input$input2build == "matrix"){
      #   cpm_output$Source_theta <- data$thetas[1:input$gene_number
      #     , 1:input$gene_number]
      # }
      sampled_from_CPMs <- NULL
      do_sampling <- input$do_sampling == "TRUE"
      if (do_sampling) {
        n_samples <- input$sample_size
        if (is.null(n_samples) | !is.numeric(n_samples) | n_samples < 100) {
          n_samples <- SHINY_DEFAULTS$cpm_samples
        }
        progress$inc(3/5, detail = paste("Running ", n_samples, " samples"))
        if (input$do_genotype_transitions) {
            sout <- c("sampled")
        }
            
        sampled_from_CPMs <-
            sample_CPMs(cpm_output, n_samples , methods,
                        out = if (input$do_genotype_transitions) { 
                                     c("sampled_genotype_counts",
                                       "obs_genotype_transitions")
                                     } else {
                                         "sampled_genotype_counts"
                                     },
                      , obs_noise = input$sample_noise)
      }
      
      progress$inc(4/5, detail = "Post processing data")
      Sys.sleep(0.5)

      orig_data <- list(data = data2run, name = data$name
        , type = input$input2build, gene_names = data$gene_names
        , thetas = data$thetas, lambdas = data$lambdas
        , dag = data$dag, dag_parent_set = data$dag_parent_set)

      tabular_data <- evamtools:::create_tabular_data(c(cpm_output, sampled_from_CPMs))

      all_evam_output <- list("cpm_output" = c(cpm_output, sampled_from_CPMs)
                             , "orig_data" = orig_data
                            , "tabular_data" = tabular_data 
                              )

      ##CPM output name
      result_index <- length(grep(sprintf("^%s", input$select_csd),
                                  names(all_cpm_out)))
      result_name <- ifelse(result_index == 0
        , input$select_csd
        , sprintf("%s__%s", input$select_csd, result_index))

      all_cpm_out[[result_name]] <- all_evam_output
      last_visited_cpm <<- result_name
      updateRadioButtons(session, "select_cpm", selected = result_name)
      progress$inc(5/5, detail = "You can see your result by going to the Results tab")
      Sys.sleep(1)
      shinyjs::enable("analysis")

      ## Maybe I do no want this
      updateTabsetPanel(session, "navbar", selected = "result_viewer")
      updateRadioButtons(session, "select_cpm", selected = result_name)

      # selected <- ifelse(is.null(input$select_cpm), last_visited_cpm, input$select_cpm)
      # if(is.null(all_cpm_out[[selected]]$Source_genotype_transitions)){
      #     updateCheckboxGroupInput(session, "cpm2show",
      #         selected = setdiff(c(input$cpm2show), "Source"))
      #     shinyjs::disable(selector = "#cpm2show input[value='Source']")
      # }else{
      #     shinyjs::enable(selector = "#cpm2show input[value='Source']")
      # }

     }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## FIXME: No default "user"
  ## cpm_out <- sample_evam_output
  ## all_cpm_out <- reactiveValues(user = cpm_out)
  all_cpm_out <- reactiveValues()

  output$sims <- renderUI({
    if ((length(names(all_cpm_out)) > 0) && (!is.null(input$select_cpm))) {
      tmp_data <- all_cpm_out[[input$select_cpm]]$cpm_output
      # column_models2show <- floor(12 / length(input$cpm2show))

      number_of_columns <- floor(12 /
        ifelse(length(input$cpm2show) <=4, 4, length(input$cpm2show)))

      lapply(input$cpm2show, function(met){
          method_data <- evamtools:::process_data(tmp_data, met,
                                                  plot_type = "trans_mat")
        output[[sprintf("plot_sims_%s", met)]] <- renderPlot({
          pl <- evamtools:::plot_method(method_data$method_info
          , method_data$parent_set
          , method_data$edges
          , met)
        })
        return(
          # column(3,
          column(number_of_columns,
            plotOutput(sprintf("plot_sims_%s", met)))
        )
      })
    }
  })

  output$sims2 <- renderUI({
    if ((length(names(all_cpm_out)) > 0)  && (!is.null(input$select_cpm))) {
      tmp_data <- all_cpm_out[[input$select_cpm]]$cpm_output
      ## Enabling donwload button
      shinyjs::enable(selector = "#download_cpm")

      ## Main display
      # column_models2show <- floor(12 / length(input$cpm2show))
      selected_plot_type <- input$data2plot

      number_of_columns <- floor(12 /
        ifelse(length(input$cpm2show) <=4, 4, length(input$cpm2show)))
      if(!(is.null(selected_plot_type))){
        lapply(input$cpm2show, function(met){
          method_data <- evamtools:::process_data(tmp_data, met, plot_type = selected_plot_type)
          output[[sprintf("plot_sims2_%s", met)]] <- renderPlot({
            pl <- evamtools:::plot_genot_fg(method_data$data2plot,
                      # observations = tmp_data$analyzed_data, # We use it to define "Observed" and "Not Observed" genotypes
                      observations = tmp_data$original_data, # We use it to define "Observed" and "Not Observed" genotypes
                      # predicted_genotypes = method_data$predicted_genotype_freqs, # To compute node sizes if sampled_counts is NULL
                      sampled_counts = method_data$sampled_genotype_counts,
                      top_paths = input$freq2label,
                      label_type = input$label2plot,
                      plot_type = selected_plot_type)
                      # , freq2label = input$freq2label)
          })
          return(
            # column(3,
            column(number_of_columns,
              plotOutput(sprintf("plot_sims2_%s", met)))
          )
        })
    } else {
      ## Disabling donwload button
      shinyjs::disable(selector = "#download_cpm")

      return(tags$h3("There are not results to show yet. Go to the input tab, select a dataset and hit the 'Run evamtools!' button"))
    }
    }
  })

  # ## Go back to input to work again with the data
  observeEvent(input$modify_data, {
    tryCatch({
      if(length(all_cpm_out) > 0){
        tmp_data <- all_cpm_out[[input$select_cpm]]$orig_data
        dataset_name <- strsplit(input$select_cpm, "__")[[1]][[1]]
        dataset_type <- tmp_data$type
        last_visited_pages[[tmp_data$type]] <<- dataset_name
        tmp_data <- datasets$all_csd[[tmp_data$type]][[dataset_name]] <- evamtools:::standarize_dataset(tmp_data)

        data <- tmp_data
        data$csd_counts <- evamtools:::get_csd(tmp_data$data)
        data$n_genes <- ncol(data$data)

        updateNumericInput(session, "gene_number", value = data$n_genes)
        updateTabsetPanel(session, "navbar",
            selected = "csd_builder")
        updateRadioButtons(session, "input2build", selected = dataset_type)
        updateRadioButtons(session, "select_csd", selected = dataset_name)
      }
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  output$customize <- renderUI({
    tags$div(class = "frame",
      tags$h3("Customize the visualization"),
      tags$div(class = "inline",
        checkboxGroupInput(inputId = "cpm2show",
          label = "CPMs to show",
          choices = c("OT", "OncoBN", "CBN", "MHN", "HESBCN", "MCCBN"),
          selected = c("OT", "OncoBN", "CBN", "MHN", "HESBCN")),

      tags$div(class = "inline",
        radioButtons(inputId = "data2plot",
          label = "Data to show",
          choiceNames =  if (input$do_genotype_transitions) {
                             c("Transition probabilities", 
                               "Transition Rate Matrix",
                               "Observed genotype transitions")
                         } else {         
                             c("Transition probabilities", 
                               "Transition Rate Matrix")
                         },
          choiceValues = if (input$do_genotype_transitions) {
                             c("trans_mat", 
                               "trans_rate_mat",
                               "obs_genotype_transitions")
                         } else {
                             c("trans_mat", 
                               "trans_rate_mat")
                         },
          selected = "trans_mat"
          )
        ),
      tags$div(class = "inline",
        radioButtons(inputId = "label2plot",
          label = "Type of label",
          choiceNames =  c("Genotype", "Last gene mutated"),
          choiceValues = c("genotype", "acquisition"),
          selected = "genotype"
          )
        ),
      ),

      tags$p(HTML("<strong>Number of most relevant paths to show</strong> "),
             "(set it to 0 to show all paths):"),
    tags$div(id="freq2label-wrap",
      sliderInput("freq2label", "", width = "500px",
        value = 3, max = 5, min = 0, step = 1)
      )
    )
  })

  output$cpm_list <- renderUI({
    all_names <- c()
    for (i in names(all_cpm_out)) {
      all_names <- c(all_names, all_cpm_out[[i]]$orig_data$name)
    }

    ## FIXME: No default "user"
    if ((length(all_names) > 0) && (last_visited_cpm != "")) {
        selected <- names(all_cpm_out)
    ## if (length(all_names) > 0) {
        ## The first selected is dead code
        ## selected <- ifelse(is.null(input$select_csd),
        ##                    "user", input$select_csd)
        ## selected <- ifelse(input$select_csd %in%
        ##                    names(all_cpm_out),
        ##                    input$select_csd, "user")
        
        
        tagList(
            radioButtons(
                inputId = "select_cpm",
                label = "",
                selected = last_visited_cpm,
                choiceNames = names(all_cpm_out),
                choiceValues = names(all_cpm_out)
            )
        )
    }
  })

  output$csd <- renderPlot({
    evamtools:::plot_genotype_counts(evamtools:::get_csd(all_cpm_out[[input$select_cpm]]$cpm_output$analyzed_data))
  })

  output$original_data <- renderUI({
    ## To see if I disable original data
    if(length(names(all_cpm_out)) > 0){
      tags$div(class="frame max_height",
        tags$h3("3. Original data"),
        plotOutput("csd"),
        tags$div(class = "download_button",
          actionButton("modify_data", "Modify data")
        )
      )
    }
  })

  output$cpm_freqs <- DT::renderDT(all_cpm_out[[input$select_cpm]]$tabular_data[[input$tabular_data2show]],
    selection = 'none', server = TRUE
    , rownames = FALSE
    , options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE)
  )

  output$tabular_data <- renderUI({
    if(length(names(all_cpm_out)) > 0){
      tags$div(class="frame max_height",
        tags$h3("4. Tabular data"),
        radioButtons(inputId = "tabular_data2show",
          label = "",
          inline = TRUE,
          choiceNames =  c( "Transition probabilities",
                            "Transition rates",
                            "Predicted genotype relative frequencies",
                            "Sampled genotype counts",
                            "Observed transitions counts"
                            ),
          choiceValues =  c("trans_mat",
                            "trans_rate_mat",
                            "predicted_genotype_freqs",
                            "sampled_genotype_counts",
                            "obs_genotype_transitions"),
          selected = "trans_mat"
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
        saveRDS(all_cpm_out[[input$select_cpm]][c("cpm_output", "tabular_data")],
                file)
    }
  )

  # ## We only want the "Source" option enable if we have the data to show it
  # observeEvent(input$select_cpm, {
  #   selected <- ifelse(is.null(input$select_cpm), last_visited_cpm, input$select_cpm)
  #   if(is.null(all_cpm_out[[selected]]$Source_genotype_transitions)){
  #     updateCheckboxGroupInput(session, "cpm2show",
  #       selected = setdiff(c(input$cpm2show), "Source"))
  #     shinyjs::disable(selector = "#cpm2show input[value='Source']")
  #   }else{
  #     shinyjs::enable(selector = "#cpm2show input[value='Source']")
  #   }
  #   if(is.na(all_cpm_out[[selected]]$MCCBN_model)){
  #     updateCheckboxGroupInput(session, "cpm2show",
  #       selected = setdiff(c(input$cpm2show), "Source"))
  #     shinyjs::disable(selector = "#cpm2show input[value='MCCBN']")
  #   }else{
  #     shinyjs::enable(selector = "#cpm2show input[value='MCCBN']")
  #   }
  # })

  ## Upload button
  # observeEvent(input$output_cpms, {
  #   cpm_out <- readRDS(input$output_cpms$datapath)
  #   ## FIXME: f_graph to trans_rate_mat?
  #   # cpm_out$MHN_f_graph <- cpm_out$MHN_trans_rate_mat
  #   if(is.null(cpm_out$name)) cpm_out$name <- "User_Data"
  #   all_cpm_out[[cpm_out$name]] <- cpm_out
  #   updateRadioButtons(session, "select_cpm", selected = cpm_out$name)
  #   ## To see if I disable original data
  #   shinyjs::disable(selector = "#variable input[value='cyl']")
  # })
}
