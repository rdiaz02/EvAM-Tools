#'  I followed this link to structure the shiny app whithin the package
#'  https://deanattali.com/2015/04/21/r-package-shiny-app/

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
  examples_csd$csd <- examples_csd$csd[1:5]
  all_csd_data <- evamtools:::standarize_all_datasets(examples_csd)
  min_genes <- SHINY_DEFAULTS$min_genes
  max_genes <- SHINY_DEFAULTS$max_genes
  default_csd_samples <- SHINY_DEFAULTS$csd_samples
  default_cpm_samples <- SHINY_DEFAULTS$cpm_samples
  more_cpms <- NULL
  keep_dataset_name <- FALSE

  last_visited_pages <- list(csd = "User", dag = "User", matrix = "User")
  last_visited_cpm <- "user"

  datasets <- reactiveValues(
    all_csd = all_csd_data
  )

  adv_options <- reactiveValues(
    do_MCCBN = FALSE,
    cpm_samples = default_cpm_samples
  )

  data <- reactiveValues(
    csd_freqs = SHINY_DEFAULTS$template_data$csd_freqs
    , data = SHINY_DEFAULTS$template_data$data
    , dag = SHINY_DEFAULTS$template_data$dag
    , dag_parent_set = SHINY_DEFAULTS$template_data$dag_parent_set
    , lambdas = SHINY_DEFAULTS$template_data$lambdas
    , thetas = SHINY_DEFAULTS$template_data$thetas
    , n_genes = SHINY_DEFAULTS$ngenes
    , gene_names = LETTERS[1: max_genes]
  )

  display_freqs <- reactive({
    evamtools:::get_display_freqs(data$csd_freqs, input$gene_number, data$gene_names)
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
    tags$div(class="inlin2",
      textInput(inputId="dataset_name", "Give your dataset a name", value = input$select_csd)
    )
  })

  listen2dataset_name <- reactive({
    list(input$dataset_name, input$save_csd_data)
  })

  ## Saving dataset
  observeEvent(input$save_csd_data,{
    tryCatch({
      ## 1 save dataset to list after user data
      if(!(input$dataset_name %in% names(datasets$all_csd[[input$input2build]]))){
        datasets$all_csd[[input$input2build]][[input$dataset_name]]$name <- input$dataset_name

        if(nrow(data$csd_freqs) > 0){
          datasets$all_csd[[input$input2build]][[input$dataset_name]]$data <- evamtools:::freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
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
      }, error=function(e){
        showModal(dataModal(e[[1]]))
      })
  })

  # observeEvent(listen2dataset_name(), {
  observeEvent(input$dataset_name, {
    
    tryCatch({
      dataset_name <- ifelse(is.null(input$dataset_name), "", input$dataset_name)
      if(dataset_name != "" &
        !(dataset_name %in% names(datasets$all_csd[[input$input2build]]))) {
        shinyjs::enable("save_csd_data")
      }else if(dataset_name %in% names(datasets$all_csd[[input$input2build]])){
        shinyjs::disable("save_csd_data")
      }
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
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
    tryCatch({

      ## Cleaning stuf
      selected <- last_visited_pages[[input$input2build]]
      tmp_data <- datasets$all_csd[[input$input2build]][[selected]]
      data$gene_names <- tmp_data$gene_names
      data$data <- tmp_data$data

      shinyjs::disable("analysis")
      if(!is.null(data$data)){
        data$csd_freqs <- get_csd(data$data)
        shinyjs::enable("analysis")
      } else{
        data$csd_freqs <- SHINY_DEFAULTS$template_data$csd_freqs
      }

      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$dag_parent_set
      data$lambdas <- tmp_data$lambdas
      data$thetas <- tmp_data$thetas
      data$name <- tmp_data$name

      if(input$input2build == "dag"){
        to_keep <- length(which(colSums(data$dag)>0 | rowSums(data$dag)>0)) - 1
        n_genes <- ifelse(to_keep < 1 , ngenes, to_keep)
      } else if(input$input2build == "matrix"){
        n_genes <- length(which(colSums(abs(data$thetas))>0
        | rowSums(abs(data$thetas))>0))
        n_genes <- ifelse(n_genes <= 0, 3, n_genes)
      } else if (input$input2build == "csd" & !is.null(data$data)){
        n_genes <- ncol(data$data)
      } else if (input$input2build == "csd" & is.null(data$data)){
        n_genes <- SHINY_DEFAULTS$ngenes
      }

      updateNumericInput(session, "gene_number", value = n_genes)
      updateNumericInput(session, "genotype_freq", value = NA)
      updateCheckboxGroupInput(session, "genotype", label = "Mutations",
        choices = lapply(1:n_genes, function(i)data$gene_names[i]), selected = NULL)
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
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
      data$csd_freqs <- new_data$csd_freqs

      datasets$all_csd[[input$input2build]][[input$select_csd]] <- new_data
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
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

  # ## Advanced option for running evamtools
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      size = "l",
      easyClose = TRUE,
      title = tags$h3("Advanced options"),
      tags$div(
        numericInput("num_steps", "Sampling steps", adv_options$cpm_samples
          , min = 0, max = 100000, step = 100, width="100%"),
        # checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("HyperTRAPS", "MCCBN"), choiceValues = c("hypertraps", "mccbn")),
        checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("MCCBN"), choiceValues = c("MCCBN"), selected=c(adv_options$do_MCCBN)),
        tags$h4("DISCLAIMER: MCCBN may take hours to run")
        # tags$h4("DISCLAIMER: Both HyperTraps and MCCBN may take hours to run")
        )
      )
    )
    ## TODO this has no effect so far
  })

  observeEvent(input$num_steps, {
    adv_options$cpm_samples <- input$num_steps
  })

  observeEvent(input$more_cpms, {
    if("MCCBN" %in% input$more_cpms) {
      adv_options$do_MCCBN <- "MCCBN"
    } else {
      adv_options$do_MCCBN <- FALSE
    }
  }, ignoreNULL = FALSE)

  # ## Define number of genes
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
            actionButton("remove_edge", "Remove edge"),
            actionButton("clear_dag", "Clear dag"),
            tags$h3("DAG table"),
            DT::DTOutput("dag_table"),
            numericInput("dag_samples", "Total genotypes to sample", value = default_csd_samples, min= 100, max = 10000, step = 100, width = "50%"),
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
              numericInput("mhn_samples", "Total genotypes to sample", value = default_csd_samples, min= 100, max= 10000, step = 100, width = "50%"),
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

  # ## ÐAG builder
  # ## Controling dag builder
  dag_data <- reactive({
    all_gene_names <- c("Root", data$gene_names)
    edges <- which(data$dag == 1, arr.ind = TRUE)
    tmp_dag_parent_set <- data$dag_parent_set
    x <- length(tmp_dag_parent_set)
    ## I have to this weird thing because using data$gene_names does not work for some unkown reason
    names(tmp_dag_parent_set) <- all_gene_names[seq(2, x + 1)]
    dag_data <- data.frame(From = all_gene_names[edges[, "row"]]
      , To = all_gene_names[edges[, "col"]]
      , Relation = tmp_dag_parent_set[edges[, "col"] - 1]
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
    tryCatch({
      tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene, operation = "add", parent_set = data$dag_parent_set)
      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$parent_set
    },error=function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Remove edge
  observeEvent(input$remove_edge, {
    from_gene <- input$dag_from
    to_gene <- input$dag_to
    tryCatch({
      tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene, operation = "remove", parent_set = data$dag_parent_set)
      data$dag <- tmp_data$dag
      data$dag_parent_set <- tmp_data$parent_set
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
      data$csd_freqs <- SHINY_DEFAULTS$template_data$csd_freqs
      data$data <- SHINY_DEFAULTS$template_data$data
      data$dag_parent_set <- tmp_data$dag_parent_set
      data$lambdas <- SHINY_DEFAULTS$template_data$lambdas
      names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
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
      tmp_data <- evamtools:::modify_lambdas_and_parent_set_from_table(dag_data(), info, data$lambdas, data$dag, data$dag_parent_set)
      data$lambdas <- tmp_data$lambdas
      data$dag_parent_set <- tmp_data$parent_set
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Building trm from dag
  observeEvent(input$resample_dag, {
    tryCatch({
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Running evamtools", value = 0)

      progress$inc(1/2, detail = "Doing sampling")

      ## FIXME all the following will be replaced
      shinyjs::disable("resample_dag")
      tmp_data <- list(edges = dag_data())
      trm <- evamtools:::cpm2tm(tmp_data)$weighted_fgraph
      samples <- evamtools:::population_sample_from_trm(trm, input$dag_samples)
      process_data <- evamtools:::process_samples(samples, input$gene_number, data$gene_names[1:input$gene_number])
      tmp_samples <- process_data$sampled_genotype_freqs
      tmp_samples <- tmp_samples[tmp_samples > 0]
      tmp_csd <- data.frame(Genotype = names(tmp_samples), Counts = tmp_samples)
      rownames(tmp_csd) <- tmp_csd$Genotype
      data$csd_freqs <- tmp_csd
      data$data <- freqs2csd(tmp_csd,data$gene_names[1:input$gene_number])
      shinyjs::enable("resample_dag")
      progress$inc(1/2, detail = "Sampling Finished")

      datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
      datasets$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
      # datasets$all_csd[[input$input2build]][[input$select_csd]]$trm <- trm
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
      title = tags$h3("How to build a DAG"),
      tags$div(
        tags$p("Select a parent node and child node and hit 'Add edge' or 'Remove edge'."),
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
    tryCatch({
      info <-input$thetas_table_cell_edit
      data$thetas[1:input$gene_number, 1:input$gene_number] <-
        DT::editData(data$thetas[1:input$gene_number, 1:input$gene_number], info, "thetas")

      datasets$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
      ## Resample based on changes
      shinyjs::click("resample_mhn")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  observeEvent(input$resample_mhn, {
    tryCatch({
      ## FIXME this will be replaced by another function
      mhn_data <-get_mhn_data(input$gene_number, input$mhn_samples,
        data$gene_names[1:input$gene_number], thetas = data$thetas[1:input$gene_number, 1:input$gene_number])
      data$csd_freqs <- mhn_data$samples
      data$data <- freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])

      datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
      datasets$all_csd[[input$input2build]][[input$select_csd]]$trm <- data$trm
      shinyjs::enable("analysis")
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
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
    tryCatch({
      genotype <- paste(input$genotype, collapse = ", ")
      genot_freq <- data$csd_freqs[, 2][data$csd_freqs[, 1] == genotype]
      updateNumericInput(session, "genotype_freq", value = genot_freq)

    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  }, ignoreNULL = FALSE)

  observeEvent(input$add_genotype, {
    tryCatch({
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
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Genotypes table
  output$csd_freqs <- DT::renderDT(display_freqs(), selection = 'none', server = TRUE, editable = list(target = "column", disable = list(columns = c(0)))
    , rownames = FALSE,
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE),
  )

  observeEvent(input$csd_freqs_cell_edit, {
    tryCatch({
      info <- input$csd_freqs_cell_edit
      info[ , "col"] <- 2
      data$csd_freqs <- DT::editData(data$csd_freqs, info, "csd_freqs")
      ## Filtering out non-positive counts
      data$csd_freqs <- data$csd_freqs[data$csd_freqs[,2] > 0,]
      data$data <- datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- evamtools:::freqs2csd(data$csd_freqs, data$gene_names[1:input$gene_number])
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  ## Plot histogram of genotypes
  output$plot <- renderPlot({
    tryCatch({
      evamtools:::plot_genotypes_freqs(display_freqs())
    }, error = function(e){})
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
      shinyjs::disable("analysis")
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "Running evamtools", value = 0)

      data2run <- evamtools:::freqs2csd(display_freqs(), data$gene_names[1:input$gene_number])
      progress$inc(1/5, detail = "Setting up data")
      Sys.sleep(0.5)
      progress$inc(2/5, detail = "Running CPMs")

      methods <- SHINY_DEFAULTS$cpms2run
      if(!is.null(input$more_cpms)){
        methods <- unique(c(methods, input$more_cpms[input$more_cpms %in% SHINY_DEFAULTS$all_cpms]))
      }
      cpm_output <- evam(data2run, methods = methods)
   
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

      n_samples <- adv_options$cpm_samples
      if(is.null(n_samples) | !is.numeric(n_samples) | n_samples < 100){
        n_samples <- SHINY_DEFAULTS$cpm_samples
      }
      progress$inc(3/5, detail = paste("Running ", n_samples, " samples"))
      n_sample <- 100
      sampled_from_CPMs <- sample_CPMs(cpm_output, n_samples
        , methods, c("sampled_genotype_freqs", "obs_genotype_transitions"))

      progress$inc(4/5, detail = "Post processing data")
      Sys.sleep(0.5)

      orig_data <- list(data = data2run, name = data$name
        , type = input$input2build, gene_names = data$gene_names
        , thetas = data$thetas, lambdas = data$lambdas
        , dag = data$dag, dag_parent_set = data$dag_parent_set)

      ## Tabular data
      tabular_data <- evamtools:::create_tabular_data(c(cpm_output, sampled_from_CPMs))

      all_evam_output <- list("cpm_output" = c(cpm_output, sampled_from_CPMs)
        , "orig_data" = orig_data
        , "tabular_data" = tabular_data 
      )

      ##CPM output name
      result_index <- length(grep(sprintf("^%s", input$select_csd), names(all_cpm_out)))
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

  cpm_out <- sample_evam_output
  all_cpm_out <- reactiveValues(user = cpm_out)

  output$sims <- renderUI({
    if(length(names(all_cpm_out)) > 0 & !is.null(input$select_cpm)){
      tmp_data <- all_cpm_out[[input$select_cpm]]$cpm_output
      column_models2show <- floor(12 / length(input$cpm2show))

      lapply(input$cpm2show, function(met){
        method_data <- evamtools:::process_data(tmp_data, met, plot_type = "trans_mat")
        output[[sprintf("plot_sims_%s", met)]] <- renderPlot({
          pl <- evamtools::plot_method(method_data$method_info
          , method_data$parent_set
          , method_data$edges
          , met)
        })
        return(
          column(3,
            plotOutput(sprintf("plot_sims_%s", met)))
        )
        })
      }
  })

  output$sims2 <- renderUI({
    if(length(names(all_cpm_out)) > 0 & !is.null(input$select_cpm)){
      tmp_data <- all_cpm_out[[input$select_cpm]]$cpm_output
      ## Enabling donwload button
      shinyjs::enable(selector = "#download_cpm")

      ## Main display
      column_models2show <- floor(12 / length(input$cpm2show))
      selected_plot_type <- input$data2plot

      if(!(is.null(selected_plot_type))){
        lapply(input$cpm2show, function(met){
          method_data <- evamtools:::process_data(tmp_data, met, plot_type = selected_plot_type)
          output[[sprintf("plot_sims2_%s", met)]] <- renderPlot({
          pl <- evamtools::plot_genot_fg(method_data$data2plot,
                      observations = tmp_data$analyzed_data, # We use it to define "Observed" and "Not Observed" genotypes
                      predicted_genotypes = method_data$predicted_genotype_freqs, # To compute node sizes if sampled_freqs is NULL
                      sampled_freqs = method_data$sampled_genotype_freqs,
                      top_paths = top_paths,
                      freq2label = input$freq2label)
          })
          return(
            column(3,
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
    }, error = function(e){
      showModal(dataModal(e[[1]]))
    })
  })

  output$customize <- renderUI({
    tags$div(class = "frame",
      tags$h3("2. Customize the visualization"),
      tags$div(class = "inline",
        checkboxGroupInput(inputId = "cpm2show",
          label = "CPMs to show",
          choices = c("OT", "OncoBN", "CBN", "MHN", "HESBCN", "MCCBN"),
          selected = c("CBN", "MHN", "HESBCN")),

      tags$div(class = "inline",
        radioButtons(inputId = "data2plot",
          label = "Data to show",
          choiceNames =  c("Probabilities", "Transitions",
                            "Transition Rate Matrix"),
          choiceValues = c("trans_mat", "obs_genotype_transitions",
                            "trans_rate_mat"),
          selected = "obs_genotype_transitions"
          )
        ),
      ),
    tags$p("Label genotypes with frequency bigger than:"),
    tags$div(id="freq2label-wrap",
      sliderInput("freq2label", "", width = "500px",
        value = 0.05, max = 1, min = 0, step = 0.05)
      )
    )
  })

  output$cpm_list <- renderUI({
    all_names <- unname(sapply(all_cpm_out, function(dataset) dataset$name))

    if(length(all_names) > 0){
      selected <- ifelse(is.null(input$select_csd), "user", input$select_csd)
      selected <- ifelse(input$select_csd %in% names(all_cpm_out),input$select_csd, "user")

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
    plot_genotypes_freqs(get_csd(all_cpm_out[[input$select_cpm]]$cpm_output$analyzed_data))
  })

  output$original_data <- renderUI({
    ## To see if I disable original data
    if(length(all_cpm_out) > 0){
      tags$div(class="frame max_height",
        tags$h3("3. The original data"),
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
                            "Predicted genotype frequencies",
                            "Observed genotype frequencies",
                            "Observed transitions counts"
                            ),
          choiceValues =  c("trans_mat",
                            "trans_rate_mat",
                            "predicted_genotype_freqs",
                            "sampled_genotype_freqs",
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
      saveRDS(all_cpm_out[[input$select_cpm]], file)
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
