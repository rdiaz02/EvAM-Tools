## Copyright 2022 Pablo Herrera-Nieto, Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.

## #'  I followed this link to structure the shiny app whithin the package
## #'  https://deanattali.com/2015/04/21/r-package-shiny-app/

dataModal <- function(error_message, type="Error: ") {
    modalDialog(
        easyClose = TRUE,
        title = tags$h3(type),
        tags$div(
                 error_message
             )
    )
}


## I have left a bunch of messages. To make it easier to dis/enable them
## mymessage <- function(...) message(...)
mymessage <- function(...) invisible(NULL)



## For DRY though gene_names and number_of_genes are also
## obtained for other purposes right before these
## are called in a couple of places. Oh well; would need
## more serious refactoring.

## We reset the number of genes back to an OK number
dag_more_genes_than_set_genes <- function(input, dag_data = dag_data(),
                                          session = session) {
    gene_names <- setdiff(unique(c(dag_data$From, dag_data$To)),
                          "Root")
    number_of_genes <- length(gene_names)
    if (!is.null(input$gene_number) &&
        (number_of_genes > input$gene_number)) {
        updateNumericInput(session, "gene_number", value = number_of_genes)
        return(TRUE)
    }
    else
        return(FALSE)
}



dag_message_more_genes_than_set_genes <- function() {
    showModal(modalDialog(
        paste("The DAG contains more genes that the number ",
              "of genes you have set under ",
              "'Set the number of genes'. ",
              "Remove edges as needed ",
              "and then decrease number of genes."),
        easyClose = TRUE))
}

## data is the big object, with a $data inside. Yeah, great naming
csd_more_genes_than_set_genes <- function(input, data,
                                          session = session) {
    number_of_genes <- sum(colSums(data$data) > 0)
    if (!is.null(input$gene_number) &&
        (number_of_genes > input$gene_number)) {
        updateNumericInput(session, "gene_number", value = number_of_genes)
        return(TRUE)
    } else return(FALSE)
}

csd_message_more_genes_than_set_genes <- function() {
    showModal(modalDialog(
        paste("The genotype data contains more genes that the number ",
              "of genes you have set under ",
              "'Set the number of genes'. ",
              "Remove genotypes as needed ",
              "and then decrease number of genes."),
        easyClose = TRUE))
}


do_gc <- function(n = 5) {
    mymessage("doing gc")
    for (i in 1:n) print(gc())
}

random_dataset_name <- function() {
    paste(c(sample(letters, 1),
            sample(c(letters, 0:9), 8, replace = TRUE)),
          collapse = "")
}



## For assumption A1_gnn
## x: the big data object's matrix
set_gene_names_after_resize <- function(x, gene_names) {
    gene_names_num <- length(gene_names)
    gene_names_in_freqs <- sort(colnames(x))
    if (length(gene_names_in_freqs) == 0) {
        ## Only WT
        return(gene_names)
    } else {
        gene_names_wo_current <- sort(setdiff(gene_names, gene_names_in_freqs))
        gene_names <-
            c(gene_names_in_freqs, gene_names_wo_current)[1:gene_names_num]
        return(gene_names)
    }
}

## More code to solve the pervasive and silly assumption that we add genes in
## order. Now this affects the DAG
## We want to obtain the actual gene names in the DAG
## And no, this does not affect MHN because it ALWAYS uses all the genes
## Somewhat similar in purpose to set_gene_names_after_resize

## ## x: the "data" object (man, what a name!)
## gene_names_from_genes_in_DAG <- function(x, gene_names) {
##     gene_names_num <- length(gene_names)
##     the_dag <- x$data
##     in_dag <- (which(colSums(the_dag) > 0))
##     if (length(in_dag) == 0) {
##         return(gene_names)
##     } else {
##         genes_in_dag <- colnames(the_dag)[in_dag]
##         gene_names_wo_current <- sort(setdiff(gene_names, genes_in_dag))
##         gene_names <-
##             c(genes_in_dag, gene_names_wo_current)[1:gene_names_num]
##         return(gene_names)
##     }
## }







server <- function(input, output, session, EVAM_MAX_ELAPSED = 1.5 * 60 * 60) {
    require(evamtools)
    require(shinyBS)
    ## Just in case
    do_gc(2)
    ## And be paranoid about making sure memory is released on disconnect
    session$onSessionEnded(
                function() {
                    mymessage("From onSessionEnded")
                    do_gc(1)}
            )
    onStop(function() {
        cat("\n Finishing")
        do_gc(1)
    })
    
    
    ## Make these deps explicit. Needed for shinytests
    data("examples_csd", package = "evamtools")
    data("SHINY_DEFAULTS", package = "evamtools")
    ## We subset, as lots of the examples are not really worth it for the web app.
    examples_csd$csd <- examples_csd$csd[1:5]
    ## The next is (one of) the functions from hell.
    ## And the "upload" component disappears; it was an empty list anyway.
    all_csd_data <- evamtools:::to_stnd_csd_all_datasets(examples_csd)
    min_genes <- .ev_SHINY_dflt$min_genes
    max_genes <- .ev_SHINY_dflt$max_genes
    default_csd_samples <- .ev_SHINY_dflt$csd_samples
    default_cpm_samples <- .ev_SHINY_dflt$cpm_samples
    default_dag_model <- .ev_SHINY_dflt$dag_model

    last_visited_pages <- list(upload="Empty", csd = "Empty", dag = "DAG_Fork_3", matrix = "MHN_all_0")

    last_visited_cpm <- ""

    datasets <- reactiveValues(
        all_csd = all_csd_data
    )

    data <- reactiveValues(
        csd_counts = .ev_SHINY_dflt$template_data$csd_counts
      , data = .ev_SHINY_dflt$template_data$data
      , dag = .ev_SHINY_dflt$template_data$dag
      , dag_parent_set = .ev_SHINY_dflt$template_data$dag_parent_set
      , lambdas = .ev_SHINY_dflt$template_data$lambdas
      , thetas = .ev_SHINY_dflt$template_data$thetas
      , n_genes = .ev_SHINY_dflt$ngenes
      , gene_names = LETTERS[1: max_genes]
    )

    display_freqs <- reactive({
        ## Remember this is called whenever changes in many places
        ## happen.
        
        ## This is often called when there is no need for it. So when you change
        ## the type of data entering (go from MHN to upload, for example) this is
        ## called again and returns the data but for nothing since what we will
        ## want to display are the new data. We now elegantly return a
        ## 0-rows data frame when nothing should be returned. Explicit and clean.

        ## For paranoia, we return always things in standard order Probably
        ## unnecessary, but just in case.
        mymessage("At display_freqs")

        if (input$input2build == "dag") {
            ## With the DAG we always return all the genotypes
            ## Other code ensures that gene number is never smaller
            ## than genes in the DAG.
            mymessage("      dag")
            dag_dataset_names <- unlist(lapply(datasets$all_csd$dag,
                                               function(x) x$name))
            if (!(data$name %in% dag_dataset_names)) {
                mymessage("       data$name not in dag_dataset_names. ",
                          "Returning a 0-rows data frame")
                return(data.frame(Genotype = character(), Counts = integer()))
            }

            ## FIXME: See comment below  A1_gnn_bisfix
            ## if (!is.null(data$n_genes) &&
            ##     input$gene_number != data$n_genes) {
            ##     mymessage("      DAG: Updating gene names in data")
            ##     new_gnames2 <- set_gene_names_after_resize(data$data,
            ##                                                data$gene_names)
            ##     data$gene_names <- new_gnames2
            ## }
            return(
                evamtools:::reorder_to_standard_order_count_df(
                                data$csd_counts[data$csd_counts$Counts > 0, ,
                                                drop = FALSE]))
        } else if (input$input2build == "upload") {
            mymessage("      upload")
            upl_dataset_names <- unlist(lapply(datasets$all_csd$upload, function(x) x$name))
            if (is.null(data$name) || !(data$name %in% upl_dataset_names)) {
                mymessage("       data$name not in upl_dataset_names.",
                          "Returning a 0-rows data frame")
                return(data.frame(Genotype = character(), Counts = integer()))
            }
            ## With upload, we do not use number of genes
            ## Return the data we have
            return(
                evamtools:::reorder_to_standard_order_count_df(
                                data$csd_counts[data$csd_counts$Counts > 0, ,
                                                drop = FALSE]))
        } else {
            mymessage("      mhn or csd")
            thisd <- input$input2build
            thisd_dataset_names <- unlist(lapply(datasets$all_csd[[thisd]],
                                                 function(x) x$name))
            if (!(data$name %in% thisd_dataset_names)) {
                mymessage("       data$name not in thisd_dataset_names. ",
                          "Returning a 0-rows data frame")
                return(data.frame(Genotype = character(), Counts = integer()))
            }

            ## FIXME: A1_gnn_bisfix The code in A1_gnn_0 is not updating
            ## data$etc. That is triggered on gene number change, but not
            ## necessarily on adding a genotype or removing a genotype, in ways
            ## that change the genes, but not the gene number.
            ## Though something I do not fully understand:
            ## It does not happen in the single observeEvent for input$gene_number
            ## and not in the updateNumericInput for input$gene_number
            ## So force it here if there have been changes in input$gene_number
            ## and if they haven't, for some weird, hard to reproduce, cases.
            ## Yes, this seem necessary to prevent BUG_Create_Rename_Click_other
            if ((input$input2build == "csd") &&
                !is.null(data$n_genes) ) {
                ## input$gene_number != data$n_genes) {
                mymessage("       CSD: Updating gene names in data")
                new_gnames2 <- set_gene_names_after_resize(data$data,
                                                           data$gene_names)
                if (identical(new_gnames2, data$gene_names)) {
                    message("A1_gnn_bisfix: identical")
                } else {
                    message("A1_gnn_bisfix: different")
                }
                data$gene_names <- new_gnames2
            }

            return(
                evamtools:::reorder_to_standard_order_count_df(
                                evamtools:::get_display_freqs(data$csd_counts,
                                                              input$gene_number,
                                                              data$gene_names,
                                                              input$input2build))
            )
        }
    })

    ## Force resample on gene number changes
    ## Can I comment the next block entirely?
    ## And the next?
    ## That removes the plotting twice issue
    ## but makes this less interactive: you must click
    ## "Sample".
    
    ## Can't make it depend on gene_number too
    ## or we get the "DAG contains more genes ..."

    ## ## Comment this out so that no resampling when renaming?
    ## ## FIXME clarify: issue_11
    ## observeEvent(input$select_csd, {
    ##     ## mymessage("at select_csd_trigger")
    ##    if (input$input2build == "dag") {
    ##         shinyjs::click("resample_dag")
    ##    } else if (input$input2build == "matrix") {
    ##        shinyjs::click("resample_mhn")
    ##    }
    ## }
    ## ## And now, display_freqs will likely be called
    ## )
    

    
    observeEvent(input$gene_number, {
        ## id: here_we_change_gene_number
        mymessage("at gene number_trigger")
        datasets$all_csd[[input$input2build]][[input$select_csd]]$n_genes <-
            input$gene_number
        if (input$input2build == "csd") {
            ## A1_gnn_0
            ## Set gene names correctly if there have been changes
            ## Counterpart to check A1_gnn in get_display_freqs
            ## FIXME: is this really necessary? Or is the similar code
            ## in display_freqs enough?
            ## Needed to prevent bug BUG_Create_Add_E_Decrease
            ## DAG does not, as it continuously monitors genes in the DAG model
            
            current_data <- datasets$all_csd[[input$input2build]][[input$select_csd]]
            new_gnames <-
                set_gene_names_after_resize(current_data$data,
                                            current_data$gene_names)
            data$gene_names <- new_gnames

            ## The next is not really necessary. And it is actually limiting what
            ## users can do. But forcing this: a) decreases likelihood of bugs;
            ## b) Makes behavior consistent with that of DAGs.
            if (csd_more_genes_than_set_genes(input, data, session)) {
                csd_message_more_genes_than_set_genes()
            }
        } else if (input$input2build == "dag") {
            if (dag_more_genes_than_set_genes(input, dag_data(), session)) {
                dag_message_more_genes_than_set_genes()
            }
        } else if (input$input2build == "matrix") {
            shinyjs::click("resample_mhn")
        }
    }
    ## And now, display_freqs will likely be called
    )


    
    ## Upload data
    observeEvent(input$csd, {
        if (grepl(".csv", input$csd$datapath)) {
            tryCatch({
                dataset_name <- input$name_uploaded
                ## repeated from obserEvent(input$save_csd_data
                if (gsub(" ", "", dataset_name, fixed = TRUE) == "") {
                    stop("Name of data cannot be an empty string")
                }
                if (grepl(" ", dataset_name, fixed = TRUE)) {
                    stop("Name of data should not contain spaces")
                }
                existing_names <- c(names(datasets$all_csd$upload),
                                    names(datasets$all_csd$csd),
                                    names(datasets$all_csd$dag),
                                    names(datasets$all_csd$matrix))
                if (dataset_name %in% existing_names) {
                    stop("That name is already in use for other data ",
                         "(if not in this input type, maybe in one ",
                         "of the others); that can be confusing. ",
                         "Use a different name for the uploaded data.")
                }
                tmp_data <- list()
                tmp_data$data <- try(read.csv(input$csd$datapath))
                if (inherits(tmp_data$data, "try-error")) {
                    if (grepl("no lines available in input", tmp_data$data)) {
                        stop("The uploaded data contains no valid input. ",
                             "Did you upload an empty file?")
                    } else {
                        stop("Error reading the uploaded data. ",
                             "Please ensure the file is of the correct format ",
                             "and it contains valid data.")
                    }
                }

                tmp_data$name <- dataset_name

                datasets$all_csd[["upload"]][[dataset_name]] <-
                    evamtools:::to_stnd_csd_dataset(tmp_data)
                datasets$all_csd[["upload"]][[dataset_name]]$name <- dataset_name
                tmp_data$gene_names <- c(
                    colnames(tmp_data$data)
                  , LETTERS[(length(colnames(tmp_data$data)) + 1):max_genes]
                )
                                        # tmp_data$gene_names <- colnames(tmp_data$data)
                tmp_data$n_genes <- ncol(tmp_data$data)
                datasets$all_csd[["upload"]][[dataset_name]] <- tmp_data
                
                last_visited_pages["upload"] <<- dataset_name
                updateRadioButtons(session, "input2build", selected = "upload")
                updateRadioButtons(session, "select_csd", selected = dataset_name)
            }, error = function(e){
                showModal(dataModal(e[[1]]))
            })
        } ## else if(grepl(".rds", input$csd$datapath, ignore.case = TRUE)){
        ## We do not accept input as rds. Nope.
        ##     tmp_data <- readRDS(input$csd$datapath)
        ##     tryCatch({
        ##         new_data <- evamtools:::to_stnd_csd_dataset(tmp_data)
        ##         datasets$all_csd[[tmp_data$type]][[new_data$name]] <- new_data
        ##         last_visited_pages[tmp_data$type] <<- tmp_data$name
        ##         updateRadioButtons(session, "input2build", selected = tmp_data$type)
        ##         updateRadioButtons(session, "select_csd", selected = tmp_data$name)
        ##     }, error = function(e){
        ##         showModal(dataModal(e[[1]]))
        ##     })
        ## }
    })

    observeEvent(input$input2build, {
        updateRadioButtons(session, "select_csd",
                           selected = last_visited_pages[[input$input2build]])
    })

    ## Define dataset name
    output$dataset_name <- renderUI({
        if(input$input2build %in% c("upload", "csd", "dag", "matrix")){
            tags$div(class = "frame inlin2",
                     tags$h3("Rename the data"),
                     tags$h5(HTML("Give the (modified) data a different name ",
                                  "that will also be used to save the CPM ",
                                  "output.")),
                     tags$div(class = "download_button",
                              ),
                     textInput(inputId = "dataset_name",
                               "Give your data a name",
                               value = random_dataset_name()
                               ),
                     actionButton("save_csd_data", "Rename the data")
                     )
        }
    })

    output$download_data <- renderUI({
        if(input$input2build %in% c("upload", "csd", "dag", "matrix")){
            tags$div(
                     class = "frame",
                     tags$h3("Download the data"),
                     tags$div(class = "download_button",
                              tags$h5(HTML("Contents of saved file: ",
                                           "the data as data frame; ",
                                           "if you built a DAG or MHN model, ",
                                           "also the model built."
                                           )),  
                              downloadButton("download_csd", "Download your data")
                              )
                 )
        }
    })


    ## Saving dataset
    observeEvent(input$save_csd_data, {
        tryCatch({
            ## 1 save dataset to list after user data
            if (gsub(" ", "", input$dataset_name, fixed = TRUE) == "") {
                stop("Name of data cannot be an empty string")
            }
            if (grepl(" ", input$dataset_name, fixed = TRUE)) {
                stop("Name of data should not contain spaces")
            }
            existing_names <- c(names(datasets$all_csd$upload),
                                names(datasets$all_csd$csd),
                                names(datasets$all_csd$dag),
                                names(datasets$all_csd$matrix))
            
            if (input$dataset_name %in% existing_names) {
                stop("That name is already in use for other data ",
                     "(if not in this input type, maybe in one ",
                     "of the others); that can be confusing. ",
                     "Use a different name.")
            }

            ## Is this ever possible?
            if ( (nrow(data$csd_counts) > 0 ) &&  (is.null(data$data))) {
                stop("This should not be possible: ",
                     "(nrow(data$csd_counts) > 0 ) &&  (is.null(data$data))")
            }
            if ( (nrow(data$csd_counts) > 0 ) && all(is.na(data$data)) ) {
                if (data$csd_counts["Genotype"] == "WT") {
                    mymessage("Data with only WT")
                }  else {
                    stop("This should not be possible",
                         "These are not data with only WT and yet ",
                         "(nrow(data$csd_counts) > 0 ) && all(is.na(data$data))"
                         ) 
                }
            }
            
            ## A couple of paranoid consistency checks
            if (nrow(data$csd_counts) > 0) {
                ddtmp <- evamtools:::genotypeCounts_to_data(data$csd_counts,
                                                            e = 0)
                ddtmp2 <- data$data[, colnames(ddtmp), drop = FALSE]
                stopifnot(all.equal(colSums(ddtmp2), colSums(ddtmp)))
                c_ddtmp2 <- evamtools:::data_to_counts(ddtmp2, "data.frame",
                                                       omit_0 = TRUE)
                rownames(c_ddtmp2) <- c_ddtmp2$Genotype
                csdc_clean <- data$csd_counts[data$csd_counts$Counts > 0, ]
                rownames(csdc_clean) <- csdc_clean$Genotype
                stopifnot(all.equal(c_ddtmp2[csdc_clean[, "Genotype"], ],
                                    csdc_clean))
                rm(ddtmp2, ddtmp)
            }

            if ((input$input2build == "upload") && is.null(data$data)) {
                stop("When no data has been uploaded, ",
                     "it makes no sense to rename the data: ",
                     "there is nothing you could do with them, ",
                     "since there are none.")
            }

            n_genes <- ifelse(input$input2build == "upload",
                       ifelse(is.null(data$data), 0, ncol(data$data)),
                       input$gene_number)
            
            tmp_data <- list(
                data = data$data
              , dag = data$dag
              , gene_names = data$gene_names
              , dag_parent_set = data$dag_parent_set
              , lambdas = data$lambdas
              , thetas = data$thetas
              , trm = data$trm
              , n_genes = n_genes
              , name = input$dataset_name)

            datasets$all_csd[[input$input2build]][[input$dataset_name]] <- tmp_data

            ## 3 update selected entry
            updateRadioButtons(session, "select_csd",
                               selected = input$dataset_name)

        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })


    ## Download csd button
    output$download_csd <- downloadHandler(
        filename = function() sprintf("%s_data.rds", input$select_csd),
        content = function(file) {
            tmp_data <- datasets$all_csd[[input$input2build]][[input$select_csd]]
            if (input$input2build == "csd") {
                saveRDS(tmp_data$data, file=file)
            } else if(input$input2build == "dag") {
                gene_names <- setdiff(unique(c(dag_data()$From, dag_data()$To)),
                                      "Root")
                number_of_genes <- length(gene_names)
                stopifnot(number_of_genes == input$gene_number)
                data2save <- list(
                    data = tmp_data$data[, 1:number_of_genes]
                  , model_edges = dag_data()
                  , model = default_dag_model
                  , dag_parent_set = tmp_data$dag_parent_set[1:number_of_genes]
                  , dag = tmp_data$dag[1:(number_of_genes + 1),
                                       1:(number_of_genes + 1)])
                saveRDS(data2save, file=file)
            } else if (input$input2build == "matrix") {
                data2save <-
                    list(data = tmp_data$data[,1:input$gene_number],
                         log_Theta_matrix = tmp_data$thetas[1:input$gene_number,
                                                            1:input$gene_number])
                saveRDS(data2save, file=file)
            }        }
    )

    observeEvent(input$select_csd, {
        tryCatch({
            ## Ufff!!! <<-
            last_visited_pages[[input$input2build]] <<- input$select_csd
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Display List of availabe CSD
    output$csd_list <- renderUI({
        if(length(names(datasets$all_csd[[input$input2build]]))>0){
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
        } else {
            tags$p("Empty until you upload a data file.")
        }
    })


    toListen <- reactive({
        list(input$input2build, input$select_csd)
    })

    ## This is probably abusing observeEvent? And mixing what would be
    ## better served with eventReactive?

    ## Notes about the logic, and the update of display_freqs.  display_freqs
    ## needs to be updated whenever we make changes that need to be replotted. It
    ## is called from output$plot. That is much more often than the events below
    ## need to be watched. That is also why the following happens:
    ## - you are at MHN
    ## - you change to upload
    ##     The first thing that gets called is display_freqs (as output$plot is called)
    ##     Then, the block below.
    ## I've left some messages (commented now), so that one can see what is happening.
    
    observeEvent(toListen(), {
        tryCatch({
            ## The next two we are observing on
            ## input$select_csd
            ## input$input2build

            mymessage("At observeEvent toListen")

            ## Cleaning stuff
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
                data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            }

            data$dag <- tmp_data$dag
            data$dag_parent_set <- tmp_data$dag_parent_set
            data$lambdas <- tmp_data$lambdas
            data$thetas <- tmp_data$thetas
            data$name <- tmp_data$name
            data$n_genes <- tmp_data$n_genes


            if (input$input2build == "dag") {
                number_of_parents <- colSums(data$dag)
                to_keep <- sum(number_of_parents > 0)
                n_genes <- ifelse(to_keep < 1, .ev_SHINY_dflt$ngenes, to_keep)
                updateRadioButtons(session, "dag_model", selected = "HESBCN")
            } else if (input$input2build == "matrix") {
                n_genes <- data$n_genes
                if (is.null(n_genes)) {
                    n_genes <- .ev_SHINY_dflt$ngenes
                } 
            } else if (input$input2build %in% c("csd", "upload")) {
                if (!is.null(data$data))  {
                    n_genes <- ncol(data$data)  
                } else { ## data$data is null
                    if (input$input2build == "csd") {
                        n_genes <- .ev_SHINY_dflt$ngenes
                    } else { ## if (input$input2build == "upload")
                        n_genes <- NULL
                    }
                }
            } else {
                stop("No appropriate input2build")
            }

            if (input$input2build %in% c("csd", "dag", "matrix")) {
                ## number of genes, in the "Set the number of genes"
                updateNumericInput(session, "gene_number", value = n_genes)
            } else if (input$input2build == "upload") {
                ## Set gene_number to NULL. Never used for upload
                updateNumericInput(session, "gene_number", value = NULL)
            }

            if (input$input2build %in% c("csd")) {
                ## Where we "Add genotypes" manually. This selects mutation (genotype)
                ##         If we only have WT, n_genes is 0, and it breaks
                ##         and we have set minimum number of genes to 2
                ## id_change_genotype_muts
                gene_options <- set_gene_names_after_resize(data$data,
                                                            data$gene_names)[1:n_genes]
                updateCheckboxGroupInput(session, "genotype", label = "Mutations",
                                         choices = lapply(1:(max(2, n_genes)),
                                                          function(i) gene_options[i]),
                                         selected = NULL)
                ## Where we "Add genotypes" manually. This sets the count
                updateNumericInput(session, "genotype_freq", value = NA)
            }

        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    observeEvent(input$display_help_change_genotype_counts, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("Changing genotype's counts"),
            tags$div(
                     tags$p("1. Double click in a Counts cell to edit it"),
                     tags$p("2. Press Tab to move to the next row"),
                     tags$p("3. Use Ctrl + Enter to save changes"),
                     tags$p("4. Set a frequency to 0 to remove a genotype"),
                     tags$p("5. Type in the Search bar to filter genotypes"),
                     tags$h4(HTML("<br/>")),
                     tags$p("Genotypes with count 0 are removed from the table. ",
                            "Thus, if you remove a genotype when editing ",
                            "genotype's counts in the DAG, MHN, or Upload data ",
                            "entries, you will need to regenerate the data ",
                            "to be able to modify those genotypes again.")
                 )
        )
        )
    })


    ## ## Updating gene names
    ## observeEvent(input$action_gene_names,{
    ##     tryCatch({
    ##         new_gene_names <-
    ##             strsplit(gsub(" ", "", input$new_gene_names), ",")[[1]]
    ##         if (isTRUE(any(duplicated(new_gene_names)))) {
    ##             stop("Duplicated new gene names.")
    ##         }
    ##         if (length(data$gene_names[1:input$gene_number]) !=
    ##             length(new_gene_names)) {
    ##             stop("Number of old and new gene names differs.")
    ##         }
    
    ##         ## Use a simple lookup-dictionary and 
    ##         ## avoid to_stnd_csd_dataset which is a function from hell.
    ##         old_gene_names <- data$gene_names
    ##         new_gene_names <- c(new_gene_names,
    ##                             LETTERS[(length(new_gene_names) + 1):max_genes]
    ##                             )
    ##         names_dict <- new_gene_names
    ##         names(names_dict) <- old_gene_names
    ##         ## For the DAG
    ##         names_dict <- c(names_dict, "Root" = "Root")
    
    ##         new_data <- list()
    ##         new_data$gene_names <- new_gene_names
    ##         new_data$name <- data$name
    ##         new_data$lambdas <- data$lambdas
    ##         new_data$dag_parent_set <- data$dag_parent_set
    ##         new_data$dag <- data$dag
    ##         new_data$thetas <- data$thetas
    ##         new_data$data <- data$data

    ##         ## To rename, use lookup
    ##         names(new_data$lambdas) <- names_dict[names(new_data$lambdas)]
    ##         names(new_data$dag_parent_set) <- names_dict[names(new_data$dag_parent_set)]
    ##         colnames(new_data$dag) <- names_dict[colnames(new_data$dag)]
    ##         rownames(new_data$dag) <- names_dict[rownames(new_data$dag)]
    ##         colnames(new_data$thetas) <- names_dict[colnames(new_data$thetas)]
    ##         rownames(new_data$thetas) <- names_dict[rownames(new_data$thetas)]
    ##         if (!is.null(new_data$data)) {
    ##             colnames(new_data$data) <- names_dict[colnames(new_data$data)]
    ##         }
    ##         ## To create
    ##         new_data$csd_counts <- get_csd(new_data$data)
    
    ##         ## Assign to the correct places
    ##         data$gene_names <- new_gene_names
    ##         data$data <- new_data$data
    ##         data$dag <- new_data$dag
    ##         data$dag_parent_set <- new_data$dag_parent_set
    ##         data$thetas <- new_data$thetas
    ##         data$lambdas <- new_data$lambdas
    ##         data$csd_counts <- new_data$csd_counts
    
    ##         datasets$all_csd[[input$input2build]][[input$select_csd]] <- new_data
    
    ##     }, error = function(e){
    ##         showModal(dataModal(e[[1]]))
    ##     })
    ## })

    
    ## Advanced option for running evamtools
    observeEvent(input$advanced_options, {
        shinyjs::toggle("all_advanced_options")
    })

    ## Define number of genes
    output$gene_number_slider <- renderUI({
        val <- ifelse(is.null(data$n_genes), 3, data$n_genes)

        if (input$input2build %in% c("csd","dag", "matrix")) {
            tags$div(class = "frame flex",
                     tags$h3("Set the number of genes"),
                     tags$h5("(Using 7 or more genes can lead ",
                             "to very long execution times for some methods ",
                             "and crowded figures.)"),
                     
                     tags$div(class="inlin",
                              tags$h3(HTML("<br/>")),
                              sliderInput("gene_number", "Number of genes",
                                          value = val, max = max_genes, min = min_genes,
                                          step = 1),
                              ## The action that takes place is
                              ## id: here_we_change_gene_number
                              ),
                     tags$h4(HTML("<br/>")),
                     ## actionButton("change_gene_names", "Change gene names"),
                     )
        } 
    })

    ## observeEvent(input$change_gene_names, {
    ##     if (input$input2build == "dag") {
    ##         gene_names_00 <- gene_names_from_genes_in_DAG(data, data$gene_names)
    ##     } else if (input$input2build == "csd") {
    ##         gene_names_00 <- set_gene_names_after_resize(data, data$gene_names)
    ##         ## Or else, it is broken in other places
    ##         data$gene_names <- gene_names_00
    ##     } else {
    ##         gene_names_00 <- data$gene_names
    ##     }

    ##     showModal(modalDialog(
    ##         title = tags$h3("Change gene names"),
    ##         tags$div(class = "inlin2",
    ##                  textInput(inputId = "new_gene_names", "Gene names",
    ##                            value = paste(gene_names_00[1:input$gene_number],
    ##                                          collapse = ", ")
    ##                            ),
    ##                  tags$h4(HTML("<br/>")),
    ##                  tags$h4("Separate you gene names with a ','. ",
    ##                          "Do no use 'WT' for any gene name. ",
    ##                          "Use only alphanumeric characters ",
    ##                          "(of course, do not use comma as part of a gene name), ",
    ##                          "and do not start ",
    ##                          "a gene name with a number; ",
    ##                          "keep gene names short (for figures)."
    ##                          ),
    ##                  tags$h4(HTML("<br/>")),
    ##                  tags$div(class = "download_button",
    ##                           tags$h4(HTML("<br/>")),
    ##                           actionButton("action_gene_names", "Change genes names"),
    ##                           )
    ##                  ),
    ##         easyClose = TRUE
    ##     ))
    ## })

    
    ## Define new genotype
    observeEvent(input$dag_model, {
        number_of_parents <- colSums(data$dag)
        if (input$dag_model == "OT" && any(number_of_parents > 1)) {
            showModal(dataModal(
                paste("This DAG cannot be transformed into a tree. ",
                      "Are there nodes with multiple parents? ",
                      "(OT cannot not have nodes with multiple parents.)")))
            updateRadioButtons(session, "dag_model", selected = "HESBCN")
        } else {
            default_dag_model <<- input$dag_model
        }
    })

    output$define_genotype <- renderUI({
        n_genes <- ifelse(is.null(input$gene_number), 3, input$gene_number)

        ## Setting it here ain't enough
        ## This can all be changed in other two lines at least too.
        ##   Search for id_change_genotype_muts
        ## gene_options <- data$gene_names[1:n_genes]
        gene_options <- set_gene_names_after_resize(data$data,
                                                    data$gene_names)[1:n_genes]


        if (input$input2build == "csd") {
            tags$div(
                     tags$h3("Add genotypes"),
                     tags$h5("WT is added by not clicking on any mutations. "),
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
                                                 choices = set_gene_names_after_resize(data$data,
                                                                                       data$gene_names)[1:n_genes])
                              ## gene_options)
                              ),
                     tags$div(id="fr",
                              numericInput(label = "Counts", value = NA, min = 0,
                                           inputId = "genotype_freq", width = NA),
                              actionButton("add_genotype", "Add genotype")
                              ),
                     )
        } else if (input$input2build == "dag") {
            ## Make sure all genes currently in the DAG are in
            ## the To and From to add/remove
            mymessage("At output$define_genotype, DAG")
            current_dag_data <- dag_data()
            ## genes_in_dag <- setdiff(unique(c(dag_data()$From, dag_data()$To)),
            ##                         "Root")
            if (is.null(current_dag_data)) mymessage("   current_dag_data is NULL")
            genes_in_dag <- setdiff(unique(c(current_dag_data$From,
                                             current_dag_data$To)),
                                    "Root")
            if (length(genes_in_dag) < n_genes) {
                genes_not_in_dag <- setdiff(gene_options, genes_in_dag)
                dag_gene_options <- c(genes_in_dag,
                                      genes_not_in_dag[1:(n_genes - length(genes_in_dag))])
            } else {
                dag_gene_options <- genes_in_dag
            }
            tags$div(                    
                     tags$div(class = "flex",
                              tags$h3("Define a Directed Acyclic Graph (DAG) ",
                                      "and generate data from it."),
                              actionButton("how2build_dag", "Help", class = "btn-info")
                              ),
                     if(!is.null(data$lambdas)){
                         tags$div(
                                  tags$h4(HTML("<br/>")),
                                  tags$h4(HTML("<u>1. Define DAG</u>")),
                                  tags$h4(HTML("<br/>")),
                                  tags$h4("Type of model"),
                                  tags$div(class = "inline",
                                           radioButtons(inputId = "dag_model",
                                                        label = "Model: ",
                                                        inline = TRUE,
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
                                                        choices =  c("Root", dag_gene_options))
                                           ),
                                  tags$div(class = "inline",
                                           radioButtons(inputId = "dag_to",
                                                        label = " To (child node)",
                                                        inline = TRUE,
                                                        choices =  dag_gene_options)
                                           ),
                                  tags$h5(HTML("<p></p>")),
                                  actionButton("add_edge", "Add edge"),
                                  actionButton("remove_edge", "Remove edge"),
                                  tags$h5(HTML("If you want to decrease the number of genes ",
                                               "first remove edges and nodes from the DAG ", 
                                               "and only then modify ",
                                               "'Set the number of genes'.",
                                               "(We cannot know which edges/nodes ",
                                               "you want to remove).")),
                                  tags$h5(HTML("If you want to increase the number of genes ",
                                               "use 'Set the number of genes' to increase ",
                                               "the available gene labels, and then ",
                                               "increase the number of nodes in the DAG")),
                                  tags$h3(HTML("<br/>DAG table")),
                                  tags$h4(HTML("Remember to hit Ctrl-Enter when you are done editing the DAG table for changes to take effect.")),
                                  DT::DTOutput("dag_table"),
                                  tags$h3(HTML("<br/>")),
                                  tags$h4(HTML("<br/>")),
                                  tags$h4(HTML("<u>2. Generate data from the DAG model</u>")),
                                  tags$h4(HTML("<br/>")),
                                  numericInput("dag_epos",
                                               HTML("epos,&epsilon;"),
                                               value = 0.01, min = 0, max = 1,
                                               step = 0.005, width = "50%"),
                                  tags$h5(HTML("For OT (epos) and OncoBN (&epsilon;) only: prob. of children "),
                                          "not allowed by model to occur. ",
                                          "(Affects predicted probabilities.) "),
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;",
                                      numericInput("dag_samples", HTML("Number of genotypes<br>to sample"),
                                                   value = default_csd_samples, min = 100, max = 10000,
                                                   step = 100, width = "70%"),
                                      ), 
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;", 
                                      numericInput("dag_noise",
                                                   HTML("Observational noise <br>(genotyping error)"),
                                                   ## HTML("Noise ",
                                                   ##      "<h5>Observational noise (genotyping error), ",
                                                   ##      "a proportion between 0 and 1.</h5>"),
                                                   value = 0.0, min = 0, max = 1,
                                                   step = 0.025, width = "70%"),
                                      )
                                  |> prompter::add_prompt(
                                                   message = paste("A proportion between 0 and 1. ", 
                                                                   "Observational noise (e.g., genotyping error) ",
                                                                   "for all models. ",
                                                                   "Added during sampling, ",
                                                                   "after predictions from model ",
                                                                   "have been obtained; ",
                                                                   "predicted probabilities are not affected.",
                                                                   "If larger than 0, this proportion of entries ",
                                                                   "in the sampled matrix will be flipped ",
                                                                   "(i.e., 0s turned to 1s and 1s turned to 0s)."
                                                                   ),
                                                   position = "right",
                                                   rounded = TRUE,
                                                   bounce = TRUE,
                                                   size = "medium"
                                               ),
                                  tags$h5(HTML("<br/>")),
                                  actionButton("resample_dag", "Generate data from DAG"),
                                  actionButton("clear_dag", HTML("Reset DAG and delete genotype data"))
                                  |> prompter::add_prompt(message = 
                                                              HTML("Resetting the DAG will replace the ",
                                                                   "contents of the named object by ",
                                                                   "those of the default one ",
                                                                   "(a three-gene fork with lambdas = 0.5)"),
                                                          position = "right",
                                                          rounded = TRUE,
                                                          bounce = TRUE,
                                                          size = "medium"
                                                          ),
                                  )
                     }
                 )
        } else if (input$input2build == "matrix") {
            tags$div(
                     tags$div(class = "flex",
                              ## tags$h3("2. Define input with a Matrix"),
                              tags$h3("Define MHN's log-Theta",
                                      HTML("matrix (log-&Theta;) ",
                                           "and generate data from it.")),
                              actionButton("how2build_matrix", "Help", class = "btn-info")
                              ),
                     if (!is.null(data$thetas)){
                         tags$div(
                                  tags$h3(HTML("<br/>")),
                                  tags$h4(HTML("<u>1. Define MHN's &theta;s</u>")),
                                  tags$h4("Entries are ",
                                          "lower case thetas, ",
                                          HTML("&theta;s, range &plusmn; &infin;"),),
                                  tags$h4(HTML("Remember to hit Ctrl-Enter when you are done editing the matrix for changes to take effect.")),
                                  DT::DTOutput("thetas_table"),
                                  tags$h3(HTML("<br/>")),
                                  tags$h4(HTML("<u>2. Generate data from the MHN model</u>")),
                                  tags$h4(HTML("<br/>")),
                                  div(style = "white-space: nowrap;", 
                                      numericInput("mhn_samples",
                                                   HTML("Number of genotypes<br>to sample"),
                                                   value = default_csd_samples, min = 100, max = 10000,
                                                   step = 100, width = "70%"),
                                      ),
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;",
                                      numericInput("mhn_noise",
                                                   HTML("Observational noise<br>(genotyping error)"),
                                                   ## HTML("Noise"),
                                                   value = 0, min = 0, max = 1,
                                                   step = 0.025, width = "70%"),
                                      )
                                  |> prompter::add_prompt(message = 
                                                              paste("A proportion between 0 and 1. ", 
                                                                    "Observational noise (e.g., genotyping error) ",
                                                                    "for all models. ",
                                                                    "Added during sampling, ",
                                                                    "after predictions from model ",
                                                                    "have been obtained; ",
                                                                    "predicted probabilities are not affected.",
                                                                    "If larger than 0, this proportion of entries ",
                                                                    "in the sampled matrix will be flipped ",
                                                                    "(i.e., 0s turned to 1s and 1s turned to 0s)."
                                                                    ),
                                                          position = "right",
                                                          rounded = TRUE,
                                                          bounce = TRUE,
                                                          size = "medium"
                                                          ),
                                  tags$h5(HTML("<br/>")),
                                  actionButton("resample_mhn", "Generate data from MHN model"),
                                  actionButton("clear_mhn",
                                               HTML("Reset log-&Theta; matrix and delete genotype data"))
                                  |> prompter::add_prompt(message = 
                                                              HTML("Resetting the log-&Theta; matrix will replace the ",
                                                                   "contents of the named object by ",
                                                                   "those of the default one ",
                                                                   "(a three-gene matrix filled with 0s)."),
                                                          position = "right",
                                                          rounded = TRUE,
                                                          bounce = TRUE,
                                                          size = "medium"
                                                         )
                                  ## An example of shinyBS::bsTooltip 
                                  ## shinyBS::bsTooltip("clear_mhn",
                                  ##                    HTML("Resetting the log-&Theta; matrix will replace the ",
                                  ##                         "contents of the named object by ",
                                  ##                         "those of the default one ",
                                  ##                         "(a three-gene matrix filled with 0s)."),
                                  ##                    "right", options = list(container = "body")
                                  ##                    )
                                  ## Prompter does not render the "log-&Theta"
                                 
                              )
                     }
                 )
        } else if (input$input2build == "upload") {
            tags$div(## class = "frame",
                     tags$h3("Upload data (CSV format)"),
                     tags$h5(HTML("If you want to give your data a specific ",
                                  "name, set it in the box below ",
                                  "before uploading the data. "
                                  )),
                     tags$div(class = "inlin3",
                              textInput(inputId = "name_uploaded",
                                        label = "Name for data",
                                        value = "Uploaded_data"
                                        )
                              ),
                     tags$h5(paste0("Format: csv ---comma separated values---,",
                                    " with first row with gene names."
                                    )),
                     tags$h5(HTML("Use only alphanumeric characters ",
                                  "for gene names, and do not start ",
                                  "a gene name with a number; ",
                                  "keep gene names short (for figures). ",
                                  "Use 0 or 1 for ",
                                  "altered/not-altered (mutated/not-mutated)."                  
                                  )),
                     tags$div(class = "upload_file",
                              fileInput("csd", "Load Data",
                                        multiple = FALSE,
                                        accept = c(
                                            "text/csv",
                                            ".csv"))),
                     tags$h5(HTML("<br/>")),
                     )
        }
    })

    ## ## Without data modification for upload.
    ## output$change_counts <- renderUI({
    ##     if (input$input2build %in% c("csd", "dag", "matrix")) {
    ##         tags$div(class = "frame",
    ##                  tags$div(class = "flex",
    ##                           tags$h3("3. Change genotype's counts"),
    ##                           actionButton("display_help_change_genotype_counts", "Help"),
    ##                           tags$h3(HTML("<br/>")),
    ##                           ),
    ##                  tags$div(id = "csd_table",
    ##                           DT::DTOutput("csd_counts")
    ##                           ),
    ##                  tags$h5(HTML("<br/>")),
    ##                  if (input$input2build %in% c("csd"))
    ##                      actionButton("clear_genotype", "Delete all genotype data")
    ##                  else if (input$input2build %in% c("matrix"))
    ##                      tags$h5(HTML("To delete genotype data, use",
    ##                                   "'Reset log-&Theta; matrix and delete genotype data'",
    ##                                   "above."))
    ##                  else if (input$input2build %in% c("dag"))
    ##                      tags$h5(HTML("To delete genotype data, use",
    ##                                   "'Reset DAG and delete genotype data'",
    ##                                   "above."))
    
    ##                  )
    ##     }
    ## })

    ## With data modification for upload
    output$change_counts <- renderUI({
        if (input$input2build %in% c("upload", "csd", "dag", "matrix")) {
            menu_num <- ifelse(input$input2build == "upload", "2", "3")
            tags$div(class = "frame",
                     tags$div(class = "flex",
                              ## tags$h3(paste0(menu_num, " . Change genotype's counts")),
                              tags$h3("Change genotype's counts"),
                              actionButton("display_help_change_genotype_counts",
                                           "Help", class = "btn-info"),
                              tags$h3(HTML("<br/>")),
                              ),
                     tags$div(id = "csd_table",
                              DT::DTOutput("csd_counts")
                              ),
                     tags$h5(HTML("<br/>")),
                     ## We could include upload and dag here, but it makes no sense
                     if (input$input2build %in% c("csd"))
                         actionButton("clear_genotype", "Delete all genotype data")
                     else if (input$input2build %in% c("matrix"))
                         tags$h5(HTML("To delete all genotype data, use",
                                      "'Reset log-&Theta; matrix and delete genotype data'",
                                      "above."))
                     else if (input$input2build %in% c("dag"))
                         tags$h5(HTML("To delete all genotype data, use ",
                                      "'Reset DAG and delete genotype data'",
                                      "above."))
                     else if (input$input2build %in% c("upload"))
                         tags$h5(HTML("To delete (or reset) all genotype data",
                                      "upload a new (or the same) data file."))
                     )
        }
    })


    ## DAG builder
    ## Controling dag builder
    dag_data <- reactive({
        if (input$input2build == "dag") {
            mymessage("At dag_data reactive call")
            input$dag_model
            input$dag_table_cell_edit
            all_gene_names <- c("Root", data$gene_names)
            edges <- which(data$dag == 1, arr.ind = TRUE)
            tmp_dag_parent_set <- data$dag_parent_set
            x <- length(tmp_dag_parent_set)
            ## The code below could fail (without sever consequences)
            ## if we have moved back and forth between DAG and upload, for example
            ## as we have not yet update the DAG data. The plot will ask for
            ## the dag data, but those are not available.
            ## Be more elegant
            ## FIXME: what should really happen is that the work that happens
            ## in "At observeEvent toListen" happened before this
            ## But this is part of the UI and the other is listening. I guess
            ## this must come first?
            
            dag_dataset_names <- unlist(lapply(datasets$all_csd$dag, function(x) x$name))

            if (!(data$name %in% dag_dataset_names)) {
                mymessage("    data$name not in dag_dataset_names. Returning a NULL")
                return(NULL)
            }
            
            ## I have to this weird thing because using data$gene_names does not work
            ## for some unkown reason. Eh??!!! What weird thing?
            names(tmp_dag_parent_set) <- all_gene_names[seq(2, x + 1)]
            dag_data <- data.frame(From = all_gene_names[edges[, "row"]]
                                 , To = all_gene_names[edges[, "col"]]
                                 , Relation = tmp_dag_parent_set[edges[, "col"] - 1]
                                 , Lambdas = data$lambdas[edges[, "col"] - 1])
            
            if ((default_dag_model %in% c("OT", "OncoBN"))
                & (any(dag_data$Lambdas < 0) | any(dag_data$Lambdas > 1))){
                showModal(dataModal("thetas/probabilities should be between 0 and 1"))
                updateRadioButtons(session, "dag_model", selected = "HESBCN")
            }

            if (default_dag_model %in% c("OT")) {
                colnames(dag_data) <- c("From", "To", "Relation", "Weight")
                dag_data$Relation <- NULL
            } else if (default_dag_model %in% c("OncoBN")) {
                if (length(unique(dag_data$Relation)) > 2) {
                    showModal(dataModal("OncoBN model must only include one type of relationship"))
                    updateRadioButtons(session, "dag_model", selected = "HESBCN")
                }
                colnames(dag_data) <- c("From", "To", "Relation", "theta")
                
            }
            return(dag_data)
        }
    })

    output$dag_table <-
        DT::renderDT(
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
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)
            shinyjs::click("resample_dag")
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
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::to_stnd_csd_dataset(data)
            if (sum(data$dag) == 0) {
                data$csd_counts <- datasets$all_csd[[input$input2build]][[input$select_csd]]$csd_counts
            } else {
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
            colnames(tmp_dag) <- rownames(tmp_dag) <- c("Root", data$gene_names)
            tmp_dag["Root", data$gene_names[1]] <- 1
            data$dag <- tmp_dag
            data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            data$data <- .ev_SHINY_dflt$template_data$data
            data$dag_parent_set <- tmp_data$dag_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::to_stnd_csd_dataset(data)
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
            tmp_data <-
                evamtools:::modify_lambdas_and_parent_set_from_table(dag_data(),
                                                                     info, data$lambdas
                                                                   , data$dag
                                                                   , data$dag_parent_set
                                                                   , dag_model = default_dag_model)
            data$lambdas <- tmp_data$lambdas
            data$dag_parent_set <- tmp_data$parent_set
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)
            shinyjs::click("resample_dag")
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Building trm from dag
    observeEvent(input$resample_dag, {
        tryCatch({
            if (sum(colSums(data$dag) > 0) < 2)
                stop("The must be at least two genes ",
                     "in the DAG.")
            the_dag_data <- dag_data()
            gene_names <- setdiff(unique(c(the_dag_data$From, the_dag_data$To)),
                                  "Root")
            tmp_dag_data <-
                evamtools:::generate_sample_from_dag(the_dag_data
                                                   , data$dag_parent_set[gene_names]
                                                   , noise = input$dag_noise
                                                   , N = input$dag_samples
                                                   , dag_model = default_dag_model
                                                   , epos = input$dag_epos)
                
            data$csd_counts <-
                tmp_dag_data$csd_counts[tmp_dag_data$csd_counts[, 2] > 0, ]
            data$data <- tmp_dag_data$data
            
            datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
            datasets$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
            datasets$all_csd[[input$input2build]][[input$select_csd]]$lambdas <- data$lambdas
            datasets$all_csd[[input$input2build]][[input$select_csd]]$dag_parent_set <- data$dag_parent_set
            shinyjs::enable("analysis")
        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    ## Help for output of downloaded CPM results
    observeEvent(input$how2downloadcpm, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("Downloading CPMs results and analyzed data"),
            tags$div(
                     tags$p(HTML("Format and contents: rds file with ",
                                 "two lists: <ul>")),
                     tags$li(" cpm_output:  ", 
                             "the concatenated output from ",
                             "evam and sample_evam. ",
                             "The analyzed data are in ",
                             " object$cpm_output$analyzed_data ."
                             ),
                     tags$li("tabular_data: the tabular data output.   "),
                     tags$p(HTML("</ul>")),
                     )
        ))
    }
    )
    ## Help for DAG building
    observeEvent(input$how2build_dag, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("How to build a DAG and generate a sample"),
            tags$div(
                     tags$p("0. Select the 'Type of model'. ",
                            "(Whenever you change the type of model, ",
                            "to avoid confusion it is best to click on ",
                            "'Reset DAG and delete genotype data', or ",
                            "you might get errors from asking for impossible ",
                            "settings ---e.g., moving from CBN with ",
                            "multiple parents to OT.)"),
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
                             "it will be assigned as descendant of Root."),
                     tags$p(HTML("</ul>")),
                     tags$p(HTML("2. To <strong>change the value of a lambda</strong> "),
                            "click on the cell, ",
                            "edit the cell's content and press Ctrl+Enter."),
                     tags$p(HTML("3. Set the value of <strong>Relation</strong> "),
                            "to one of 'Single' (single parent), ",
                            "AND, OR, XOR.",
                            "OT only accepts 'Single' as each node ",
                            "has a single parent. ",
                            "OncoBN accepts 'Single', 'AND', 'OR' ",
                            "or combinations of either Single and AND or ",
                            "Single and OR ",
                            "(and terms not among Single, AND, OR ",
                            "will be converted to ORs).",
                            "'CBN/H-ESBCN' models can be specified with ",
                            "AND, OR, XOR, Single, or combinations of the above ",
                            "(and terms not among Single, AND, OR, XOR ",
                            "will be converted to ANDs)[1]. ",
                            "Edit the cell's content and press Ctrl+Enter. ",
                            "All incoming edges to a node must have the same ",
                            "Relation (the program will force this)."),
                     tags$p(HTML("4. Modify, if you want, the <strong>size of the sample</strong> "),
                            "('Number of genotypes to sample') and ",
                            HTML("the <strong>Observational noise</strong> (genotyping error) "),
                            HTML("and click on <strong>'Generate data from DAG'</strong> to generate a sample. ")),
                     tags$p("After the sample is generated for the first time, ",
                            "the sample should be generated again automatically ",
                            "whenever you change the model ",
                            "(add or remove edges, change lambdas, etc) ",
                            "if you hit Ctrl-Enter after you are done editing ",
                            "the DAG table."),
                     tags$h3(HTML("<br/>")),
                     tags$p(HTML("<h5>[1] The models fitted by CBN contain only ANDs, but the ",
                                 "specification of a model for data simulation from H-ESBCN and CBN is ",
                                 "the same, except for the type of relationship; ",
                                 "both H-ESBCN and CBN can contain ANDs, but ",
                                 "only H-ESBCN can contain ORs and XORs. ",
                                 "In both cases, the DAG specification is used to ",
                                 "obtain the transition rate matrix ",
                                 "that respects the encoded AND/OR/XOR.",
                                 "The data are generated from this ",
                                 "transition rate matrix.</h5>"))
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
            tmp_thetas <- .ev_SHINY_dflt$template_data$thetas
            colnames(tmp_thetas) <- rownames(tmp_thetas) <- data$gene_names
            data$thetas <- tmp_thetas
            data$dag <- NULL
            data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            data$data <- .ev_SHINY_dflt$template_data$data
            data$dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
            shinyjs::disable("analysis")
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::to_stnd_csd_dataset(data)
        }, error=function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    observeEvent(input$how2build_matrix, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3(HTML("How to input &theta;s and generate a sample")),
            tags$div(
                     tags$p("0. Select the number of genes with the slider, above."
                            ),
                     tags$p(HTML("1. &Theta;<sub>i,j</sub> ",
                                 "(i.e., <em>e<sup>&theta;<sub>i,j</sub></sup></em>) ",
                                 "is the multiplicative ",
                                 "effect of gene in column <em>j</em> on ",
                                 "gene in row <em>i</em>. ",
                                 "&Theta;<sub>i,i</sub> is the baseline hazard rate ",
                                 "of event <em>i</em>. ",
                                 "The entries you provide in the matrix are ",
                                 "the &theta;<sub>i,j</sub>, ",
                                 "the lower-case theta (not the multiplicative ",
                                 "effects themselves, but their log)."## ,
                                 )),
                     tags$p("2. Double click in a cell to edit it."),
                     tags$p("3. Press Tab to move to the next row."),
                     tags$p("4. Use Ctrl + Enter to save changes. ",
                            HTML("You <strong>must</strong> save the changes.")),
                     tags$p("5. Modify, if you want, the size of the sample ",
                            "('Number of genotypes to sample') and ",
                            HTML("the <strong>Observational noise</strong> (genotyping error) "),
                            HTML("and click on <strong>'Generate data from MHN model'</strong> to generate a sample. "),
                            "The sample is also updated as soon as you save an entry ",
                            "in the matrix or change the number of genes."),
                     tags$p(HTML("You can make sure <b>the &theta;s have been updated</b> "),
                            "by checking the figure of the matrix on the right.")
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
            tmp_dag <- .ev_SHINY_dflt$template_data$dag
            colnames(tmp_dag) <- rownames(tmp_dag) <- c("WT", data$gene_names)
            tmp_dag["WT", data$gene_names[1]] <- 1
                                        # data$dag <- NULL
            data$dag <- tmp_dag
            data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            data$data <- .ev_SHINY_dflt$template_data$data
            data$dag_parent_set <- .ev_SHINY_dflt$template_data$dag_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$dag_parent_set) <- data$gene_names
            shinyjs::disable("analysis")
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::to_stnd_csd_dataset(data)
        }, error=function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    observeEvent(input$add_genotype, {
        tryCatch({
            genotype <- paste(evamtools:::evam_string_sort(input$genotype), collapse = ", ")
            genotype <- ifelse(genotype == "", "WT", genotype)
            genot_count <- ifelse(is.na(input$genotype_freq), -1,
                                  input$genotype_freq)

            if (genot_count > 0) {
                data$csd_counts[genotype, ] <- c(genotype, genot_count)
                rownames(data$csd_counts) <- data$csd_counts$Genotype
                data$csd_counts[, 2] <- as.numeric(data$csd_counts[, 2])
                ## Filtering out non-positive counts
                data$csd_counts <-
                    evamtools:::reorder_to_standard_order_count_df(
                                    data$csd_counts[data$csd_counts[, 2] > 0, ])
                data$data <-
                    datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
                        evamtools:::genotypeCounts_to_data(data$csd_counts, e = 0)
                
                
                shinyjs::enable("analysis")
            } else {
                showModal(modalDialog(paste("Counts <= 0 present. ",
                                            "They will be removed.")))
            }
            updateNumericInput(session, "genotype_freq", value = NA)

            ## id_change_genotype_muts
            gene_options <- set_gene_names_after_resize(data$data,
                                                        data$gene_names)[1:input$gene_number]
            updateCheckboxGroupInput(session, "genotype", label = "Mutations",
                                     choices = lapply(1:input$gene_number,
                                                      function(i) gene_options[i]),
                                     selected = NULL)
        ## updateCheckboxGroupInput(session, "genotype", label = "Mutations",
            ##                          choices = lapply(1:input$gene_number,
            ##                                           function(i) data$gene_names[i]),
            ##                          selected = NULL)
        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    ## Genotypes table
    output$csd_counts <- DT::renderDT( {
        d1 <- display_freqs()
        if (nrow(d1)) {
            d1 <- data.frame(Index = 1:nrow(d1),
                             d1)
        } else {
            d1 <- data.frame(Index = integer(),
                             Genotype = character(),
                             Counts = integer())
        }
        d1
    }
   ,
    selection = 'none', server = TRUE,
    editable = list(target = "column",
                    disable = list(columns = c(0, 1)))
  , rownames = FALSE,
    options = list(
        columnDefs =
            list(list(
                orderable = TRUE,
                className = 'dt-center',
                targets = "_all")),
        info = FALSE, paginate= FALSE),
    )

    observeEvent(input$csd_counts_cell_edit, {
        tryCatch({
            info <- input$csd_counts_cell_edit
            info[ , "col"] <- 2
            data$csd_counts <- DT::editData(data$csd_counts, info, "csd_counts")

            ## Filtering out non-positive counts
            if (any(data$csd_counts[, 2] < 0))
                        showModal(modalDialog(paste("Counts < 0 present. ",
                                                    "They will be removed.")))
            ## We want to purge 0 entries
            data$csd_counts <- data$csd_counts[data$csd_counts[, 2] > 0, ]
            data$data <-
                datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
                    evamtools:::genotypeCounts_to_data(data$csd_counts, e = 0)
            ##}
        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    ## ## Plot histogram of genotypes
    output$plot <- plotly::renderPlotly({
        tryCatch({
            mymessage("At output$plot")
            evamtools:::plot_genotype_counts_plly(display_freqs())
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Plot dag or MHN of dataset
    output$dag_plot <- renderPlot({
        data2plot <- NULL
        edges <- NULL

        if (input$input2build %in% c("dag")
            && sum(data$dag) > 0
            && !is.null(input$gene_number)
            ) {
            data2plot <- igraph::graph_from_adjacency_matrix(data$dag)
            data2plot <- igraph::decompose(data2plot)[[1]]
            edges <- igraph::as_data_frame(data2plot)
            colnames(edges) <- c("From", "To")
            if(!is.null(data$dag_parent_set)) edges$Relation <- data$dag_parent_set[edges$To]
        } else if (input$input2build %in% c("matrix") 
                   && !is.null(data$thetas)
                   && !is.null(input$gene_number)
                   && length(data$thetas[1:input$gene_number, 1:input$gene_number])>0
                   ) {
            data2plot <- data$thetas[1:input$gene_number, 1:input$gene_number]
        }
        evamtools:::plot_method(data2plot, data$dag_parent_set, edges)
    })

    ## Run CPMs
    observeEvent(input$analysis, {
        ## Calculate TRM for DAG and for matrices

        tryCatch({

            if (input$gene_number >= 7) {
                showModal(
                    dataModal("Beware! You are analyzing data ",
                              "with 7 or more genes. ",
                              "This can take longer than usual ",
                              "and plots may be crowded. ",
                              "We recommend using top_paths options in ",
                              "the Results' tab.",
                              type = "Warning: "))
            }

            if (is.null(input$cpm_methods) ||
                (length(input$cpm_methods) ==1 && is.na(input$cpm_methods)))
                stop("You must use at least one method ",
                     "(check 'CPMs to use' under 'Advanced options ",
                     "and CPMs to use).")

            
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

            cbn_opts <- list(init_poset = input$CBN_init_poset,
                             omp_threads = input$CBN_omp_threads)
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

            data2run <- evamtools:::genotypeCounts_to_data(display_freqs(),
                                                           e = 0)
            if (ncol(data2run) < 1) {
                stop("Your data contains less than one column (genes, events). ",
                     "Do you have a single genotype? Maybe only WT?")
                
            }

            
            progress$inc(1/5, detail = "Setting up data")
            Sys.sleep(0.5)
            progress$inc(2/5, detail = "Running CPMs")

            ## methods <- .ev_SHINY_dflt$cpms2run
            if (!is.null(input$cpm_methods)) {
                methods <- unique(input$cpm_methods)
            }
            
            
            cpm_output <- R.utils::withTimeout({evam(data2run
                                                   , methods = methods
                                                   , paths_max = input$return_paths_max
                                                   , mhn_opts = mhn_opts
                                                   , ot_opts = ot_opts
                                                   , cbn_opts = cbn_opts
                                                   , hesbcn_opts = hesbcn_opts
                                                   , oncobn_opts = oncobn_opts
                                                   , mccbn_opts = mccbn_opts)},
                                               elapsed = EVAM_MAX_ELAPSED, 
                                               timeout = EVAM_MAX_ELAPSED, 
                                               cpu = Inf,
                                               onTimeout = "silent")

            if (is.null(cpm_output)) stop("Error running evam. ",
                                          "Most likely you exceeded maximum ",
                                          "allowed time (EVAM_MAX_ELAPSED).")

            sampled_from_CPMs <- NULL
            do_sampling <- input$do_sampling == "TRUE"
            if (do_sampling) {
                n_samples <- input$sample_size
                if ((is.null(n_samples)) ||
                    (!is.numeric(n_samples)) ||
                    (n_samples < 100)) {
                    n_samples <- .ev_SHINY_dflt$cpm_samples
                }
                progress$inc(3/5, detail = paste("Running ", n_samples, " samples"))
                ## if (input$do_genotype_transitions) {
                ##     ## disabled when removal_note_sogt_1
                ##     sout <- c("sampled")
                ## }
                
                sampled_from_CPMs <-
                    sample_evam(cpm_output, N = n_samples, methods = methods,
                                output = "sampled_genotype_counts",
                                ## ## disabled when removal_note_sogt_1
                                ## if (input$do_genotype_transitions) { 
                                ##     c("sampled_genotype_counts",
                                ##       "obs_genotype_transitions")
                                ##                                } else {
                                ##           "sampled_genotype_counts"
                                ##       },
                                obs_noise = input$sample_noise)
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
                                  , "do_sampling" = do_sampling
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

            updateTabsetPanel(session, "navbar", selected = "result_viewer")
            updateRadioButtons(session, "select_cpm", selected = result_name)

        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
            })

            all_cpm_out <- reactiveValues()

            ## To increase lag in the redrawing of output, change the number
            ## Values less than 500 can break the app. Larger values might break
            ## with complex plots. Deactivate this completely commenting the debounce
            plot2show <- debounce(reactive({
                input$cpm2show
            }),  900)

            ## ## No delay showing plots. 
            ## plot2show <- reactive({
            ##     input$cpm2show
            ## })

            
            output$sims <- renderUI({
                if ((length(names(all_cpm_out)) > 0) && (!is.null(input$select_cpm))) {
                    tmp_data <- all_cpm_out[[input$select_cpm]]$cpm_output

                    number_of_columns <- floor(12 /
                                               ifelse(length(plot2show()) <=0, 1, length(plot2show())))

                    lapply(plot2show(), function(met) {
                        method_data <- evamtools:::process_data(tmp_data, met,
                                                                plot_type = "trans_mat")
                        output[[sprintf("plot_sims_%s", met)]] <- renderPlot({
                            pl <- evamtools:::plot_method(method_data$method_info
                                                        , method_data$parent_set
                                                        , method_data$edges
                                                        , met)
                        })
                        return(
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
            selected_plot_type <- input$data2plot

            number_of_columns <- floor(12 /
                                       ifelse(length(plot2show()) <=0, 1, length(plot2show())))
            if(!(is.null(selected_plot_type))){
                if(selected_plot_type %in% c("trans_mat", "trans_rate_mat")){
                    lapply(plot2show(), function(met){
                        method_data <- evamtools:::process_data(tmp_data, met,
                                                                plot_type = selected_plot_type)
                        output[[sprintf("plot_sims2_%s", met)]] <- renderPlot({
                            ## This occasionally fails for no reason
                            ## so try up to max_plot_tries
                            max_plot_tries <- 5
                            for (i in 1:max_plot_tries) {
                                pl <-
                                    try(evamtools:::plot_genot_fg(method_data$data2plot,
                                                                  ## We use it to define "Observed" and "Not Observed" genotypes
                                                                  observations = tmp_data$original_data, 
                                                                  sampled_counts = method_data$sampled_genotype_counts,
                                                                  top_paths = input$freq2label,
                                                                  label_type = input$label2plot,
                                                                  plot_type = selected_plot_type))
                                if (!inherits(pl, "try-error")) break
                                if (i >= max_plot_tries) stop("Could not produce requested plot ",
                                                              "after ", max_plot_tries, " attempts. ",
                                                              "A problem with the data? ",
                                                              "Can any of the other plot types be produced?")
                            }

                        })
                        return(
                            column(number_of_columns,
                                   plotOutput(sprintf("plot_sims2_%s", met)))
                        )
                    })
                } else if(selected_plot_type %in% c("predicted_genotype_freqs", "sampled_genotype_counts")){
                    lapply(plot2show(), function(met){
                        method_data <- evamtools:::process_data(tmp_data, met,
                                                                plot_type = selected_plot_type)$data2plot

                        if(selected_plot_type %in% c("predicted_genotype_freqs")){
                            data2plot <- data.frame("Genotype" = names(method_data), 
                                                    "Freq" = as.vector(method_data))

                        }
                        if(selected_plot_type %in% c("sampled_genotype_counts")){
                            data2plot <- data.frame("Genotype" = names(method_data), 
                                                    "Counts" = as.vector(method_data))
                        }
                        output[[sprintf("plot_sims2_%s", met)]] <- renderPlot(
                            pl <- evamtools:::plot_genotype_counts(data2plot)
                        )
                        return(
                            column(number_of_columns,
                                   plotOutput(sprintf("plot_sims2_%s", met)))
                        )
                        ## Trying to use plotly. Does not work. I give up
                        ## Anyway, not really needed: would be slower, and the
                        ## numbers can be easily seen in the table.
                        ## output[[sprintf("plot_sims2_%s", met)]] <-
                        ##     plotly::renderPlotly(
                        ##                 pl <- evamtools:::plot_genotype_counts_plly(data2plot)
                        ##             )
                        ## return(
                        ##     column(number_of_columns,
                        ##            plotly::plotlyOutput(sprintf("plot_sims2_%s", met)))
                        ## )

                    })
                } 
            } else {
                ## Disabling donwload button
                shinyjs::disable(selector = "#download_cpm")

                return(tags$h3("There are not results to show yet. ",
                               "Go to the input tab, select ",
                               "the data to analyze and ",
                               "hit the 'Run evamtools' button."))
            }
        }
            })

    ## ## Go back to input to work again with the data
    ## observeEvent(input$modify_data, {
    ##     tryCatch({
    
    ##         if (length(all_cpm_out) > 0){
    ##             tmp_data <- all_cpm_out[[input$select_cpm]]$orig_data
    ##             dataset_name <- strsplit(input$select_cpm, "__")[[1]][[1]]
    ##             dataset_type <- tmp_data$type
    ##             last_visited_pages[[tmp_data$type]] <<- dataset_name

    ##             tmp_data <- datasets$all_csd[[tmp_data$type]][[dataset_name]] <-
    ##                 evamtools:::to_stnd_csd_dataset(tmp_data)

    ##             data <- tmp_data
    ##             data$csd_counts <- evamtools:::get_csd(tmp_data$data)
    ##             data$n_genes <- ncol(data$data)
    ##             updateNumericInput(session, "gene_number", value = data$n_genes)
    ##             updateTabsetPanel(session, "navbar",
    ##                               selected = "csd_builder")
    ##             updateRadioButtons(session, "input2build", selected = dataset_type)
    ##             updateRadioButtons(session, "select_csd", selected = dataset_name)
    ##         }
    ##     }, error = function(e){
    ##         showModal(dataModal(e[[1]]))
    ##     })
    ## })


    output$customize <- renderUI({
        do_sampling <- tryCatch({
            sampling <- ifelse(
                is.null(all_cpm_out[[input$select_cpm]]$do_sampling), FALSE, 
                all_cpm_out[[input$select_cpm]]$do_sampling)
            sampling
        }, error = function(e){
            return(FALSE)
        })


        tagList(
            tags$div(class = "frame",
                     tags$h3("Customize the visualization"),
                     tags$div(class = "inline",
                              checkboxGroupInput(inputId = "cpm2show",
                                                 label = HTML("CPMs to show"),
                                                 ## "<h5>(Some or all of those used to analyze the data; ",
                                                 ##              "use 'Modify data' ---below--- to go back ",
                                                 ##              "and click on 'Advanced options' if you",
                                                 ##              "want to use other methods)</h5>"), 
                                                 choices = input$cpm_methods,
                                                 ## c("OT", "OncoBN", "CBN", "MHN", "HESBCN", "MCCBN"),
                                                 selected = input$cpm_methods)
                              |> prompter::add_prompt(message = 
                                                          HTML("Show graphical output of the CPMs used to analyze the data.  "
                                                             , "Go back to \"User input\" "
                                                             , "and click on \"Advanced options\" if you"
                                                             , "want to use other methods."
                                                               ),
                                                      position = "right",
                                                      rounded = TRUE,
                                                      bounce = TRUE,
                                                      size = "medium"
                                                      ),
                              tags$h4(HTML("<hr style=\"height:1px; width:80%; background-color:black;text-align:left\">")),
                              tags$h4(HTML("<br/>")),
                              tags$div(class = "inline",
                                       radioButtons(inputId = "data2plot",
                                                    label = HTML("Predictions from models to display"
                                                                 ## , "<h5>(This output is also displayed in tabular form on the bottom right.)</h5>"
                                                                 ),
                                                    choiceNames = 
                                                        if (do_sampling) {
                                                            c("Transition probabilities", 
                                                              "Transition rates",
                                                              "Predicted genotype relative frequencies",
                                                              "Sampled genotype counts")
                                                        } else {         
                                                            c("Transition probabilities", 
                                                              "Transition rates",
                                                              "Predicted genotype relative frequencies")
                                                        }
                                                   ,
                                                    choiceValues = 
                                                        if (do_sampling) {
                                                            c("trans_mat", 
                                                              "trans_rate_mat",
                                                              "predicted_genotype_freqs",
                                                              "sampled_genotype_counts")
                                                        } else {
                                                            c("trans_mat", 
                                                              "trans_rate_mat",
                                                              "predicted_genotype_freqs")
                                                        }
                                                   ,
                                                    selected = "trans_mat"
                                                    )
                                       |> prompter::add_prompt(message = 
                                                                   HTML("This output is also displayed in tabular form on the bottom right."
                                                                        ),
                                                               position = "right",
                                                               rounded = TRUE,
                                                               bounce = TRUE,
                                                               size = "medium"
                                                               ),
                                       ),
                              tags$h4(HTML("<hr style=\"height:1px; width:80%; background-color:black;text-align:left\">")),
                              tags$h4(HTML("<br/>")),
                              tags$div(class = "inline",
                                       radioButtons(inputId = "label2plot",
                                                    label = HTML("Type of label"), ## <h5>(for transition [rate] plots).</h5>"),
                                                    choiceNames =  c("Genotype", "Last gene mutated"),
                                                    choiceValues = c("genotype", "acquisition"),
                                                    selected = "genotype"
                                                    )
                                       |> prompter::add_prompt(message = 
                                                                   HTML("Type of label for transition [rate] plots."
                                                                        ),
                                                               position = "right",
                                                               rounded = TRUE,
                                                               bounce = TRUE,
                                                               size = "medium"
                                                               ),
                                       ),
                              ),
                     tags$h4(HTML("<hr style=\"height:1px; width:80%; background-color:black;text-align:left\">")),
                     tags$h4(HTML("<br/>")),
                     ## tags$p(HTML("<strong>Number of most relevant paths to show</strong> "),
                     ##        ## "(set it to 0 to show all paths or ",
                     ##        ## "all genotype labels):"
                     ##        ),
                     tags$div(id="freq2label-wrap",
                              tags$p(HTML("<strong>Number of most relevant paths to show</strong> ")),
                              sliderInput("freq2label",
                                          "",
                                          width = "500px",
                                          value = 5, max = 10, min = 0, step = 1)
                              |> prompter::add_prompt(message = 
                                                          HTML("Set it to 0 to show all paths or all genotype labels"
                                                               ),
                                                      position = "right",
                                                      rounded = TRUE,
                                                      bounce = TRUE,
                                                      size = "medium"                                                 ),
                              )
                     )
        )
    })

    output$cpm_list <- renderUI({
        all_names <- c()
        for (i in names(all_cpm_out)) {
            all_names <- c(all_names, all_cpm_out[[i]]$orig_data$name)
        }

        if ((length(all_names) > 0) && (last_visited_cpm != "")) {
            selected <- names(all_cpm_out)
            
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

    output$original_data <- renderUI({
        ## To see if I disable original data        
        if (length(names(all_cpm_out)) > 0) {
            tags$div(class="frame max_height", tags$h3("Original data"),
                     plotly::renderPlotly(
                                 evamtools:::plot_genotype_counts_plly(
                                                 evamtools:::get_csd(all_cpm_out[[input$select_cpm]]$cpm_output$analyzed_data))
                             ),
                     ## Disable because it led to inconsistent behavior such as this
                     ## Sample from a DAG, say Fork, and analyze
                     ## Go back, and sample, but now sample only 10 cases. Analyze
                     ## Go to first analysis, click on "Modify data" and .. you
                     ## are shown the 10 samples.
                     ## Even worse, you can modify the model in between.
                     
                     ## tags$div(class = "download_button",
                     ##          actionButton("modify_data", "Modify data")
                     ##          )
                     )
        }
    })
    
    output$cpm_freqs <- DT::renderDT(all_cpm_out[[input$select_cpm]]$tabular_data[[input$data2plot]],
                                     selection = 'none', server = TRUE
                                   , rownames = FALSE
                                   , options = list(
                                         columnDefs = list(list(className = 'dt-center', targets = "_all")), info = FALSE, paginate= FALSE)
                                     )

    output$tabular_data <- renderUI({
        if (length(names(all_cpm_out)) > 0) {
            tags$div(class="frame max_height",
                     id = "table_out1",
                     tags$div(class=" max_height",
                              ## A hack to get the tooltip to only show on hover
                              ## over heading of table
                              id = "table_out3", 
                              tags$h3(paste("Tabular output of predictions from models: ",
                                            switch(ifelse(is.null(input$data2plot),
                                                          "not_valid_or_not_yet_existent",
                                                          input$data2plot),
                                                   "trans_mat" = "Transition probabilities.",
                                                   "trans_rate_mat" = "Transition rates.",
                                                   "predicted_genotype_freqs" = "Predicted genotype relative frequencies.",
                                                   "sampled_genotype_counts" = "Sampled genotype counts.",
                                                   "Not a valid input$data2plot"
                                                   )))
                              ) , ## prompter does not work here
                     shinyBS::bsTooltip("table_out3",
                                        ## message = 
                                        HTML("This output is also displayed as the second row of figures. "
                                           , "Choose the output to display from the left radio buttons "
                                           , "\"Predictions from models to display\"."
                                             ),
                                        "bottom", options = list(container = "body")
                                        ## position = "bottom",
                                        ## rounded = TRUE,
                                        ## bounce = TRUE,
                                        ## size = "medium"
                                        ),
                     ## tags$h4("(This output is also displayed as the second row of figures. ",
            ##         "Choose the output to display from the left radio buttons ",
            ##         "'Predictions from models to display')"),
            tags$div(id = "table_out2",
                     DT::DTOutput("cpm_freqs"),
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
}



## What must be checked when dealing with the reactive
## MHN
##  - changing data set resamples and plots
##  - changing number of genes resamples and plots

## DAG:
##  - as for MHN
##  - as for MHN
##  - Can go from linear to User to Linear to User and same with AND_OR_XOR
##    without error message of number of genes.


## Both:
##   -plot in the very first data (before any changes)


## Things that could avoid repeating, but I give up
## With DAGs, if we do a resample trigger, no need for a gene number trigger.
##  but it is cheap, since no resampling when we do number trigger.
##  Problem is it creates another figure.

## In many DAGs, calling the plot function twice or even thrice.
