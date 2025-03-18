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


dataModal <- function(error_message, type = "Error: ") {
    if (type == "Error: ") {
        type <- HTML("<font color ='red'>", type, "</font color>")
    }

    error_message_htmlized <- gsub("\n", "<br>", error_message)
    modalDialog(
        easyClose = TRUE,
        title = tags$h3(type),
        tags$div(HTML(error_message_htmlized))
    )
}

options(warnPartialMatchDollar = TRUE)





######################################################################
######
######  In lieu of proper refactoring, this is some info
######  about some atrocious names
######
######################################################################

## data$data:  the subjects by genes matrix (or data frame) of 0/1
## csd_counts: the Genotype and Counts data frame
## data$dag:   the DAG as an adjacency matrix
## dag_data:   object and functions for the dag model as a data frame with
##              From, To, Weights/Lambdas/etc
##



## I have left a bunch of messages. To make it easier to dis/enable them
## by turning them to a no-op (https://stackoverflow.com/a/10933721)
## mymessage <- function(...) message(...)
mymessage <- function(...) invisible(NULL)




sanity_file_name <- function(x) {
    if (gsub(" ", "", x, fixed = TRUE) == "") {
        stop("Name of data or file name cannot be an empty string.")
    }
    if (grepl(" ", x, fixed = TRUE)) {
        stop("Name of data or file name cannot contain spaces.")
    }

    gn_space <- !(stringi::stri_count_regex(x, "^[a-zA-Z].*"))
    if (any(gn_space))
        stop("All file and data names should start with a letter. ",
             "Yours don't; that is not allowed.")

    gn_space <- stringi::stri_count_regex(x,  "[^a-zA-z0-9_-]+")
    
    if (any(gn_space))
        stop("Use only letters, numbers, the hyphen - and the underscore _ .",
             "The name you provided contains other characters.")
}


sanity_new_gene_names <- function(x) {
    ## From similar code in evam-wrapper
    gn_backslash <- stringi::stri_count_fixed(x, "\\") 
    if (any(gn_backslash))
        stop("At least one of your new gene names has a backslash. That is not allowed.")

    gn_space <- stringi::stri_count_regex(x, "[\\s]") 
    if (any(gn_space))
        stop("At least one of your new gene names has a space. That is not allowed.")

    gn_space <- stringi::stri_count_regex(x, "[-]") 
    if (any(gn_space))
        stop("At least one of your new gene names has a hyphen. That is not allowed.")

    gn_space <- stringi::stri_count_regex(x, "^\\d") 
    if (any(gn_space))
        stop("At least one of your new gene names starts with a number. That is not allowed.")

    gn_space <- !(stringi::stri_count_regex(x, "^[a-zA-Z].*"))
    if (any(gn_space))
        stop("All gene names should start with a letter. Yours don't; that is not allowed.")

    gn_space <- stringi::stri_count_regex(x,  "[^a-zA-z0-9_]+")
    if (any(gn_space))
        stop("Use only letters, numbers, and the underscore _ .",
             "The gene names you provided contain other characters.")
    
    if (any(x == "WT"))
        stop("One of your new gene names is called WT. That is not allowed.")
}
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
    } else {
        return(FALSE)
    }
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
    number_of_genes <- ifelse(is.null(data$data),
                              0, sum(colSums(data$data) > 0))
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
    gene_names_in_freqs <- evamtools:::evam_string_sort(colnames(x))
    if (length(gene_names_in_freqs) == 0) {
        ## Only WT
        return(evamtools:::evam_string_sort(gene_names))
    } else {
        gene_names_wo_current <-
            evamtools:::evam_string_sort(setdiff(gene_names, gene_names_in_freqs))
        gene_names <-
            c(gene_names_in_freqs, gene_names_wo_current)[1:gene_names_num]
        return(gene_names)
    }
}

server <- function(input, output, session, EVAM_MAX_ELAPSED = 1.5 * 60 * 60) {
    require(evamtools)
    require(tippy)
    ## require(shinyBS)
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

    ## FIXME: a lot of this defaults logic sucks: values depend on 
    ## the object .ev_SHINY_dflt, which would need to be regenerated
    ## separately, which is a PITA. Those defaults should be here.
    ## I make these deps explicit. Needed for shinytests
    data("examples_csd", package = "evamtools")
    data("SHINY_DEFAULTS", package = "evamtools")
    ## We subset, as lots of the examples are not really worth it for the web app.
    examples_csd$csd <- examples_csd$csd[1:5]
    ## The next is (one of) the functions from hell.
    ## And the "upload" component disappears; it was an empty list anyway.
    all_csd_data <- evamtools:::to_stnd_csd_all_datasets(examples_csd)
    min_genes <- .ev_SHINY_dflt$min_genes
    max_genes <- .ev_SHINY_dflt$max_genes


    default_dag_model <- .ev_SHINY_dflt$dag_model
    default_number_genes <- 4
    max_allowed_num_samples <- 100000

    default_num_samples <- 100
    default_obs_noise <- 0
    default_epos <- 0
    
    ## We could hack this, and use a global variable, as we did in the past,
    ## and assign via "<<-" but this is arguably a reactive value:
    ## we do/did last_visited_pages["upload"] <<- dataset_name
    ## where dataset_name is input$name_uploaded
    ## According to my understanding of
    ## https://shiny.rstudio.com/articles/scoping.html
    ## these should be reactive.
    ## If the variable is defined inside this server block,
    ## no bad consequences: it is not changed in all of the
    ## same user's session.
    ## It would if the variable was defined outside the block.

    ## This: https://stackoverflow.com/questions/20333399/are-there-global-variables-in-r-shiny
    ## does not give me much more insight.

    ## But we do not really use the reactivity functionality anyway. So...
    ## might as well just use a variable, and use "<<-"
    reactive_last_visited_pages <- list( 
        upload = "Empty",
        csd    = "Empty",
        dag    = "DAG_Fork_4",
        matrix = "MHN_all_0"
    )

    reactive_last_visited_cpm <- list( 
        ## As above. We assign to last_visited_cpm <- result_name
        ## where result_name is a function of an input, input$select_csd
        the_last_visited_cpm = ""
    )

    ## This is "system-wide"
    ## and different from "data$this_d_dag_model"
    the_dag_model <- list( ## reactiveValues(
        stored_dag_model = default_dag_model
    )

    datasets <- reactiveValues(
        all_csd = all_csd_data
    )

    ## I've not had a problem using identical ids, etc
    ## for DAG and MHN. But just in case
    generate_data <- reactiveValues(
        dag_num_samples = default_num_samples,
        dag_obs_noise   = default_obs_noise,
        mhn_num_samples = default_num_samples,
        mhn_obs_noise   = default_obs_noise,
        epos            = default_epos 
    )

    resample_trigger_from_data_change <- function() {
        if (nrow(data$csd_counts) || !is.null(data$data)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
    
    data <- reactiveValues(
        csd_counts = .ev_SHINY_dflt$template_data$csd_counts
      , data = .ev_SHINY_dflt$template_data$data
      , dag = .ev_SHINY_dflt$template_data$dag
      , DAG_parent_set = .ev_SHINY_dflt$template_data$DAG_parent_set
      , lambdas = .ev_SHINY_dflt$template_data$lambdas
      , thetas = .ev_SHINY_dflt$template_data$thetas
      , n_genes = default_number_genes ## .ev_SHINY_dflt$ngenes
      , gene_names = LETTERS[1: max_genes]
      , name = NULL
      , this_d_dag_model = default_dag_model
    )

    ## The value, true or false, does not really matter; only matters that it
    ## changes.
    changed_dag_model <- reactiveValues(
        invalidate = FALSE
    )


    
    ## Things that could avoid repeating, but I give up
    ## With DAGs, if we do a resample trigger, no need for a gene number trigger.
    ##  but it is cheap, since no resampling when we do number trigger.
    ##  Problem is it creates another figure.
    ## In many DAGs, calling the plot function twice or even thrice.
    
    
    display_freqs <- reactive({
        ## Remember this is called whenever changes in many places
        ## happen.
        
        ## This is often called when there is no need for it. So when you change
        ## the type of data entering (go from MHN to upload, for example) this is
        ## called again and returns the data but for nothing since what we will
        ## want to display are the new data. We now return a
        ## 0-rows data frame when nothing should be returned to make it explicit.

        ## For paranoia, we return always things in standard order. Probably
        ## unnecessary, but just in case.
        
        mymessage("At display_freqs")

        ## provide_gene_names is being enabled somewhere I can't locate
        ## So make sure we catch it right on the redisplay
        ## There are a bunch of calls like this.
        if ((!is.null(data$data) ||
             (nrow(data$csd_counts) > 0))) {
            mymessage("    disabled provide_gene_names under display_freqs")
            shinyjs::disable("provide_gene_names")
        }  

        ## If no data to display, return empty data frame
        thisd <- input$input2build

        if (is.null(data$name) ) {
            mymessage("       NULL data ",
                      "Returning a 0-rows data frame")
            return(data.frame(Genotype = character(), Counts = integer()))
        } else {
            thisd_dataset_names <- unlist(lapply(datasets$all_csd[[thisd]],
                                                 function(x) x$name))
            if (!(data$name %in% thisd_dataset_names)) {
                mymessage("       data$name not in ", thisd_dataset_names, ". ",
                          "Returning a 0-rows data frame")
                return(data.frame(Genotype = character(), Counts = integer()))
            }
        }
        return(
            evamtools:::reorder_to_standard_order_count_df(
                            data$csd_counts[data$csd_counts$Counts > 0, ,
                                            drop = FALSE]))
    })


    observeEvent(input$gene_number, {
        ## id: here_we_change_gene_number
        mymessage("at gene_number trigger")
        
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
            ## b) Makes behavior consistent with that of DAGs.  And we simplify a
            ## lot the display_freqs code: A21_gnn_numfix (also in shiny-utils.R)
            ## BEWARE: if we were not to enforce the next lines, we would need to
            ## go back to using the commented code in A21_gnn_numfix, in
            ## display_freqs.
            if (csd_more_genes_than_set_genes(input, data, session)) {
                csd_message_more_genes_than_set_genes()
            }
        } else if (input$input2build == "dag") {
            the_dag_data <- dag_data()
            if (dag_more_genes_than_set_genes(input, the_dag_data, session)) {
                dag_message_more_genes_than_set_genes()
            }
        } else if (input$input2build == "matrix") {
            if (resample_trigger_from_data_change()) shinyjs::click("resample_mhn")
        }
    }
    ## And now, display_freqs will likely be called
    )


    ## Upload data
    observeEvent(input$csd, {
        if (grepl(".csv", input$csd$datapath)) {
            tryCatch({
                dataset_name <- input$name_uploaded
                sanity_file_name(dataset_name)
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
                tmp_data$data <- try(read.csv(input$csd$datapath),
                                     silent = TRUE)
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

                tmp_data$n_genes <- ncol(tmp_data$data)
                datasets$all_csd[["upload"]][[dataset_name]] <- tmp_data
                
                ## last_visited_pages["upload"] <<- dataset_name
                reactive_last_visited_pages$upload <<- dataset_name
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
                           selected = reactive_last_visited_pages[[input$input2build]])
    })

    ## Define dataset name
    output$dataset_name <- renderUI({
        if(input$input2build %in% c("upload", "csd", "dag", "matrix")){
            tags$div(class = "frame inlin2",
                     tags$h3("Rename the data"),
                     tags$h5(HTML("Give the (modified) data a different name ",
                                  "that will also be used to save the CPM ",
                                  "output.",
                                  "Names should start with a letter, ",
                                  "and can contain only letters, numbers, ",
                                  "hyphen, and underscore, but no ",
                                  "other characters (no periods, spaces, etc.)"
                                  )),
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
                     tags$div(class = "flex",
                              tags$h3(HTML("Download the data")),
                              actionButton("how2downloaddata", "Help", class = "btn-info"),
                              ),
                     tags$div(class = "download_button",
                              ## tags$h5(HTML("Contents of saved file: ",
                              ##              "the data as data frame; ",
                              ##              "if you built a DAG or MHN model, ",
                              ##              "also the model built."
                              ##              )),  
                              downloadButton("download_csd", "Download your data")
                              )
                 )
        }
    })


    ## Saving dataset
    observeEvent(input$save_csd_data, {
        tryCatch({
            sanity_file_name(input$dataset_name)
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
                stop("When no data have been uploaded, ",
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
              , DAG_parent_set = data$DAG_parent_set
              , lambdas = data$lambdas
              , thetas = data$thetas
              , trm = data$trm
              , n_genes = n_genes
              , name = input$dataset_name
              , this_d_dag_model = data$this_d_dag_model)

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
            if (input$input2build %in%  c("csd", "upload")) {
                saveRDS(tmp_data$data, file=file)
            } else if(input$input2build == "dag") {
                gene_names <- setdiff(unique(c(dag_data()$From, dag_data()$To)),
                                      "Root")
                number_of_genes <- length(gene_names)
                stopifnot(number_of_genes == input$gene_number)
                data2save <- list(
                    data = tmp_data$data[, 1:number_of_genes]
                  , model_edges = dag_data()
                  , model = data$this_d_dag_model
                  , DAG_parent_set = tmp_data$DAG_parent_set[1:number_of_genes]
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
            reactive_last_visited_pages[[input$input2build]] <<- input$select_csd
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Display List of availabe CSD
    output$csd_list <- renderUI({
        if (length(names(datasets$all_csd[[input$input2build]]))>0){
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
                    ## selected = last_visited_pages[[input$input2build]],
                    selected = reactive_last_visited_pages[[input$input2build]],
                    choiceNames = all_choice_names,
                    choiceValues = all_names
                )
            )
        } else {
            tags$p("Empty until you upload a data file.")
        }
    })

    ## This is probably abusing observeEvent? And mixing what would be
    ## better served with eventReactive?
    ## But we are only using it for its side effects, not assigning,
    ## so observeEvent seems right.

    ## Notes about the logic, and the update of display_freqs.  display_freqs
    ## needs to be updated whenever we make changes that need to be replotted. It
    ## is called from output$plot. That is much more often than the events below
    ## need to be watched. That is also why the following happens:
    ## - you are at MHN
    ## - you change to upload
    ##     The first thing that gets called is display_freqs (as output$plot is called)
    ##     Then, the block below.
    ## I've left some messages (commented now), so that one can see what is happening.

    toListen <- reactive({
        list(input$input2build, input$select_csd)
    })
    
    observeEvent(toListen(), {
        tryCatch({
            
            mymessage("At observeEvent toListen")

            ## Cleaning stuff
            ## selected <- last_visited_pages[[input$input2build]]
            selected <- reactive_last_visited_pages[[input$input2build]]
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

            mymessage("At observeEvent toListen, b")
            mymessage("    disabled provide_gene_names under toListen")
            shinyjs::disable("provide_gene_names")
            
            data$dag <- tmp_data$dag
            data$DAG_parent_set <- tmp_data$DAG_parent_set
            data$lambdas <- tmp_data$lambdas
            data$thetas <- tmp_data$thetas
            data$name <- tmp_data$name
            data$n_genes <- tmp_data$n_genes
            data$this_d_dag_model <- tmp_data$this_d_dag_model

            if (input$input2build == "dag") {
                number_of_parents <- colSums(data$dag)
                to_keep <- sum(number_of_parents > 0)
                n_genes <- ifelse(to_keep < 1, default_number_genes, to_keep)
            } else if (input$input2build == "matrix") {
                n_genes <- data$n_genes
                if (is.null(n_genes)) {
                    n_genes <- default_number_genes
                } 
            } else if (input$input2build %in% c("csd", "upload")) {
                if (!is.null(data$data))  {
                    n_genes <- ncol(data$data)  
                } else { ## data$data is null
                    if (input$input2build == "csd") {
                        n_genes <- default_number_genes
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
                                                            data$gene_names)[1:max(2, n_genes)]

                updateCheckboxGroupInput(session, "genotype", label = "Mutations",
                                         choices = lapply(1:(max(2, n_genes)),
                                                          function(i) gene_options[i]),
                                         selected = NULL)
                ## Where we "Add genotypes" manually. This sets the count
                updateNumericInput(session, "genotype_freq", value = NA)
            }

            if (is.null(data$data) || (nrow(data$csd_counts) == 0)) {
                mymessage("    enabled provide_gene_names at end of toListen")
                shinyjs::enable("provide_gene_names")
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
                     tags$p("Genotypes are always shown with gene names ",
                            "sorted alphabetically. "),
                     tags$p("Genotypes with count 0 are removed from the table. ",
                            "Thus, if you remove a genotype when editing ",
                            "genotype's counts in the DAG, MHN, or Upload data ",
                            "entries, you will need to regenerate the data ",
                            "to be able to modify those genotypes again."),
                     tags$p("The first column, 'Index' allows us to sort ",
                            "by 'standard order': by number of mutations first ",
                            "and then by alphabetical order of ",
                            "genotype names, where genotypes themselves ",
                            "have genes sorted alphabetically ",
                            "(actually, by what is often called 'natural order').")
                 )
        )
        )
    })


    observeEvent(input$display_help_upload_data, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("Format for uploading data"),
            tags$div(
                     tags$h5(paste0("Format: csv ---comma separated values---,",
                                    "with subjects in rows and gens in columns. ",
                                    "Use 0 or 1 for ",
                                    "altered/not-altered (mutated/not-mutated). ", 
                                    "The first row should contain ",
                                    " the gene names ",
                                    "and there should be no subject names ",
                                    "(i.e., no column of subject names)."
                                    )),
                     tags$h5(HTML("Use only alphanumeric characters ",
                                  "for gene names, and do not start ",
                                  "a gene name with a number. <br> <br> ",
                                  "Keep gene names short (for figures). <br>"
                                  )),
                     tags$h5(HTML("<br><br>A short <b>example file</b> is available here<br><br> ",
                                  "<a href=\"https://raw.githubusercontent.com/rdiaz02/EvAM-Tools/main/examples_for_upload/tinydata.csv\">",
                                  "https://raw.githubusercontent.com/rdiaz02/EvAM-Tools/main/examples_for_upload/tinydata.csv </a>",
                                  "<br><br>",
                                  "and there are additional examples in that ",
                                  "directory ",
                                  "<a href=\"https://github.com/rdiaz02/EvAM-Tools/tree/main/examples_for_upload\">",
                                  "(https://github.com/rdiaz02/EvAM-Tools/tree/main/examples_for_upload)</a> ---make sure to view the files as raw or ",
                                  "download them and open them with a ",
                                  "text editor in your computer.")),
                     )
        )
        )
    })

    
    
    ## Advanced option for running evamtools
    observeEvent(input$advanced_options, {
        shinyjs::toggle("all_advanced_options")
    })

    ## Define number of genes
    output$gene_number_slider <- renderUI({
        val <- ifelse(is.null(data$n_genes),
                      default_number_genes,
                      data$n_genes)
        ## BEWARE!!! Never, ever, add this without isolate.
        ## Without isolate, it forces the slider
        ## of Number of genes back. And creates a mess.
        ## A simple example belos: even that code below screws things up. 
        ## Yes, this is the gene number slider. 
        if ((!is.null(isolate(data$data)) ||
             (nrow(isolate(data$csd_counts)) > 0))) {
            mymessage("    disabled provide_gene_names renderUI")
            shinyjs::disable("provide_gene_names")
        }
        
        ## This breaks things. 
        ## if ((!is.null(data$data)) ||
        ##     (nrow(data$csd_counts) > 0)) {
        ##     uu <- 3 + 2
        ##     message("    just a message from renderUI, no isolate")
        ##     ## shinyjs::disable("provide_gene_names")
        ## }

        
        if (input$input2build %in% c("csd","dag", "matrix")) {
            tags$div(class = "frame flex",
                     tags$h3("Set the number of genes"),
                     tags$h5("(Using 7 or more genes can lead ",
                             "to very long execution times for some methods ",
                             "and crowded figures.)"),
                     tags$div(class="inlin",
                              tags$h3(HTML("<br/>")),
                              sliderInput("gene_number", "Number of genes",
                                          value = val,
                                          max = max_genes, min = min_genes,
                                          step = 1),
                              ## The action that takes place is
                              ## id: here_we_change_gene_number
                              
                              ),
                     tags$h4(HTML("<br/>")),
                     tags$div(class="inlin",
                              tags$div(id="dummy_for_tooltip1",
                                       actionButton("provide_gene_names", "Use different gene names"),
                                       tippy::tippy_this(element = "dummy_for_tooltip1", ## "provide_gene_names",
                                                         tooltip = paste("<span style='font-size:1.5em; text-align:left;'>",
                                                                         "<p><b>Use different gene names</b></p>",
                                                                         "<p>Create new models/new data using gene names you provide.</p> ",
                                                                         "<br>",
                                                                         "<p><b>Can only be used for models/data that are empty. </b>",
                                                                         "Therefore, you might need to click on ",
                                                                         "'Delete all genotype data'",
                                                                         "'Reset DAG and delete genotype data'",
                                                                         "or 'Reset log-&Theta; matrix and delete genotype data' ",
                                                                         "before being able to use this.</p>",
                                                                         "<br>",
                                                                         "<p><b>IMPORTANT:</b> Do not think about this option as a way ",
                                                                         "to rename genes in existing data or models. ",
                                                                         "That becomes confusing very quickly. ",
                                                                         "Instead, think of this as a way to build ",
                                                                         "<b>new models from scratch</b> ",
                                                                         "with the names that you want.</p><span>"
                                                                         ),
                                                         arrow = TRUE, animation = "shift-toward"
                                                         )
                                       )
                              ),
                     )
        }
    })

    observeEvent(input$provide_gene_names, {
        
        showModal(modalDialog(
            title = tags$h3("Use different gene names"),
            tags$div(class = "inlin2",
                     textInput(inputId = "new_gene_names", "New gene names",
                               value = paste(LETTERS[1:max_genes],
                                             collapse = ", ")
                               ),
                     tags$h4(HTML("<br/>")),
                     tags$h4(paste("Provide up to ", max_genes,
                                   " gene names that you then can use to build ",
                                   "models or create data sets from scratch.")),
                     tags$h4(HTML("<br/>")),
                     tags$h4("Separate you gene names with a ','. ",
                             "Do no use 'WT' for any gene name. ",
                             "Use only letters, numbers, and underscore ",
                             "but no other characters ",
                             "(of course, do not use comma as part of a gene name). ",
                             " Start gene names with letters ",
                             "(i.e., do not start them with numbers ",
                             "or other characters such as ._-, etc.). ",
                             "Try to keep gene names short (for figures)."
                             ),
                     tags$h4(HTML("<br/>")),
                     tags$h4(HTML("<b>IMPORTANT:</b> Again, do not think about this option as a way ",
                                  "to rename genes in existing data or models. ",
                                  "That becomes confusing very quickly. ",
                                  "Instead, think of this as a way to build ",
                                  "<b>new models from scratch</b> ",
                                  "with the names that you want. ")),
                     tags$h4(HTML("<br/>")),
                     tags$h4(HTML("Likewise, enter all new names in one go ",
                                  "(to try to enforce this, the button \"Use these gene names\" ",
                                  "will be disabled/closed as soon as you click it). ",
                                  "Repeatedly editing the new gene names ",
                                  "can easily lead to mismatch errors ",
                                  "(that will be displayed in red, with text ",
                                  "that contains \"all(gene_names_in_freqs %in% valid_gene_names \"). ",
                                  "This is inconsequential; delete all genotype data, and start again.")),
                     tags$div(class = "download_button",
                              tags$h4(HTML("<br/>")),
                              actionButton("action_provide_gene_names",
                                           "Use these gene names"),
                              ),
                     ),
            easyClose = TRUE
        ))  
    })


    ## Updating gene names
    observeEvent(input$action_provide_gene_names,{
        tryCatch({
            mymessage("At action_provide_gene_names")
            ## The disabling of the option does not work well sometimes.
            ## So instead of filling up the code with disables,
            ## catch it here.
            ## I still leave a bunch of them, as it would be ideal not to
            ## end here.
            if ((!is.null(data$data) ||
                 (nrow(data$csd_counts) > 0))) {
                mymessage("    disabled provide_gene_names action_provide_gene_names")
                shinyjs::disable("provide_gene_names")
                stop("As the tooltip and box text explained, ",
                     "you should ONLY try to ",
                     "change gene names on empty models/data. ",
                     "This data already has values; ",
                     "you need to delete all data first. ",
                     "We will abort this operation.")
            }
            
            old_gene_names <- data$gene_names
            new_gene_names <-
                strsplit(gsub(" ", "", input$new_gene_names), ",")[[1]]
            if (isTRUE(any(duplicated(new_gene_names)))) {
                stop("Duplicated new gene names.")
            }
            sanity_new_gene_names(new_gene_names)
            if (length(new_gene_names) < max_genes) {
                new_gene_names <- c(new_gene_names,
                                    LETTERS[(length(new_gene_names) + 1):max_genes]
                                    )
                if (isTRUE(any(duplicated(new_gene_names)))) {
                    stop("Duplicated gene names between new and old entries.")
                }
            }
            data$gene_names <- new_gene_names
            
            ## Use a simple lookup-dictionary and 
            ## avoid to_stnd_csd_dataset which is a function from hell.

            names_dict <- new_gene_names
            names(names_dict) <- old_gene_names
            ## For the DAG
            names_dict <- c(names_dict, "Root" = "Root")

            new_data <- list()
            new_data$gene_names <- new_gene_names
            new_data$name <- data$name
            new_data$lambdas <- data$lambdas
            new_data$DAG_parent_set <- data$DAG_parent_set
            new_data$this_d_dag_model <- data$this_d_dag_model
            
            ## BEWARE! If we do not do this, new_data$dag,
            ## because of partial matching, gets DAG_parent_set
            ## Not anymore, as I rewrote as DAG_parent_set. Anyway.
            if (!is.null(data$dag)) {
                new_data[["dag"]] <- data[["dag"]]
            } else {
                new_data$dag <- NULL
            }
            if (!is.null(data$thetas)) {
                new_data$thetas <- data$thetas
            } else {
                new_data$thetas <- NULL
            }
            new_data$data <- data$data

            mymessage("        At action_provide_gene_names: 2")
            ## To rename, use lookup
            names(new_data$lambdas) <- names_dict[names(new_data$lambdas)]
            names(new_data$DAG_parent_set) <- names_dict[names(new_data$DAG_parent_set)]

            if (!is.null(new_data[["dag"]])) {
                colnames(new_data$dag) <- names_dict[colnames(new_data$dag)]
                rownames(new_data$dag) <- names_dict[rownames(new_data$dag)]
            }
            if (!is.null(new_data$thetas)) {
                colnames(new_data$thetas) <- names_dict[colnames(new_data$thetas)]
                rownames(new_data$thetas) <- names_dict[rownames(new_data$thetas)]
            }
            if (!is.null(new_data$data)) {
                colnames(new_data$data) <- names_dict[colnames(new_data$data)]
            }
            ## To create
            new_data$csd_counts <- evamtools:::get_csd(new_data$data)

            ## Assign to the correct places: update the info
            data$gene_names <- new_gene_names
            data$data <- new_data$data
            data$dag <- new_data$dag
            data$DAG_parent_set <- new_data$DAG_parent_set
            data$thetas <- new_data$thetas
            data$lambdas <- new_data$lambdas
            data$csd_counts <- new_data$csd_counts
            ## data$this_d_dag_model has not changed
            mymessage("        At action_provide_gene_names: 3")
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- new_data
            ## Disable as soon as clicked on "Use these genes"
            shinyjs::disable("action_provide_gene_names")
            removeModal()
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    

    
    output$define_genotype <- renderUI({
        n_genes <- ifelse(is.null(input$gene_number), default_number_genes,
                          input$gene_number)

        gene_options <- set_gene_names_after_resize(data$data,
                                                    data$gene_names)[1:n_genes]
        if (input$input2build == "csd") {
            tags$div(
                     tags$h3("Add genotypes"),
                     tags$h5("WT is added by not clicking on any mutations. "),
                     tags$h5("Any gene without mutations is excluded from the data, ",
                             "regardless of the setting for number of genes. "),
                     tags$h5("For the CPM analysis, if any gene is always observed mutated ",
                             "(i.e., has a constant value of 1 for all observations), ",
                             "one observation with no genes mutated is added ",
                             "to the sample before the analysis."),
                     tags$h5("For the CPM analysis, genes that have identical patterns",
                             "(i.e., that lead to identical columns in the data matrix), ",
                             "are fused into",
                             "a single gene."),
                     tags$h5(" "),
                     tags$div(class = "inline",
                              checkboxGroupInput(inputId = "genotype",
                                                 label = "Mutations",
                                                 choices = set_gene_names_after_resize(data$data,
                                                                                       data$gene_names)[1:n_genes]),
                              tippy::tippy_this("genotype",
                                                HTML("<span style='font-size:1.5em; text-align:left;'><p>The list of genes next to \"Mutations\" is kept sorted ",
                                                     "(\"natural order\") showing first the genes ",
                                                     "in existing genotypes, and then other genes up to ",
                                                     "\"Number of genes\". ",
                                                     "The list will be resorted   ",
                                                     "when new genes are added to genotypes. </p>",
                                                     "<br> <p>\"Mutations\" is just a list of candidate gene names. You can ",
                                                     "see more (or fewer, up to the number of genes in your genotypes) ",
                                                     "by moving the slider of \"Number of genes\".</p>",
                                                     "<br>",
                                                     "<p>On renaming of the data, we trigger a counting of the ",
                                                     "genes used, reset the slider ",
                                                     "of \"Number of genes\", and reorder the genes next to  ",
                                                     "\"Mutations\".</p><span>"),
                                                arrow = TRUE, animation = "shift-toward"
                                                ),
                              ),
                     tags$div(id="fr",
                              numericInput(label = "Counts", value = NA, min = 0,
                                           inputId = "genotype_freq", width = NA),
                              actionButton("add_genotype", "Add genotype")
                              ),
                     ## shinyBS::removeTooltip(session, "genotype"),
                     )
        } else if (input$input2build == "dag") {
            ## Make sure all genes currently in the DAG are in
            ## the To and From to add/remove
            mymessage("At output$define_genotype, DAG")

            ## ## Otherwise, it gets enabled again occasionally.
            if ((!is.null(data$data) ||
                 (nrow(data$csd_counts) > 0))) {
                mymessage("    disabled provide_gene_names under select_csd")
                shinyjs::disable("provide_gene_names")
            }
            
            current_dag_data <- dag_data()
            if (is.null(current_dag_data)) mymessage("   current_dag_data is NULL")
            genes_in_dag <- setdiff(unique(c(current_dag_data$From,
                                             current_dag_data$To)),
                                    "Root")
            if (length(genes_in_dag) < n_genes) {
                genes_not_in_dag <- setdiff(gene_options, genes_in_dag)
                dag_gene_options <-
                    evamtools:::evam_string_sort(c(genes_in_dag,
                                                   genes_not_in_dag[1:(n_genes - length(genes_in_dag))]))
            } else {
                dag_gene_options <- evamtools:::evam_string_sort(genes_in_dag)
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
                                                        choiceNames = list("OT", "OncoBN", "CBN/H-ESBCN"),
                                                        choiceValues = list("OT", "OncoBN", "HESBCN"),
                                                        selected = data$this_d_dag_model)
                                           ),
                                  tags$h4("New edge"),
                                  tags$h5(HTML("<p></p>")),
                                  tags$div(class = "inline",
                                           radioButtons(inputId = "dag_from",
                                                        label = "From (parent node)",
                                                        inline = TRUE,
                                                        choices =  c("Root", dag_gene_options)),
                                           tippy::tippy_this("dag_from",
                                                             HTML("<span style='font-size:1.5em; text-align:left;'>",
                                                                  "The list of genes next to 'From' and 'To' is kept sorted ",
                                                                  "(alphabetically, often callen 'natural order'), with 'Root' first. ",
                                                                  " You can ",
                                                                  "see more genes (or fewer, up to the number of genes in your genotypes) ",
                                                                  "by moving the slider of 'Number of genes'.<span>"),
                                                             arrow = TRUE, animation = "shift-toward", placement = "right")
                                           ),
                                  tags$div(class = "inline",
                                           radioButtons(inputId = "dag_to",
                                                        label = " To (child node)",
                                                        inline = TRUE,
                                                        choices =  dag_gene_options),
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
                                               "increase the number of nodes in the DAG.")),
                                  tags$h3(HTML("<br/>DAG table")),
                                  tags$h4(HTML("Remember to hit Ctrl-Enter when you are done editing the DAG table for changes to take effect.")),
                                  DT::DTOutput("dag_table"),
                                  tags$h3(HTML("<br/>")),
                                  tags$h4(HTML("<br/>")),
                                  tags$h4(HTML("<u>2. Generate data from the DAG model</u>")),
                                  tags$h4(HTML("<br/>")),
                                  numericInput("epos",
                                               HTML("epos,&epsilon;"),
                                               value = generate_data$epos,
                                               min = 0, max = 0.995,
                                               step = 0.005, width = "12em"),
                                  tags$h5(HTML("For <b>OT (epos) and OncoBN (&epsilon;) only</b>: ",
                                               "probability that children nodes "),
                                          "not allowed by the model (the DAG) occur. ",
                                          "Accepted values: [0, 1). ",
                                          "This setting affects predicted probabilities. "),
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;",
                                      numericInput("dag_num_samples",
                                                   HTML("Number of genotypes<br>to sample"),
                                                   value = generate_data$dag_num_samples,
                                                   min = 100,
                                                   max = 10000,
                                                   step = 100,
                                                   width = "22em"),
                                      ),
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;", 
                                      numericInput("dag_obs_noise",
                                                   HTML("Observational noise <br>(genotyping error)"),
                                                   value = generate_data$dag_obs_noise,
                                                   min = 0, max = 0.95,
                                                   step = 0.02500, width = "18em"),
                                      ), 
                                  tippy::tippy_this("dag_obs_noise",
                                                    paste("<span style='font-size:1.5em; text-align:left;'>",
                                                          "A proportion between 0 and 1 ",
                                                          "(open interval on the right, so 1 is not allowed).",
                                                          "Observational noise (e.g., genotyping error) ",
                                                          "for all models. ",
                                                          "Added during sampling, ",
                                                          "after predictions from model ",
                                                          "have been obtained; ",
                                                          "predicted probabilities are not affected.",
                                                          "If larger than 0, this proportion of entries ",
                                                          "in the sampled matrix will be flipped ",
                                                          "(i.e., 0s turned to 1s and 1s turned to 0s).",
                                                          "<span>"
                                                          ),
                                                    arrow = TRUE, animation = "shift-toward", placement = "right"),
                                  tags$h5(HTML("<br/>")),
                                  actionButton("resample_dag", "Generate data from DAG"),
                                  actionButton("clear_dag", HTML("Reset DAG and delete genotype data")) ,
                                  tippy::tippy_this("clear_dag",
                                                    paste("<span style='font-size:1.5em; text-align:left;'>",
                                                          "Resetting the DAG will replace the ",
                                                          "contents of the named object by ",
                                                          "those of the default one ",
                                                          "(a three-gene fork with lambdas = 0.5), ",
                                                          "and will remove the generated genotype ",
                                                          "data.",
                                                          "<span>"),
                                                    arrow = TRUE, animation = "shift-toward", placement = "right"),
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
                                      numericInput("mhn_num_samples",
                                                   HTML("Number of genotypes<br>to sample"),
                                                   value = generate_data$mhn_num_samples,
                                                   min = 100,
                                                   max = 10000,
                                                   step = 100, width = "22em"),
                                      ),
                                  tags$h3(HTML("<br/>")),
                                  div(style = "white-space: nowrap;",
                                      numericInput("mhn_obs_noise",
                                                   HTML("Observational noise<br>(genotyping error)"),
                                                   value = generate_data$mhn_obs_noise,
                                                   min = 0, max = 1,
                                                   step = 0.025, width = "18em"),
                                      ),
                                  tippy::tippy_this("mhn_obs_noise",
                                                    paste("<span style='font-size:1.5em; text-align:left;'>",
                                                          "A proportion between 0 and 1 ",
                                                          "(open interval on the right, so 1 is not allowed).",
                                                          "Observational noise (e.g., genotyping error) ",
                                                          "for all models. ",
                                                          "Added during sampling, ",
                                                          "after predictions from model ",
                                                          "have been obtained; ",
                                                          "predicted probabilities are not affected.",
                                                          "If larger than 0, this proportion of entries ",
                                                          "in the sampled matrix will be flipped ",
                                                          "(i.e., 0s turned to 1s and 1s turned to 0s).",
                                                          "<span>"
                                                          ),
                                                    arrow = TRUE, animation = "shift-toward"
                                                  , placement = "right"),
                                  tags$h5(HTML("<br/>")),
                                  actionButton("resample_mhn", "Generate data from MHN model"),
                                  actionButton("clear_mhn",
                                               HTML("Reset log-&Theta; matrix and delete genotype data")),
                                  tippy::tippy_this("clear_mhn",
                                                    paste("<span style='font-size:1.5em; text-align:left;'>",
                                                          "Resetting the log-&Theta; matrix will replace the ",
                                                          "contents of the named object by ",
                                                          "those of the default one ",
                                                          "(a three-gene matrix filled with 0s), ",
                                                          "and will remove the generated genotype ",
                                                          "data."),
                                                    arrow = TRUE, animation = "shift-toward"
                                                  , placement = "right"
                                                    )
                              )
                     }
                 )
        } else if (input$input2build == "upload") {
            tags$div(## class = "frame",
                     tags$h3("Upload data (CSV format)"),
                     actionButton("display_help_upload_data",
                                  "Help", class = "btn-info"),
                     tags$h5(HTML("If you want to give your data a specific ",
                                  "name, set it in the box below ",
                                  "before uploading the data. ",
                                  "File names should start with a letter, ",
                                  "and can contain only letters, numbers, ",
                                  "hyphen, and underscore, but no ",
                                  "other characters (no periods, spaces, etc)."
                                  )),
                     tags$div(class = "inlin3",
                              textInput(inputId = "name_uploaded",
                                        label = "Name for data",
                                        value = "Uploaded_data"
                                        )
                              ),
                     
                     ## tags$h5(paste0("Format: csv ---comma separated values---,",
                     ##                "with subjects in rows and gens in columns. ",
                     ##                "Use 0 or 1 for ",
                     ##                "altered/not-altered (mutated/not-mutated).", 
                     ##                "The first row should contain ",
                     ##                " the gene names ",
                     ##                "and there should be no subject names. "
                     ##                )),
                     ## tags$h5(HTML("Use only alphanumeric characters ",
                     ##              "for gene names, and do not start ",
                     ##              "a gene name with a number. <br> ",
                     ##              "Keep gene names short (for figures). <br>"
                     ##              )),
                     ## tags$h5(HTML("A short example file is available here ",
                     ##              "https://raw.githubusercontent.com/rdiaz02/EvAM-Tools/main/examples_for_upload/tinydata.csv" ,

                     ##              "and there are additional examples in that ",
                     ##              "directory ",
                     ##              "(https://github.com/rdiaz02/EvAM-Tools/tree/main/examples_for_upload ---make sure to view the files as raw)")),
                     tags$div(class = "upload_file",
                              fileInput("csd", "Load data",
                                        multiple = FALSE,
                                        accept = c(
                                            "text/csv",
                                            ".csv"))),
                     tags$h5(HTML("<br/>")),
                     )
        }
    })
    
    ## With data modification for upload
    output$change_counts <- renderUI({
        ## if (input$input2build %in% c("upload", "csd", "dag", "matrix")) {
        menu_num <- ifelse(input$input2build == "upload", "2", "3")
        tags$div(class = "frame",
                 tags$div(class = "flex",
                          ## tags$h3(paste0(menu_num, " . Change genotype's counts")),
                          tags$h3("Change genotype's counts"),
                          actionButton("display_help_change_genotype_counts",
                                       "Help", class = "btn-info"),
                          tags$h3(HTML("<br/>")),
                          ## shinyBS::removeTooltip(session, "genotype"),
                          ),
                 tags$div(id = "csd_table",
                          DT::DTOutput("csd_counts")
                          ),
                 tags$h5(HTML("<br/>")),
                 ## We could include upload and dag here, but it makes no sense
                 if (input$input2build %in% c("csd")) {
                     actionButton("clear_genotype", "Delete all genotype data")
                 } else if (input$input2build %in% c("matrix")) {
                     
                     tags$h5(HTML("To delete all genotype data, use",
                                  "'Reset log-&Theta; matrix and delete genotype data'",
                                  "above. ",
                                  "That will also reset the model. ",
                                  "If you just want new genotype data from ",
                                  "the same model, click on ",
                                  "'Generate data from MHN model'."))
                     ## If you uncomment below, that is what takes
                     ## effect, and the h5 message is not shown.
                     ## actionButton("clear_genotype", "Delete all genotype data")
                 } else if (input$input2build %in% c("dag")) {
                     tags$h5(HTML("To delete all genotype data, use ",
                                  "'Reset DAG and delete genotype data'",
                                  "above. ",
                                  "That will also reset the model. ",
                                  "If you just want new genotype data from ",
                                  "the same model, click on ",
                                  "'Generate data from DAG'."))
                     ## To remove just the genotype
                     ## it makes little sense. For new DAGs, with new names
                     ## we want people to start from scratch.
                     ## And clear_genotype clears everything for now.
                     ## actionButton("clear_genotype", "Delete all genotype data")
                 } else if (input$input2build %in% c("upload")) {
                     tags$h5(HTML("To delete (or reset) all genotype data",
                                  "upload a new (or the same) data file."))
                 }
                 )
        ##  }
    })

    observeEvent(input$dag_model, {
        tryCatch({
            mymessage("At observeEvent input$dag_model")
            former_dag_model <- data$this_d_dag_model
            can_change_dag_model <- TRUE

            if ((input$dag_model %in% c("OT", "OncoBN")) &&
                (!is.null(data$lambdas)) &&
                (any(data$lambdas > 0.99999999))
                ) {
                can_change_dag_model <- FALSE
                updateRadioButtons(session, "dag_model", selected = former_dag_model)
                stop("thetas/probabilities should be between 0 and 1 ",
                     "(actually, for numerical reasons, 0.99999999).")
            }
            if (input$dag_model == "OncoBN") {
                if (any(data$DAG_parent_set == "XOR")) {
                    can_change_dag_model <- FALSE
                    updateRadioButtons(session, "dag_model", selected = former_dag_model)
                    stop("The OncoBN model cannot include ",
                         "XOR relationships.")
                } else if (length(unique(data$DAG_parent_set)) > 2) {
                    can_change_dag_model <- FALSE
                    changed_dag_model$invalidate <- FALSE
                    updateRadioButtons(session, "dag_model", selected = former_dag_model)
                    stop("The OncoBN model can only include ",
                         "one type of relationship",
                         "(conjunctive ---AND--- or disjunctive ---OR---, ",
                         "as specified in \"Advanced options\").")
                }
            } else if (input$dag_model == "OT") {
                number_of_parents <- colSums(data$dag)
                if (any(number_of_parents > 1)) {
                    can_change_dag_model <- FALSE
                    updateRadioButtons(session, "dag_model", selected = former_dag_model)
                    stop("This DAG has nodes with multiple parents. ",
                         "OT can only use trees ",
                         "(i.e., no node can have with multiple parents.).")
                } else if (length(unique(data$DAG_parent_set)) > 2) {
                    can_change_dag_model <- FALSE
                    updateRadioButtons(session, "dag_model", selected = former_dag_model)
                    stop(HTML("The OT model is only for trees."))
                }
            }
            
            if (can_change_dag_model) {
                changed_dag_model$invalidate <- !changed_dag_model$invalidate
                data$this_d_dag_model <- input$dag_model
                datasets$all_csd[["dag"]][[input$select_csd]]$this_d_dag_model <- input$dag_model
                if (resample_trigger_from_data_change())  shinyjs::click("resample_dag")
            } else {
                message("Are we ever here??? Value is ", can_change_dag_model)
                changed_dag_model$invalidate <- FALSE
                updateRadioButtons(session, "dag_model", selected = former_dag_model)
            }
        }, 
        error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    
    ## DAG builder
    ## Controling dag builder
    toListen2 <- reactive({
        list(## input$dag_model,
            input$dag_table_cell_edit,
            data$name,
            ## because add_edge and remove_edge modify
            ## DAG_parent_set and the dag, leaving
            ## this in a possibly inconsistent state
            ## unless we recreate it 
            data$DAG_parent_set,
            data$dag,
            data$lambdas,
            ## dag_data will be triggered on successful dag_model change
            changed_dag_model$invalidate
        )
    })

    ## With the dag model itself being part of data,
    ## we can take for granted that anything that is in data must be correct
    ## so no need to check that weights are in (0, 1), etc.
    dag_data <- eventReactive(toListen2(), {
        if (isolate(input$input2build) != "dag") {
            mymessage("dag_data reactive: Why are you calling ",
                      "dag_data if not dealing with DAGs?")
            return()
        }
        mymessage("At dag_data reactive call")

        all_gene_names <- c("Root", data$gene_names)
        edges <- which(data$dag == 1, arr.ind = TRUE)
        tmp_DAG_parent_set <- data$DAG_parent_set
        x <- length(tmp_DAG_parent_set)
        
        dag_dataset_names <- unlist(lapply(datasets$all_csd$dag,
                                           function(x) x$name))

        if (!(data$name %in% dag_dataset_names)) {
            mymessage("    data$name not in dag_dataset_names. Returning a NULL")
            return(NULL)
        }
        mymessage("    dag_data_reactive, position 3")
        
        names(tmp_DAG_parent_set) <- all_gene_names[seq(2, x + 1)]
        dag_data <- data.frame(From = all_gene_names[edges[, "row"]]
                             , To = all_gene_names[edges[, "col"]]
                             , Relation = tmp_DAG_parent_set[edges[, "col"] - 1]
                             , Lambdas = data$lambdas[edges[, "col"] - 1])

        if (data$this_d_dag_model %in% c("OT")) {
            colnames(dag_data) <- c("From", "To", "Relation", "Weight")
            dag_data$Relation <- NULL
        } else if (data$this_d_dag_model %in% c("OncoBN")) {
            colnames(dag_data) <- c("From", "To", "Relation", "theta")
        }

        ## isolate(changed_dag_model$invalidate <- FALSE)
        return(dag_data)
    })


    output$dag_table <-
        DT::renderDT(
                dag_data(),
                escape = FALSE,
                selection = 'none',
                server = FALSE,
                rownames = FALSE,
                ## DO NOT set target to "columns". Yeah, it'd be easier for
                ## moving around with tab, but if you do, editing
                ## breaks (function modify_lambdas_and_parent_set_from_table )
                ## and leads to very confusing bugs. So, JUST DON'T.
                editable = list(target = "all",
                                disable = list(columns = c(0, 1))),
                options = list(dom = 't', paging = FALSE, ordering = FALSE,
                               columnDefs = list(list(className = 'dt-center',
                                                      targets = "_all")))
            )

    ## Adding new edge
    observeEvent(input$add_edge, {
        mymessage("At input$add_edge 1")
        from_gene <- input$dag_from
        to_gene <- input$dag_to
        tryCatch({
            tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene,
                                               operation = "add",
                                               parent_set = data$DAG_parent_set,
                                               dag_model = data$this_d_dag_model)
            data$dag <- tmp_data$dag
            data$DAG_parent_set <- tmp_data$parent_set
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)
            if (resample_trigger_from_data_change()) shinyjs::click("resample_dag")
            shinyjs::disable("provide_gene_names")
            mymessage("At input$add_edge 1")
        },error=function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Remove edge
    observeEvent(input$remove_edge, {
        mymessage("At input$remove_edge 2")
        from_gene <- input$dag_from
        to_gene <- input$dag_to
        tryCatch({
            tmp_data <- evamtools:::modify_dag(data$dag, from_gene, to_gene,
                                               operation = "remove",
                                               parent_set = data$DAG_parent_set,
                                               dag_model = data$this_d_dag_model)
            data$dag <- tmp_data$dag
            data$DAG_parent_set <- tmp_data$parent_set
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)
            if (resample_trigger_from_data_change())  shinyjs::click("resample_dag")
            mymessage("At input$remove_edge 2")
        },error=function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Clear DAG
    observeEvent(input$clear_dag, {
        tryCatch({
            mymessage("At clear_dag")
            tmp_data <- evamtools:::modify_dag(data$dag, NULL, NULL, operation = "clear")
            tmp_dag <- tmp_data$dag
            colnames(tmp_dag) <- rownames(tmp_dag) <- c("Root", data$gene_names)
            tmp_dag["Root", data$gene_names[1]] <- 1
            data$dag <- tmp_dag
            data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            data$data <- .ev_SHINY_dflt$template_data$data
            data$DAG_parent_set <- tmp_data$DAG_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$DAG_parent_set) <- data$gene_names
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)
            shinyjs::disable("analysis")
            mymessage("enable provide_gene_names under clear_dag")
            shinyjs::enable("provide_gene_names")
            mymessage("       At clear_dag: exit")
        }, error=function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    observeEvent(input$dag_table_cell_edit, {
        tryCatch({
            names(data$DAG_parent_set) <- data$gene_names[1:length(data$DAG_parent_set)]
            names(data$lambdas) <- data$gene_names[1:length(data$DAG_parent_set)]
            info <- input$dag_table_cell_edit

            ## Allowing the setting of a value of 0 to remove an edge
            ## FIXME: this block should be a function
            ##  For being kinda of a side issue makes the function
            ## very complicated.
            tryCatch({
                
                if (any(info$value == 0L)) {
                    mymessage("At dag_table_cell_edit, passed a 0 value; ",
                              "trying to remove the edge")
                    info_from <- info[info["col"] == 0, "value"]
                    info_to <- info[info["col"] == 1, "value"]
                    if (data$this_d_dag_model == "OT") {
                        lambda_col <- 2
                    } else {
                        lambda_col <- 3
                    }
                    lambda_val <-
                        suppressWarnings(as.numeric(info[info["col"] == lambda_col,
                                                         "value"]))
                    info_rows_to_rm <- info[info$value == 0L, "row"]
                    lambda_val_0 <- which(lambda_val == 0L)
                    if ( (length(info_rows_to_rm) != length(lambda_val_0)) ||
                         !(all(info_rows_to_rm == lambda_val_0))) {
                        stop("This delete operation using 0s cannot be performed. ",
                             "Try removing setting 0s one by one and/or ",
                             "using 'Remove edge' ")
                    }

                    ## For removing entries from info, and because more likely
                    ## to be terminal nodes
                    lambda_val_0 <- sort(lambda_val_0, decreasing = TRUE)

                    for (nd in lambda_val_0) {
                        from_gene <- info_from[nd]
                        to_gene   <- info_to[nd]
                        
                        post_delete <- evamtools:::modify_dag(data$dag, from_gene, to_gene,
                                                              operation = "remove",
                                                              parent_set = data$DAG_parent_set,
                                                              dag_model = data$this_d_dag_model)
                        data$dag <- post_delete$dag
                        data$DAG_parent_set <- post_delete$parent_set
                        info <- info[!(info[, "row"] == nd), ]
                    }
                }
            }, error = function(e) {
                e0 <- paste("Deleting one or more edge(s) by ",
                            "setting it/them to 0 did not succeed ",
                            "when you tried to remove edge ",
                            "<b>", from_gene, " -> ", to_gene, "</b> ",
                            "because: ",
                            "\n-----------------------\n")
                e2 <- paste("\n-----------------------\n",
                            "Some of the operations requested ",
                            "might have been performed successfully. ",
                            "Check the DAG. ")
                stop(paste(e0, "\n\n", e[[1]], "\n\n", e2))
            })

            mymessage("At dag_table_cell_edit, 1")
            tmp_data <-
                evamtools:::modify_lambdas_and_parent_set_from_table(dag_data(),
                                                                     info,
                                                                     data$lambdas
                                                                   , data$dag
                                                                   , data$DAG_parent_set
                                                                   , dag_model = data$this_d_dag_model)
            data$lambdas <- tmp_data$lambdas
            data$DAG_parent_set <- tmp_data$parent_set
            datasets$all_csd[[input$input2build]][[input$select_csd]] <-
                evamtools:::to_stnd_csd_dataset(data)

            mymessage("At dag_table_cell_edit, 2")
            shinyjs::click("resample_dag")
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    ## Building trm from dag
    observeEvent(input$resample_dag, {
        tryCatch({
            mymessage("At resample_dag")
            if (sum(colSums(data$dag) > 0) < 2) {
                data$csd_counts <- data.frame(Genotype = character(),
                                              Counts = integer())
                data$data <- NULL
                shinyjs::disable("analysis")
                stop("There must be at least two genes ",
                     "in the DAG.\n\n",
                     "<b>Any frequency data is suspect until ",
                     "you solve this problem.</b> ",
                     "so will try to delete all genotype data. ")
            }
            the_dag_data <- dag_data()
            gene_names <- setdiff(unique(c(the_dag_data$From, the_dag_data$To)),
                                  "Root")
            if ((input$dag_num_samples < 1) ||
                (input$dag_num_samples >  max_allowed_num_samples)
                ) stop("Generate data: number of ",
                       "genotypes to sample cannot be ",
                       "less than 1 or greater than ",
                       max_allowed_num_samples,
                       ".")
            if ((input$dag_obs_noise < 0) ||
                (input$dag_obs_noise >= 0.99999999))
                stop(HTML("Generate data: observational noise ",
                          "cannot be ",
                          "less than 0 or greater than (or equal to) 1 ",
                          "(to prevent numerical problems, no larger than 0.99999999)."))

            if ((input$epos < 0) ||
                (input$epos >= 0.99999999))
                stop("Generate data: epos,e  ",
                     "cannot be ",
                     "less than 0 or greater than (or equal to) 1 ",
                     "(to prevent numerical problems, no larger than 0.99999999).")
            
            tmp_dag_data <-
                evamtools:::generate_sample_from_dag(the_dag_data
                                                   , data$DAG_parent_set[gene_names]
                                                   , noise = input$dag_obs_noise
                                                   , N = input$dag_num_samples
                                                   , dag_model = data$this_d_dag_model
                                                   , epos = input$epos)
            
            data$csd_counts <-
                tmp_dag_data$csd_counts[tmp_dag_data$csd_counts[, 2] > 0, ]
            data$data <- tmp_dag_data$data

            ## FIXME: most of this should not be needed
            ## FIXME: and csd counts are not being assigned?
            ## nope: they are derived from data$data
            datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- data$data
            datasets$all_csd[[input$input2build]][[input$select_csd]]$dag <- data$dag
            datasets$all_csd[[input$input2build]][[input$select_csd]]$lambdas <- data$lambdas
            datasets$all_csd[[input$input2build]][[input$select_csd]]$DAG_parent_set <- data$DAG_parent_set
            datasets$all_csd[[input$input2build]][[input$select_csd]]$this_d_dag_model <- data$this_d_dag_model
            
            generate_data$dag_num_samples <- input$dag_num_samples
            generate_data$dag_obs_noise   <- input$dag_obs_noise
            generate_data$mhn_num_samples <- input$dag_num_samples
            generate_data$mhn_obs_noise   <- input$dag_obs_noise
            generate_data$epos        <- input$epos
            
            shinyjs::enable("analysis")
            ## The next is not really necessary, but we do it for consistency
            shinyjs::disable("provide_gene_names")
            
        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })



    ## Help for output of downloaded before results
    observeEvent(input$how2downloaddata, {
        showModal(modalDialog(
            easyClose = TRUE,
            title = tags$h3("Download data"),
            tags$div(
                     tags$h5(HTML("Contents of saved file: ",
                                  "the data as an R data frame in an RDS file.  ",
                                  "If you built a DAG or MHN model, ",
                                  "also the model built."
                                  )),
                     )
        ))
    }
    )
    ## Help for output of downloaded CPM results
    observeEvent(input$how2downloadcpm, {
        ## sometimes the tooltips misbehave
        ## shinyBS::removeTooltip(session, "table_out3")
        ## shinyBS::removeTooltip(session, "freq2label")
        
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
                     tags$p("0. Select the 'Type of model'. "
                            ##, "(Whenever you change the type of model, ",
                            ## "to avoid confusion it is best to click on ",
                            ## "'Reset DAG and delete genotype data', or ",
                            ## "you might get errors from asking for impossible ",
                            ## "settings ---e.g., moving from CBN with ",
                            ## "multiple parents to OT.)"
                            ),
                     tags$p("1. Select 'From' (parent node) and 'To' (child node) ",
                            HTML("and hit 'Add edge' or 'Remove edge'.<ul>")),
                     tags$li(HTML("An edge wont be allowed if: <ul>")),
                     tags$li("it is already present;"),
                     tags$li("it introduces cycles."),
                     tags$p(HTML("</ul>")),
                     tags$li(HTML("To remove an edge you can also "),
                             "set the lambda/weight/theta of the relationship to 0."),
                     tags$li(HTML("Removing edges might restructure the DAG."),
                             "If a node has no parent, " ,
                             "it will be assigned as descendant of Root."),
                     
                     tags$p(HTML("</ul>")),
                     tags$p(HTML("2. To <strong>change the value of a lambda/theta/Weight</strong> ",
                                 "click on the cell, ",
                                 "edit the cell's content, and press Ctrl+Enter.",
                                 "<ul>",
                                 "<li>theta (OncoBN) and Weight (OT) denote conditional probabilities ",
                                 "of an event occurring given the parent conditions are satisfied; ",
                                 "thus they must be between 0 and 1 ",
                                 "(open interval, so neither 0 nor 1 are allowed).</li>",
                                 "<li>lambdas (CBN, H-ESBCN) denote rates (time to occurrence of an event, ",
                                 "given its parents are satisfied, is modeled as an exponential ",
                                 "process with this rate). Thus, lambdas must be larger than 0.</li>",
                                 "<li>(Thus, we use setting weights/lambdas/thetas to 0 as one way of ",
                                 "removing edges from the DAG.)</li>",
                                 "</ul>")),
                     tags$p(HTML("3. Set the value of <strong>Relation</strong> ",
                                 "to one of 'Single' (single parent), ",
                                 "AND, OR, XOR.",
                                 "<ul>",
                                 "<li>OT only accepts 'Single' as each node ",
                                 "has a single parent. </li>",
                                 "<li>OncoBN accepts 'Single', 'AND', 'OR' ",
                                 "or combinations of either Single and AND or ",
                                 "Single and OR ",
                                 "(and terms not among Single, AND, OR ",
                                 "will be converted to ORs).</li>",
                                 "<li>'CBN/H-ESBCN' models can be specified with ",
                                 "AND, OR, XOR, Single, or combinations of the above ",
                                 "(and terms not among Single, AND, OR, XOR ",
                                 "will be converted to ANDs)[1]. </li>",
                                 "<li>All incoming edges to a node must have the same ",
                                 "Relation (the program will force this). </li>",
                                 "<li>Edit the cell's content and press Ctrl+Enter. </li>",
                                 "</ul>"
                                 )),
                     tags$p(HTML("4. Modify, if you want, ",
                                 "the <strong>'Number of genotypes to sample'</strong> ",
                                 "(the size of the sample) and ",
                                 "the <strong>Observational noise</strong> (genotyping error) ",
                                 "and, for OT and OncoBN, the <strong>epos</strong>, ",   
                                 "and click on <strong>'Generate data from DAG'</strong> to generate a sample. ")),
                     tags$p(HTML("<br>")),
                     tags$p(HTML("You can create, and simulate from, ",
                                 "<b>non-transitively reduced DAGs</b> ",
                                 "(not for OT, of course) ", 
                                 "but CBN only estimates transitively ",
                                 "reduced graphs. ",
                                 "More discussion about this issue is ",
                                 "available in the 'Additional documentation'")),
                     tags$p(HTML("<br>")),
                     tags$p("After the sample is generated for the first time, ",
                            "the sample should be generated again automatically ",
                            "whenever you change the model ",
                            "(add or remove edges, change lambdas, etc) ",
                            "if you hit Ctrl-Enter after you are done editing ",
                            "the DAG table."),
                     tags$p(HTML("<br>")),
                     tags$p("Whenever edges are removed/added or method changed ",
                            "in a way that could result in an inconsistent ",
                            "state of the generated data, we force that ",
                            "new data be generated."),
                     tags$p(HTML("<br>")),
                     tags$p("Values in the input boxes have arbitrary default values, ",
                            "but your last used values are preserved after generating data.",
                            "(Last includes also MHN, so it is easier  ",
                            "to use the same set of values across models)."),
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
                                                                                      targets = "_all")),
                                            info = FALSE, paginate= FALSE),
                                        )

    observeEvent(input$thetas_table_cell_edit, {
        tryCatch({
            info <-input$thetas_table_cell_edit
            data$thetas[1:input$gene_number, 1:input$gene_number] <-
                DT::editData(data$thetas[1:input$gene_number, 1:input$gene_number], info, "thetas")

            datasets$all_csd[[input$input2build]][[input$select_csd]]$thetas <- data$thetas
            datasets$all_csd[[input$input2build]][[input$select_csd]]$n_genes <- input$gene_number
            ## Resample based on changes
            if (resample_trigger_from_data_change()) shinyjs::click("resample_mhn")
            shinyjs::disable("provide_gene_names")            
        }, error = function(e){
            showModal(dataModal(e[[1]]))
        })
    })

    observeEvent(input$resample_mhn, {
        tryCatch({
            if ((input$mhn_num_samples < 1) ||
                (input$mhn_num_samples >  max_allowed_num_samples)
                ) stop("Generate data: number of ",
                       "genotypes to sample cannot be ",
                       "less than 1 or greater than ",
                       max_allowed_num_samples,
                       ".")
            if ((input$mhn_obs_noise < 0) ||
                (input$mhn_obs_noise >= 0.99999999)) stop("Generate data: observational noise ",
                                                          "cannot be ",
                                                          "less than 0 or greater than 1",
                                                          "(for numerical issues, no larger than 0.99999999).")
            
            mhn_data <- evamtools:::get_mhn_data(data$thetas[1:input$gene_number
                                                           , 1:input$gene_number]
                                               , noise = input$mhn_obs_noise 
                                               , N = input$mhn_num_samples)
            data$csd_counts <- mhn_data$csd_counts
            data$data <- mhn_data$data
            datasets$all_csd[[input$input2build]][[input$select_csd]]$data <- mhn_data$data

            generate_data$mhn_num_samples <- input$mhn_num_samples
            generate_data$mhn_obs_noise   <- input$mhn_obs_noise
            generate_data$dag_num_samples <- input$mhn_num_samples
            generate_data$dag_obs_noise   <- input$mhn_obs_noise
            
            shinyjs::enable("analysis")
            ## The next is not really necessary, but we do it for consistency
            shinyjs::disable("provide_gene_names")
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
            data$DAG_parent_set <- .ev_SHINY_dflt$template_data$DAG_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$DAG_parent_set) <- data$gene_names
            shinyjs::disable("analysis")
            datasets$all_csd[[input$input2build]][[input$select_csd]] <- evamtools:::to_stnd_csd_dataset(data)
            mymessage("    enabled provide_gene_names under clear_mhn")
            shinyjs::enable("provide_gene_names")
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
                     tags$p(HTML("5. Modify, if you want, ",
                                 "the <strong>'Number of genotypes to sample'</strong> ",
                                 "(the size of the sample) and ",
                                 "the <strong>Observational noise</strong> (genotyping error) ",
                                 "and click on <strong>'Generate data from MHN model'</strong> to generate a sample. ",
                                 "The sample is also updated as soon as you save an entry ",
                                 "in the matrix or change the number of genes.")),
                     tags$p(HTML("<br>")),
                     tags$p(HTML("You can make sure <b>the &theta;s have been updated</b> "),
                            "by checking the figure of the matrix on the right."),
                     tags$p(HTML("<br>")),
                     tags$p("Values in the input boxes have arbitrary default values, ",
                            "but your last used values are preserved after generating data.",
                            "(Last includes also the DAG models, so it is easier  ",
                            "to use the same set of values across models)."),
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
            data$dag <- tmp_dag
            data$csd_counts <- .ev_SHINY_dflt$template_data$csd_counts
            data$data <- .ev_SHINY_dflt$template_data$data
            data$DAG_parent_set <- .ev_SHINY_dflt$template_data$DAG_parent_set
            data$lambdas <- .ev_SHINY_dflt$template_data$lambdas
            names(data$lambdas) <- names(data$DAG_parent_set) <- data$gene_names
            shinyjs::disable("analysis")
            mymessage("enabled provide_gene_names uder clear_genotype")
            shinyjs::enable("provide_gene_names")
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
                shinyjs::disable("provide_gene_names")
            } else {
                showModal(modalDialog(paste("Counts <= 0 ",
                                            " (or non-numeric values) present. ",
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
  , selection = 'none', server = TRUE,
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

            if (nrow(data$csd_counts) == 0) {
                data$data <- NULL
            } else {
                data$data <-
                    datasets$all_csd[[input$input2build]][[input$select_csd]]$data <-
                        evamtools:::genotypeCounts_to_data(data$csd_counts, e = 0)
            }
            shinyjs::disable("provide_gene_names")
        }, error = function(e) {
            showModal(dataModal(e[[1]]))
        })
    })

    ## ## Plot histogram of genotypes
    output$plot <- plotly::renderPlotly({
        tryCatch({
            mymessage("At output$plot")
            ## Not relevant anymore. Remove eventually.
            ## provide_gene_names is being enabled somewhere I can't locate
            ## So make sure we catch it right on the redisplay
            ## if ((!is.null(data$data) ||
            ##      (nrow(data$csd_counts) > 0))) {
            ##     mymessage("    disabled provide_gene_names under output_plot")
            ##     shinyjs::disable("provide_gene_names")
            ## }
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
            if (!is.null(data$DAG_parent_set)) edges$Relation <- data$DAG_parent_set[edges$To]
        } else if (input$input2build %in% c("matrix") 
                   && !is.null(data$thetas)
                   && !is.null(input$gene_number)
                   && length(data$thetas[1:input$gene_number, 1:input$gene_number])>0
                   ) {
            data2plot <- data$thetas[1:input$gene_number, 1:input$gene_number]
        }
        evamtools:::plot_method(data2plot, data$DAG_parent_set, edges)
    })

    ## Run CPMs
    observeEvent(input$analysis, {
        ## Calculate TRM for DAG and for matrices
        tryCatch({

            if (ncol(data$data) >= 7) {
                showModal(
                    dataModal(paste("Beware! You are analyzing data ",
                                    "with 7 or more genes. ",
                                    "This can take longer than usual ",
                                    "and plots may be crowded. "),
                              ## "We recommend using top_paths options in ",
                              ## "the Results' tab.",
                              type = "Warning: "))
            }

            if (is.null(input$cpm_methods) ||
                (length(input$cpm_methods) ==1 && is.na(input$cpm_methods)))
                stop("You must use at least one method ",
                     "(check 'CPMs to use' under 'Advanced options ",
                     "and CPMs to use').")

                
                shinyjs::disable("analysis")

                progress <- shiny::Progress$new()
                ## Make sure it closes when we exit this reactive, even if there's an error
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
                    MCMC_iter = input$HESBCN_MCMC_iter,
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
                                                       , mccbn_opts = mccbn_opts
                                                         ## FIXME: remove next?
                                                         ## This is just in case shiny's code
                                                         ## depends on NA entries
                                                       , only_used_methods = FALSE
                                                         )},
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
                                , dag = data$dag, DAG_parent_set = data$DAG_parent_set)

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
                ## last_visited_cpm <<- result_name
                reactive_last_visited_cpm$the_last_visited_cpm <<- result_name
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
        gsub("H-ESBCN", "HESBCN", input$cpm2show, fixed = TRUE)
    }),  900)

    ## ## No delay showing plots. 
    ## plot2show <- reactive({
    ##     gsub("H-ESBCN", "HESBCN", input$cpm2show, fixed = TRUE)
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
                    lapply(plot2show(), function(met) {
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
                                                 choices = gsub("HESBCN", "H-ESBCN",
                                                                input$cpm_methods, fixed = TRUE),
                                                 selected = gsub("HESBCN", "H-ESBCN",
                                                                 input$cpm_methods, fixed = TRUE)
                                                 ),
                              tippy::tippy_this("cpm2show",
                                                HTML("<span style='font-size:1.5em; text-align:left;'>",
                                                     "<p>Show graphical output of the CPMs used to analyze the data.  "
                                                   , "Go back to \"User input\" "
                                                   , "and click on \"Advanced options\" if you"
                                                   , "want to use other methods.</p>"
                                                   , "<p>When re-displaying previous results, "
                                                   , "if subsequent analyses used different methods, "
                                                   , "you will need to go back to "
                                                   , "\"User input\" , \"Advanced options\", "
                                                   , "and click on the methods you used. "
                                                   , "But do not re-run; click on the methods "
                                                   , "and come back to the 'Results' tab."
                                                   , "You can click on all methods: this is innocuous "
                                                   , "though the figure will contain empty wholes "
                                                   , "corresponding to methods clicked but not used in the run."
                                                   , "(Remove the empty holes by unclicking those methods here.)"
                                                  ,  "<span>"
                                                    ),
                                                arrow = TRUE, animation = "shift-toward"
                                              , placement = "right"),
                              tags$h4(HTML("<hr style=\"height:1px; width:80%; background-color:black;text-align:left\">")),
                              tags$h4(HTML("<br/>")),
                              tags$div(class = "inline",
                                       radioButtons(inputId = "data2plot",
                                                    label = HTML("Predictions from fitted models to display"
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
                                                    ), ## Prompter does not allow including HTML formatting.
                                       tippy::tippy_this("data2plot",
                                                         HTML("<span style='font-size:1.5em; text-align:left'>",
                                                              "<p>This output is also displayed in tabular form on the bottom right.</p>",
                                                              "<br><p><u>\"Sampled genotype counts\"</u> is only available if you selected ",
                                                              "\"Sample genotypes\" under \"Advanced options\". </p>",
                                                              "",
                                                              "<p>For <u>\"Predicted genotype relative frequencies\"</u> ",
                                                              "and <u>\"Sampled genotype counts\"</u>, the histograms only ",
                                                              "show the 20 most frequent genotypes ", ## 20: argument to  plot_genotype_counts
                                                              "because of figure size and legibility of genotype labels .",
                                                              "The table shows all the genotypes. </p>",
                                                              "",
                                                              "<p>For <u>\"Predicted genotype relative frequencies\"</u> ",
                                                              " in the table we add the relative frequencies of genotypes in the original data ",
                                                              "to make it easier to visually asses how close predictions are to observed data. ",
                                                              ## "(It would not be sensible to do this with \"Sampled genotype counts\" as ",
                                                              ## "those include additional sampling noise.) For easier comparison, both genotypes ",
                                                              ## "with no observed counts in the original data ",
                                                              ## "and genotypes with predicted frequency exactly 0 are left as empty (not as 0) ",
                                                              ## "in the displayed table.</p>",
                                                              "</p><span>"
                                                              ),
                                                         arrow = TRUE, animation = "shift-toward"
                                                       , placement = "right"
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
                                                    ),
                                       tippy::tippy_this("label2plot",
                                                         paste("<span style='font-size:1.5em; text-align:left;'>",
                                                               "Type of label for transition [rate] plots."
                                                               ),
                                                         arrow = TRUE, animation = "shift-toward"
                                                       , placement = "right"
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
                                          value = 5, max = 10, min = 0, step = 1),
                              tippy::tippy_this("freq2label-wrap",
                                                HTML("<span style='font-size:1.5em; text-align:left;'><p>Set it to 0 to show all paths or all genotype labels. <p>",
                                                     "<p>You probably <u>do NOT want to set it to 0</u> when there are many ",
                                                     "genes, specially if you are also displaying the output ",
                                                     "from MHN (where any genotype is connected to each of its ",
                                                     "descendants) or from models for which the DAG is close to a star. </p><span>"
                                                     ),
                                                arrow = TRUE, animation = "shift-toward"
                                              , placement = "right")
                              )
                     ),
            
            )
    })

    output$cpm_list <- renderUI({
        all_names <- c()
        for (i in names(all_cpm_out)) {
            all_names <- c(all_names, all_cpm_out[[i]]$orig_data$name)
        }
        
        ## if ((length(all_names) > 0) && (last_visited_cpm != "")) {
        if ((length(all_names) > 0) &&
            (reactive_last_visited_cpm$the_last_visited_cpm != "")) {
            selected <- names(all_cpm_out)
            
            tagList(
                radioButtons(
                    inputId = "select_cpm",
                    label = "",
                    selected = reactive_last_visited_cpm$the_last_visited_cpm,
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
                     ## Disabled because it led to inconsistent behavior such as this
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

    output$cpm_freqs <- DT::renderDT({
        d1 <- all_cpm_out[[input$select_cpm]]$tabular_data[[input$data2plot]]

        if (input$data2plot == "predicted_genotype_freqs") {
            d2 <-
                evamtools:::get_csd(all_cpm_out[[input$select_cpm]]$cpm_output$analyzed_data)
            d2$Counts <- round(d2$Counts/sum(d2$Counts), 3)
            colnames(d2)[2] <- "Original_data"
            ## Yes, set to NA for easier visualization.
            ## They will also be NA in the Original data if not present
            d1[d1 == 0] <- NA
            d3 <- dplyr::full_join(d1, d2, by = "Genotype")
            d3 <- evamtools:::reorder_to_standard_order_arbitrary_df(d3)
            d1 <- d3
        }
        d1
    },
    selection = 'none', server = TRUE
  , rownames = FALSE
  , options = list(
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        info = FALSE, paginate= FALSE)
    )

    output$tabular_data <- renderUI({
        if (length(names(all_cpm_out)) > 0) {
            tags$div(class="frame max_height",
                     id = "table_out1",
                     tags$div(class=" max_height",
                              ## A hack to get the tooltip to only show on hover
                              ## over heading of table
                              id = "table_out3", 
                              tags$h3(paste("Tabular output of predictions from fitted models: ",
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
                     tippy::tippy_this("table_out3",
                                       ## message = 
                                       HTML("<span style='font-size:1.5em; text-align:left;'><p>This output is also displayed as the second row of figures. ",
                                            "Choose the output to display from the left radio buttons ",
                                            "\"Predictions from models to display\".</p>",
                                            "<br><p><u>\"Sampled genotype counts\"</u> is only available if you selected ",
                                            "\"Sample genotypes\" under \"Advanced options\". </p>",
                                            "<br>",
                                            "<p>For <u>\"Predicted genotype relative frequencies\"</u> ",
                                            "and <u>\"Sampled genotype counts\"</u>, the histograms only ",
                                            "show the 20 most frequent genotypes ", ## 20: argument to  plot_genotype_counts
                                            "for reasons of figure size and legibility of genotype labels .",
                                            "The table shows all the genotypes. </p></span>",
                                            )
                                     , arrow = TRUE, animation = "shift-toward"),
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
        filename = function() sprintf("%s_cpm.rds", input$select_cpm),
        content = function(file) {
            saveRDS(all_cpm_out[[input$select_cpm]][c("cpm_output", "tabular_data")],
                    file)
        }
    )
}










######################################################################
######################################################################
#########
#########    Old code  
#########
######################################################################
######################################################################

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

## Old consistency checks in dag_table_cell_edit
## but catching them here is a little bit late.
## And this problem can no longer exist here, I think.
## 
## cnames <- colnames(dag_data())
## if (data$this_d_dag_model == "OT") {
##     if ( (any(c("Lambdas", "theta", "Relation") %in% cnames)) ||
##          !("Weight" %in% cnames ))
##         stop("OT model with impossible DAG table columns. ",
##              "Did you change the model while editing? ",
##              "Reset the DAG and expect possible app errors.")
## } else if (data$this_d_dag_model == "OncoBN") {
##     if ( (any(c("Lambdas", "Weight") %in% cnames)) ||
##          !("theta" %in% cnames ) ||
##          !("Relation" %in% cnames))
##         stop("OncoBN model with impossible DAG table columns. ",
##              "Did you change the model while editing? ",
##              "Reset the DAG and expect possible app errors.")
## } else if (data$this_d_dag_model == "HESBCN") {
##     if ( (any(c("theta", "Weight") %in% cnames)) ||
##          !("Lambdas" %in% cnames ) ||
##          !("Relation" %in% cnames))
##         stop("HESBCN model with impossible DAG table columns. ",
##              "Did you change the model while editing? ",
##              "Reset the DAG and expect possible app errors.")
## }



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
##         new_data$DAG_parent_set <- data$DAG_parent_set
##         new_data$dag <- data$dag
##         new_data$thetas <- data$thetas
##         new_data$data <- data$data

##         ## To rename, use lookup
##         names(new_data$lambdas) <- names_dict[names(new_data$lambdas)]
##         names(new_data$DAG_parent_set) <- names_dict[names(new_data$DAG_parent_set)]
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
##         data$DAG_parent_set <- new_data$DAG_parent_set
##         data$thetas <- new_data$thetas
##         data$lambdas <- new_data$lambdas
##         data$csd_counts <- new_data$csd_counts

##         datasets$all_csd[[input$input2build]][[input$select_csd]] <- new_data

##     }, error = function(e){
##         showModal(dataModal(e[[1]]))
##     })
## })



## From the older implementation. Will remove this commented code.
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




##    observeEvent(input$dag_model, {
##     ## FIXME Unnecessary, as caught at a much more sensible place
##     ## but leave it here anyway, just in case, until much more testing done.
##     ## The other message is "This DAG has nodes with multiple parents. "
##     ## way below.
##     ## But here, it is often caught even earlier and we avoid the
##     ## flickering screen that happens when we left the error handling below.

##     number_of_parents <- colSums(data$dag)

##     if (input$dag_model == "OncoBN") {
##         if (any(data$DAG_parent_set == "XOR")) {
##             updateRadioButtons(session, "dag_model", selected = "HESBCN")
##             showModal(dataModal(HTML("The OncoBN model cannot include ",
##                                      "XOR relationships.")))
##         } else if (length(unique(data$DAG_parent_set)) > 2) {
##             updateRadioButtons(session, "dag_model", selected = "HESBCN")
##             showModal(dataModal(HTML("The OncoBN model can only include ",
##                                      "one type of relationship",
##                                      "(conjunctive or disjunctive, ",
##                                      "as specified in \"Advanced options\").")))
##         } else {
##             ## default_dag_model <<- input$dag_model
##             the_dag_model$stored_dag_model <- input$dag_model
##         }
##     } else if (input$dag_model == "OT") {
##         if (any(number_of_parents > 1)) {
##             updateRadioButtons(session, "dag_model", selected = "HESBCN")
##             showModal(dataModal(
##                 paste("This DAG has nodes with multiple parents. ",
##                       "OT can only use trees ",
##                       "(i.e., no node can have with multiple parents.)")))
##             ## paste("This DAG cannot be transformed into a tree. ",
##             ##   "Are there nodes with multiple parents? ",
##             ##   "(OT cannot not have nodes with multiple parents.)")))
##         } else if (length(unique(data$DAG_parent_set)) > 2) {
##             updateRadioButtons(session, "dag_model", selected = "HESBCN")
##             showModal(dataModal(HTML("The OT model  ",
##                                      "is only for trees. ")))
##         } else {
##             ## default_dag_model <<- input$dag_model
##             the_dag_model$stored_dag_model <- input$dag_model
##         }
##     } else {
##         ## default_dag_model <<- input$dag_model
##         the_dag_model$stored_dag_model <- input$dag_model 
##     }
##     the_dag_model$stored_dag_model <- input$dag_model 
## })



## Older version. Contains comments about A21_gnn_numfix
## display_freqs <- reactive({        
##     ## Remember this is called whenever changes in many places
##     ## happen.

##     ## This is often called when there is no need for it. So when you change
##     ## the type of data entering (go from MHN to upload, for example) this is
##     ## called again and returns the data but for nothing since what we will
##     ## want to display are the new data. We now elegantly return a
##     ## 0-rows data frame when nothing should be returned. Explicit and clean.

##     ## For paranoia, we return always things in standard order Probably
##     ## unnecessary, but just in case.
##     mymessage("At display_freqs")

##     ## provide_gene_names is being enabled somewhere I can't locate
##     ## So make sure we catch it right on the redisplay
##     ## There are a bunch of calls like this.
##     if ((!is.null(data$data) ||
##          (nrow(data$csd_counts) > 0))) {
##         mymessage("    disabled provide_gene_names under display_freqs")
##         shinyjs::disable("provide_gene_names")
##     }  

##     ## If no data to display, return empty data frame
##     thisd <- input$input2build

##     if (is.null(data$name) ) {
##         mymessage("       NULL data ",
##                   "Returning a 0-rows data frame")
##         return(data.frame(Genotype = character(), Counts = integer()))
##     } else {
##         thisd_dataset_names <- unlist(lapply(datasets$all_csd[[thisd]],
##                                              function(x) x$name))
##         if (!(data$name %in% thisd_dataset_names)) {
##             mymessage("       data$name not in ", thisd_dataset_names, ". ",
##                       "Returning a 0-rows data frame")
##             return(data.frame(Genotype = character(), Counts = integer()))
##         }
##     } 

##     if (input$input2build == "dag") {
##         ## With the DAG we always return all the genotypes
##         ## Other code ensures that gene number is never smaller
##         ## than genes in the DAG.
##         mymessage("      dag")

##         ## FIXME: See comment below  A1_gnn_bisfix
##         ## if (!is.null(data$n_genes) &&
##         ##     input$gene_number != data$n_genes) {
##         ##     mymessage("      DAG: Updating gene names in data")
##         ##     new_gnames2 <- set_gene_names_after_resize(data$data,
##         ##                                                data$gene_names)
##         ##     data$gene_names <- new_gnames2
##         ## }
##         return(
##             evamtools:::reorder_to_standard_order_count_df(
##                             data$csd_counts[data$csd_counts$Counts > 0, ,
##                                             drop = FALSE]))
##     } else if (input$input2build == "matrix") {
##         ## With the MHN we always return all the genotypes
##         ## The genotypes are given by the number of genes directly.
##         mymessage("      matrix (i.e., mhn)")
##         return(
##             evamtools:::reorder_to_standard_order_count_df(
##                             data$csd_counts[data$csd_counts$Counts > 0, ,
##                                             drop = FALSE]))
##     } else if (input$input2build == "upload") {
##         mymessage("      upload")
##         ## With upload, we do not use number of genes
##         ## Return the data we have
##         return(
##             evamtools:::reorder_to_standard_order_count_df(
##                             data$csd_counts[data$csd_counts$Counts > 0, ,
##                                             drop = FALSE]))
##     } else if (input$input2build == "csd") {
##         mymessage("      csd")

##         return(
##             evamtools:::reorder_to_standard_order_count_df(
##                             data$csd_counts[data$csd_counts$Counts > 0, ,
##                                             drop = FALSE]))

##         ## I think what follows is a relic from the past, when we could change
##         ## gene names on data sets with data already.
##         ## Also would play a role if we allowed decreasing number of genes
##         ## below those in use. Not anymore. Search for
##         ## csd_more_genes_than_set_genes
##         ## A21_gnn_numfix
##         ## I do not allow that
##         ## anymore. So all of this could just be the same as above, which
##         ## shows we always do the same thing


##         ## ## FIXME: A1_gnn_bisfix The code in A1_gnn_0 is not updating
##         ## ## data$etc. That is triggered on gene number change, but not
##         ## ## necessarily on adding a genotype or removing a genotype, in ways
##         ## ## that change the genes, but not the gene number.
##         ## ## Though something I do not fully understand:
##         ## ## It does not happen in the single observeEvent for input$gene_number
##         ## ## and not in the updateNumericInput for input$gene_number
##         ## ## So force it here if there have been changes in input$gene_number
##         ## ## and if they haven't, for some weird, hard to reproduce, cases.
##         ## ## Yes, this seem necessary to prevent BUG_Create_Rename_Click_other

##         ## if ((input$input2build == "csd") &&
##         ##     !is.null(data$n_genes) ) {
##         ##     ## input$gene_number != data$n_genes) {
##         ##     mymessage("       CSD: Updating gene names in data")
##         ##     new_gnames2 <- set_gene_names_after_resize(data$data,
##         ##                                                data$gene_names)
##         ##     if (identical(new_gnames2, data$gene_names)) {
##         ##         mymessage("A1_gnn_bisfix: identical")
##         ##     } else {
##         ##         mymessage("A1_gnn_bisfix: different")
##         ##     }
##         ##     data$gene_names <- new_gnames2
##         ## }

##         ## return(
##         ##     evamtools:::reorder_to_standard_order_count_df(
##         ##                     evamtools:::get_display_freqs(data$csd_counts,
##         ##                                                   input$gene_number,
##         ##                                                   data$gene_names,
##         ##                                                   input$input2build))
##         ## )
##     }
## })


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




## Disable because it led to inconsistent behavior such as this
## Sample from a DAG, say Fork, and analyze
## Go back, and sample, but now sample only 10 cases. Analyze
## Go to first analysis, click on "Modify data" and .. you
## are shown the 10 samples.
## Even worse, you can modify the model in between.
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



## Examples with prompter and shinyBS



## |> prompter::add_prompt(
##                  message =
##                      paste("The list of genes next to 'From' and 'To' is kept sorted ",
##                            "(alphabetically, often callen 'natural order'), with 'Root' first. ",
##                            " You can ",
##                            "see more genes (or fewer, up to the number of genes in your genotypes) ",
##                            "by moving the slider of 'Number of genes'."),
##                  position = "right",
##                  rounded = TRUE,
##                  bounce = TRUE,
##                  size = "medium")
## shinyBS::bsTooltip("dag_from",
##                    HTML("<p>The list of genes next to \"From\" and \"To\" is kept sorted ",
##                         "(alphabetically, often callen \"natural order\"), with \"Root\" first. ",
##                         " You can ",
##                         "see more (or fewer, up to the number of genes in your genotypes) ",
                                           ##                         "by moving the slider of \"Number of genes\".</p>")
                                           ##                  , "right", options = list(container = "body")
                                           ##                    ),
## Prompter is not opaque. Changing opacity possible?
## https://github.com/etiennebacher/prompter/issues/3
## But need to edit the CSS. PITA



## shinyBS::bsTooltip("data2plot",
##                    HTML("<p>This output is also displayed in tabular form on the bottom right.</p>",
##                         "<br><p><u>\"Sampled genotype counts\"</u> is only available if you selected ",
##                         "\"Sample genotypes\" under \"Advanced options\". </p>",
##                         "<br>",
##                         "<p>For <u>\"Predicted genotype relative frequencies\"</u> ",
##                         "and <u>\"Sampled genotype counts\"</u>, the histograms only ",
##                         "show the 20 most frequent genotypes ", ## 20: argument to  plot_genotype_counts
##                         "for reasons of figure size and legibility of genotype labels .",
##                         "The table shows all the genotypes. </p>",
##                         "<br>",
##                         "<p>For <u>\"Predicted genotype relative frequencies\"</u> ",
##                         "we add, to the table, the relative frequencies of genotypes in the original data ",
##                         "to make it easier to visually asses how close predictions are to observed data. ",
##                         "(It would not be sensible to do this with \"Sampled genotype counts\" as ",
##                         "those include additional sampling noise.) For easier comparison, both genotypes ",
##                         "with no observed counts in the original data ",
##                         "and genotypes with predicted frequency exactly 0 are left as empty (not as 0) ",
##                         "in the displayed table.</p>"
##                         ),
##                    "right", options = list(container = "body")
##                    )
## |> prompter::add_prompt(message = 
##                             HTML("This output is also displayed in tabular form on the bottom right.",
##                                  "                   ",
##                                  paste("<p>\"Sampled genotype counts\" is only available if you selected ",
##                                        "\"Sample genotypes\" under \"Advanced options\". </p>"),
##                                  "                   ",
##                                  paste("For \"Predicted genotype relative frequencies\" ",
##                                        "and \"Sampled genotype counts\", the histograms only ",
##                                        "show the 20 most frequent genotypes ", ## 20: argument to  plot_genotype_counts
##                                        "for reasons of figure size and legibility of genotype labels .",
##                                        "The table shows all the genotypes. "),
##                                  sep = "\n"
##                                  ),
##                         position = "right",
##                         rounded = TRUE,
##                         bounce = TRUE,
##                         size = "large"
##                         )



