## Copyright 2022 Pablo Herrera Nieto, Ramon Diaz-Uriarte

cpm_info <- function(){
  tags$div(id = "background",
    tags$head(
      tags$style(HTML("
        html, body{
          margin: 0;
        }
        #background{
          background-color: #f1f5f8;
          margin: 0;
          padding: 0;
          width: 100%;
          # height: 100vh;
          top: 0;
          left: 0;
          position: absolute;
        }

        #cpm_info{
          background: white;
          padding: 10px 50px;
          width: 75%;
          margin: auto;
          margin-top: 100px;
          border: 5px solid white;
          border-radius: 5px;
          box-shadow: 0 0 10px 5px rgba(160,160,160, 0.5);
        }

        #cpm_info hr{
          border-top: 3px solid #0892d0;
        }

        #cpm_info table{
          width: 50%;
          text-align: center;
        }

        #cpm_info h1{
          color: #0892d0;
          margin-top: 0;
          margin-bottom: 0;
        }

        #cpm_info td, #cpm_info th {
          border: 1px solid #ddd;
          padding: 8px;
        }

        #cpm_info tr:nth-child(even){background-color: #08B5FF;}

        #cpm_info th {
          padding-top: 12px;
          padding-bottom: 12px;
          text-align: center;
          background-color: #08B5FF;
          color: white;
        }
        .container-fluid{
          padding:0;
        }
        "
    ))),
    tags$div(id = "cpm_info",
      includeMarkdown("assets/landing_page.md")
    )
  )
}

results_simple <- function(){
    fluidPage(
        tags$head(
            tags$style(
                HTML("
                    .row{
                      width: 98%;
                      margin: auto;
                        # background-color:red;
                        # height: 45vh;
                        # overflow: auto;
                    }

                    .max_height{
                      max-height: 70vh;
                      overflow: auto;
                    }
                    #sims > .col-sm-3{
                      padding: 0;
                    }

                    #noprogress input.form-control{
                      visibility: hidden;
                    }

                    .form-control{
        height: 42px !important;
      }

      #freq2label-label{
        width: 0%;
      }

      #freq2label-wrap .form-group{
        display: block !important;
        margin-bottom: 0;
      }

       #freq2label-wrap .irs-single, .irs-min, .irs-max{
        visibility: visible !important;
      }

       #freq2label-wrap .irs-grid-text{
        visibility: hidden !important;
      }

      .col-sm-12{
        padding-left:5px;
        padding-right:5px;
      }
      .col-sm-10{
        padding-left:0px;
        padding-right:0px;
      }

      .col-sm-11{
        padding-right:0px;
      }

      .col-sm-2{
        padding-left:0px;
        padding-right:0px;
      }

      .col-sm-4{
        padding-left:0px;
        padding-right:0px;
      }

      .col-sm-7{
        padding-left:0px;
        padding-right:0px;
      }

      .col-sm-1{
        padding-left:0px;
      }

      .row{
        margin-right:0px;
      }

                   
                #select_cpm div.radio{
                  background-color: rgba(200,200,200, 0.5);
                  text-align: center;
                  border: 2px solid gray;
                  border-radius: 3px;
                  margin-top: 5px;
                  # margin-left: -20px;
                  white-space: nowrap;
                  overflow: hidden;
                  text-overflow: ellipsis;
                }

                @media only screen and (min-width: 1400px) {
                  #select_cpm div.radio{
                    max-width: 200px;
                  }
                  #select_cpm div.radio{
                    max-width: 200px;
                  }
                }

                @media only screen and (min-width: 1900px) {
                  #select_cpm div.radio{
                    max-width: 250px;
                  }
                  #select_cpm div.radio{
                    max-width: 250px;
                  }

                }

                #select_cpm div.radio:hover {
                  background-color: rgba(160,160,160, 0.5);
                  box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
                  cursor: pointer;
                }

                #select_cpm div.radio{
                  background-color: rgba(200,200,200, 0.5);
                  text-align: center;
                  border: 2px solid gray;
                  border-radius: 3px;
                  margin-top: 5px;
                  white-space: nowrap;
                  overflow: hidden;
                  text-overflow: ellipsis;

                    }
                    "
                )
            )
        ),
        tags$div(class = "row",
          # column(1,
          # column(11,
            column(2,
              tags$div(class = "frame",
                       tags$h3("Active data set"),
                uiOutput("cpm_list")
              ),
              uiOutput("customize"),
              tags$div(class = "frame",
                       tags$h3(HTML("Download CPM results <br> and analyzed data")),
                       tags$h5(HTML("Format and contents: rds file with ",
                                    "two lists: 1. cpm_output, ",
                                    "the concatenated output from ",
                                    "evam and sample_evam; 2. the tabular data ",
                                    "Analyzed data in ",
                                    " object$cpm_output$analyzed_data")),
                tags$div(class = "download_button",
                  downloadButton("download_cpm", "Download")
                )
                ),
              tags$div(
                       tags$h5(paste("evamtools version: ",
                          packageVersion("evamtools")))
                       )
              ),
            column(10,
              column(12, uiOutput("sims")),
              column(12, uiOutput("sims2"))
            )
            ,
            column(4,
                   ## FIXME zzply
                   uiOutput("original_data")
            ),
            column(6,
              uiOutput("tabular_data")
            )

          # )
        )
    )

}

user_input <- function() {
    fluidPage(
        prompter::use_prompt(),
        shinyjs::useShinyjs(),
        tags$head(
                 tags$style(HTML("
      body{
        font-size: 15px;
      }
      
      .irs-grid-pol.small{
        height: 0px !important; 
      }

      #define_genotype>*>*>label{
        margin-right: 10px;
        text-align: right;
        margin-bottom: 0;
      }
      input[type='checkbox']{
        zoom: 1.5;
        margin-right: 4px;
      }

      .checkbox{
        vertical-align: middle;
        margin-top: 10px !important; 
      }

      #genotype  input[type=checkbox] {
        zoom: 0;
      }

      #genotype>.shiny-options-group{
        display: flex;
        margin-right: 15px;
      }

      label{
        margin-right: 10px;
      }
      .frame .form-group{
        display: flex;
      }

      .form-control{
        height: 42px !important;
      }

      .form-group{
        display: flex;
        margin-bottom: 0;
      }
      .upload_file .form-group{
        width: 500px !important; 
      }
      input.formc-control{
        height: 40px !important;
      }

      #genotype>label{
        margin-right: 20px;
      }
      #add_genotype{
        width: 35%;
        height: 40px;
      }
      #fr{
        display: flex;
        justify-content: space-between;
      }

      #genotype_freq-label{
        margin-right: 20px;
      }

      div.inlin label { 
        width: 15%; 
      }

     div.inlin3 label { 
        width: 100%;
        text-align: left;
        vertical-align: middle; 
      }
      div.inlin2 label { 
        width: 100%; 
      }
      .inlin2>div{ 
        width: 100% !important; 
      }

      input{
        z-index: 100;
      }
      .inlin label{ 
        display: table-cell; 
        text-align: left; 
        vertical-align: middle; 
      }

      .inlin>.form-group>.irs--shiny.irs-with-grid{
        margin-left: 10%;
        width: 75%;
      }

      .btn, input.form-control{
        font-size: 15px;
      }
      .inlin .form-group { 
        display: table-row;
        font-size: 15px;
      }

      .irs > span{
        font-size: 15px;
      }

      span [class*='irs'] { 
        font-size: 15px !important; 
      }

      .irs-single, .irs-min, .irs-max{
        visibility: hidden !important;
      }

      .download_button{
        heigth: auto;
        display: block;
        display:flex;
        # height: 80px;
        margin-top: 15px;
        justify-content: center;
      }

      .submit_button>button{
        font-size: 30px !important;
        padding: 15px 30px;
        background: #FAA0A0;
      }

      .submit_button>button:hover{
        font-size: 30px !important;
        padding: 15px 30px;
        background: #F20D0D;
      }

      .shiny-download-link{
        margin: auto;
        # font-size: 18px;
      }
      .not_show{
        color: rgba(0,0,0,0);
      }

      .frame{
        margin-bottom:15px;
        padding: 15px;
        border: 3px solid rgba(100, 100, 100, 0.5);
        border-radius: 5px;
      }
      .upload_file{
        width: 100%;
        display: table-cell;
        # justify_content: center
      }
      h3{
        margin-top:0;
      }

      # #csd_table{
      #   height: 30vh;
      #   overflow: auto;
      # }
      .row{
        padding-top: 20px;
      }

      .flex>*{
        display: inline;
      }
      .flex>h3{
        width: 50%;
        margin-right: 0;
      }

      #move div.radio:hover {
        background-color: rgba(160,160,160, 0.5);
        box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
        cursor: pointer;
      }

      #select_csd div.radio{
        background-color: rgba(200,200,200, 0.5);
        text-align: center;
        border: 2px solid gray;
        border-radius: 3px;
        margin-left: -20px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }

      @media only screen and (min-width: 1200px) {
        #input2build div.radio{
          width: 140px;
        }
        #select_csd div.radio{
          width: 140px;
        }
        #select_cpm div.radio{
          width: 170px;
        }
      }

      @media only screen and (min-width: 1900px) {
        #input2build div.radio{
          width: 150px;
        }
        #select_csd div.radio{
          width: 150px;
        }
        #select_cpm div.radio{
          width: 140px;
        }
      }

      #input2build div.radio:hover {
        background-color: rgba(160,160,160, 0.5);
        box-shadow: 0 0 2px 2px rgba(160,160,160, 0.5);
        cursor: pointer;
      }

      #input2build div.radio{
        background-color: rgba(200,200,200, 0.5);
        text-align: center;
        border: 2px solid gray;
        border-radius: 3px;
        margin-left: -20px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      
      # #input2build input[type=radio]{
      #   visibility:hidden;
      # }

      # #csd_list input[type=radio]{
      #   visibility:hidden;
      # }

      # #cpm_list input[type=radio]{
      #   visibility:hidden;
      # }

      #input2build div label{
        padding-left: 20px;
      }
      #all_advanced_options{
        background-color:  #f1f5f8;
        margin-top: 20px;
        border: 5px solid white;
        border-radius: 5px;
        box-shadow: 0 0 10px 5px rgba(160,160,160, 0.5);
      }
      #all_advanced_options label{
        width: 50%;
      }

        ") # end HTML
     ) # end tags$style
     ),
     column(width = 12,
            
            sidebarLayout(
                column(width = 1,
                       tags$div(
                                tags$h3("Input to build"),
                                tagList(
                                    radioButtons(inputId = "input2build", label = "",
                                                 choiceNames = list(
                                                     HTML("Cross-sectional <br> data"),
                                                     HTML("DAG and <br> rates/probs."),
                                                     "MHN thetas"),
                                                 choiceValues = list("csd", "dag", "matrix"),
                                                 selected = "dag"
                                                 )
                                    |> prompter::add_prompt(
                                                     position = "bottom",
                                                   bounce = "TRUE",
                                                                     message = "a text for input2build from add_prompt")
                                    ## If we wanted shinyBS, comment above and use something like
                                    ## shinyBS::bsTooltip("input2build", "Working example of a tooltip on server.R",
                                    ##                    "right", options = list(container = "body")),
                                ),
                                tags$footer(tags$script(HTML("
              tmp_label = document.createElement('p');
              tmp_label.innerHTML = 'Raw data';
              document.querySelector('#input2build div.radio').before(tmp_label)
              tmp_label = document.createElement('p');
              tmp_label.innerHTML = 'CPM types';
              document.querySelector('#input2build div.radio').after(tmp_label)
              "))),
              tags$h4("Examples and user's data"),
              ## How many to show is controlled from server function
              ## in server.R, examples_csd$dag, etc
              uiOutput("csd_list"),

              tags$h3(HTML("<br/>")),
              tags$h5(paste("evamtools version: ",
                            packageVersion("evamtools")))
              
              )
              ## do it with a render UI
              ),
              

              
              column(width = 11,
                     titlePanel(HTML("&ensp; Cross-sectional data input")),
                     column(width = 6,
                            column(width = 12,
                                   ## Upload
                                   tags$div(class = "frame",
                                            tags$h3(" Define your data interactively"),
                                            ## Save/Download/Rename/Use
                                            tags$div(class = "frame",
                                                     tags$h3("(Re)name the data"),
                                                     tags$h5(HTML("Give the modified data a name ",
                                                                  "that will also be used to save the CPM ",
                                                                  "output.")),
                                                     tags$div(class = "download_button",
                                                              ),
                                                     uiOutput("dataset_name"),
                                                     actionButton("save_csd_data", "Use this name"),
                                                     ),
                                            
                                            tags$div(class = "frame",
                                                     uiOutput("upload_data"),
                                                     tags$div(class = "flex",
                                                              tags$h3("1. Set the number of genes"),
                                                              tags$h5("(Using 7 or more genes can lead ",
                                                                      "to very long execution times for some methods ",
                                                                      "and crowded figures.)"),
                                                              actionButton("change_gene_names", "Change gene names")
                                                              
                                                              ),
                                                     uiOutput("gene_number")),
                                            tags$div(class = "frame",
                                                     uiOutput("define_genotype"),
                                                     ),
                                            uiOutput("change_counts"),

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
                                            )

                                   )
                            ), 
                     column(width = 6,
                            
                            tags$div(class = "download_button submit_button",
                                     actionButton("analysis", "Run evamtools")
                                     ),
                            tags$div(class = "download_button",
                                     actionButton("advanced_options", "Advanced options and CPMs to use")
                                     ),
                            tags$div(id="all_advanced_options", 
                                     ## title = tags$h3("Advanced options"),
                                     tags$div(class="inlin",
                                              tags$h5(HTML("(See additional details for all options ",
                                                           "in the help of the <tt>evam</tt> and <tt>sample_evam</tt> functions ",
                                                           "available from the 'Package evamtools' ",
                                                           "help files in the ",
                                                           "<a href='https://rdiaz02.github.io/EvAM-Tools/pdfs/Additional_doc_all.pdf'>",
                                                           "Additional documentation</a>).")),
                                              tags$hr(style="border-color: darkgrey;"),
                                              checkboxGroupInput("cpm_methods",
                                                                 "CPMs to use", 
                                                                 width = "100%",
                                                                 choiceNames = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN",
                                                                     "MCCBN",
                                                                     "H-ESBCN"),
                                                                 choiceValues = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN",
                                                                     "MCCBN",
                                                                     "HESBCN"
                                                                 ),
                                                                 selected = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN"
                                                                 ),
                                                                 inline = FALSE),
                                              tags$h5("Beware: MCCBN may take hours to run. ",
                                                      "H-ESBCN often takes much longer than the ",
                                                      " remaining methods (often > 20 seconds)."),
                                              tags$hr(style="border-color: darkgrey;"),

                                              selectInput("return_paths_max",
                                                          "Return paths to maximum(a)",
                                                          c(TRUE, FALSE),
                                                          selected = FALSE),
                                              tags$h5("(Paths to the maximum/maxima and their ",
                                                      " probabilities. ",
                                                      "These are not part of the tabular output ",
                                                      "(because of their possibly huge number) ",
                                                      "but if requested are included in the result object ",
                                                      "you can download)"),
                                              tags$hr(style="border-color: darkgrey;"),
                                              selectInput("do_sampling", "Sample genotypes: ",
                                                          c(TRUE, FALSE),
                                                          selected = FALSE),
                                              tags$h5("Generate a finite sample of genotypes ",
                                                      "according to the predicted frequencies of ",
                                                      "the model."),
                                              ## selectInput("do_genotype_transitions",
                                              ##             "Sample for observed genotype transitions",
                                              ##             c("True" = TRUE, "False" = FALSE),
                                              ##             selected = FALSE
                                              ##             ),
                                              ## tags$h5("Obtain observed genotype transitions? ",
                                              ##         "Requires simulating sampling from the  ",
                                              ##         "continuous time Markov chain and is ",
                                              ##         "slower than simply obtaining a sample of genotypes. ",
                                              ##         "For this to have an effect, set ",
                                              ##         "'Sample genotypes' to true."),
                                              ##  This is removal_note_sogt_1
                                              ## The functionality is still present in the R package itself
                                              numericInput("sample_size", "Number of samples",
                                                           10000
                                                           ## next fails with shinytest, so give the number
                                                           ## .ev_SHINY_dflt$cpm_samples
                                                         , min = 0, max = 100000, step = 100,
                                                           width = "100%"),
                                              tags$h5("Number of genotypes to generate ",
                                                      "when generating a finite sample of genotypes ",
                                                      "according to the predicted frequencies of ",
                                                      "from model."),                     
                                              numericInput("sample_noise", "Observation noise", 0
                                                         , min = 0, max = 1, step = 0.1, width="100%"),
                                              tags$h5("If > 0, the proportion of observations ",
                                                      "in the sampled matrix with error ",
                                                      "(for instance, genotyping error). ",
                                                      "This proportion of observations will have 0s flipped ",
                                                      "to 1s, and 1s flipped to 0s."),                     
                                              
                                              tags$hr(style="border-color: darkgrey;"),
                                              
                                              tags$h4("MHN options"),
                                              numericInput("MHN_lambda", "Lambdas: ", NULL, min=0),
                                              tags$h5("Penalty term; default: 1/number of rows of data set. ",
                                                      "(Do not enter anything, unless you want to use a value ",
                                                      "different from the default)."),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("OT options"),
                                              selectInput("OT_with_error", "Return errors: ", c("True" = TRUE, 
                                                                                                "False" = FALSE), selected = "True"),
                                              tags$h5("For large models this may take quite some time"),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("CBN options"),
                                              selectInput("CBN_init_poset", "Initial poset: ",
                                                          c("OT" = "OT", "Linear" = "linear"), selected="OT"),
                                              numericInput("CBN_omp_threads", "Number CBN OMP threads: ",
                                                           1, min = 1, max = parallel::detectCores()),
                                              tags$h5("Number of OMP threads; large numbers do ",
                                                      "necessarily lead to faster computations."),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("HESBCN options"),
                                              numericInput("HESBCN_steps", "Steps: ", 100000, min=100, max=10000000),
                                              numericInput("HESBCN_seed", "Seed: ", NULL, min=0),
                                              selectInput("HESBCN_reg", "Regularization: ",
                                                          c("BIC" = "bic", "AIC" = "aic", "Loglik" = "loglik"),
                                                          selected = "BIC"),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("OncoBN options"),
                                              selectInput("OncoBN_model", "Model: ",
                                                          c("Disjunctive (DBN)" = "DBN",
                                                            "Conjunctive (CBN)" = "CBN"),
                                                          selected = "DBN"
                                                          ),
                                              selectInput("OncoBN_algorithm", "Algorithm: ",
                                                          c("Dynamic programming (DP)" = "DP",
                                                            "Genetic algorithm (GA)" = "GA"),
                                                          selected = "DP"),
                                              numericInput("OncoBN_k",
                                                           "k: In-degree bound on the estimated network: ",
                                                           3, min = 0),
                                              numericInput("OncoBN_epsilon", "Epsilon: ", NULL, min=0),
                                              tags$h5("Penalty term for mutations not conformig ",
                                                      "to estimated network.",
                                                      "Default is min(colMeans(data)/2). ",
                                                      "(Do not enter anything, unless you want to a value ",
                                                      "different from the default)."),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("MCCBN options"),
                                              selectInput("MCCBN_model", "Model: ",
                                                          c("OT-CBN" = "OT-CBN", "H-CBN2" = "H-CBN2"),
                                                          selected = "OT-CBN"),
                                              numericInput("MCCBN_L",
                                                           "L: Number of samples to be drawn from the proposal in the E-step: ", 100, min=0),
                                              selectInput("MCCBN_sampling", "Sampling: ", c("forward", "add-remove", "backward", "bernoulli", "pool"), selected = "forward"),
                                              numericInput("MCCBN_max_iter",
                                                           "max.iter: Maximum number of EM iterations: ", 100, min=0, max=100000),
                                              numericInput("MCCBN_update_step_size",
                                                           "update.step.size: Number of EM steps after which the number of samples, ‘L’, is doubled: ", 20L, min=0, max=100000),
                                              numericInput("MCCBN_tol",
                                                           "tol: Convergence tolerance: ", 0.001, min=0, max=100),
                                              numericInput("MCCBN_max_lambda_val",
                                                           "max.lambda.val: Upper bound on the value of the rate parameters: ", 1e6, min=0, max=1e9),
                                              numericInput("MCCBN_T0",
                                                           "T0: Initial value of the temperature: ", 50, min=0, max=1000),
                                              numericInput("MCCBN_adapt_rate",
                                                           "adap.rate: Constant adaptation rate: ", 0.3, min=0, max=1000),
                                              numericInput("MCCBN_acceptance_rate",
                                                           "acceptance.rate: Desirable acceptance rate: (if NULL, defaults to 1/number of mutations)",
                                                           NULL, min=0, max=1000),
                                              numericInput("MCCBN_step_size",
                                                           "step.size: Number of iterations after which the temperature should be updated: (if NULL, defaults to 50)", NULL, min=0, max=1000),
                                              numericInput("MCCBN_max_iter_asa", "max.iter.asa: Maximun number of iterations: ", 10000L, min=0, max=1000000L),
                                              numericInput("MCCBN_neighborhood_dist", "neighborhood.dist: Hamming distance between the observation and the samples generated by backward sampling: ", 1L, min=0, max=1000000L),
                                              selectInput("MCCBN_adaptive", "adaptive: Use an adaptive
                    annealing schedule?: ", c(TRUE, FALSE), selected = TRUE),
                    numericInput("MCCBN_seed", "Seed: ", NULL, min=0, width="50%")
                                        # tags$h4("DISCLAIMER: Both HyperTraps and MCCBN may take hours to run")
                                        # )
                    )
                                        # )
                                        # ),
                    )
                    ## FIXME zzply
                  , plotly::plotlyOutput("plot") ## and this calls output$plot
                    ## ,  plotOutput("plot")
                  , plotOutput("dag_plot")

                    )
                    )
            )
            )
  )
}
ui <- 
  navbarPage( 
    "", ## "Evamtools",
    id = "navbar",
    header = tags$head(
      tags$style(
        HTML(
          "
        h2{
          margin: 0;
          padding: 0;
        }

        .navbar{
          padding-left: 20px;
          margin-bottom: 0px;
          position: fixed;
          width: 100%;
        }

        .dropdown-menu > li > a{
          font-size: 15px;
          padding: 10px 10px;
        }

        body .tab-content{
          margin-top: 50px;
        }
        "
      ))),
    tabPanel("About EVAM-tools",
      cpm_info()
    ),
    tabPanel("User input", 
      value = "csd_builder",
      user_input()
    ),
    tabPanel("Results",
      value = "result_viewer",
      results_simple()
    )
)

