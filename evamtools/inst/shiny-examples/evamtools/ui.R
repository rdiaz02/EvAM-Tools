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
                tags$h3("Outputs"),
                uiOutput("cpm_list")
              ),
              uiOutput("customize"),
              tags$div(class = "frame",
                tags$h3("Download"),
                # tags$div(id = "noprogress",
                # fileInput("output_cpms", "Load your results"
                #   , multiple = FALSE,
                #   accept = c(".Rdata", ".rds", ".RDS"))
                # ),
                tags$div(class = "download_button",
                  downloadButton("download_cpm", "Download")
                )
              )
            ),
            column(10,
              column(12, uiOutput("sims")),
              column(12, uiOutput("sims2"))
            )
            ,
            column(4,
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

      #csd_table{
        height: 40vh;
        overflow: auto;
      }
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
    column(width=12,
      
      sidebarLayout(
        column(width = 1,
          tags$div(
            tags$h3("Input to build"),
            tagList(
              radioButtons(inputId = "input2build", label = "", 
                choiceNames = c("Cross sectional data", "DAG builder", "Matrix Builder"),
                choiceValues = c("csd", "dag", "matrix"),
                selected = "csd"
              )
            ),
            tags$h3("Some examples"),
            uiOutput("csd_list")
          )
          ## do it with a render UI
        ),

        column(width = 11,
          titlePanel("Cross sectional data input"),
          column(width = 6,
              column(width = 12,
            tags$div(class = "frame",
                tags$div(class = "flex",
                    tags$h3("1. Set the number of genes"),
                    actionButton("change_gene_names", "Change gene names")
                  ),
                uiOutput("genes_number")),
            tags$div(class = "frame",
              uiOutput("define_genotype"),
            ),
            uiOutput("change_freqs"),

            tags$div(class = "frame",
                     tags$h3("Upload your own data"),
                     tags$h4(paste0("Format: csv ---comma separated values---,",
                                    " with first row with gene names."
                                    )),
                     tags$h5(paste0("Use only alphanumeric characters ",
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
                                            ## ".rds", ".RDS",
                                            "text/csv",
                                            ".csv")) 
            )),

            tags$div(class = "frame",
              tags$h3("Save & Download data"),
              uiOutput("dataset_name"),
              tags$div(
                tags$div(class = "download_button",
                  shinyjs::disabled(actionButton("save_csd_data", "Save Data")),
                ),
                tags$div(class = "download_button",
                  downloadButton("download_csd", "Download your data")
                )
              )
            )
            )
          ),
          column(width = 6,
          
            tags$div(class = "download_button submit_button",
              actionButton("analysis", "Run evamtools!")
            ),
            tags$div(class = "download_button",
              actionButton("advanced_options", "Advanced Options")
            ),
            tags$div(id="all_advanced_options", 
              # title = tags$h3("Advanced options"),
              # tags$div(
                numericInput("num_steps", "Sampling steps", SHINY_DEFAULTS$cpm_samples
                  , min = 0, max = 100000, step = 100, width="100%"),
                numericInput("sample_noise", "Sampling noise", 0
                  , min = 0, max = 1, step = 0.1, width="100%"),
                # checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("HyperTRAPS", "MCCBN"), choiceValues = c("hypertraps", "mccbn")),
                checkboxGroupInput("more_cpms", "Additional CPMs", width = "100%", choiceNames = c("MCCBN"), choiceValues = c("MCCBN")),
                tags$h4("DISCLAIMER: MCCBN may take hours to run"),
                tags$hr(),
                tags$div(class="inlin",
                tags$h4("MHN options"),
                numericInput("MHN_lambda", "Lambdas: ", NULL, min=0),
                tags$hr(),
                tags$h4("OT options"),
                selectInput("OT_with_error", "Return errors: ", c("True" = TRUE, 
                  "False" = FALSE), selected = "True"),
                tags$h5("For large models this may take quite some time"),
                tags$hr(),
                tags$h4("CBN options"),
                selectInput("CBN_init_poset", "Initial poset: ", c("OT" = "OT", "Linear" = "linear"), selected="OT"),
                tags$hr(),
                tags$h4("HESBCN options"),
                numericInput("HESBCN_steps", "Steps: ", 100000, min=100, max=10000000),
                numericInput("HESBCN_seed", "Seed: ", NULL, min=0),
                selectInput("HESBCN_reg", "Regularization: ", c("BIC" = "bic", "AIC" = "aic", "Loglik" = "loglik"), selected = "BIC"),
                tags$hr(),
                tags$h4("OncoBN options"),
                selectInput("OncoBN_model", "Model: ", c("Conjunctive" = "CBN", "Disjunctive" = "DBN"), selected = "Disjunctive"),
                selectInput("OncoBN_algorithm", "Algorithm: ", c("Dynamic programming" = "DP", "Genetic algorithm" = "GA"), selected = "Dynamic programming"),
                numericInput("OncoBN_k", "In-degree bound on the estimated network: ", 3, min=0),
                numericInput("OncoBN_epsilon", "Penalty: ", NULL, min=0),
                tags$hr(),
                  tags$h4("MCCBN options"),
                  selectInput("MCCBN_model", "Model: ", c("OT-CBN" = "OT-CBN", "H-CBN2" = "H-CBN2"), selected = "OT-CBN"),
                  numericInput("MCCBN_L", "Number of samples to be drawn from the proposal in the E-step: ", 100, min=0),
                  selectInput("MCCBN_sampling", "Sampling: ", c("forward", "add-remove", "backward", "bernoulli", "pool"), selected = "forward"),
                  numericInput("MCCBN_max_iter", "Maximum number of EM iterations: ", 100, min=0, max=100000),
                  numericInput("MCCBN_update_step_size", "Number of EM steps after which the number of samples, ‘L’, is doubled: ", 20L, min=0, max=100000),
                  numericInput("MCCBN_tol", "Convergence tolerance: ", 0.001, min=0, max=100),
                  numericInput("MCCBN_max_lambda_val", "Upper bound on the value of the rate parameters: ", 1e6, min=0, max=1e9),
                  numericInput("MCCBN_T0", "Initial value of the temperature: ", 50, min=0, max=1000),
                  numericInput("MCCBN_adapt_rate", "Constant adaptation rate: ", 0.3, min=0, max=1000),
                  numericInput("MCCBN_acceptance_rate", "Desirable acceptance rate: ", NULL, min=0, max=1000),
                  numericInput("MCCBN_step_size", "Number of iterations after which the temperature should be updated: ", NULL, min=0, max=1000),
                  numericInput("MCCBN_max_iter_asa", "Maximun number of iterations: ", 10000L, min=0, max=1000000L),
                  numericInput("MCCBN_neighborhood_dist", "Hamming distance between the observation and the samples generated by backward sampling: ", 1L, min=0, max=1000000L),
                  selectInput("MCCBN_adaptive", "Use an adaptive
                    annealing schedule?: ", c(TRUE, FALSE), selected = TRUE),
                  numericInput("MCCBN_seed", "Seed: ", NULL, min=0, width="50%")
                  # tags$h4("DISCLAIMER: Both HyperTraps and MCCBN may take hours to run")
                  # )
                )
              # )
            # ),
            ),
            plotOutput("plot")
            ,
            plotOutput("dag_plot")

          )
        )
      )
      )
  )
}
ui <- 
  navbarPage( 
    "Evamtools",
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

