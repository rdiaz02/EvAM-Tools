## Copyright 2022, 2025 Pablo Herrera-Nieto, Ramon Diaz-Uriarte,
## Javier Pérez de Lema Díez

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


tutorial <- function(){
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
                 includeMarkdown("assets/tutorial.md")
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
                 #select_cpm .radio label {
                     white-space: normal !important;
                     overflow-wrap: break-word !important;
                     text-align: left !important;
                    }
                    "
      )
      )
      ),
      tags$div(class = "row",
               column(2,
                      tags$div(id="dummy_data_name",
                               class = "frame",
                               tags$h3("Data name"),
                               uiOutput("cpm_list"),
                               tippy::tippy_this("dummy_data_name",
                                                 paste("<span style='font-size:1.5em; text-align:left;'>",
                                                       "Analyzed data. ",
                                                       "Beware: we cannot guarantee you have ",
                                                       "not made additional changes ",
                                                       "to the data ",
                                                       "(model changes, genotype ",
                                                       "frequencies changes) ",
                                                       "in the corresponding 'User input' ",
                                                       " screen since ",
                                                       "you analyzed them.",
                                                       "<span>"),
                                                 arrow = TRUE, animation = "shift-toward"
                                                 ## , placement = "right"
                                                 ),
                               ),
                      uiOutput("customize"),
                      tags$div(class = "frame",
                               tags$h3(HTML("Download EvAM results <br> and analyzed data")),
                               actionButton("how2downloadcpm", "Help", class = "btn-info"),
                               tags$div(class = "download_button",
                                        downloadButton("download_cpm", "Download")
                                        )
                               ),
                      tags$div(
                               tags$h5(paste("evamtools R package version: ",
                                             packageVersion("evamtools"))) ,
                               ## tags$h5(paste("commit: ",
                               ##               this_string_to_be_replaced_by_git_hash
                               ##               ))
                               ## substr(system("git rev-parse HEAD", intern=TRUE), 1, 7)))

                           )
                      ),
               column(10,
                      column(12, uiOutput("sims")),
                      column(12, uiOutput("sims2"))
                      ),
                      column(10,
                          uiOutput("HyperTraPSSummary")
                        ),
                      column(10,
                             uiOutput("BML_bootstrap")
                        ),
               column(4,
                      ## FIXME zzply
                      uiOutput("original_data")
                      ),
               column(6,
                      uiOutput("tabular_data")
                      ),
               )
    )
}

## Avoid repeating the same structure. This combines the input function
## with helper text not in bold. But we need to deal with the annoying
## arrows on top of short labels.
## inputWithHelper <- function(inputFunc, inputId, label,  helperText,
##                             ...) {
##     helperTag = tags$h5
##     helperStyle = "margin-top: 0.5em; margin-bottom: 1em;"

##     tagList(
##         inputFunc(inputId, label, ...),
##         helperTag(helperText) ##, style = helperStyle)
##     )
## }

## This is both a sensible function and a kludge. Avoid repeating the
## inputWhatever and h5, (that is good). But selectInput boxes
## can be tiny if label is tiny, so increase label artificially (kludge)
inputWithHelper <- function(inputFunc, inputId, label,
                            helperText, ...) {

    ## Add non-breaking spaces to ensure minimum label width
    ## This forces the input control to have adequate space
    if (deparse(substitute(inputFunc)) == "selectInput") {
        min_length <- 15
        label_length <- nchar(label)
        if (label_length < min_length) {
            padding_spaces <- paste(rep("\u00A0", min_length - label_length),
                                    collapse = "")
            label <- paste0(label, padding_spaces)
        }
    }

    ## Create the input with padded label
    input_el <- inputFunc(inputId, label, ...)

    ## Return both elements
    tagList(
        input_el,
        tags$h5(helperText)
    )
}


inputWithHelperUnbold <- function(inputFunc, inputId, label,
                                  label_unbold = "",
                                  helperText = "", ...) {

    ## Add non-breaking spaces to ensure minimum label width
    ## This forces the input control to have adequate space
    if (deparse(substitute(inputFunc)) == "selectInput") {
        min_length <- 15
        label_length <- nchar(label) + nchar(label_unbold)
        if (label_length < min_length) {
            padding_spaces <- paste(rep("\u00A0", min_length - label_length),
                                    collapse = "")
            label_unbold <- paste0(label_unbold, padding_spaces)
        }
    }

    mixedLabel <- HTML(paste0(
        label, " ",
        "<span style='font-weight: normal;'>", label_unbold, "\u00A0 </span>"))


    ## Create the input with padded label
    input_el <- inputFunc(inputId, mixedLabel,  ...)

    width_hack <- TRUE
    if (!width_hack) {
        ## Return both elements
        tagList(
            input_el,
            tags$h5(helperText)
        )
    } else {
        uniqueID <- paste0(inputId, "-container")
        ## Add custom CSS for this specific input

        ## Superkludge: for the few cases where the above ain't enough
        if (inputId %in% c("MCCBN_model", "MCCBN_sampling", "MCCBN_adaptive")) {
            select_box_width <- 110
        } else {
            select_box_width <- 60
        }

        customCSS <- tags$style(HTML(paste0(
                              "#", uniqueID, " .selectize-input, #", uniqueID, " select {",
                              "  min-width: ", select_box_width, "px !important;",
                              "  width: ", select_box_width,"px !important;",
                              "}"
                          )))

        tags$div(
                 id = uniqueID,
                 customCSS,
                 input_el,
                 tags$h5(helperText)
             )
    }
}


user_input <- function() {
    fluidPage(
        ## require(shinyBS),
        ## prompter::use_prompt(),
        shinyjs::useShinyjs(),
        ## tippy::use_tippy(),
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
                                tags$h3(HTML("<h3 style=\"margin-left:-5px\">Cross-sectional data. Upload, create, generate, modify: </h3>")),
                                ## tags$h3(HTML('<hr style="height:4px;background-color:black;width:100%; text-align:left;margin-left:0">')),
                                tags$h4(HTML("<hr style=\"height:1px; width:70px; background-color:black;text-align:left;margin-left:5px\">")),
                                tags$h4(HTML("<h4 style=\"margin-left:-20px\"> Enter<br>cross-sectional data: </h4>")),
                                tagList(
                                    radioButtons(inputId = "input2build", label = "",
                                                 choiceNames = list(
                                                     HTML("Upload file"),
                                                     HTML("Enter<br>genotype <br> frequencies <br> manually"),
                                                     HTML("DAG and<br>rates/cond.<br> probs."),
                                                     HTML("MHN <br>log-&Theta;<br> matrix")
                                                 ),
                                                 choiceValues = list("upload", "csd", "dag", "matrix"),
                                                 selected = "dag"
                                                 )
                                    ## |> prompter::add_prompt(
                                    ##                  position = "bottom",
                                    ##                bounce = "TRUE",
                                    ##                                  message = "a text for input2build from add_prompt")
                                            ## If we wanted shinyBS, comment above and use something like
                                            ## shinyBS::bsTooltip("input2build", "Working example of a tooltip on server.R",
                                            ##                    "right", options = list(container = "body")),
                                        ),

                                        tags$footer(tags$script(HTML("
              tmp_label = document.createElement('p');
              tmp_label.innerHTML = '<hr style=\"height:1px; width:70px; background-color:black;text-align:left;margin-left:-10px\"> <h4 style=\"margin-left:-30px\">Generate<br>cross-sectional data from EvAM models:</h4>';
              document.querySelectorAll('#input2build div.radio')[1].after(tmp_label)
              document.querySelectorAll('#input2build div.radio')[1].after(tmp_label)
              "))),
              tags$h3(HTML('<hr style="height:4px;background-color:black;width:100%;margin-left:-20px">')),
              tags$h4(HTML("<br/>")),
              tags$h4(HTML("Examples and user's data:")),
              tags$h5(HTML("<br/>")),
              ## How many to show is controlled from server function
              ## in server.R, examples_csd$dag, etc
              uiOutput("csd_list"),

              tags$h3(HTML("<br/>")),
              tags$h3(HTML('<hr style="height:1px;background-color:black;margin-left:-20px">')),
              tags$h5(paste("evamtools R package version: ",
                            packageVersion("evamtools"))),
              ## tags$h5(paste("commit: ",
              ##               this_string_to_be_replaced_by_git_hash
              ##               ))
              ## substr(system("git rev-parse HEAD", intern=TRUE), 1, 7)))

              )
              ## do it with a render UI
              ),



              column(width = 11,
                                        #  titlePanel(HTML("&ensp; Cross-sectional data input")),
                     column(width = 6,
                            column(width = 12,
                                   uiOutput("gene_number_slider"),
                                   tags$div(class = "frame",
                                            uiOutput("define_genotype"),
                                            ),
                                   uiOutput("change_counts"),
                                   uiOutput("dataset_name"),
                                   uiOutput("download_data")
                                   )
                            ),
                     column(width = 6,

                            tags$div(class = "download_button submit_button",
                                     actionButton("analysis", "Run evamtools")
                                     ),
                            tippy::tippy_this("analysis",
                                              paste("Inactive unless data available. ",
                                                    "Even if the button is active, ",
                                                    "the run will be aborted if the data ",
                                                    "contain fewer than two genes or two genotypes."),
                                              arrow = TRUE, animation = "shift-toward"
                                              ),
                            ## |> prompter::add_prompt(message = paste("Inactive unless data available. ",
                            ##                                         "Even if the button is active, ",
                            ##                                         "the run will be aborted if the data ",
                            ##                                         "contain fewer than two genes or two genotypes."),
                            ##                        ,position = "bottom",
                            ##                  rounded = TRUE,
                            ##                  bounce = TRUE,
                            ##                  size = "medium"),
                            tags$div(class = "download_button",
                                     actionButton("advanced_options", "Advanced options and EvAMs to use")
                                     ),
                            tags$div(id="all_advanced_options",
                                     tags$div(class="inlin",
                                              tags$h5(HTML("(See additional details for all options ",
                                                           "in the help of the <tt>evam</tt> and <tt>sample_evam</tt> functions ",
                                                           "available from the 'Package evamtools' ",
                                                           "help files in the ",
                                                           "<a href='https://rdiaz02.github.io/EvAM-Tools/pdfs/Additional_doc_all.pdf'>",
                                                           "Additional documentation</a>).")),
                                              tags$hr(style="border-color: darkgrey;"),
                                              checkboxGroupInput("cpm_methods",
                                                                 "EvAMs to use",
                                                                 width = "100%",
                                                                 choiceNames = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN",
                                                                     "MCCBN",
                                                                     "H-ESBCN",
                                                                     "HyperTraPS",
                                                                     "BML"),
                                                                 choiceValues = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN",
                                                                     "MCCBN",
                                                                     "HESBCN",
                                                                     "HyperTraPS",
                                                                     "BML"
                                                                 ),
                                                                 selected = c(
                                                                     "CBN",
                                                                     "OT",
                                                                     "OncoBN",
                                                                     "MHN",
                                                                     "HyperTraPS"
                                                                 ),
                                                                 inline = FALSE),
                                              tags$h5("Beware: MCCBN may take hours to run. ",
                                                      "H-ESBCN often takes longer than the ",
                                                      "remaining methods (except MCCBN) for small ",
                                                      "numbers of genes (5 or less).",
                                                      "For 7 or more genes, CBN can be much slower ",
                                                      "than OT, OncoBN, or MHN (e.g., data analyzed in < 1 second ",
                                                      "by those three methods can take 45 with CBN),  ",
                                                      "and also often slower than H-ESBCN. ",
                                                      "BML's execution time increases with sample ",
                                                      "size; for real, you will want to ",
                                                      "run bootstrap analysis, but it is disabled ",
                                                      "(number of bootstrap replicates set to 0) ",
                                                      "because it can be very slow with large sample sizes; ",
                                                      "experiment before launching many bootstrap iterations. ",
                                                      "HyperTraPS-CT, especially depending on the options, ",
                                                      "can also take a long time to finish; experiment ",
                                                      "especially with arguments 'Number of walkers' ",
                                                      "'Inference chain length' and 'Model structure'."),
                                              tags$hr(style="border-color: darkgrey;"),

                                           inputWithHelper(numericInput,
                                                           "evam_run_num_cores", "Number of cores",
                                                           paste("Number of cores used for the evam call. Using ",
                                                                 "more than the number of methods selected is pointless. ",
                                                                 "Setting it to 1 will run methods sequentially instead of ",
                                                                 "in parallel."),
                                                           8, min = 1, max = 10, step = 1),
                                           ## numericInput("evam_run_num_cores",
                                           ##              "Number of cores", 8, min = 1, max = 10,
                                           ##              step = 1),
                                           ## tags$h5("Number of cores used for the evam call. Using ",
                                           ##         "more than the number of methods selected is pointless. ",
                                           ##          "Setting it to 1 will run methods sequentially instead of ",
                                           ##         "in parallel."),
                                           tags$hr(style="border-color: darkgrey;"),
                                              selectInput("return_paths_max",
                                                          "Return paths to maximum(a)",
                                                          c(TRUE, FALSE),
                                                          selected = FALSE),
                                              tags$h5("(Paths to the maximum/maxima and their ",
                                                      " probabilities. ",
                                                      "These are not part of the tabular or graphical output ",
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
                                              numericInput("MHN_lambda", "Lambda: ", NULL, min=0),
                                              tags$h5("Lambda: penalty term in fitting algorithm. Default = 1/number of rows of data set. ",
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
                                              tags$h4("H-ESBCN options"),
                                              numericInput("HESBCN_MCMC_iter", "Number of MCMC iterations: ",
                                                           200000, min=100, max=10000000),
                                              tags$h5("Number of MCMC iterations: ",
                                                      "Argument '-n | --number_samples' in the H-ESBCN C code. ",
                                                      "EvAM's web app default is 200000, larger than the ",
                                                      "original default of 100000. ",
                                                      "You might want to increase it to 500000 or 1000000 ",
                                                      "being aware that this will result in longer running times. "),
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
                                              tags$h5("Epsilon: Penalty term for mutations not conformig ",
                                                      "to estimated network.",
                                                      "Default is min(colMeans(data)/2). ",
                                                      "(Do not enter anything, unless you want to a value ",
                                                      "different from the default)."),
                                              tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("MCCBN options"),
                                           inputWithHelperUnbold(selectInput, "MCCBN_model", "Model: ", "", "",
                                                          c("OT-CBN" = "OT-CBN", "H-CBN2" = "H-CBN2"),
                                                          selected = "OT-CBN"),
                                           inputWithHelperUnbold(numericInput, "MCCBN_L",
                                                                 "L", " (number of samples to be drawn from the proposal in the E-step: )", "", 100, min=0),
                                           inputWithHelperUnbold(selectInput, "MCCBN_sampling", "Sampling: ", "", "", c("forward", "add-remove", "backward", "bernoulli", "pool"), selected = "forward"),
                                           inputWithHelperUnbold(numericInput, "MCCBN_max_iter",
                                                                 "max.iter", " (maximum number of EM iterations)", "", 100, min=0, max=100000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_update_step_size",
                                                                 "update.step.size ", "(number of EM steps after which the number of samples, ‘L’, is doubled): ", "", 20L, min=0, max=100000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_tol",
                                                                 "tol", " (convergence tolerance):", "", 0.001, min=0, max=100),
                                           inputWithHelperUnbold(numericInput, "MCCBN_max_lambda_val",
                                                                 "max.lambda.val", " (upper bound on the value of the rate parameters):", "", 1e6, min=0, max=1e9),
                                           inputWithHelperUnbold(numericInput, "MCCBN_T0",
                                                                 "T0", " (initial value of the temperature):", "", 50, min=0, max=1000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_adapt_rate",
                                                                 "adap.rate", " (constant adaptation rate):", "", 0.3, min=0, max=1000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_acceptance_rate",
                                                                 "acceptance.rate", " (desirable acceptance rate: (if NULL, defaults to 1/number of mutations):", "",
                                                           NULL, min=0, max=1000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_step_size", "step.size", " (number of iterations after which the temperature should be updated; if NULL, defaults to 50):", "", NULL, min=0, max=1000),
                                           inputWithHelperUnbold(numericInput, "MCCBN_max_iter_asa", "max.iter.asa", " (maximun number of iterations):", "", 10000L, min=0, max=1000000L),
                                           inputWithHelperUnbold(numericInput, "MCCBN_neighborhood_dist", "neighborhood.dist", " (hamming distance between the observation and the samples generated by backward sampling):", "", 1L, min=0, max=1000000L),
                                           inputWithHelperUnbold(selectInput, "MCCBN_adaptive", "adaptive", " (use an adaptive annealing schedule?):", "", c(TRUE, FALSE), selected = TRUE),
                                                numericInput("MCCBN_seed", "Seed: ", NULL, min=0, width="50%"),
                                              tags$hr(style="border-color: darkgrey;"),

                                              tags$h4("HyperTraPS options"),

                                           inputWithHelperUnbold(selectInput, "HyperTraPS_model", "model structure",
                                                                 "(-1: full ---all edges---; 1: main effects, no interaction; 2: pairwise; 3: 3-way; 4: 4-way):", "", choices = c(-1, 1:4), selected = 2),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_walkers", "walkers ",
                                                                 "(number of walkers): ", "", 200, min=0),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_length", "length ",
                                                                 "(inference chain length): ", "", 3, min=0),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_kernel", "kernel ",
                                                                 "(perturbation kernel): ", "", 5, min=0),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_seed", "seed ",
                                                                 "(random number seed; a value of -1 means use a random seed): ",
                                                                 "", -1, min=-1),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_losses", "losses ",
                                                                 "(gains = 0, losses = 1): ", "",
                                                                 choices = c(0, 1), selected = 0),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_apm", "apm ",
                                                                 "(use APM; No = 0, Yes = 1): ", "",
                                                                 c(0, 1), selected = 0),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_sa", "sa    ",
                                                                 "(use SA: 0/1):", "", c(0, 1), selected = 0),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_sgd", "sgd ",
                                                                 "(use SGD 0/1):", "", c(0, 1), selected = 0),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_pli", "pli ",
                                                                 "(use PLI 0/1):", "", c(0, 1), selected = 0),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_samplegap", "samplegap ",
                                                                 "(gap between samples):", "", 10, min = -1),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_penalty", "penalty ",
                                                                 "(penalty for regularisation by penalised likelihood):", "", 0,
                                                                 min=0, step = 0.01),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_LASSO", "lasso ",
                                                                 "(regularisation by LASSO):", "", 0, min=0,
                                                                 step = 0.01),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_regul", "regularise ",
                                                                 "(stepwise regularise model after best parameterisation; enabling it can considerably increase running time):", "", c(0, 1), selected = 0),
                                           inputWithHelperUnbold(numericInput, "HyperTraPS_samprow", "samples_per_row ",
                                                                 "(number of simulations per parameter sample):", "", 10, min=1),
                                           inputWithHelperUnbold(selectInput, "HyperTraPS_outtrans", "output_transitions ",
                                                                 "(output exact transitions):", "", c(0, 1), selected = 1),
                                           ## inputWithHelperUnbold(numericInput, "HyperTraPS_nsampl", "nsample ",
                                           ##                       "(number of samples from an exponential distribution to estimate the predicted genotype frequencies; make it larger if you want more accurate estimates):", "", 1000, min=10),
                                           ## inputWithHelperUnbold(numericInput, "HyperTraPS_cores", "cores ",
                                           ##                       "(number of cores to use when estimating the predicted genotype frequencies):", "", 10, min=1, max = 15, step = 1),
                                           tags$hr(style="border-color: darkgrey;"),
                                              tags$h4("BML options"),
                                           inputWithHelper(numericInput,
                                                           "BML_ntree", "ntree: ",
                                                           "Number of random restarts for searching the tree space.",
                                                           100, min=0),
                                           inputWithHelper(numericInput, "BML_threshold", "threshold:      ",
                                                           "Threshold for inferring paths.", 0.3, min=0, max=1, step = 0.1),
                                           inputWithHelper(numericInput,
                                                           "BML_rep", "Bootstrap replicates: ",
                                                           "Number of bootstrap replicates, if nrep = 0 no bootstrap will be performed. The default, 10, is way too small for real use.",
                                                           10, min=0),
                    ),

                    )
                  , plotly::plotlyOutput("plot") ## and this calls output$plot
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
        tabPanel("About EvAM-Tools",
                 cpm_info()
                 ),
        tabPanel("User input",
                 value = "csd_builder",
                 user_input()
                 ),
        tabPanel("Results",
                 value = "result_viewer",
                 results_simple()
                 ),
        ## tabPanel("Tutorial",
        ##   tutorial()
        ## )
        )

## |> prompter::add_prompt(message =
##                             paste("Analyzed data. ",
##                                   "Beware: we cannot guarantee you have ",
##                                   "not made additional changes ",
##                                   "to the data ",
##                                   "(model changes, genotype ",
##                                   "frequencies changes) ",
##                                   "in the corresponding \'User input\' ",
##                                   " screen since ",
##                                   "you analyzed them."),
##                         position = "bottom",
##                         rounded = TRUE,
##                         bounce = TRUE,
##                         size = "medium")
