user_input <- function(){
  fluidPage(
    tags$head(
      tags$style(HTML("
      body{
        font-size: 20px;
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

      #genotype>.shiny-options-group>*{
        margin-right: 15px;
      }

      #genotype>.shiny-options-group{
        display: flex;
        margin-right: 15px;
      }

      label{
        margin-right: 10px;
      }
      # #genotype{
      #   display: flex;
      # }
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

      div#inlin label { 
        width: 15%; 
      }
      
      #inlin label{ 
        display: table-cell; 
        text-align: left; 
        vertical-align: middle; 
      }

      #inlin>.form-group>.irs--shiny.irs-with-grid{
        margin-left: 10%;
        width: 75%;
      }

      .btn, input.form-control{
        font-size: 20px;
      }
      #inlin .form-group { 
        display: table-row;
        font-size: 20px;
      }

      .irs > span{
        font-size: 20px;
      }

      span [class*='irs'] { 
        font-size: 20px !important; 
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

      .shiny-download-link{
        margin: auto;
        # font-size: 18px;
      }
      .not_show{
        color: rgba(0,0,0,0);
      }

      .frame{
        margin-bottom:15px;
        padding: 20px;
        border: 5px solid rgba(100, 100, 100, 0.5);
        border-radius: 5px;
        # max-height: 40vh;
        # overflow: auto;
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
        ") # end HTML
      ) # end tags$style
    ),
    column(width=12,
      
      sidebarLayout(
        column(width = 1,
          "load datasets",
          "WIP"
          ## do it with a render UI
        ),

        column(width = 11,
          titlePanel("Cross sectional data input"),
          column(width = 6,
              column(width = 12,
            tags$div(class = "frame",
                tags$h3("1. Set the number of genes"),
                uiOutput("genes_number")),
            tags$div(class = "frame",
              tags$h3("2. Add new genotypes"),
              uiOutput("define_genotype"),
            ),
            tags$div(class = "frame",
              tags$div(class = "flex",
                tags$h3("3. Change frequencies"),
                actionButton("display_help", "Help"),
              ),
              tags$div(id = "csd_table", 
                DTOutput("csd_freqs")
                )
              )
            )
          ),
          column(width = 6,
          tags$div(class = "frame",
              tags$h3("... or load your own data"),
              tags$div(class = "upload_file",
                  fileInput("csd", "Load a CSV",
                    multiple = FALSE,
                    accept = c("text/csv",
                                ".csv")) 
            )),
            plotOutput("plot"),
            tags$div(
              tags$label(class="not_show", "Download your data"), 
              tags$div(class = "download_button submit_button",
                actionButton("analysis", "Run guloMAM!")
              ),
              tags$div(class = "download_button",
                downloadButton("download", "Download your data")
              )
            )
          )
        )
      )
      )
  )
}