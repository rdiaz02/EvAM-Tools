

## Read Random Data

results_simple <- function(){
    fluidPage(
        # tags$head(
        #     tags$style(
        #         HTML("
        #             .row{
        #                 background-color:red;
                        
        #             }
        #             .download_button{
        #                 heigth: auto;
        #                 # display: block;
        #                 display:flex;
        #                 # height: 80px;
        #                 margin-top: 15px;
        #                 justify-content: center;
        #             }
        #             #modify_data{
        #                 margin: auto;
        #             }
        #             "
        #         )
        #     )
        # ),
        tags$div(class = "row",
        column(4,
            "Here I put the model"
        ),
        column(4,
            "Here comes the plot"
        ),
        column(4,
            plotOutput("csd"),
            tags$div(class = "download_button submit_button",
                actionButton("modify_data", "Modify data")
            )
        )
        )
    )

}