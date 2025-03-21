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