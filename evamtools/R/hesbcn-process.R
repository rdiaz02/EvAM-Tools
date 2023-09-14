## Copyright 2021, 2022 Pablo Herrera Nieto, Ramon Diaz-Uriarte

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


## Testing HESBCN programs available
.._EvamTools_test.hesbcn <- Sys.which("h-esbcn") == ""
if(.._EvamTools_test.hesbcn) {
    warning(paste(
        "\n\n",
        "\n******************************************************",
        "\n******************************************************\n",
        "\n  WARNING \n",
        "The external program h-esbcn was not found.",
        "You will not be able to use HESBCN.",
        "You can it from https://github.com/danro9685/HESBCN",
        "See further instructions in the README.",
        "\n******************************************************",
        "\n******************************************************",
        "\n\n"
        )
            )
}
rm(.._EvamTools_test.hesbcn)


do_HESBCN <- function(data,
                      MCMC_iter = 100000,
                      reg = c("bic", "aic", "loglik"),
    seed = NULL,
    tmp_dir = NULL,
    addname = NULL,
    silent = TRUE) {

    reg <- match.arg(reg)
    # Setting tmp folder
    if (is.null(tmp_dir)) {
        tmp_dir <- tempfile()
        dirname0 <- NULL
        if (!is.null(addname)) {
            dirname0 <- tmp_dir
            tmp_dir <- paste0(tmp_dir, "/",
                              "_hesbcn_", addname)
        }
        if (!silent)
            message(paste("\n Using dir", tmp_dir))
        if (dir.exists(tmp_dir)) {
            stop("dirname ", tmp_dir, "exists")
        }
        dir.create(tmp_dir, recursive = TRUE)
    }

    ## We need to change column names because of a weird behavior of H-ESBCN
    ## See: https://github.com/danro9685/HESBCN/issues/3
    orig_gene_names <- colnames(data)
    colnames(data) <- LETTERS[1:ncol(data)]

    write.csv(data, file = paste0(tmp_dir, "/input.txt"),
              row.names = FALSE, quote = FALSE)

    # Launching
    ## message("Running HESBCN")
    command <- paste0("h-esbcn -d ",
                      tmp_dir, "/input.txt -o ",
                      tmp_dir, "/output.txt -n ",
                      format(round(MCMC_iter), scientific = FALSE),
                      " --reg ", reg)
    if (!is.null(seed))
        command <- paste0(command,
                          " -s ",
                          format(round(seed), scientific = FALSE)
                          )
    
    if (!silent) message("HESBCN command: ", command)
    system(command, ignore.stdout = TRUE)

    # Reading output
    model_info <- import.hesbcn(paste0(tmp_dir, "/output.txt"),
                                genes = orig_gene_names)
 
    indexes <- which(model_info$lambdas_matrix > 0, arr.ind = TRUE)
    lambdas <- model_info$lambdas[indexes[, "col"] - 1]
    from <- rownames(model_info$lambdas_matrix)[indexes[, "row"]]
    to <- colnames(model_info$lambdas_matrix)[indexes[, "col"]]

    ## Check lambdas
    ## Remember their R code divides the lambdas when placing them
    ## in the matrix, but the lambdas, as rates, are the original, undivided ones.
    stopifnot(isTRUE(all.equal(colSums(model_info$lambdas_matrix)[-1],
                               model_info$lambdas,
                               check.attributes = FALSE)))
    
    model_info$edges <- data.frame(From = from,
                                   To = to,
                                   Edge = paste(from, to, sep = " -> "),
                                   Lambdas = lambdas)

    ## Check we are not in the strange case of AND when hanging from Root
    for (ic in seq_along(model_info$parent_set)) {
        if (model_info$parent_set[ic] != "AND") next
        name_child <- names(model_info$parent_set)[ic]
        parents <- model_info$edges[model_info$edges$To == name_child, "From"]
        if ((length(parents) == 1) && (parents == "Root")) {
            warning("The strange case of AND when hanging from root happened. ",
                    "Setting 'AND' to 'Single'")
            model_info$parent_set[ic] <- "Single"
        }
    }
    
    ## Create a "Relation" column, so the $edges component is easy to understand
    ## That is also used when translating to OncoSimulR
    model_info$edges$Relation <- vapply(
        model_info$edges$To,
        function(x) model_info$parent_set[x],
        "some_string"
    )
    model_info["command"] <- command
    return(model_info)
}
