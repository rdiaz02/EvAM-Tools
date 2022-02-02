## Copyright 2021, 2022 Pablo Herrera Nieto

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
        "You can download them from https://github.com/danro9685/HESBCN",
        "For proper compilation you should modify", 
        "line 23 in *Makefile* to include *LDLIBS = -lgsl -lm -lgslcblas*.", 
        "You can also change compilation from gcc-5 to something a bit more up to date, like gcc-10.", 
        "Finally, add this folder to your $PATH.",
        "\n******************************************************",
        "\n******************************************************",
        "\n\n"
        )
            )
}
rm(.._EvamTools_test.hesbcn)


#' @title Run HESBCN
#' 
#' @param data Cross secitonal data. Matrix of genes (columns)
#' and individuals (rows)
#' @param n_steps Number of steps to run. Default: 100000
#' @param tmp_folder Folder name where the oput is located. 
#' It will be place under /tmp/HESBCN/tmp_folder
#' @param seed Seed to run the experiment
#' @param clean_dir Whether to delete the folder upon completion
#' 
#' @return A list with the adjacency matrix, the lambdas, the parent set
#' and a data.frame with From-To edges and associated lambdas.
do_HESBCN <- function(data, n_steps=100000, 
    tmp_folder="", seed=NULL, clean_dir=FALSE){
    orig_folder <- getwd()

    # Setting tmp folder
    if(is.null(tmp_folder)) {
        tmp_folder <- tempfile()
        dirname0 <- NULL
        if(!is.null(addname)) {
            dirname0 <- tmp_folder
            tmp_folder <- paste0(tmp_folder, "/",
                              "_cbn_", init.poset, "_",
                              addname)
        }
        if(!silent)
            message(paste("\n Using dir", tmp_folder))
        if(dir.exists(tmp_folder)) {
            stop("dirname ", tmp_folder, "exists")
        }
        dir.create(tmp_folder, recursive = TRUE)
    }
    setwd(tmp_folder)

    orig_gene_names <- colnames(data)
    colnames(data) <- LETTERS[1:ncol(data)]

    write.csv(data, "input.txt", row.names = FALSE, quote = FALSE)
 
    # Launching
    print("Running HESBCN")
    if (is.null(seed)) {
        command <- sprintf("h-esbcn -d input.txt -o output.txt -n %1.f", 
            n_steps)
    } else if (is.numeric(seed) & seed > 0) {
        command <- sprintf("h-esbcn -d input.txt -o output.txt -n %1.f -s %1.f", 
            n_steps, seed)
    }
    system(command, ignore.stdout = TRUE)

    # Reading output
    model_info <- import.hesbcn("output.txt", genes = orig_gene_names)

    indexes_array <- data.frame(which(model_info$lambdas_matrix > 0, arr.ind = TRUE))
    indexes_list <- which(model_info$lambdas_matrix > 0, arr.ind = TRUE)
    lambdas <- model_info$lambdas_matrix[indexes_list]
    from <- rownames(model_info$lambdas_matrix)[indexes_array$row]
    to <- colnames(model_info$lambdas_matrix)[indexes_array$col]
    edges <- paste(from, to, sep = "->")
    adjacency_matrix_2 <- data.frame(From = from, To = to, Edge = edges, Lambdas = lambdas)
    model_info$edges <- adjacency_matrix_2

    # Housekeeping
    setwd(orig_folder)
    if (clean_dir) {
        unlink(tmp_folder, recursive = TRUE)
    }

    return(model_info)
}


