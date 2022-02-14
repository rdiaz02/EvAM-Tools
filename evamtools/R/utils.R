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


## TODO add support for custom genes names in the str conversions

#' @title Integer to binary
#' 
#' @description  Transform a genotype from integer nomenclature to binary coding
#' 
#' @param int_state Ingeter >=0
#' @param n Integer. Number of digits to return
#' 
#' @return vector with 0 and 1 with the binary coding
int2binary <- function(int_state, n = NULL){
    int_state <- as.integer(int_state)

    if (is.na(int_state)){
        stop("Invalid input")
    }

    if(int_state < 0){
        stop("Only integers >= 0 are valid")
    }

    num_digits <- NULL
    tmp_num_digits <- 0
    while(is.null(num_digits)){
        tmp_num_digits <- tmp_num_digits + 1
        if(2**tmp_num_digits > int_state){
            num_digits <- tmp_num_digits 
        }
    }
    if(is.null(n) || n < 1 ){
        n <- num_digits
    }

    base_binary <- as.integer(base::intToBits(int_state))[1 : num_digits]
    remaining_digits <- n - num_digits
    if(remaining_digits < 0) remaining_digits <- 0

    return(c(base_binary[1:num_digits], rep(0, remaining_digits)))
}

#' @title Integer to string
#' 
#' @description Transform a genotype from integer nomenclature to string binary coding
#' 
#' @param int_state Ingeter >=0
#' @param sep String. Separator between genes letters
#' @param wt String. How to define the wild type
#' 
#' @return string with mutated genes
int2str <- function(int_state, sep = ", ", wt = "WT"){
    if(int_state < 0){
        stop("Only integers >= 0 are valid")
    }

    if(int_state == 0) return(wt)

    binary_state <- int2binary(int_state)
    return(binary2str(binary_state, sep = sep))
}

#' @title Binary to Integer
#' 
#' @description  Transform a genotype from binary nomenclature to integer value
#' 
#' @param binary_state vector with 0 and 1
#' 
#' @return Integer values of the genotype
binary2int <- function(binary_state){
    if(length(binary_state) == 0) stop("Empty state")
    
    if(sum(!(binary_state %in% c(0, 1)))){
        stop("Binary state should only contain 0 and 1")
    }
    powers.of.two <- 2^(0:(length(binary_state) - 1))
    return(as.integer(binary_state %*% powers.of.two))
}

#' @title Binary to string
#' 
#' @description  Transform a genotype from binary nomenclature to integer string
#' 
#' @param binary_state vector with 0 and 1
#' @param sep String. Separator between genes letters
#' @param wt String. How to define the wild type
#' 
#' @return string with mutated genes
binary2str <- function(binary_state, sep = ", ", wt = "WT"){
    if(length(binary_state) == 0) stop("Empty state")

    if(sum(!(binary_state %in% c(0, 1)))){
        stop("Binary state should only contain 0 and 1")
    }

    if(sum(binary_state) == 0) return(wt)
    return(paste(LETTERS[
        which(binary_state == 1)
        ], collapse = sep))
}

#' @title String to binary
#' 
#' @description Transform a genotype from string nomenclature to binary coding
#' 
#' @param str_state String
#' @param sep String. Separator between genes letters
#' @param wt String. How to define the wild type
#' @param n Integer. Number of digits to return
#' 
#' @return vector with 0 and 1 with the binary coding
str2binary <- function(str_state, sep =", ", wt = "WT", n = NULL){
    if(str_state == wt && !(is.null(n)) && n > 0) return(rep(0, n))
    else if (str_state == wt) return(c(0))
    
    str_state <- sort(str_split(str_state, sep)[[1]])
    num_digits <- which(LETTERS == str_state[length(str_state)])

    if(is.null(n) || n < 0){
        n <- num_digits
    }

    remaining_digits <- n - num_digits
    if(remaining_digits < 0) remaining_digits <- 0
    return(c(as.integer(LETTERS %in% str_state)[1:num_digits]
        , rep(0, remaining_digits)))
} 

#' @title String to Integer
#' 
#' @description Transform a genotype from string nomenclature to binary coding
#' 
#' @param str_state String
#' @param sep String. Separator between genes letters
#' @param wt String. How to define the wild type
#' @param n Integer. Number of digits to return
#' 
#' @return vector with 0 and 1 with the binary coding
str2int <- function(str_state, sep =", ", wt = "WT", n = NULL){
    binary_state <- str2binary(str_state, sep = sep, wt = wt, n = n)
    return(binary2int(binary_state))
}



## #' @title List of sorted genotypes
## #' 
## #' @description Returns genotypes sorted for a given number of genes
## #'    with sorting the same as given by sample_to_pD_order
## #'    gene_names is always sorted inside the function to
## #'    ensure results do not differ by gene_names order
## #' 
## #' @param n_genes Integer >0
## #' @param gene_names Vector of strings
## #'      If NULL, defaults gene_names will the firt n_genes letters of the alphabet
## #' 
## #' @return Vector with the sorted genotypes
## #' 
generate_sorted_genotypes <- function(n_genes, gene_names = NULL,
                                      sort_gene_names = TRUE){
    if(n_genes <= 0) stop("Number of genes should be > 0")
    n_states <- 2**n_genes
    if(is.null(gene_names)) gene_names <- LETTERS[seq_len(n_genes)]
    stopifnot(n_genes == length(gene_names))
    if(sort_gene_names)
        gene_names <- sort(gene_names)
    
    sorted_genotypes <- vapply(0:(n_states - 1), function(x) {
        tmp_genotype <- paste(gene_names[int2binary(x, n_genes) == 1]
            , collapse = ", ")
        tmp_genotype <- ifelse(tmp_genotype == "", "WT", tmp_genotype)
        return(tmp_genotype)
    }, character(1))

    return(sorted_genotypes)
}





## RDU: FIXME: Esto está roto!!!
##  Por ejemplo, con n_genes = 0 no funciona
##  Y qué sentido tiene dar algo con 0 genes?
##  Y LETTERS[1:0] devuelve A
##  No se debe añadir funcionalidad que no se usa a menos que se teste!!!!
## #' @title List of sorted genotypes
## #' 
## #' @description Returns all sorted genotypes for a given number of genes
## #' 
## #' @param n_genes Ingeter >=0
## #' @param gene_names Vector of strings
## #'      If NULL, defaults gene_names will the firt n_genes letters of the alphabet
## #' 
## #' @return Vector with the sorted genotypes
## generate_sorted_genotypes <- function(n_genes, gene_names = NULL){
##     if (is.null(gene_names)) gene_names <- LETTERS[1:n_genes]
##     # if(n_genes < 0) stop("Number of genes should be >= 0")
##     if(n_genes == 0) {states <- c()}
##     else {n_states <- 2**n_genes}

##     if(is.null(gene_names)) gene_names <- LETTERS[1:n_genes]

##     sorted_genotypes <- vapply(0:(n_states - 1), function(x){
##         tmp_genotype <- paste(gene_names[int2binary(x, n_genes) == 1]
##             , collapse = ", ")
##         tmp_genotype <- ifelse(tmp_genotype == "", "WT", tmp_genotype)
##         return(tmp_genotype)
##     }, character(1))

##     return(sorted_genotypes)
## }

