library(stringr)
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
    # print(n)
    # print(num_digits)

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
#' @param int_state Ingeter >=0
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
#' @param int_state Ingeter >=0
#' @param sep String. Separator between genes letters
#' @param wt String. How to define the wild type
#' @param n Integer. Number of digits to return
#' 
#' @return vector with 0 and 1 with the binary coding
str2int <- function(str_state, sep =", ", wt = "WT", n = NULL){
    binary_state <- str2binary(str_state, sep = sep, wt = wt, n = n)
    return(binary2int(binary_state))
}