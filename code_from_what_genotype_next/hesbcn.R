

do_hesbcn <- function(data, n_steps=100000, tmp_folder="", seed=NULL, clean_dir=FALSE){
    # Setting tmp folder
    dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    random_letters <- paste(c("_", tmp_folder, "_", LETTERS[floor(runif(4, min=1, max=26))]), collapse="")
    tmp_folder <- paste(c(dateTime, random_letters), collapse="")

    tmp_folder <- file.path("/", "tmp", "HESBCN", tmp_folder)
    dir.create(tmp_folder, recursive=TRUE)
    orig_folder <- getwd()

    setwd(tmp_folder)

    write.csv(data, "input.txt", row.names=FALSE, quote=FALSE)
 
    # Launching
    print("Running HESBCN")
    if (is.null(seed)){
        command = sprintf("h-esbcn -d input.txt -o output.txt -n %1.f", n_steps)
    } else if (is.numeric(seed) & seed > 0){
        command = sprintf("h-esbcn -d input.txt -o output.txt -n %1.f -s %1.f", n_steps, seed)
    }
    system(command, ignore.stdout = TRUE)

    # Reading output
    # model_info <- cleaning_output("output.txt")
    model_info <- import.hesbcn("output.txt")
    # browser()

    # Transforming data from model

    # Housekeeping
    setwd(orig_folder)
    if(clean_dir){
        unlink(tmp_folder, recursive = TRUE)
    }

    # browser()

    return(list(
            # edges = dbn_out,
            # weighted_fgraph = weighted_fgraph,
            # trans_mat_genots = trans_mat_genots,
            # likelihood = out$score
            ))
}

import.hesbcn <- function( file, genes = NULL ) {
    ## Copied from https://github.com/BIMIB-DISCo/PMCE/blob/main/R/utils.R

    # read results from file
    results = read.table(file=file,header=FALSE,sep="\n",stringsAsFactors=FALSE,check.names=FALSE)

    # build adjacency matrix
    start = 2
    end = grep("accepted",results[,1])-1
    adjacency_matrix = NULL
    for(entry in results[start:end,1]) {
        adjacency_matrix = rbind(adjacency_matrix,as.numeric(strsplit(entry,split=",")[[1]]))
    }
    tmp_genes = LETTERS[1:nrow(adjacency_matrix)]
    # rownames(adjacency_matrix) = paste0("Gene_",1:nrow(adjacency_matrix))
    # colnames(adjacency_matrix) = paste0("Gene_",1:ncol(adjacency_matrix))

    rownames(adjacency_matrix) = paste0(tmp_genes)
    colnames(adjacency_matrix) = paste0(tmp_genes)

    # get lambda values and parent set types
    lambdas = results[(grep("Best Lambdas",results[,1])+1),1]
    lambdas = strsplit(lambdas,split=" ")[[1]][1:(length(strsplit(lambdas,split=" ")[[1]])-1)]
    lambdas = as.numeric(lambdas)
    parent_set = strsplit(results[(grep("best theta types",results[,1])+1),1],split=",")[[1]]
    for(i in 1:length(parent_set)) {
        if(parent_set[i]=="-1") {
            parent_set[i] = "Single"
        }
        else if(parent_set[i]=="0") {
            parent_set[i] = "AND"
        }
        else if(parent_set[i]=="1") {
            parent_set[i] = "OR"
        }
        else if(parent_set[i]=="2") {
            parent_set[i] = "XOR"
        }
    }
    names(parent_set) = paste0("Gene_",1:nrow(adjacency_matrix))

    # make final results and normalize lambda values
    adjacency_matrix = cbind(rep(0,nrow(adjacency_matrix)),adjacency_matrix)
    adjacency_matrix = rbind(rep(0,ncol(adjacency_matrix)),adjacency_matrix)
    rownames(adjacency_matrix)[1] = "Root"
    colnames(adjacency_matrix)[1] = "Root"
    adjacency_matrix["Root",as.numeric(which(colSums(adjacency_matrix)==0))[-1]] = 1
    lambdas_matrix = adjacency_matrix
    for(i in 2:nrow(lambdas_matrix)) {
        curr_in_lambda = lambdas[(i-1)]
        # we assume equally distributed lambdas among components of a single formula
        lambdas_matrix[as.numeric(which(adjacency_matrix[,i]==1)),i] = curr_in_lambda / length(which(adjacency_matrix[,i]==1))
    }

    # set genes names
    if(!is.null(genes)) {
        rownames(adjacency_matrix) = c("Root",genes)
        colnames(adjacency_matrix) = c("Root",genes)
        rownames(lambdas_matrix) = c("Root",genes)
        colnames(lambdas_matrix) = c("Root",genes)
        names(parent_set) = genes
    }

    # adjacency_matrix_2
    indexes_array = data.frame(which(lambdas_matrix>0, arr.ind=TRUE))
    indexes_list = which(lambdas_matrix>0, arr.ind=TRUE)
    lambdas = lambdas_matrix[indexes_list]
    from = rownames(lambdas_matrix)[indexes_array$row]
    to = colnames(lambdas_matrix)[indexes_array$col]
    edges = paste(from, to, sep="->")
    adjacency_matrix_2 <- data.frame(cbind(from, to, edges, lambdas))
    colnames(adjacency_matrix_2) <- c("From", "To", "Edge", "Lambdas")

    # estimated epsilon
    epsilon = as.numeric(gsub("Best epsilon = ","",results[nrow(results),1]))
    
    # return results
    hesbcn = list(adjacency_matrix=adjacency_matrix,lambdas_matrix=lambdas_matrix,parent_set=parent_set,epsilon=epsilon,hesbcn_out=adjacency_matrix_2)
    return(hesbcn)

}
