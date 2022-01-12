## Code from https://github.com/BIMIB-DISCo/PMCE/blob/main/Utilities/R/utils.R
## Commit 5578c79 from 2021-09-29

## License: Apache License 2.0, which is compatible with the GPL 3 used by
## the rest of this project.

## Author of code: from commit history, most likely D. Ramazotti (danro9685)

## Authors of project: F. Angaroni, K. Chen, C. Damiani, G. Caravagna,
## A. Graudenzi, D. Ramazotti



"import.hesbcn" <- function( file, genes = NULL ) {

    # read results from file
    results = read.table(file=file,header=FALSE,sep="\n",stringsAsFactors=FALSE,check.names=FALSE)

    # build adjacency matrix
    start = 2
    end = grep("accepted",results[,1])-1
    adjacency_matrix = NULL
    for(entry in results[start:end,1]) {
        adjacency_matrix = rbind(adjacency_matrix,as.numeric(strsplit(entry,split=",")[[1]]))
    }
    
    rownames(adjacency_matrix) = paste0("Gene_",1:nrow(adjacency_matrix))
    colnames(adjacency_matrix) = paste0("Gene_",1:ncol(adjacency_matrix))

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

    # estimated epsilon
    epsilon = as.numeric(gsub("Best epsilon = ","",results[nrow(results),1]))
    
    # return results
    hesbcn = list(adjacency_matrix=adjacency_matrix,lambdas_matrix=lambdas_matrix,parent_set=parent_set,epsilon=epsilon)
    return(hesbcn)

}
