

do_HESBCN <- function(data, n_steps=100000, tmp_folder=""){
    dateTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    random_letters <- paste(c("_", tmp_folder, "_", LETTERS[floor(runif(4, min=1, max=26))]), collapse="")
    tmp_folder <- paste(c(dateTime, random_letters), collapse="")

    tmp_folder <- file.path("/", "tmp", "HESBCN", tmp_folder)
    dir.create(tmp_folder, recursive=TRUE)
    orig_folder <- getwd()

    setwd(tmp_folder)

    write.csv(data, "input.txt", row.names=FALSE, quote=FALSE)

    #Launching
    print("Running HESBCN")
    time_hesbcn <- system.time(
        system(sprintf("h-esbcn -d input.txt -o output.txt -n %1.f -s 42", n_steps), intern=FALSE))["elapsed"]

    # Reading output
    output <- cleaning_output("output.txt")
    
    setwd(orig_folder)

    return(list(
            # edges = dbn_out,
            # weighted_fgraph = weighted_fgraph,
            # trans_mat_genots = trans_mat_genots,
            # likelihood = out$score
            time = time_hesbcn
            ))
}

cleaning_output <- function(file_name){
    connection <- file(file_name)
    open(connection)
    is_graph_completed <- FALSE
    line <- 1
    while(length(line) > 0){
        line <- readLines(connection, n=1)
        print(line)
        if (!is.null(line) & line == "best poset: "){
            print("create Graph")
            vector <- readLines(connection, n=1)
            first_vector <- as.numeric(unlist(strsplit(vector, ",")[1]))
            remaining <- readLines(connection, n=(length(first_vector) - 1))

            poset <- c(vector, remaining)
            poset <- sapply(poset, function(x){
                as.numeric(unlist(strsplit(x, ",")[1]))})
            poset <- t(poset)

            colnames(poset) <- LETTERS[1: length(first_vector)]
            rownames(poset) <- LETTERS[1: length(first_vector)]

        }
    
        if(line == FALSE & is_graph_completed == FALSE){
            ## DO GRAPH
        }

        if(startsWith(line, "best theta types")){
            theta_types <- readLines(connection, n=1)
            theta_types <- as.numeric(unlist(strsplit(theta_types, ",")[1]))
        }

        if(line == "Best Lambdas:"){
            lambdas <- readLines(connection, n=1)
            lambdas <- as.numeric(unlist(strsplit(lambdas, " ")[1]))
        }

        if(startsWith(line, "Best epsilon")){
            epsilon <- as.double(unlist(strsplit(line, "=")[2]))
            ## This is the last thing we want to get
            break
        }
    }

    close(connection)

    return(list(poset = poset,
                theta_types = theta_types,
                lambdas = lambdas,
                epsion = epsilon))
}
N <- 100
na <- N
nc <- N + round( 10 * runif(1))
nab <- 1.6 * N + round( 10 * runif(1))
ncd <- 1.5 * N + round( 10 * runif(1))
n00 <- round( 10 * runif(1))
dB <- matrix(
  c(
    rep(c(1, 0, 0, 0), na) 
    , rep(c(0, 0, 1, 0), nc)
    , rep(c(1, 1, 0, 0), nab)
    , rep(c(0, 0, 1, 1), ncd)        
    , rep(c(0, 0, 0, 0), n00)
  ), ncol = 4, byrow = TRUE
)

colnames(dB) <- LETTERS[1:4]

# do_HESBCN(dB)
cleaning_output("/home/pablo/HESBCN/output3.txt")