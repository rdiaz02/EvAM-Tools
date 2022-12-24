do_HyperHMM <- function(xoriginal,
                        precursors = NA, # precursors (e.g. ancestors) -- blank for cross-sectional
                        nboot = 1, # number of boostrap resamples
                        random.walkers = 0, # run random walkers for each resample? 0 no, 1 yes
                        label = "label", # label for file I/O
                        simulate = TRUE, # actually run HyperHMM? If not, try and pull precomputed output using "label"
                        fork = FALSE # if we are running it, fork to a new process, or keep linear?
                        ){
  # get number and names of features
  L <- ncol(xoriginal)
  features <- colnames(xoriginal)
  
  # deal with non-binary values
  if(any(xoriginal != 0 & xoriginal != 1)) {
    message("Non-binary observations found. Interpreting all nonzero entries as presence")
    xoriginal[xoriginal != 0 & xoriginal != 1] <- 1
  }
  
  #compile into strings
  obs.sums <- rowSums(xoriginal)
  obs.rows <- apply(xoriginal, 1, paste, collapse="")
  final.rows <- obs.rows
  
  # check to see if precursor states have been provided; set cross-sectional flag accordingly
  if(!is.matrix(precursors)) { 
    message("Couldn't interpret precursors observations as matrix!")
    message("No precursor observations; assuming cross-sectional setup.")
    cross.sectional = 1
  } else {
    cross.sectional = 0
    # deal with issues in precursor set
    if(nrow(precursors) != nrow(xoriginal)) { message("Number of observations and precursors didn't match!"); return()}
    if(ncol(precursors) != ncol(xoriginal)) { message("Number of features in observations and precursors didn't match!"); return }
    if(any(precursors != 0 & precursors != 1)) {
      message("No-binary precursors found. Interpreting all nonzero entries as presence")
      precursors[precursors != 0 & precursors != 1] <- 1
    }
    # check to make sure all precursors states come before their descendants
    precursor.sums = rowSums(precursors)
    if(any(precursor.sums > obs.sums)) {
      message("Precursors found more advanced than observations!")
      return()
    }
    # compile precursor-observations pairs into strings
    precursor.rows = apply(precursors, 1, paste, collaspse="")
    for(i in 1:length(obs.rows)) {
      final.rows[i] <- paste(c(precursor.rows[i], obs.rows[i]), collapse = "")
    }
  }
  
  if(any(obs.sums == 0)) {
    message("Dropping O^L observations")
    final.rows = final.rows[obs.sums]
  }
  
  # create input file for HyperHMM
  filename <- paste(c("HyperHMM-in-", label, ".txt"), collapse="")
  write(final.rows, filename)
  
  # create commandname for HyperHMM
  # possible names of executable
  cmds <- c("hyperhmm.ce", "hyperhmm.exe", "hyperhmm")
  #see if any are here
  for(cmd in cmds) {
    if(cmd %in% list.files()) {
      commandname <- cmd
    }
  }
  # create call to HyperHMM
  if(simulate == TRUE) {
    if(fork == TRUE) {
      syscommand <- paste(c("./", commandname, " ", filename, " ", L, " ", nboot, " ", label, " ", cross.sectional, " ", random.walkers, " &"), collapse = "")
      message(paste(c("Attempting to externally execute ", syscommand), collapse = ""))
      system(syscommand)
      message("Forked: not retrieving output")
      return(NULL)
    } else {
      syscommand <- paste(c("./", commandname, " ", filename, " ", L, " ", nboot, " ", label, " ", cross.sectional, " ", random.walkers), collapse = "")
      message(paste(c("Attempting to externally execute ", syscommand), collapse = ""))
      system(syscommand)
    }
  }
  
  # attempt to read the output of the run
  mean.filename = paste(c("mean_", label, ".txt"), collapse="")
  sd.filename = paste(c("sd_", label, ".txt"), collapse="")
  transitions.filename = paste(c("transitions_", label, ".txt"), collapse="")
  viz.filename = paste(c("graph_viz_", label, ".txt"), collapse="")
  
  message(paste(c("Attempting to read output... e.g. ", mean.filename), collapse=""))
  if(!(mean.filename %in% list.files()) | !(sd.filename %in% list.files()) | !(transitions.filename %in% list.files()) | !(viz.filename %in% list.files())) {
    message(paste(c("Couldn't find file output of externally executed HyperHMM!"), collapse=""))
    return()
  } else {
    mean.table = read.table(mean.filename)
    sd.table = read.table(sd.filename)
    transitions = read.csv(transitions.filename, sep=" ")
    viz.tl = readLines(viz.filename)
    
    message("All expected files found.")
    L <- nrow(mean.table)
    # pull the wide output data into long format
    stats.df <- data.frame(feature=rep(0, L*L), order = rep(0, L*L),
                           mean = rep(0, L*L), sd = rep(0, L*L))
    index <- 1
    for(i in 1:L) {
      for(j in 1:L) {
        stats.df$feature[index] <- L-i+1
        stats.df$order[index] <- j
        stats.df$mean[index] <- mean.table[i,j]
        stats.df$sd[index] <- sd.table[i,j]
        index <- index+1
      }
    }
    fitted <- list(stats.df, transitions, features, viz.tl)
    return(fitted)
    #create a list of useful objects and return it
  }
}
