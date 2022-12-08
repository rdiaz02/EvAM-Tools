## Code from file common.r from https://github.com/cbg-ethz/MC-CBN
## Commit a6eceb8 from 2022-02-05

## License: GNU General Public License v2.0
## GNU GPL be combined with software under the
## AGPL 3 as used by the rest of this project.

## Author of code: from comments and commit history (e.g., d2d8dfd, f08add1),
## most likely Hesam Montazeri

## Authors in DESCRIPTION file of R package:
##  Hesam Montazeri, Susana Posada-Cespedes


## RDU: I add the "evam_" to make it clear which one
## we are using. I also use the relations namespace
## explicitly in calls to as.relation, transitive_reduction, and
## relation_incidence.

## Why don't we use the same approach as in code in
## non-exported-from-others.R? Because we do not want to depend
## on having mccbn installed.

evam_random_poset <- function(p, graph_density=0.15, trans_reduced = TRUE){
    evam_random_posets(nr_pos=1, nr_muts=p , ldenses=graph_density, trans_reduced = trans_reduced)[[1]]
}

evam_random_posets <- function(nr_pos, nr_muts , ldenses, trans_reduced = TRUE){

    if(length(ldenses) == 1 ) {
        ldenses = rep(ldenses, nr_pos)
    } else if(length(ldenses) < nr_pos ) {
        stop( paste("Invalid value for ldenses ", ldenses) )
    }
    
    if(length(nr_muts) == 1 ) {
        nr_muts = rep(nr_muts, nr_pos)
    } else if(length(nr_muts) < nr_pos ) {
        stop( paste("Invalid value for nr_muts ", nr_edges) )
    }
    
    nr_edges = ceiling(choose(nr_muts, 2)*ldenses)
    
                                        # Initialisation
    all_posets = rep(0, nr_pos)
    pos = 1
    
    while (pos <= nr_pos){
        
        poset <- matrix(0, nr_muts[pos], nr_muts[pos])
        poset[upper.tri(poset)][sample(choose(nr_muts[pos], 2), nr_edges[pos])] <- 1

                                        # added by Hesam. Add the transitive reduced poset to the list
        if(trans_reduced) {
            poset = evam_trans_reduction(poset)  
        }
        
        all_posets[pos] = list(poset)
        pos = pos + 1  
    }
    all_posets
}

evam_trans_reduction <- function(A) {
    
    old_names = dimnames(A)
    
    colnames(A) =  rownames(A) = paste("n", 1:nrow(A), sep='')
    R <- relations::as.relation(A)
    
    RT = relations::transitive_reduction(R)
    res = relations::relation_incidence(RT)
    res = res[rownames(A),colnames(A)]
    
    dimnames(res) = old_names
    res
}
