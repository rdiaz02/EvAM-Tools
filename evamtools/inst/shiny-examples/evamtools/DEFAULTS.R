## Source this file whenever you make changes to it

max_genes <- 10
all_gene_names <- LETTERS[1: max_genes]
template_dag <- matrix(0, ncol = max_genes + 1, nrow = max_genes + 1)
rownames(template_dag) <- colnames(template_dag) <- c("Root", all_gene_names)
template_dag[1, 2] <- 1
template_dag[1, 3] <- 1
template_dag[2, 4] <- 1
template_parent_set <- rep("Single", max_genes)
names(template_parent_set) <- all_gene_names
template_lambdas <- rep(0.5, max_genes)
names(template_lambdas) <- all_gene_names
template_thetas <- matrix(0, ncol = max_genes, nrow = max_genes)
rownames(template_thetas) <- colnames(template_thetas) <- all_gene_names
template_csd_counts <- data.frame(Genotype = character(), Counts = integer())
template_csd_data <- matrix(0, ncol = 4, nrow = 0)

.ev_SHINY_dflt <- list(
  max_genes = 10,
  min_genes = 2,
  cpm_samples = 10000,
  ngenes = 4,
  csd_samples = 1000,
  dag_model = "HESBCN",
  all_cpms = c("OT", "CBN", "OncoBN", "MHN", "MCCBN", "HESBCN"),
  template_data = list(
      csd_counts =  template_csd_counts ## data frame of Genotypes and Counts.
    , data = NULL ## Yeah, what is data? The data matrix of 0/1 with subjects
      ## as rows and genes as columns. 
    , dag = template_dag ## the dag as adjacency matrix.
    , DAG_parent_set = template_parent_set
    , lambdas = template_lambdas ## the lambdas, but I guess also the thetas/Weights
    , thetas = template_thetas ## For MHN.
    , gene_names = LETTERS[1: max_genes]
    , name = "New_CSD"
  )
)

## This expects to be run from the evamtools directory. O.w. adjust path
## Note how the object itself is called (to minimize risk of overwritting)
save(.ev_SHINY_dflt, file = "./data/SHINY_DEFAULTS.RData")
