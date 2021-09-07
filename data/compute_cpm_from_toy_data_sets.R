pwd2 <- getwd()
source("toy_datasets.R")
setwd("../examples/")
source("access_genots_from_oncosimul.R")
setwd("../code_from_what_genotype_next/")
source("code-all-methods-minimal.R")
setwd(pwd2)
rm(pwd2)

for (i in names(all_examples)){
    print(sprintf("Dataset %s", i))
    tmp_dataset <- all_examples[[i]]
    out <- all_methods_2_trans_mat(tmp_dataset)
    saveRDS(out, sprintf("./toy_datasets_cpms/%s.rds", i))
}