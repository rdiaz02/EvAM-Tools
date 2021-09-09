# source("../data/toy_datasets.R")
# source("access_genots_from_oncosimul.R")

for (i in names(all_examples)){
    print(sprintf("Dataset %s", i))
    tmp_dataset <- all_examples[[i]]
    out <- all_methods_2_trans_mat(tmp_dataset)
    saveRDS(out, sprintf("../data/toy_datasets_cpms/%s.rds", i))
}