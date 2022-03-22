## Avoid use of ":::". Instead, use approach in
## https://stackoverflow.com/a/46098814
## And note that this is all code that I (RDU) control as is from OncoSimulR
## except for GA_Likelihood

evam_GA_Likelihood <- utils::getFromNamespace("GA_Likelihood",
                                              "OncoBN")

evam_allGenotypes_to_matrix <-
    utils::getFromNamespace("allGenotypes_to_matrix",
                            "OncoSimulR")

evam_wrap_accessibleGenotypes <-
    utils::getFromNamespace("wrap_accessibleGenotypes",
                            "OncoSimulR")

evam_checkProperOTAdjMat <- utils::getFromNamespace("checkProperOTAdjMat",
                                                    "OncoSimulR")

evam_checkProperMinimalAdjMat <-
    utils::getFromNamespace("checkProperMinimalAdjMat",
                            "OncoSimulR")


evam_OTtoPoset <- utils::getFromNamespace("OTtoPoset",
                                          "OncoSimulR")

evam_posetToGraph <- utils::getFromNamespace("posetToGraph",
                                             "OncoSimulR")

evam_shannonI <- utils::getFromNamespace("shannonI",
                                         "OncoSimulR")

