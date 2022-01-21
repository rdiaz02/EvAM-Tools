## Avoid use of ":::". Instead, use approach in
## https://stackoverflow.com/a/46098814

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

## Commented out until we can use OncoBN
## evam_inferTheta <- utils::getFromNamespace("inferTheta", "OncoBN")




