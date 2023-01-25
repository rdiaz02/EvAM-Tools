\name{every_which_way_data}
\alias{every_which_way_data}
\docType{data}
\title{
  Cancer data sets from Diaz-Uriarte and Vasallo, 2019; also used in
  Diaz-Colunga and Diaz-Uriarte, 2021.
}

\description{
  Twenty two cancer data sets used in Diaz-Uriarte and Vasallo, 2019, as
  well as Diaz-Colunga and Diaz-Uriarte, 2021. The data cover six
  different cancer types (breast, glioblastoma, lung, ovarian,
  colorectal, and pancreatic cancer), use different types of features
  (nonsynonymous somatic mutations, copy number alterations, or both)
  were analyzed in terms of pathways, functional modules, genes, gene
  events, or mutations (yielding from 3 to 192 different features), and
  have samples sizes from 27 to 594. 

  The original sources are listed below. Most of these data sets have
  been used before in CPM research.

  
  Complete details about sources, processing, and former use in CPM papers are available from S5 Text (https://doi.org/10.1371/journal.pcbi.1007246.s007) of Diaz-Uriarte and Vasallo, 2019. 
    

}
\usage{data("every_which_way_data")}

\format{
  A list of length 22. Each element of the list is a data set, with
  subjects in rows and genes/probes in columns.
}

\references{

  Bamford S, Dawson E, Forbes S, Clements J, Pettett R, Dogan A, et al. The COSMIC (Catalogue of Somatic Mutations in Cancer) Database and Website. Br J Cancer. 2004;91(2):355–358. 

  Cancer Genome Atlas Research Network. Comprehensive Genomic Characterization Defines Human Glioblastoma Genes and Core Pathways. Nature. 2008;455(7216):1061–1068.

  Cancer Genome Atlas Research Network. Comprehensive Genomic Characterization Defines Human Glioblastoma Genes and Core Pathways. Nature. 2008;455(7216):1061–1068.

  Jones S, Zhang X, Parsons DW, Lin JCH, Leary RJ, Angenendt P, et al. Core Signaling Pathways in Human Pancreatic Cancers Revealed by Global Genomic Analyses. Science (New York, NY). 2008;321(5897):1801–6. 

  Parsons DW, Jones S, Zhang X, Lin JCH, Leary RJ, Angenendt P, et al. An Integrated Genomic Analysis of Human Glioblastoma Multiforme. Science. 2008;321(5897):1807–1812. 

  Wood LD, Parsons DW, Jones S, Lin J, Sjoblom T, Leary RJ, et al. The Genomic Landscapes of Human Breast and Colorectal Cancers. Science. 2007;318(5853):1108–1113. 

  Brennan CW, Verhaak RGW, McKenna A, Campos B, Noushmehr H, Salama SR, et al. The Somatic Genomic Landscape of Glioblastoma. Cell. 2013;155(2):462–477. 

  Ding L, Getz G, Wheeler DA, Mardis ER, McLellan MD, Cibulskis K, et al. Somatic Mutations Affect Key Pathways in Lung Adenocarcinoma. Nature. 2008;455(7216):1069–1075. 

  Cancer Genome Atlas Research Network. Integrated Genomic Analyses of Ovarian Carcinoma. Nature. 2011;474(7353):609–615.

  Knutsen T, Gobu V, Knaus R, Padilla-Nash H, Augustud M, Strausberg RL, et al. The Interactive Online SKY/M-FISH & CGH Database and the Entrez Cancer Chromosomes Search Database: Linkage of Chromosomal Aberrations with the Genome Sequence. Genes, Chromosomes and Cancer. 2005;44(1):52–64. 

  Piazza R, Valletta S, Winkelmann N, Redaelli S, Spinelli R, Pirola A, et al. Recurrent SETBP1 Mutations in Atypical Chronic Myeloid Leukemia. Nature Genetics. 2013;45(1):18–24. 

  Cancer Genome Atlas Research Network. Comprehensive Molecular Characterization of Human Colon and Rectal Cancer. Nature. 2012;487(7407):330–337. 

  Diaz-Uriarte, R., & Vasallo, C. (2019). Every which way? On predicting tumor evolution using cancer progression models. PLOS Computational Biology, 15(8), 1007246. http://dx.doi.org/10.1371/journal.pcbi.1007246 
 
 Diaz-Colunga, J., & Diaz-Uriarte, R. (2021). Conditional prediction of
 consecutive tumor evolution using cancer progression models: What
 genotype comes next? PLoS Computational Biology, 17(12):
 e1009055. https://doi.org/10.1371/journal.pcbi.1009055 
 
}


\examples{
data(every_which_way_data)
lapply(every_which_way_data, colnames)
lapply(every_which_way_data, dim)

## Run on a piece of one of the above data sets
Dat1 <- every_which_way_data[[16]][1:40, 2:6]
out <- evam(Dat1,
            methods = c("OT", "OncoBN",
                        "MHN"))
}
\keyword{datasets}
