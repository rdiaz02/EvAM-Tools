#!/bin/bash

## We need to use "--shell-escape" for latexgit.
## (usual caveat: do not use shell-escape with documents you do not trust)

texi2pdf --shell-escape intro_additional_docs.tex

texi2pdf --shell-escape evamtools_methods_details_faq.tex

texi2pdf --shell-escape ../evamtools/inst/miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex

## ./runKnitr-2.sh  evamtools_examples.Rnw
Rscript --vanilla -e 'library(knitr); knit("evamtools_examples.Rnw")'
texi2pdf --shell-escape evamtools_examples.tex


pdftk intro_additional_docs.pdf evamtools_methods_details_faq.pdf evamtools-examples.pdf ../evamtools.Rcheck/evamtools-manual.pdf Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.pdf cat output Additional_doc_all.pdf

cp Additional_doc_all.pdf ../../evamtools-gh-pages/pdfs/.

cp evamtools_methods_details_faq.pdf ../../evamtools-gh-pages/pdfs/.
cp evamtools_examples.pdf ../../evamtools-gh-pages/pdfs/.
cp ../evamtools.Rcheck/evamtools-manual.pdf ../../evamtools-gh-pages/pdfs/.
cp Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.pdf ../../evamtools-gh-pages/pdfs/.

