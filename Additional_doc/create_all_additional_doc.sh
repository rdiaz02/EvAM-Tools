#!/bin/bash

## We need to use "--shell-escape" for latexgit.
## (usual caveat: do not use shell-escape with documents you do not trust)

texi2pdf --shell-escape intro_additional_docs.tex

texi2pdf --shell-escape Additional_tech_doc.tex

texi2pdf --shell-escape ../evamtools/inst/miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex

./runKnitr-2.sh  evamtools-examples.Rnw

pdftk intro_additional_docs.pdf Additional_tech_doc.pdf evamtools-examples.pdf ../evamtools.Rcheck/evamtools-manual.pdf Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.pdf cat output Additional_doc_all.pdf

cp Additional_doc_all.pdf ../../evamtools-gh-pages/pdfs/.

cp Additional_tech_doc.pdf ../../evamtools-gh-pages/pdfs/.
cp evamtools-examples.pdf ../../evamtools-gh-pages/pdfs/.
cp ../evamtools.Rcheck/evamtools-manual.pdf ../../evamtools-gh-pages/pdfs/.
cp Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.pdf ../../evamtools-gh-pages/pdfs/.

