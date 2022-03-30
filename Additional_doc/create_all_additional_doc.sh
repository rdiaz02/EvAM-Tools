#!/bin/bash

texi2pdf intro_additional_docs.tex
texi2pdf Additional_doc.tex
texi2pdf ../evamtools/inst/miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex
pdftk intro_additional_docs.pdf Additional_doc.pdf ../evamtools.Rcheck/evamtools-manual.pdf Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.pdf cat output Additional_doc_all.pdf

cp Additional_doc_all.pdf ../../evamtools-gh-pages/pdfs/.
