#!/bin/bash

texi2pdf Additional_doc.tex

pdftk Additional_doc.pdf ../evamtools.Rcheck/evamtools-manual.pdf cat output Addtional_doc_all.pdf
