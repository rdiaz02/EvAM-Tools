#!/bin/bash
## if second argument is purl, it purls
## if "bib" it runs texi2pdf to try to, well, get all refs and bibs done
## if "both" it purls and runs texi2pdf

## Will often fail if run from a non-X ssh

## Using
## https://github.com/EmmanuelCharpentier/patchSynctex


## https://github.com/jan-glx/patchKnitrSynctex

FILE=$1
BASENAME=$(basename $FILE .Rnw)
SECOND=$2

## RSCRIPT="/usr/bin/Rscript"
RSCRIPT="/home/ramon/bin/Rscript --vanilla"


$RSCRIPT -e 'library(knitr); library(patchSynctex); opts_knit$set("concordance" = TRUE); knit2pdf("'$1'"); patchSynctex("'$BASENAME'")'

if [[ $SECOND == "purl" ]]; then
    $RSCRIPT -e 'library(knitr); purl("'$1'")'
fi

if [[ $SECOND == "bib" ]]; then
    texi2pdf $BASENAME.tex
fi

if [[ $SECOND == "both" ]]; then
    $RSCRIPT -e 'library(knitr); purl("'$1'")'
    texi2pdf $BASENAME.tex
fi



rm $BASENAME.Rnw.synctex.gz
rm $BASENAME.Rnw.pdf
ln -s $BASENAME.synctex.gz $BASENAME.Rnw.synctex.gz
ln -s $BASENAME.pdf $BASENAME.Rnw.pdf

## Without the links, you can go from the PDF to
## Emacs, but not the other way around because emacs complaints that it
## cannot find a *.Rnw.pdf. And not, I cannot get past that.
## If you create that, via symlink, then when you
## try to go from Emacs to the PDF, it opens a new viewer. And from that
## one you cannot move back to Emacs.
## I just need to link $.synctex.gz to $.Rnw.synctex.gz
## and .pdf to Rnw.pdf and C-c C-v from emacs. That's it.
## From PDF to emacs by shift-left click (always as "browse tool" with mouse)
## Of course, needs emacs to have the server started

## if run from shell, can jump to error, but goes to .tex, not .Rnw

## It would be great to have this work from Emacs with C-c C-c
## maybe with start-process or similar. I have it now linked to M-n r

## (defun run-knit-synctex-on-buffer () 
## (interactive)
## (shell-command (format "runKnitr-2.sh %s" (shell-quote-argument buffer-file-name)))
## )
## (define-key ess-noweb-minor-mode-map "\M-nr" 'run-knit-synctex-on-buffer)
