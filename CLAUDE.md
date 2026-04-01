# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Scientific context for the project

This repository contains an R package (evamtools directory) and shiny R application. The main code directory is the evamtools directory. Other files, including under Additional_doc (LaTeX files provided) give additional information and tutorials; evamtools/inst/miscell contains additional documentation, examples, and tests.


## R conventions

- Stick to standard R as much as possible and minimize use of the tidyverse. Specifically:

- Do not use dplyr. If data.frame is slow and we need an alternative, use data.table.

- Never use pipes (from tidyverse or the new ones in R).

- Do not use the stringr package. Use stringi.

- ggplot2 is OK, but respecting the above rules.

- Avoid (unless essential) any non-standard evaluation. We want clean, standard, non-surprises R code that will run 5 years from now.

- Try not to go beyond 80-columns (with common sense).

- Do not use Roxygen2 or similar for documentation (you'll see remnants are commented out). I write Rd files directly and use comments in the code.

- Never put 2 or more statements per line separated by a ";".

- Use spaces around "+", "-", "=", "*", etc.

- Please, DO NOT ask for permission to execute R code that tests the files you create. You have my permission.
