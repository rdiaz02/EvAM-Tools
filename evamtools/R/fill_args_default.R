## Copyright 2022 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.


## list, list -> list

## If any entry en defaults is not in inputl, add to inputl
## If inputl has an name component not in defaults, give warning.

## I use it to add defaults (possibly not explicitly passed) to the
##      list inputl, if the inputl has no value

fill_args_default <- function(inputl, defaults) {

    not_valid <- which(!(names(inputl) %in% names(defaults)))

    if(length(not_valid)) {
        name_input <- deparse(substitute(inputl))
        not_valid_opt <- names(inputl)[not_valid]
        warning("Option(s) ", paste(not_valid_opt, sep = ", ", collapse = ", "),
                " in argument ", name_input, " is(are) invalid.",
                " They will be ignored.")
        inputl <- inputl[-not_valid]
    }
        
    not_passed <- vapply(names(defaults),
                         function(x) !exists(x, where = inputl),
                         TRUE)
    inputl <- c(inputl, defaults[not_passed])
    return(inputl)
}
