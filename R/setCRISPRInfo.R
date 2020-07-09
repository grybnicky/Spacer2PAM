#' Set the Name of the CRISPR System
#'
#' This function allows you to set the genus, species, and strain of the organism that encodes the CRISPR system. It also sets a unique identifier for the CRISPR system as there can be multiple per organism. These values will be used by other functions to name files.
#' @param genus A string value that indicates the genus name of the organism encoding the CRISPR array. Make sure to capitalize the first letter of the genus name.
#' @param species A string value that indicates the species name of the organism encoding the CRISPR array. Make sure not to capitalize the first letter of the species name.
#' @param strain A string value that indicates the strain name of the organism encoding the CRISPR array.
#' @param crisprSystemNumber A string value that will act as an identifier to distinguish multiple CRISPR systems in the same organism.
#' @export
#' @examples
#' Given a CRISPR system in Bacillus halodurans C-125 that you designate as "1", you would enter:
#' setCRISPRInfo(genus = "Bacillus", species = "halodurans", strain = "C-125", crisprSystemNumber = "1")
setCRISPRInfo = function(genus,species,strain,crisprSystemNumber){
  assign("genus", genus, envir = .GlobalEnv)
  assign("species", species, envir = .GlobalEnv)
  assign("strain", strain, envir = .GlobalEnv)
  assign("crisprSystemNumber", crisprSystemNumber, envir = .GlobalEnv)
  assign("systemName", sprintf("%s %s %s System %s", genus, species, strain, crisprSystemNumber), envir = .GlobalEnv)
}

setCRISPRInfoInput = function(){
  assign("genus", readline(prompt = "What genus is the organism that encodes this CRISPR system? : "), envir = .GlobalEnv)
  assign("species", readline(prompt = "What species is the organism that encodes this CRISPR system? : "), envir = .GlobalEnv)
  assign("strain", readline(prompt = "What strain is the organism that encodes this CRISPR system? : "), envir = .GlobalEnv)
  assign("crisprSystemNumber", readline(prompt = sprintf("Please provide an ID for this CRISPR system within %s %s : ", genus, species)), envir = .GlobalEnv)
  assign("systemName", sprintf("%s %s %s System %s", genus, species, strain, crisprSystemNumber), envir = .GlobalEnv)
}
