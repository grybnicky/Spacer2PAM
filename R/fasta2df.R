#' Convert a FASTA of Spacers to a Dataframe
#'
#' This function allows you to convert a FASTA containing CRISPR array spacers to a dataframe compatible with the rest of the package.
#' @param fastaFile The location of the FASTA file containing the CRISPR array. Spacers will be numbered in the order that they appear in the FASTA file.
#' @param arrayNumbers A vector of the numbers to denote which CRISPR array the spacers came from.
#' @param arrayLengths A vector of the number of spacers in each CRISPR array.
#' @param arrayOrientations A vector indicating the direction of each CRISPR array. Can input either "Forward" or "Reverse".
#' @param arrayRepeats A vector indicating the repeat sequence for each CRISPR array.
#' @param spacerDataFrameName The name of the dataframe output.
#' @export
#' @examples
#' Given a set of 3 CRISPR arrays of lengths 20, 45, and 100 spacers, all in the forward orientation with the sequence repeat "ATCGATCGATCGATCG", you would enter:
#' df2fasta(fastaFile = "yourFileLocation", arrayNumbers = c(1,2,3), arrayLengths = c(20,45,100), arrayOrientations = c("Forward", "Forward", "Forward"), arrayRepeats = c("ATCGATCGATCGATCG","ATCGATCGATCGATCG","ATCGATCGATCGATCG"), spacerDataFrameName = "yourDataframeName")
fasta2df = function(fastaFile,arrayNumbers,arrayLengths,arrayOrientations,arrayRepeats, spacerDataFrameName){
  readFASTA = seqinr::read.fasta(file=fastaFile, as.string = TRUE)
  realArrayNumbers = c()
  realArraySpacer = c()
  realArrayOrientations = c()
  realArrayRepeats = c()
  for (i in 1:length(arrayNumbers)){
    realArrayNumbers = append(realArrayNumbers, rep(arrayNumbers[i], times = arrayLengths[i]), after = length(realArrayNumbers))
    realArrayOrientations = append(realArrayOrientations, rep(arrayOrientations[i], times = arrayLengths[i]), after = length(realArrayOrientations))
    realArrayRepeats = append(realArrayRepeats, rep(arrayRepeats[i], times = arrayLengths[i]), after = length(realArrayRepeats))
    realArraySpacer = append(realArraySpacer, 1:arrayLengths[i], after = length(realArraySpacer))
  }

  Strain = c()
  Spacers = c()
  Array.Orientation = c()
  Repeat = c()
  Array = c()
  Spacer = c()
  for (i in 1:length(readFASTA)){
    Strain[i] = strain
    Spacers[i] = unlist(readFASTA[[i]])
    Array.Orientation[i] = realArrayOrientations[i]
    Repeat[i] =realArrayRepeats[i]
    Array[i] = realArrayNumbers[i]
    Spacer[i] = realArraySpacer[i]
  }

  assign(spacerDataFrameName, data.frame(Strain, Spacers, Array.Orientation, Repeat, Array, Spacer, stringsAsFactors = FALSE), envir = .GlobalEnv)
}

