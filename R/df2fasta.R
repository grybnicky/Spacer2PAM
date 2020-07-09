#' Convert a Spacer Dataframe to a FASTA file
#'
#' This function allows you to convert a dataframe containing the spacer data to a FASTA file with annotation nomenclature compatible with the rest of this package. The FASTA is written to the current file location with file name "_spacers.fasta" where the blank is filled with the systemName as defined by the setCRISPRInfo function.
#' @param spacerDataFrame The appropriately formatted data frame containing the CRISPR array spacers and associated data.
#' @export
#' @examples
#' df2fasta()
df2fasta = function(spacerDataFrame){
  name = c()
  for (i in 1:nrow(spacerDataFrame)){
    name[i] = sprintf("%sA%sS%s",strain,spacerDataFrame$Array[i],spacerDataFrame$Spacer[i])
  }

  #Write FASTA file
  seqinr::write.fasta(as.list(spacerDataFrame$Spacers),name, file.out=sprintf("%s spacers.fasta", systemName), open="w", as.string= FALSE)

}
