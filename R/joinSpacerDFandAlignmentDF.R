#' Join Spacer and Alignment Dataframes
#'
#' This function allows you to join a dataframe containing the CRISPR array spacer information with a dataframe containing the results of the BLAST alignment. Produces a dataframe with the same number of rows as the alignment dataframe and the spacer information added as additional columns. This function also takes the accession numbers for the alignment and adds a column with the Taxa IDs and organism name. Taxa IDs
#' @param alignmentDF The dataframe containing the CRISPR array spacer information.
#' @param spacerDF The dataframe containing the spacer alignment information.
#' @param accessionDatabaseLocation The location of your Taxonomizer SQL database.
#' @export
#' @examples
#' joinSpacerDFandAlignmentDF(alignmentDF = yourAlignmentDataframe, spacerDF = yourSpacerDataframe, accessionDatabaseLocation = "theLocationOfYourTaxonomizerDatabase")
joinSpacerDFandAlignmentDF = function(alignmentDF, spacerDF, accessionDatabaseLocation){

  #Mutate spacerDF to have query.acc.ver column
  mutspacerDF = spacerDF %>%
    dplyr::mutate(query.acc.ver = sprintf("%sA%sS%s",strain, spacerDF$Array, spacerDF$Spacer))

  #Add Protospacer data to alignment table
  joineddata = alignmentDF %>%
    dplyr::left_join(mutspacerDF, by="query.acc.ver")

  #Add Taxa information and omit NAs
  joineddata = joineddata%>%
    dplyr::mutate("TaxaID"= taxonomizr::accessionToTaxa(as.character(joineddata$subject.acc.ver), accessionDatabaseLocation))%>%
    dplyr::mutate("Species"= as.character(taxonomizr::getTaxonomy(TaxaID, accessionDatabaseLocation)[,7]))%>%
    na.omit()

  assign("joineddata", joineddata, envir = .GlobalEnv)
}
