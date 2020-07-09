#' Submits Alignment Hits for Prophage Prediction
#'
#' This function allows you to submit the accession numbers of organisms identified by BLAST alignment to Phaster for prophage prediction. This process can take a bit of time to run and depends on the number of sequences being submitted.
#' @param joinedData A dataframe containing the joined spacer and alignment information. Defaults to joineddata, the output of the joinSpacerDFandAlignmentDF function.
#' @param startPosition A value to indicate where in the list of accession numbers to start. Defaults to 1.
#' @export
#' @examples
#' submit2Phaster()
submit2Phaster = function(joinedData = joineddata, startPosition = 1){
  totalAccList = unique(joinedData$subject.acc.ver)

  for (i in startPosition:length(totalAccList)){
    if(i %% 8==0){
      print("60 sec pause to prevent Phaster Server Overload")
      print(sprintf("%s / %s Submitted", i, length(totalAccList)))
      Sys.sleep(60)
    }

    if(i %% 32==0){
      print("additional 75 sec pause to prevent Phaster Server Overload")
      Sys.sleep(75)
    }

    if(i %% 64==0){
      print("additional 60 sec pause to prevent Phaster Server Overload")
      Sys.sleep(60)
    }
    #submit accession number to phaster.ca using API
    httr::GET(url = sprintf("http://phaster.ca/phaster_api?acc=%s", totalAccList[i]))
  }
}
