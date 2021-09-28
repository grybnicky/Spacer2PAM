#' Submit spacer FASTA to BLAST and format the results
#'
#' This function allows you to programatically submit your spacer FASTA to BLAST and retrieve the results. The results are then parsed and formatted for use with the joinSpacerDFandAlignmentDF function.
#' @param FASTAfile The file location of the FASTA you want to submit to BLAST.
#' @export
#' @examples
#' FASTA2AlignmentDF()
FASTA2AlignmentDF = function(FASTAfile){
fastaFile = read_file(FASTAfile)

urlEncodedFASTA = URLencode(fastaFile)

spacersRequest = httr::PUT(sprintf("https://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY=%s&DATABASE=nt&PROGRAM=blastn&CMD=put", urlEncodedFASTA))

if(spacersRequest$status_code == 414){
  print("URL is too long. Please resubmit with fewer spacers or submit BLAST request through the web interface manually. Manually BLAST results can be integrated using the function alignmentCSV2DF.")
  stop()
}

requestContent = content(spacersRequest, as="text")

RID = substr(requestContent,(regexpr("QBlastInfoBegin\n    RID = ", requestContent))+26, ((regexpr("QBlastInfoBegin\n    RID = ", requestContent))+36))

repeat{
  print(sprintf("Checking status of BLAST Search RID: %s", RID))
  requestStatus = httr::GET(url = sprintf("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=get&FORMAT_OBJECT=SearchInfo&RID=%s", RID))
  requestStatusContent = content(requestStatus, as="text")
    
  Status = substr(requestStatusContent,(regexpr("QBlastInfoBegin\n\t                Status=",requestStatusContent))+40,(regexpr("QBlastInfoBegin\n\t                Status=",requestStatusContent))+44)
    
  if(Status == "READY")
    {print("BLAST Results are ready") 
    break}
    
  print("Will check status again in one minute")
  Sys.sleep(60)}

resultGET = httr::GET(url = sprintf("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=get&FORMAT_TYPE=Tabular&RID=%s", RID))

resultGETContent = rawToChar(resultGET$content)

alignmentDF = readr::read_tsv(substr(resultGETContent,(regexpr("hits found",resultGETContent)+11), nchar(resultGETContent)-8), col_names = FALSE)

colnames(alignmentDF)<-c("query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")

alignmentDF = alignmentDF[complete.cases(alignmentDF),]

assign("alignmentDF", alignmentDF, envir = .GlobalEnv)
}
