#' Predict Protospacer Adjacent Motifs from Joined Spacer and Alignment Data
#'
#' This function allows you to filter alignment data and then predict protospacer adjacent motifs. Predicted motifs are scored based on number of alignments and the amount of information the predicted PAM encodes relative to the maximum possible information encoded.
#' @param joinedData The dataframe containing both CRISPR array spacer data and alignment data. Oftentimes will be the output of the joinSpacerDFandAlignmentDF function.
#' @param uniqueAlignsRange A vector setting the range of values determining unique alignment filter criteria. The abundance of certain sequences can be biased by representation in the database. Given a TRUE input, this filter will only keep one instance of alignments that are from the same organism and align to the same spacer. Values can be TRUE or FALSE. Defaults to TRUE.
#' @param excludeSelfRange A vector setting the range of values determining whether the organism encoding the CRISPR array is excluded. The most abundant alignment will likely be to the organism encoding the CRISPR array analyzed. Given a TRUE input, this filter will remove all alignmens to the organism that encodes the CRISPR array being analyzed. Values can be TRUE of FALSE. Defaults to TRUE.
#' @param numGapsRange A vector setting the range of values determining the maximum number of alignment gaps permissable in an alignment. Alignments with equal or fewer number of alignment gaps as the input will pass through the filter. Values can range from 0 to the length of the CRISPR array. Defaults to 0.
#' @param e.valueRange A vector setting the range of values determining the maximum e-value permissable for an alignment. Alingments with an equal or lower e-value will pass through the filter. Values must be positive. Defaults to 0.05.
#' @param nucleotidesShorterThanProtospacerRange A vector setting the range of values determining the maximum length that an alignment can be shorter than a spacer. Alignments with a length less than the spacer length minus the value will be filtered out. Values can range from 0 to the length of the spacer. Defaults to 0.
#' @param queryStartRange A vector setting the range of values determining the maximum nucleotide position in ths spacer that the alignment can start. Alignments that start 3' of the indicated position will be filtered out. Values can range from 1 to the length of the spacer. Defaults to 1.
#' @param prophageOnlyRange A vector setting the range of values determining whether prophage content is used as a filter criteria. Given a TRUE input, only alignments that are located within predicted prophage regions will pass through the filter. Values can be TRUE or FALSE. Defaluts to FALSE.
#' @param flankLength A number indicating the number of length of DNA sequence to be searched for a PAM. This value influences the PAM score calculation. Defaults to 10.
#' @param RangeStart A number indicating what filter combination in the filter criteria space the function should start at. Defaults to 1.
#' @param saveLogo A value that determines whether the sequence logos generated are saved. Given a TRUE input, a PDF of the sequence logo is saved to the working directory. Defaults to TRUE.
#' @param savePAMSeqs A value that determines whether the vector of sequences used to build the sequence logo are saved. Given a TRUE input, the vectors of upstream and downstream sequences are assigned to the global environment as upstreamPAMSeqs and downstreamPAMSeqs, respectively. Defaults to FALSE.
#' @param removeFASTA A value that determines whether the temporary FASTA file from eFetch is deleted from the working directory.  Given a TRUE input, the FASTA file is deleted. Defaults to TRUE.
#' @export
#' @examples
#' If you wanted to search the prediction space over the e-value cut-offs 0.01, 0.05, and 0.1 keeping the rest of the values at default and resulting in 3 total predictions, you would enter:
#' join2PAM(joinedData = nameOfJoinedDataframe, e.valueRange = c(0.01, 0.05, 0.1))
#'
#' If you wanted to search the prediction space with and without prophage prediction keeping the rest of the values at default and resulting in 2 total predictions, you would enter:
#' join2PAM(joinedData = nameOfJoinedDataframe, prophageOnlyRange = c(TRUE, FALSE))
#'
#' If you wanted to search the prediction space combine the searhes above, resulting in 6 total predictions, you would enter:
#' join2PAM(joinedData = nameOfJoinedDataframe, e.valueRange = c(0.01, 0.05, 0.1), prophageOnlyRange = c(TRUE, FALSE))
join2PAM = function(joinedData,
                    uniqueAlignsRange = T,
                    excludeSelfRange = T,
                    numGapsRange = 0,
                    e.valueRange = 0.05,
                    nucleotidesShorterThanProtospacerRange = 0,
                    queryStartRange = 1,
                    prophageOnlyRange = F,
                    flankLength = 10,
                    RangeStart = 1,
                    saveLogo = T,
                    savePAMSeqs = F,
                    removeFASTA = T
                    ){
  nConditions = length(uniqueAlignsRange)*length(excludeSelfRange)*length(numGapsRange)*length(e.valueRange)*length(nucleotidesShorterThanProtospacerRange)*length(queryStartRange)*length(prophageOnlyRange)

  collectionFrame = stats::setNames(data.frame(matrix(nrow = nConditions, ncol = 18)), c("uniqueAligns","excludeSelf", "numGaps", "e.value", "nucleotidesShorterThanProtospacer", "queryStart", "prophageOnly", "Filter0","Filter1", "Filter2", "Filter3", "Filter4", "Filter5", "Filter6","upPAM","upScore","downPAM","downScore"))

  counter = RangeStart

  uniqueAlignsVec = c()
  excludeSelfVec = c()
  numGapsVec = c()
  e.valueVec = c()
  nucleotidesShorterThanProtospacerVec = c()
  queryStartVec = c()
  prophageOnlyVec = c()

  for(j in 1:length(uniqueAlignsRange)){

    for(o in 1:length(excludeSelfRange)){

      for(p in 1:length(numGapsRange)){

        for(k in 1:length(e.valueRange)){

          for(l in 1:length(nucleotidesShorterThanProtospacerRange)){

            for(m in 1:length(queryStartRange)){

             for(n in 1:length(prophageOnlyRange)){

                uniqueAlignsVec = append(uniqueAlignsVec, uniqueAlignsRange[j], after = length(uniqueAlignsVec))
                excludeSelfVec = append(excludeSelfVec, excludeSelfRange[o], after = length(excludeSelfVec))
                numGapsVec = append(numGapsVec, numGapsRange[p], after = length(numGapsVec))
                e.valueVec = append(e.valueVec, e.valueRange[k], after = length(e.valueVec))
                nucleotidesShorterThanProtospacerVec = append(nucleotidesShorterThanProtospacerVec, nucleotidesShorterThanProtospacerRange[l], after = length(nucleotidesShorterThanProtospacerVec))
                queryStartVec = append(queryStartVec, queryStartRange[m], after = length(queryStartVec))
                prophageOnlyVec = append(prophageOnlyVec, prophageOnlyRange[n], after = length(prophageOnlyVec))

            }
          }
        }
      }
    }
  }
}

  upstreamLabels = c()
  downstreamLabels = c()
  for( i in 1:flankLength){
    upstreamLabels = append(upstreamLabels, sprintf("-%s", flankLength+1-i), after = length(upstreamLabels))
    downstreamLabels = append(downstreamLabels, sprintf("+%s", i), after = length(downstreamLabels))
  }

  for (n in counter:nConditions){

              uniqueAligns = uniqueAlignsVec[n]
              excludeSelf = excludeSelfVec[n]
              numGaps = numGapsVec[n]
              e.value = e.valueVec[n]
              nucleotidesShorterThanProtospacer = nucleotidesShorterThanProtospacerVec[n]
              queryStart = queryStartVec[n]
              prophageOnly = prophageOnlyVec[n]

              #File names generated based on systemName and value
              upstreamLogoOutput = sprintf("%s Protospacers Upstream %s.pdf", systemName,counter)
              downstreamLogoOutput = sprintf("%s Protospacers Downstream %s.pdf", systemName,counter)


              #Unique Aligns
              if (uniqueAligns == TRUE){
                uniqueSet = joinedData%>%
                  dplyr::distinct(query.acc.ver, TaxaID, .keep_all = T)
              }
              else{
                uniqueSet = joinedData
              }

              #Filter Step 1, Remove NAs and alignments to genome being analyzed
              if (excludeSelf == TRUE){
                withoutself = uniqueSet%>%
                  dplyr::filter(grepl(paste0(genus," ",species), as.character(Species))==FALSE | grepl(" phage", as.character(Species)) == TRUE | grepl("plasmid", as.character(Species))==TRUE)
              }
              else{
                withoutself = uniqueSet
              }

              #Filter step 2, filter based on alignment gaps
              nogap = withoutself %>%
                dplyr::filter(as.integer(gap.opens)<=numGaps)

              #Filter step 3, filter based on alignment E Value
              lowEValue = nogap%>%
                dplyr::filter(as.numeric(evalue)<=e.value)

              #Filter step 4, filter based on alignment length
              noshort = lowEValue %>%
                dplyr::mutate(PS.length = nchar(as.vector(Spacers))) %>%
                dplyr::filter(alignment.length >= (PS.length - nucleotidesShorterThanProtospacer))

              #Filter step 5, filter based on query start
              seedonly = noshort %>%
                dplyr::filter(q..start <= queryStart)

              #Retrieve prophage data from PHASTER.ca
              if (prophageOnly == TRUE){
                #Determine if alignment is to prophage in bacterial genome
                #Build empty dataframe with proper column names
                vecTitle = c("subject.acc.ver","REGION","REGION_LENGTH","COMPLETENESS(score)","SPECIFIC_KEYWORD","REGION_POSITION","TRNA_NUM","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM","HYPOTHETICAL_PROTEIN_NUM","PHAGE+HYPO_PROTEIN_PERCENTAGE","BACTERIAL_PROTEIN_NUM","ATT_SITE_SHOWUP","PHAGE_SPECIES_NUM","MOST_COMMON_PHAGE_NAME(hit_genes_count)","FIRST_MOST_COMMON_PHAGE_NUM","FIRST_MOST_COMMON_PHAGE_PERCENTAGE","GC_PERCENTAGE")
                cumResultsFrame = stats::setNames(data.frame(matrix(ncol = 18, nrow = 0)), vecTitle)

                #function to figure out which row has the dashes denoting the begining of the results
                recheck = function(){while(length(vertSplitPhaster)<5){
                  assign("phasterGet", httr::GET(url = sprintf("http://phaster.ca/phaster_api?acc=%s", conAccList[i])), envir = .GlobalEnv)
                  assign("rawPhaster", rawToChar(phasterGet$content), envir = .GlobalEnv)
                  assign("jsonPhaster", jsonlite::fromJSON(rawPhaster), envir = .GlobalEnv)
                  #split by \n
                  assign("vertSplitPhaster", strsplit(as.character(jsonPhaster), "\n"), envir = .GlobalEnv)
                  print(content(phasterGet)[[1]])
                  print(content(phasterGet)[[2]])
                  if (phasterGet$status_code != 200){goAhead()}
                  else {print("will check again in 2.5 minutes")
                    Sys.sleep(150)
                  }
                }}

                #function to proceed and add individual results to cumulative results
                goAhead = function(){
                  if(content(phasterGet)[[2]] =="No prophage region detected!\n"){
                    print(sprintf("%s has no prophages", conAccList[i]))
                  }
                  else{
                    for (j in (grep("---",vertSplitPhaster[[5]])+1):length(vertSplitPhaster[[5]])){
                      if (grep("---",vertSplitPhaster[[5]])+1 > length(vertSplitPhaster[[5]])){
                        print(sprintf("%s has no prophages", conAccList[i]))
                        break
                      }
                      #split by spaces
                      assign("horSplitPhasterEntry", strsplit(as.character(vertSplitPhaster[[5]][j]), " "), envir = .GlobalEnv)
                      #flatten list
                      assign("unlistedEntry", unlist(horSplitPhasterEntry), envir = .GlobalEnv)
                      #remove spaces in vector
                      assign("vecEntry", c(vertSplitPhaster[[1]],unlistedEntry[unlistedEntry!=""]), envir = .GlobalEnv)
                      #create dataframe with individual result
                      assign("phasterTable", stats::setNames(data.frame(t(data.frame(vecEntry))), vecTitle), envir = .GlobalEnv)
                      #bind individual result to all results for that genome
                      assign("indResultsFrame", rbind(indResultsFrame,phasterTable), envir = .GlobalEnv)
                    }
                    #bind individual result to cumulative result
                    assign("cumResultsFrame", rbind(indResultsFrame,cumResultsFrame), envir = .GlobalEnv)
                    #print confirmation of adding results
                    print(sprintf("result for %s complete", conAccList[i]))
                  }
                }

                #function to check if full result or still running
                check = function(){
                  if(phasterGet$status_code != 200){
                    recheck()
                  }
                  else {
                    goAhead()
                  }
                }

                #make consolidated list of accession numbers to push to Phaster
                conAccList = unique(seedonly$subject.acc.ver)

                #start for loop to go access phaster results from each phage
                for (i in 1:length(conAccList)){
                  #retreive data from phaster.ca using API
                  phasterGet = httr::GET(url = sprintf("http://phaster.ca/phaster_api?acc=%s", conAccList[i]))

                  while(phasterGet$status_code != 200){
                    print(sprintf("status code %s", phasterGet$status_code))
                    print("60 second pause to let server reset")
                    Sys.sleep(60)
                    phasterGet = httr::GET(url = sprintf("http://phaster.ca/phaster_api?acc=%s", conAccList[i]))
                  }

                  rawPhaster = rawToChar(phasterGet$content)
                  jsonPhaster = jsonlite::fromJSON(rawPhaster)
                  #split by \n
                  vertSplitPhaster = strsplit(as.character(jsonPhaster), "\n")

                  if(i %% 8==0){
                    print("60 sec pause to prevent Phaster Server Overload")
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

                  #initialize individual
                  indResultsFrame = stats::setNames(data.frame(matrix(ncol = 18, nrow = 0)), vecTitle)

                  #check if still running
                  check()

                }

                #split cumResultsFrame$REGION_POSTITION
                posSplit = cumResultsFrame %>%
                  dplyr::mutate(Region_Start = substring(cumResultsFrame$REGION_POSITION,1,regexpr("-", cumResultsFrame$REGION_POSITION)-1))%>%
                  dplyr::mutate(Region_Stop = substring(cumResultsFrame$REGION_POSITION,regexpr("-", cumResultsFrame$REGION_POSITION)+1,nchar(as.character(cumResultsFrame$REGION_POSITION))))

                #create dataframe with only subject.acc.ver and prophage regions
                finalPhaster = data.frame(posSplit$subject.acc.ver, posSplit$Region_Start, posSplit$Region_Stop)%>%
                  dplyr::rename(subject.acc.ver = posSplit.subject.acc.ver,
                                Region_Start = posSplit.Region_Start,
                                Region_stop = posSplit.Region_Stop)

                #join seedOnly and finalPhaster
                withProphageRegions = dplyr::left_join(seedonly,finalPhaster, by="subject.acc.ver")

                #Filter step 6, filter based on prophage identity
                isProphage = withProphageRegions%>%
                  dplyr::mutate(prophage = (as.numeric(s..start) >= as.numeric(sprintf("%s",Region_Start)))
                                & (as.numeric(s..start) <= as.numeric(sprintf("%s",Region_stop)))
                                & (as.numeric(s..end) >= as.numeric(sprintf("%s",Region_Start)))
                                & (as.numeric(s..end) <= as.numeric(sprintf("%s",Region_stop))))%>%
                  dplyr::filter(prophage == TRUE | (grepl("plasmid", withProphageRegions$Species)==TRUE))%>%
                  dplyr::distinct()
              }
              else{
                isProphage = seedonly
              }

              nrow(isProphage)

              #Retrieve genome sequences from NCBI

              finalAccList = unique(isProphage$subject.acc.ver)
              requestString = ""
              for (i in 1:length(finalAccList)){
                requestString = sprintf("%s,%s",requestString,finalAccList[i])
              }

              if (nrow(isProphage)>= 1){
                requestString = substr(requestString, 2, base::nchar(requestString))

                eFetchGet = httr::GET(url = sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta", requestString), add_headers("http_version" = "http_version$1"))
                raweFetch = rawToChar(eFetchGet$content)
                writeLines(raweFetch, "eFetch FASTA.fasta")

                #Read in Entrez generated FASTA
                trimFASTA = seqinr::read.fasta(file = "eFetch FASTA.fasta", as.string = TRUE)
              }
              else {seqinr::write.fasta(sequences = c(),names = c(), file.out = "eFetch FASTA.fasta")
                trimFASTA = seqinr::read.fasta(file = "eFetch FASTA.fasta", as.string = TRUE)}

               #Take data from seqinr list and put into dataframe
              subject.acc.ver = c()
              for (i in 1:length(trimFASTA)){
                subject.acc.ver[i] = attr(trimFASTA[[i]], "name")
              }

              genomeSequence = c()
              for (i in 1:length(trimFASTA)){
                genomeSequence[i] = unlist(trimFASTA[[i]])
              }

              annotation = c()
              for (i in 1:length(trimFASTA)){
                annotation[i] = attr(trimFASTA[[i]], "Annot")
              }

              accSeq = data.frame(subject.acc.ver, annotation, genomeSequence)

              #Join isProphage dataframe and accSeq dataframe
              withGenome = isProphage %>%
                dplyr::left_join(accSeq, by="subject.acc.ver")

              #Determine if alignment is on plus or minus strand (plus is TRUE, minus is FALSE)
              whichStrand = function (start, end)
              {
                ifelse(end - start > 0, TRUE, FALSE)
              }

              stranded = withGenome %>%
                dplyr::mutate(strand = whichStrand(s..start, s..end))

              #Extract flankLength nucleotides upstream and downstream of alignment from genome sequence based on direction of alignment

              upStrandAccount = function (seq, strand, start, length)
              {
                ifelse(strand == TRUE, substr(seq, start-flankLength,start-1), substr(seq, start-length-flankLength+1,start-length))
              }

              downStrandAccount = function (seq, strand, start, length)
              {
                ifelse(strand == TRUE, substr(seq, start+length, start+length+flankLength-1), substr(seq, start+1, start+flankLength))
              }

              revcomplement = function (strand, flank){
                ifelse(strand == FALSE, sapply(lapply(strsplit(chartr("atgc","tacg",as.character(flank)), NULL), rev), paste, collapse=""), sapply(lapply(strsplit(chartr("atgc","tacg",as.character(flank)), NULL), rev), paste, collapse=""))
              }

              upnDown = stranded %>%
                dplyr::mutate(upstream = upStrandAccount(genomeSequence, strand, s..start, PS.length))%>%
                dplyr::mutate(downstream = downStrandAccount(genomeSequence, strand, s..start, PS.length))%>%
                dplyr::mutate(upstream.rev = revcomplement(strand, upstream))%>%
                dplyr::mutate(downstream.rev = revcomplement(strand, downstream))

              #Generate vector of sequences for alignment
              whichFlank = function(orientation, strand, flank, revcompflank){
                ifelse((orientation == "Forward") == strand, flank, revcompflank)
              }

              alignmentUp = c()
              if(nrow(upnDown)>=1){
                for (i in 1:nrow(upnDown)){
                  alignmentUp[i] = as.character(whichFlank(upnDown$Array.Orientation[i],upnDown$strand[i],upnDown$upstream[i],upnDown$downstream.rev[i]))
                }
              }
              else{alignmentUp = c("")}

              alignmentDown = c()
              if(nrow(upnDown)>=1){
                for (i in 1:nrow(upnDown)){
                  alignmentDown[i] = as.character(whichFlank(upnDown$Array.Orientation[i],upnDown$strand[i],upnDown$downstream[i],upnDown$upstream.rev[i]))
                }
              }
              else{alignmentDown = c("")}

              #Calculate PAM Scores
              #Generate frequency dataframe
              upfreqFrame = as.data.frame(matrix(nrow = flankLength, ncol = 5))
              names(upfreqFrame) = c("Position", "fa", "ft", "fc", "fg")
              upfreqFrame$Position = c(1:flankLength)

              for (i in 1:flankLength){
                nucs = c()
                for (j in 1:length(alignmentUp)){
                  nucs = append(nucs, substr(alignmentUp[j], i, i), after = length(nucs))
                }
                tempFrame = as.data.frame(table(nucs))
                names(tempFrame) = c("posNucs", "freq")
                posNucs = c("a", "t", "c","g")
                posNucsFrame = data.frame(posNucs)
                uptotalNucsFrame = dplyr::left_join(posNucsFrame, tempFrame, by = "posNucs")
                uptotalNucsFrame[is.na(uptotalNucsFrame)] = 0
                upfreqFrame$fa[i] = uptotalNucsFrame$freq[1]/ length(alignmentUp)
                upfreqFrame$ft[i] = uptotalNucsFrame$freq[2]/ length(alignmentUp)
                upfreqFrame$fc[i] = uptotalNucsFrame$freq[3]/ length(alignmentUp)
                upfreqFrame$fg[i] = uptotalNucsFrame$freq[4]/ length(alignmentUp)
              }

              #Calculate entropy and R score for each position
              upwithIndEntropy = upfreqFrame %>%
                dplyr::mutate(Hai = (upfreqFrame$fa * log2(upfreqFrame$fa)))%>%
                dplyr::mutate(Hti = (upfreqFrame$ft*log2(upfreqFrame$ft)))%>%
                dplyr::mutate(Hci = (upfreqFrame$fc*log2(upfreqFrame$fc)))%>%
                dplyr::mutate(Hgi =(upfreqFrame$fg*log2(upfreqFrame$fg)))

              is.nan.data.frame <- function(x)
                do.call(cbind, lapply(x, is.nan))

              upwithIndEntropy[is.nan(upwithIndEntropy)] = 0

              upwithTotalEntropy = upwithIndEntropy%>%
                dplyr::mutate(Hi = -1*(Hai+Hti+Hci+Hgi))

              upwithRScore = upwithTotalEntropy%>%
                dplyr::mutate(R = log2(4)-(Hi+((1/log(2))*((4-1)/(2*length(alignmentUp))))))

              upmeanRScore = mean(upwithRScore$R)
              upsdRScore = stats::sd(upwithRScore$R)

              #Significant position as defined by position R Score >= mean R score + 1/2 standard deviation of R scores
              upisSigPos = upwithRScore%>%
                dplyr::mutate(Significant = (R >=upmeanRScore+(0.5*upsdRScore)))

              ###CONTINUE PAM SCORE CALCULATIONS, NEED Hi AVERAGE and Havgdev
              upHavg = sum(upisSigPos$Hi)/flankLength

              upsigHi = c()
              for(i in 1:flankLength){
                if (upisSigPos$Significant[i] == TRUE){
                  upsigHi = append(upsigHi, upisSigPos$Hi[i], after= length(upsigHi))
                }
              }

              upHavgdev = sum(abs((upsigHi-upHavg)))/sqrt(length(upsigHi)+1)
              #uppamScore = sqrt(length(alignmentUp))*(1-upHavgdev)
              #experimantal PAM Score
              uppamfreqVec = c()
              uppamRVec = c()
              for (i in 1:flankLength){
                if (upisSigPos$Significant[i] == TRUE){
                  uppamfreqVec = append(uppamfreqVec, paste(upisSigPos[i,apply(upisSigPos[i,2:5],1,function(x) which(x>0.25))+1]), after= length(uppamfreqVec))
                  for(j in 1:length(paste(upisSigPos[i,apply(upisSigPos[i,2:5],1,function(x) which(x>0.25))+1]))){
                    uppamRVec = append(uppamRVec, upisSigPos$R[i], after= length(uppamRVec))
                  }
                }
              }
              uppamScore = length(alignmentUp)*(sum(as.numeric(uppamfreqVec)*as.numeric(uppamRVec))/(log2(4)*length(upsigHi)))


              uppamSeq = c()
              for (i in 1:flankLength){
                if (upisSigPos$Significant[i] == TRUE){
                  uppamSeq = append(uppamSeq, sprintf("%s", paste(toupper(substr(colnames(upisSigPos)[apply(upisSigPos[i,2:5],1,function(x) which(x>0.25))+1],2,2)), collapse = "/"), after= length(uppamSeq)))
                }
                else {
                  uppamSeq = append(uppamSeq, "N", after= length(uppamSeq))

                }
              }

              print(sprintf("%s %s with a PAM Score of: %s", paste(c("The upstream consensus PAM #",counter), collapse=""),paste(c( "is: ", uppamSeq), collapse=""), uppamScore))

              downfreqFrame = as.data.frame(matrix(nrow = flankLength, ncol = 5))
              names(downfreqFrame) = c("Position", "fa", "ft", "fc", "fg")
              downfreqFrame$Position = c(1:flankLength)

              for (i in 1:flankLength){
                nucs = c()
                for (j in 1:length(alignmentDown)){
                  nucs = append(nucs, substr(alignmentDown[j], i, i), after = length(nucs))
                }
                tempFrame = as.data.frame(table(nucs))
                names(tempFrame) = c("posNucs", "freq")
                posNucs = c("a", "t", "c","g")
                posNucsFrame = data.frame(posNucs)
                downtotalNucsFrame = dplyr::left_join(posNucsFrame, tempFrame, by = "posNucs")
                downtotalNucsFrame[is.na(downtotalNucsFrame)] = 0
                downfreqFrame$fa[i] = downtotalNucsFrame$freq[1]/ length(alignmentDown)
                downfreqFrame$ft[i] = downtotalNucsFrame$freq[2]/ length(alignmentDown)
                downfreqFrame$fc[i] = downtotalNucsFrame$freq[3]/ length(alignmentDown)
                downfreqFrame$fg[i] = downtotalNucsFrame$freq[4]/ length(alignmentDown)
              }

              #Calculate entropy and R score for each position
              downwithIndEntropy = downfreqFrame %>%
                dplyr::mutate(Hai = (downfreqFrame$fa * log2(downfreqFrame$fa)))%>%
                dplyr::mutate(Hti = (downfreqFrame$ft*log2(downfreqFrame$ft)))%>%
                dplyr::mutate(Hci = (downfreqFrame$fc*log2(downfreqFrame$fc)))%>%
                dplyr::mutate(Hgi =(downfreqFrame$fg*log2(downfreqFrame$fg)))

              is.nan.data.frame <- function(x)
                do.call(cbind, lapply(x, is.nan))

              downwithIndEntropy[is.nan(downwithIndEntropy)] = 0

              downwithTotalEntropy = downwithIndEntropy%>%
                dplyr::mutate(Hi = -1*(Hai+Hti+Hci+Hgi))

              downwithRScore = downwithTotalEntropy%>%
                dplyr::mutate(R = log2(4)-(Hi+((1/log(2))*((4-1)/(2*length(alignmentDown))))))

              downmeanRScore = mean(downwithRScore$R)
              downsdRScore = stats::sd(downwithRScore$R)

              #Significant position as defined by position R Score >= mean R score + 1/2 standard deviation of R scores
              downisSigPos = downwithRScore%>%
                dplyr::mutate(Significant = (R >=downmeanRScore+(0.5*downsdRScore)))

              ###CONTINUE PAM SCORE CALCULATIONS, NEED Hi AVERAGE and Havgdev
              downHavg = sum(downisSigPos$Hi)/flankLength

              downsigHi = c()
              for(i in 1:flankLength){
                if (downisSigPos$Significant[i] == TRUE){
                  downsigHi = append(downsigHi, downisSigPos$Hi[i], after= length(downsigHi))
                }
              }

              downHavgdev = sum(abs((downsigHi-downHavg)))/sqrt(length(downsigHi)+1)
              #downpamScore = sqrt(length(alignmentDown))*(1-downHavgdev)
              #experimantal PAM Score
              downpamfreqVec = c()
              downpamRVec = c()
              for (i in 1:flankLength){
                if (downisSigPos$Significant[i] == TRUE){
                  downpamfreqVec = append(downpamfreqVec, paste(downisSigPos[i,apply(downisSigPos[i,2:5],1,function(x) which(x>0.25))+1]), after= length(downpamfreqVec))
                  for(j in 1:length(paste(downisSigPos[i,apply(downisSigPos[i,2:5],1,function(x) which(x>0.25))+1]))){
                    downpamRVec = append(downpamRVec, downisSigPos$R[i], after= length(downpamRVec))
                  }
                }
              }
              downpamScore = length(alignmentDown)*(sum(as.numeric(downpamfreqVec)*as.numeric(downpamRVec))/(log2(4)*length(downsigHi)))

              downpamSeq = c()
              for (i in 1:flankLength){
                if (downisSigPos$Significant[i] == TRUE){
                  downpamSeq = append(downpamSeq, sprintf("%s", paste(toupper(substr(colnames(downisSigPos)[apply(downisSigPos[i,2:5],1,function(x) which(x>0.25))+1],2,2)), collapse = "/"), after= length(downpamSeq)))
                }
                else {
                  downpamSeq = append(downpamSeq, "N", after= length(downpamSeq))

                }
              }

              print(sprintf("%s %s with a PAM Score of: %s", paste(c("The downstream consensus PAM #",counter), collapse=""),paste(c( "is: ", downpamSeq), collapse=""), downpamScore))

              #Plot upstream WebLogo and save
              if(nrow(upnDown)>=1 && saveLogo == T){
                ggplot2::ggplot()+ggseqlogo::geom_logo( as.character(toupper(alignmentUp)), seq_typ="DNA")+ ggplot2::annotate('text', x=ceiling(flankLength/2), y=2.1, label=sprintf("Consensus: %s Score: %s", paste(uppamSeq, collapse=""),round(uppamScore, digits = 2)), hjust = 0.5, size = 12)+ggseqlogo::theme_logo()+ggplot2::scale_x_discrete(name ="Position",limits=upstreamLabels)+ggplot2::theme(axis.text.x = element_text(size=32),axis.text.y = element_text(size=32), axis.title=element_text(size=48) ,axis.ticks = element_line(size = 1.5), axis.ticks.length = unit(20, "pt"))
                ggplot2::ggsave(upstreamLogoOutput, width = 10, height = 7, units = "in")



                #Plot downstream WebLogo
                ggplot2::ggplot()+ggseqlogo::geom_logo( as.character(toupper(alignmentDown)), seq_typ="DNA")+ ggplot2::annotate('text', x=ceiling(flankLength/2), y=2, label=sprintf("Consensus: %s Score: %s", paste(downpamSeq, collapse=""),round(downpamScore, digits = 2)), hjust = 0.5, size = 12)+ggseqlogo::theme_logo()+ggplot2::scale_x_discrete(name ="Position",limits=downstreamLabels)+ggplot2::theme(axis.text.x = element_text(size=32),axis.text.y = element_text(size=32), axis.title=element_text(size=48) ,axis.ticks = element_line(size = 1.5), axis.ticks.length = unit(20, "pt"))
                ggplot2::ggsave(downstreamLogoOutput, width = 10, height = 7, units = "in")
              }


              #Add filter data to data frame
              collectionFrame$uniqueAligns[counter] = uniqueAligns
              collectionFrame$excludeSelf[counter] = excludeSelf
              collectionFrame$numGaps[counter] = numGaps
              collectionFrame$e.value[counter] = e.value
              collectionFrame$nucleotidesShorterThanProtospacer[counter] = nucleotidesShorterThanProtospacer
              collectionFrame$queryStart[counter] = queryStart
              collectionFrame$prophageOnly[counter] = prophageOnly
              collectionFrame$Filter0[counter] = nrow(uniqueSet)
              collectionFrame$Filter1[counter] = nrow(withoutself)
              collectionFrame$Filter2[counter] = nrow(nogap)
              collectionFrame$Filter3[counter] = nrow(lowEValue)
              collectionFrame$Filter4[counter] = nrow(noshort)
              collectionFrame$Filter5[counter] = nrow(seedonly)
              collectionFrame$Filter6[counter] = nrow(isProphage)
              collectionFrame$upPAM[counter] = paste(c(uppamSeq), collapse="")
              collectionFrame$upScore[counter] = uppamScore
              collectionFrame$downPAM[counter] = paste(c(downpamSeq), collapse="")
              collectionFrame$downScore[counter] = downpamScore

              counter = counter+1
              Sys.sleep(3)

            }
  assign("collectionFrame", collectionFrame, envir= .GlobalEnv)

  if(removeFASTA == T){
  file.remove("eFetch FASTA.fasta")
  }

  if(savePAMSeqs == T){
    assign("upstreamPAMSeqs", alignmentUp, envir= .GlobalEnv)
    assign("downstreamPAMSeqs", alignmentDown, envir= .GlobalEnv)
  }
}
