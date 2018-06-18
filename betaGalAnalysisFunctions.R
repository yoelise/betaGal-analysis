meanExpressionOneReplicate <- function(oneFileName) {

  OD420data <- read_excel(oneFileName, sheet = 2)
  OD600data <- read_excel(oneFileName, sheet = 3)
  
  meanCellDensities <-  filter(OD600data, strain != "blank") %>%
    mutate(trueDensity = `Blank 600`*as.numeric(dilutionFactor)) %>%
    group_by(strain, RpoS) %>% 
    summarise(meanDensity = mean(trueDensity)) 
  
  combinedData <-left_join(OD420data, meanCellDensities, by=c("strain", "RpoS")) %>%
    select(`Blank 420`, experimenter, date, strain, RpoS, dilutionFactor, volumeAssayed, reactionTime, meanDensity) %>%
    filter(strain != "blank") %>%
    rename(A420 = `Blank 420`) %>%
    type_convert() %>%
    mutate(`Miller Units` = 1000*A420*dilutionFactor/(meanDensity*volumeAssayed*reactionTime))
  
  expressionLevels <- group_by(combinedData, strain, RpoS, experimenter, date) %>% 
    summarise(meanExpression = mean(`Miller Units`))
  expressionLevels
}


calculateSingleSens <- function(RpoS, meanExpression) {
  RpoSLevels <- sort(unique(RpoS))
  maxRpoS <- RpoSLevels[3]
  midRpoS <- RpoSLevels[2]
  minRpoS <- RpoSLevels[1]
  
  rise <- meanExpression[RpoS == maxRpoS] - meanExpression[RpoS == minRpoS]
  run <- maxRpoS - minRpoS
  slope <- rise/run
  predicted <- slope*midRpoS + meanExpression[RpoS == minRpoS] #This is the expected value on the line
  
  dif <- meanExpression[RpoS == midRpoS] - predicted #This is the difference between the observed - predicted
  
  sens <- dif/rise #This is A/B, the metric for sensitivity.
  
  sens
}
