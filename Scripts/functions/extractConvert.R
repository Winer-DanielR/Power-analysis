# I wrote all of these functions
# Function converts data from list of lists to data frame
# and then converts to long form for graphing
# takes as input:
# data = extracted dataset 
# output is long form data frame
convert.long <- function(data){
  #convert to dataframe
  data <- as.data.frame(do.call(rbind, lapply(data, as.vector)))
  # move replicates from rows to columns
  data <- rownames_to_column(data, var = "replicate")
  # gather all effect data into two columns (effects and eval values)
  data <- gather(data, `5`, `10`, `15`, `20`, `25`, key="effect", value="eval")
  #transform eval column into numeric
  data <- transform(data, eval=as.numeric(eval))
  #transform effect column into factor with ordered levels
  data$effect <- factor(data$effect, levels=c('5','10','15','20','25'))
  return(data)
  }

# Function extracts data for each evaluation paramater 
# (True positives, specificity, etc)
# and converts to proper long form dataset 
# takes as input:
# data = results from model
# replicates = list of different numbers of replicates modelled
# effects = list of different numbers of effect sizes modelled
# eval = the evaluation paramter to extract
# output is long form dataset 
extract.data <- function(data, replicates, effects, eval){
  #c reate a list with the replicates as names 
  data <-setNames(replicate(length(replicates), list()), replicates)
  # create a list of indexes for every replicate
  indexes = c(1:length(replicates))
  # go through each index in index list
  for (index in indexes){
    #set a second index
    index2 = 1
    # go through every effect
    for (effect in effects){
      # add the data for each effect at the right index
      # by searching the original dataset using the indexes 
      # and the name of the parameter you are extracting 
      data[[index]][as.character(effect)] = model.res[[index]][[index2]][eval]
      index2 <- index+1
    }
  }
  # convert to long data frame and return
  return(convert.long(data))
}