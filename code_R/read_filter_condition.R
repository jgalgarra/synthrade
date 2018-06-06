filtered_file  <-  read.table("filtered_condition.txt",header=TRUE)
if (filtered_file$Condition == 1){
  filtered_string <- "_FILT" # if set to "_FILT" uses filtered matrix
} else
  filtered_string <- "" # if set to "" uses unfiltered matrix