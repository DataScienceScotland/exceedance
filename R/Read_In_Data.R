
if (!is.element('readxl', installed.packages()[,1]))
  install.packages(p, dep = TRUE)
require('readxl', character.only = TRUE)


read<-function(filename){
  cat(c("Loading data from ", filename, ".\n"))
  col_types = c("date", rep("numeric", 34))
  Positive = read_xlsx(filename, sheet = "Positive Local Auth", col_types = col_types)
  dP = dim(Positive)[1]
  Positive = as.data.frame(Positive)[,-35]
  names(Positive)[2] = 'Null'

  Negative = read_xlsx(filename, sheet = "Negative Local Auth", col_types = col_types)
  dN = dim(Negative)[1]
  Negative = as.data.frame(Negative)[,-35]
  names(Negative)[2] = 'Null'

  Place_names = names(Positive)[-1]
  Both = merge(Negative, Positive, by = 'SpecimenDate', all = TRUE,  suffixes = c(".Negative", ".Positive"))
  SpecimenDate = Both$SpecimenDate

  Dates = data.frame(
    seq(min(SpecimenDate[-1]), max(SpecimenDate), "days")
  )
  names(Dates) = 'SpecimenDate'
  Both = merge(Both, Dates,
    by = 'SpecimenDate', all = TRUE)

  last_date = min(max(Negative$`SpecimenDate`), max(Positive$`SpecimenDate`))
  cat("Using data until ", as.character(last_date), ".\n")
  idx = as.Date(Both$`SpecimenDate`) <= as.Date(last_date)

  Both = Both[idx, ]



  rownames(Both) = Both$`SpecimenDate`
  Both$`SpecimenDate` = NULL



  Both[is.na(Both)] = 0
  col.idx.pos = sapply(Place_names, paste0, '.Positive')
  col.idx.neg = sapply(Place_names, paste0, '.Negative')

  Negatives = Both[,col.idx.neg]
  names(Negatives) = Place_names

  Positives = Both[,col.idx.pos]
  names(Positives) = Place_names

  save("Negatives","Positives", "SpecimenDate", "Place_names", "filename", file = 'Scottish_Data.RData')
  cat("Data from ", filename, " saved in `Scottish_Data.RData`.\n")


}
