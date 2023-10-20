library(readxl)
primers <- read_xlsx("ncbi_ribo_test_primerAnnotations.xlsx")
hist(primers$Minimum, breaks = 100)
median(primers$Minimum)
