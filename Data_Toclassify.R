library("readxl")
TOclassify <- read_excel("TOclassify.xlsx", col_names = T) 

gumunifdata <- read_excel("gumunifdata.xlsx", col_names = T)
norunifdata <- read_excel("norunifdata.xlsx", col_names = T)
tunifdata <- read_excel("tunifdata.xlsx", col_names = T)

lstX <- list(gumunifdata, norunifdata, tunifdata)
lstX <- lapply(lstX, as.matrix)
TOclassify <- as.matrix(TOclassify)
mul.discr.eigen.reg(lstX, TOclassify, w = rep(1, 3)) 

