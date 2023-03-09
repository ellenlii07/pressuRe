source("R/pressuRe_functions.R")

filepath1 <- "example_data/emed test.lst"
filepath2 <- "C:/Users/telfe/OneDrive - UW/CoRE Projects/Insoles/test subjects/S01/emed files/test__0912_1.lst"

f1 <- load_emed(filepath1)
f2 <- load_emed2(filepath1)
f3 <- load_emed2(filepath2)

identical(f1, f2)
