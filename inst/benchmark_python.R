library(tidync)
pythondir <- '/Net/Groups/BGI/people/twutz/devpy/ecosystem-transpiration'
ecf <- tidync(file.path(pythondir,"ec_tea_tutorial.nc"))
ecf
ec <- ecf %>% hyper_tibble()

df <- tea_preprocess(FIHyy)

summary(ec$CSWI - df$CSWI)
all(abs(ec$CSWI - df$CSWI) < 1e-8)

#plot(ec$CSWI ~ df$CSWI)
