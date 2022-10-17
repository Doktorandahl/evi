library(devtools)
library(usethis)

use_package('dplyr')
use_package('maxLik')

load_all()
document()

hks <- readRDS("/Users/davra945/Documents/PhD/pl_reg/hks_jvdr.rds")%>% mutate(brv_AllLag = log(brv_AllLag + 1),
                                                                              troopLag_log = log(troopLag+1),
                                                                              epdur_log = log(epduration),
                                                                              policeLag_log = log(policeLag+1),
                                                                              militaryobserversLag_log = log(militaryobserversLag +1))

use_data(hks,overwrite = T)

f1 <- osvAll ~ troopLag + policeLag + militaryobserversLag + brv_AllLag + osvAllLagDum + incomp + epduration + lntpop
f2 <- osvAll ~ troopLag_log + epdur_log + brv_AllLag + lntpop
tt <- run_evzinb(formula_nb = f1,
           formula_zi = f1,
           formula_evi = f1,
           formula_pareto = f2,
           data=hks)





