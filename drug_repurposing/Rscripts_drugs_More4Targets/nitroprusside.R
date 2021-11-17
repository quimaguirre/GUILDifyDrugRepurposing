library(guildifyR)
library(tidyr)
drug_databases = c("drugbank", "dgidb", "drugcentral", "chembl")
result.table = query("nitroprusside", species="9606", tissue="all", network.source = "BIANA")
sep.result.table = separate_rows(result.table, source, sep=", ")
filt.result.table = sep.result.table[ which( sep.result.table$source %in% drug_databases & sep.result.table$in.network==1 ), ]
submit.job(filt.result.table, species="9606", tissue="all", network.source="BIANA", scoring.options = list(netscore=T, repetitionSelector=3, iterationSelector=2))
