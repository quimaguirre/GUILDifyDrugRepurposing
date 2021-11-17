#!/usr/bin/env Rscript
library(guildifyR)
library("optparse")

option_list = list(
	make_option(c("-k", "--keywords"), type="character", default=NULL, 
				help="disease name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
            	help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$keywords)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

species="9606"
tissue="All"
result.table = query("alzheimer", species, tissue)
job.id = submit.job(result.table, species, tissue, list(netscore=T, repetitionSelector=1, iterationSelector=1))
result = retrieve.job(job.id)
result = retrieve.job(job.id, n.top=120, fetch.files=T, output.dir="./")
result.table = query("lung neoplasms", species, tissue)
job.id2 = submit.job(result.table, species, tissue, list(netscore=T, repetitionSelector=1, iterationSelector=1))
result = retrieve.overlap(job.id, job.id2)
keywords = "Sfn;Krt6b;Krt12;Krt21"
result.table = query(keywords, species, tissue)
job.id = submit.job(result.table, species, tissue, list(netcombo=T))
keywords = "nonsensemakingtestquery"
result.table = query(keywords, species, tissue)

job.id = "4a720f37-de37-434a-8334-0dee3b702d9c"
result = retrieve.job(job.id)
result = retrieve.overlap("format% error", "lung neoplasms")
result = retrieve.overlap("alzheimer", "lung neoplasms", fetch.files=T)

species="10090"
tissue="All"
result.table = query("oxidative", species, tissue)
job.id = submit.job(result.table, species, tissue, list(netscore=T, repetitionSelector=3, iterationSelector=2))
result = retrieve.job(job.id)
result = retrieve.job(job.id, fetch.files=T, output.dir="test")
