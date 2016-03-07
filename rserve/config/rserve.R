## This script will be loaded when Rserve starts
## and before any client connect.
## Use it to pre-load packages and any data you want
## as well as define any global variables you want all
## scripts to see

## Today, pretty much everyone speaks UTF-8, it makes the life easier
Sys.setlocale(,"en_US.UTF-8")

## it is useful to have access to the root of your
## installation from R scripts
root <- Sys.getenv("ROOT")
if (is.null(root) || nchar(root) == 0) root <- "/Users/tjc29/"

## run the server in the "tmp" directory of the root in
## case some files need to be created
#setwd(paste(root,"tmp",sep='/'))

## if you have multiple servers it's good to know which machine this is
host <- tolower(system("hostname -s", TRUE))
cat("Starting Rserve on", host,"\n")

## This is jsut a friendly way to load package and report success/failure
## You will definiteily need FastRWeb, others are optional
.libPaths( c( .libPaths(), "/home/tjc29/R/x86_64-pc-linux-gnu-library/3.1/") )
pkgs <- c("vcf2ld", "Matrix", "snpStats", "rjson", "data.table")

cat("Loading packages...\n")
for (pkg in pkgs) cat(pkg, ": ",require(pkg, quietly=FALSE, character.only=TRUE),"\n",sep='')

# RDS data
cat("Loading 1000 genome data...\n")
KG_DIR='/usr/share/rserve/data/'
CURRENT_BUILD='grch38'
snp_data_src <- c("EUR")
for (src in snp_data_src) {
  cat(src, "\n")
  data.files <- list.files(path=paste0(KG_DIR,'/',CURRENT_BUILD,'/',src), pattern='*.rds$', full.names = TRUE)
  tmp <- lapply(data.files, function(x){readRDS(x)})
  names(tmp) <- gsub('.*(chr[^\\.]+).*', '\\1', basename(data.files))
  assign(src, tmp)
  rm(tmp)
}
cat("Loaded data\n")

# read population and samples
tmp <- read.table(paste0(KG_DIR,'/',CURRENT_BUILD,'/20131219.populations.tsv'),header=TRUE,sep = "\t",colClasses = c(rep("NULL", 1), rep("character", 2), rep("NULL", 6)))
tmp[tmp==""]<-NA
tmp<-na.omit(tmp)
POPS<-data.frame(tmp[,-1], row.names=tmp[,1])
rm(tmp)
SAMPLES <- read.table(paste0(KG_DIR,'/',CURRENT_BUILD,'/integrated_call_samples_v3.20130502.ALL.panel'),row.names=1,header=TRUE)

cat("Done\n")

## init() is a special function that will be called from
## each script. Do what you want here - it is usually a good idea
## to have a "common" script that is loaded on each request
## so you don't need re-start Rserve for global code changes
init <- function() {
    set.seed(Sys.getpid()) # we want different seeds so we get different file names
    
    ## get a temporary file name for this session
    tmpfile<<-paste('tmp-',paste(sprintf('%x',as.integer(runif(4)*65536)),collapse=''),'.tmp',sep='')

    ## if there is a common script, source it first
    common <- paste(root,"/web.R/common.R",sep='')
    if (isTRUE(file.exists(common))) source(paste(root,"/web.R/common.R",sep=''))
}

