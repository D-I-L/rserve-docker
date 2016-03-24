#' Get markers in LD.
#' library(data.table)
#' library(snpStats)
#' library(rjson)
#' KG_DIR='/Users/tjc29/1000genome/'
#' CURRENT_BUILD='grch38'
#' ld_run('EUR',1,'rs2476601')
#' 
#' @param chromosome, marker chromosomal location
#' @param dataset to use for calculation (super or sub-population e.g. EUR or CEU)
#' @param marker1 marker to identify markers in LD for
#' @param marker2 optional marker to calculate LD for marker1, Defaults to NULL.
#' @param window_size window size to look for markers in LD, Defaults to 1000000.
#' @param dprime dprime to use, Defaults to 0
#' @param rsq R square threshold, Defaults to 0.8
#' @param maf if TRUE report the MAF in the result
#' @param position if TRUE report the position in the result
#' @param m1_pos only required when reading from VCF
#' run()

ld_run <- function(dataset, chromosome, marker1, marker2=NULL,
		           window_size=1000000, dprime=0, rsq=0.8, maf=FALSE,
				   position=FALSE, read_from_rds=TRUE, m1_pos=0,
				   build_version=CURRENT_BUILD) {
  # ptm <- proc.time()

  snp_ids = NULL
  
  # find if dataset is a super- or sub-population
  super_pop<-as.character(POPS[dataset,])
  sub_pop<-NA
  if(is.na(super_pop)){
    super_pop<-dataset # calculation over super-population
  }else{
    sub_pop<-dataset   # calculation over sub-population
  }
  # look for pre-loaded snpMatrix
  sm <- try(eval(parse(text = super_pop))[[paste('chr',chromosome,sep="")]])
  
  if(length(sm) == 1 || build_version != CURRENT_BUILD) {
    rds_file <- list.files(paste0(KG_DIR,'/',build_version,'/',super_pop), paste0("ALL\\.",'chr',chromosome,'\\..*',super_pop,'.*.rds$'), full.names=TRUE)
    if(length(rds_file) > 0 && file.exists(rds_file) && read_from_rds) {
       print('READ RDS')
       sm<-readRDS(rds_file)
       m1_pos = sm$map[sm$map$ID == marker1,]$POS
       snp_ids <-sm$map[sm$map[,'POS'] > (m1_pos-window_size) & sm$map[,'POS'] < (m1_pos+window_size), ]$ID
    } else {
       print(paste('READ VCF ', marker1, m1_pos))
       startPos = m1_pos - window_size
       endPos = m1_pos + window_size
       region = paste(startPos,endPos,sep='-')
       region<-paste(chromosome, region, sep=':')
       sm <- vcf2sm(region,super_pop)
    }
  } else {
    print('READ FROM PRE-LOADED snpMatrix')
    m1_pos = sm$map[sm$map$ID == marker1,]$POS
    snp_ids <-sm$map[sm$map[,'POS'] > (m1_pos-window_size) & sm$map[,'POS'] < (m1_pos+window_size), ]$ID
  }

  pos=sm$map[sm$map$ID==marker1,]$POS

  if(length(pos) == 0){
    err <- paste("Marker ",marker1," was not found in this dataset",sep="")
  } else {
  	if(length(marker2) != 0) {
  		snp_ids <-c(marker1,marker2)
  		snp_data_subset <-sm$gt[,snp_ids]
  	} else if(length(snp_ids) != 0) {
  	  snp_data_subset <-sm$gt[,snp_ids]
  	} else {
  		snp_data_subset <-sm$gt
  	}
	
	# filter by sub-population if requested
    if(!is.na(sub_pop)){
      sub_pop_samples <- row.names(subset(SAMPLES, pop==sub_pop))
      snp_data_subset<-snp_data_subset[sub_pop_samples,]
    }
    #Filter based on some quality scores
    snp_stats       <-col.summary(snp_data_subset)
    snp_data_subset <-snp_data_subset[, snp_stats$MAF >= 0.01 & snp_stats$Call.rate >= 0.90 & snp_stats$z.HWE^2 < 25]

    #Check that snp makes it through quality control
    if (marker1 %in% colnames(snp_data_subset)) {
      ld_data <- ld(snp_data_subset, snp_data_subset[,marker1], stats=c("D.prime", "R.squared"))
      ld_data_formatted <- data.frame(  rownames(ld_data$D.prime), ld_data )
      colnames(ld_data_formatted) = c("marker2", names(ld_data))
      #remove the instance where marker is compared with self
      ld_data_formatted <- ld_data_formatted[ marker1!=rownames(ld_data_formatted), ]
      ld_data_formatted <- na.omit(ld_data_formatted)
      #format R squared before filtering
      ld_data_formatted$R.squared <- round(ld_data_formatted$R.squared, digits=2)
      ld_data_formatted <- ld_data_formatted[ ld_data_formatted$R.squared>=rsq, ]
  
      #print(ld_data_formatted)
      msg <- ld_data_formatted[ ld_data_formatted$D.prime>=dprime, ]
      result_markers <- as.character(msg$marker2)
      if(maf) {
        #print(snp_stats[result_markers,])
        #print(sm$map[sm$map$ID == result_markers,])
        
        # 
        # RAF: The "risk" allele (allele B) frequency. So if 
        # RAF < 0.5 allele B is the minor allele otherwise A.
        #for (i in 1:NROW(snp_stats[result_markers,])) {
          #print(sm$map[sm$map$ID == result_markers,]$ALT[i])
          #print(sm$map[sm$map$ID == result_markers,]$REF[i])
          
          #if(snp_stats[result_markers,]$RAF[i] < 0.5) {
          #  msg$minor[i] <- sm$map[sm$map$ID == result_markers,]$ALT[i]
          #  msg$major[i] <- sm$map[sm$map$ID == result_markers,]$REF[i]
          #} else {
          #  msg$minor[i] <- sm$map[sm$map$ID == result_markers,]$REF[i]
          #  msg$major[i] <- sm$map[sm$map$ID == result_markers,]$ALT[i]
          #}
        #}
        msg$MAF <- c(round(snp_stats[result_markers,]$MAF, digits=7))
      }
     if(position) {
       msg$position <- c(sm$map[sm$map$ID %in% result_markers,]$POS)
     }
    } else {
      err <- paste("Marker ",marker1," did not pass QC",sep="")
    }
  }
  # print(proc.time() - ptm)
  if(exists("err")) {
    toJSON(list(error=err))
  } else {
    # make a list of lists
    lol <- do.call(Map, c(list, msg))
    toJSON(list(ld=lol))
  }
}


get_pop <- function(dataset, chromosome, marker1, build_version=CURRENT_BUILD) {
  # find if dataset is a super- or sub-population
  super_pop<-as.character(POPS[dataset,])
  sub_pop<-NA
  if(is.na(super_pop)){
    super_pop<-dataset # calculation over super-population
  }else{
    sub_pop<-dataset   # calculation over sub-population
  }
  sub_pops <- row.names(subset(POPS, tmp....1. == dataset))
  sm <- try(eval(parse(text = super_pop))[[paste('chr',chromosome,sep="")]])

  m1_pos = sm$map[sm$map$ID == marker1,]$POS
  pos=sm$map[sm$map$ID==marker1,]$POS
  
  if(length(pos) == 0){
    err <- paste("Marker ",marker1," was not found in this dataset",sep="")
  } else {
    snp_data_subset <-sm$gt[,marker1]
    msg <- list()
    msg[[ super_pop ]] <- pops(marker1, snp_data_subset, sm, NA)
    for (sub_pop in sub_pops){ 
      msg[[ sub_pop ]] <- pops(marker1, snp_data_subset, sm, sub_pop)
    }
  }

  if(exists("err")) {
    toJSON(list(error=err))
  } else {
    toJSON(msg)
  }
}

pops<-function(marker1, snp_data_subset, sm, sub_pop){
  # filter by sub-population if requested
  if(!is.na(sub_pop)){
    sub_pop_samples <- row.names(subset(SAMPLES, pop==sub_pop))
    snp_data_subset<-snp_data_subset[sub_pop_samples,]
  }
  snp_stats <-col.summary(snp_data_subset)
  msg <- {}
  
  #Check that snp makes it through quality control
  if (marker1 %in% colnames(snp_data_subset)) {
    #print(snp_stats[marker1,])
    #print(sm$map[sm$map$ID == marker1,])

    # 
    # RAF: The "risk" allele (allele B) frequency. So if 
    # RAF < 0.5 allele B is the minor allele otherwise A.
    if(snp_stats[marker1,]$RAF < 0.5) {
      msg$minor <- sm$map[sm$map$ID == marker1,]$ALT
      msg$major <- sm$map[sm$map$ID == marker1,]$REF
    } else {
      msg$minor <- sm$map[sm$map$ID == marker1,]$REF
      msg$major <- sm$map[sm$map$ID == marker1,]$ALT
    }
    msg$MAF <- c(round(snp_stats[marker1,]$MAF, digits=7))
  } else {
    msg <- paste("Marker ",marker1," did not pass QC",sep="")
  }
  msg
}

tabix_bin="/usr/local/bin/tabix"

## convert VCF format to snpMatrix format
vcf2sm<-function(region,dataset){
  chrom <- gsub("^([^:]+).*","\\1",region)
  kg_file <- list.files(paste0(KG_DIR,'/',dataset), paste0("ALL\\.",'chr',chrom,'\\..*',dataset,'.*vcf.gz$'), full.names=TRUE)
  ## get the header
  my.pipe<-pipe(paste(tabix_bin,'-H',kg_file,region))
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
  tmp<-as.data.frame(fread(paste(tabix_bin,kg_file,region,' | grep PASS'),sep="\t",header=FALSE,stringsAsFactors=FALSE))
  colnames(tmp)<-cnames
  return(convert2sm(tmp))
}

convert2sm<-function(vcfrows){
  ## remove SNPs which have more than 2 alleles as we cannot process these
  idx<-which(nchar(vcfrows$ALT)>=8 | nchar(vcfrows$REF) >=8)
  if(length(idx)>0)
    vcfrows<-vcfrows[-idx,]
  gt<-vcfrows[,10:ncol(vcfrows)]
  if(nrow(gt)==0)
    return(NA)
  info<-vcfrows[,1:9]
  
  sm<-apply(gt,1,gtcode)
  sm<-t(apply(sm,1,as.raw))
  
  colnames(sm)<-make.names(info$ID,unique=TRUE)
  rownames(sm)<-colnames(gt)
  sm<-new("SnpMatrix", sm)
  return(list(map=info,gt=sm))
}

## genotype data coded as 0, 1, 2, or 3
gtcode<-function(x) {
  x=sub("0|0","1",x,fixed=TRUE)
  x=sub("0|1","2",x,fixed=TRUE)
  x=sub("1|0","2",x,fixed=TRUE)
  x=sub("1|1","3",x,fixed=TRUE)
  ## set anything else to a missing value
  sub("[0-9]\\|[0-9]","0",x)
}


# convert VCF files to snpMatrix and save
save_rdata<-function(dataset,build_version) {
  data.files <- list.files(path=paste0(KG_DIR,'/',build_version,'/',dataset,'/'), pattern='*.vcf.gz$', full.names = TRUE)

  for (df in data.files) {
    if(!file.exists(paste(df, 'rds', sep='.'))) {
      my.pipe<-pipe(paste(tabix_bin,'-H',df))
      header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
      close(my.pipe)
      cnames<-unlist(strsplit(header,"\t"))
      
      fn<-paste(tabix_bin,df)
      chr<-gsub('.*(chr[^\\.]+).*', '\\1', basename(df))
      print(paste('Reading',fn, '...'))
      tmp<-as.data.frame(fread(paste('gzcat',df, ' | grep -v "^#"',' | grep PASS'),sep="\t",header=FALSE,stringsAsFactors=FALSE))
      colnames(tmp)<-cnames
      print(paste('Create snpMatrix for',chr, '...'))
      sm = convert2sm(tmp)
      print(paste('Saving',chr, '...'))
      saveRDS(sm, file = paste(df, 'rds', sep='.'))
    }
  }
}
