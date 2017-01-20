#source function definitions or have just one document
require(doParallel)
require(foreach)

## load bam file using the which argument to ScanBamParam
bamfile3 <- "STAR/11647X5_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfile2 <- "STAR/11647X3_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfile1 <- "STAR/11647X1_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem1 <- "STAR/11647X2_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem2 <- "STAR/11647X4_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"
bamfilem3 <- "STAR/11647X6_150723_D00550_0277_BC7TFTANXX_1/Aligned.sortedByCoord.out.bam"

wt_list <- BamFileList(c(bamfile1))#, bamfile2, bamfile3))
mut_list <- BamFileList(c(bamfilem1))#, bamfilem2, bamfilem3))
bf_list <- BamFileList(c(wt_list, mut_list))

source('~/MMAPPR2/read_bam.r')

myrange <- as(seqinfo(bf_list), "GRanges")
# myrange <- GRanges(seqnames = 'chr4:1-77000000')
#cut to chromosomes and make short for testing
myrange <- GenomeInfoDb::keepStandardChromosomes(myrange)
# shorten range for faster testing (shortens each chromosome to just 20000 bp)
# width(ranges(myrange)) <- 100000



loess_fit <- function(chrRange){
  dist_df = getDistanceDf(chrRange)
  if(class(dist_df) == 'try-error') stop('empty dataframe')
  
  aicc_start <- proc.time(dist_df)
  aicc_df <- aicc_opt() #loess, store values; return best aicc
  message("aicc for chr ", names(chrRange), " took: ", proc.time() - aicc_start)
  
  #loess fit(df, aiccParam) #keep all and pick best of list? or takes too much memory
  #  recursive: (min, max, interval parameters)
  # each chromosome has list: $fit(pos, unfitted, fitted), $aicc(vector of aicc data)
  return(list(fit = dist_df, aicc = aicc_df)) #actually return list of fit and aicc
}

chr_list <- list()
for (i in 1:length(myrange)){
  chr_list[[toString(seqnames(myrange[i]))]] <- myrange[i]
}


start <- proc.time()


#automatic core number generation (only on Linux)
#mem_req is memory required per base per replicate (per two files) in MB
#I calculated it myself in a test using only the memory taken up by 
#the apply_pileup and dividing it by reps and sequence length
#in other words, it doesn't take into account all the memory used by 
#each core as the loess_fit (hence the getChrDf) function is run
core_calc <- function(mem_req = 33.3){
  tryCatch({
    #this gets free RAM (on Linux)
    free_ram <- system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE) #gets mem info in MB
    #get numeric, in bytes (not MB)
    free_ram <- 1000 * as.numeric(free_ram)
    #get longest chromosome, calculate needed memory, and divide into free RAM
    num_cores <- ceiling(free_ram / (max(seqlengths(myrange)) * mem_req * length(wt_list)))
    message("num_cores = ", num_cores)
    return(num_cores)
    },
  warning = function(w) {
    message('Warning: arbitrary number of clusters used (are you not on Linux?)')
    return(ceiling(detectCores() / 3))
  }
  )
}

#cluster generation
cl <- makeCluster(core_calc(), type = "SOCK", outfile = "")


# register the cluster
registerDoParallel(cl)

# insert parallel computation here
try({
  loess_fit_list <- foreach(i = chr_list,
                            .packages = c("Rsamtools", "tidyr", "dplyr")) %dopar% loess_fit(i)
  names(loess_fit_list) <- names(chr_list)
}
)
# loess_fit_list <- lapply(chr_list, loess_fit)

time <- proc.time() - start
print(time)

stopCluster(cl)
# insert serial backend, otherwise error in repetetive tasks
registerDoSEQ()

# clean up a bit.
invisible(gc)





