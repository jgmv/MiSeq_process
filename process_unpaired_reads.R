process_unpaired_reads <- function(data_folder, output_string = "output",
  sample_names) {

  require(dada2)


  ## input files
  fwd_reads <- sort(list.files(data_folder, pattern = "_R1.fastq.bz2",
    full.names = T))
  #names(fwd_reads) <- gsub('^.*cox2R-\\s*|\\s*_R1.*$', '', fwd_reads)
  names(fwd_reads) <- sample_names

  rev_reads <- sort(list.files(data_folder, pattern = "_R2.fastq.bz2",
    full.names = T))
  #names(rev_reads) <- gsub('^.*cox2R-\\s*|\\s*_R2.*$', '', rev_reads)
  names(rev_reads) <- sample_names


  ## create output folders
  output_main <- paste(output_string, Sys.Date(), sep = "_")
  if (!dir.exists(output_main)) dir.create(output_main)

  output_QC <- paste(output_main, "QC", sep = "/")
  if (!dir.exists(output_QC)) dir.create(output_QC)

  output_derep <- paste(output_main, "dereplicated_reads", sep = "/")
  if (!dir.exists(output_derep)) dir.create(output_derep)

  output_err <- paste(output_main, "learned_errors", sep = "/")
  if (!dir.exists(output_err)) dir.create(output_err)

  output_dada <- paste(output_main, "ASVs", sep = "/")
  if (!dir.exists(output_dada)) dir.create(output_dada)


  ## create log table with summary data
  log_tab <- data.frame(row.names = names(fwd_reads))
  log_tab$filename_fwd <- rep(NA, nrow(log_tab))
  log_tab$filename_rev <- rep(NA, nrow(log_tab))
  log_tab$timestamp_start <- as.POSIXlt(rep(NA, nrow(log_tab)))
  log_tab$timestamp_end <- as.POSIXlt(rep(NA, nrow(log_tab)))
  log_tab$reads_fwd_filt_in <- rep(NA, nrow(log_tab))
  log_tab$reads_fwd_filt_out <- rep(NA, nrow(log_tab))
  log_tab$reads_rev_filt_in <- rep(NA, nrow(log_tab))
  log_tab$reads_rev_filt_out <- rep(NA, nrow(log_tab))
  log_tab$reads_fwd_derep <- rep(NA, nrow(log_tab))
  log_tab$reads_rev_derep <- rep(NA, nrow(log_tab))
  log_tab$fwd_asv <- rep(NA, nrow(log_tab))
  log_tab$rev_asv <- rep(NA, nrow(log_tab))


  ## individual sample processing in loop
  for(i in names(fwd_reads)) {
    cat("Processing:", i, "\n")

    # record time for analysis end
    log_tab[i, "timestamp_start"] <- Sys.time()

    tryCatch({

      # quality-filter sequences
      fwd_filt <- tempfile(fileext = ".fastq.gz")  
      rev_filt <- tempfile(fileext = ".fastq.gz")
      temp <- fastqFilter(fwd_reads[i], fwd_filt, truncQ = 2, truncLen = 0,
        trimLeft = 0, maxN = 0, minQ = 0, maxEE = Inf, rm.phix = F, n = 1e+06,
        compress = T, verbose = T)
      log_tab[i , "filename_fwd"] <- basename(fwd_reads[i])
      log_tab[i , "reads_fwd_filt_in"] <- temp[1]
      log_tab[i , "reads_fwd_filt_out"] <- temp[2]

      temp <- fastqFilter(rev_reads[i], rev_filt, truncQ = 2, truncLen = 0,
        trimLeft = 0, maxN = 0, minQ = 0, maxEE = Inf, rm.phix = F, n = 1e+06,
        compress = T, verbose = T)
      log_tab[i , "filename_rev"] <- basename(rev_reads[i])
      log_tab[i , "reads_rev_filt_in"] <- temp[1]
      log_tab[i , "reads_rev_filt_out"] <- temp[2]

      # plot QC plots
      pdf(paste0(output_QC, "/", i, "_fwd.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(fwd_reads[i])
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_rev.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(rev_reads[i])
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_fwd_filt.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(fwd_filt)
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_rev_filt.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(rev_filt)
      print(x)
      dev.off()

      # learn error rates
      fwd_err <- learnErrors(fwd_filt, multithread = F)
      saveRDS(fwd_err, file = paste0(output_err, "/", i, "_fwd_err.rds"))
      rev_err <- learnErrors(rev_filt, multithread = F)
      saveRDS(rev_err, file = paste0(output_err, "/", i, "_rev_err.rds"))

      pdf(paste0(output_err, "/", i, "_fwd_err.pdf"), w = 4, h = 4)
      x <- plotErrors(fwd_err, nominalQ = T)
      print(x)
      dev.off()

      pdf(paste0(output_err, "/", i, "_rev_err.pdf"), w = 4, h = 4)
      x <- plotErrors(rev_err, nominalQ = T)
      print(x)
      dev.off()
  
      # dereplicate sequences
      fwd_derep <- derepFastq(fwd_filt, verbose = T)
      log_tab[i , "reads_fwd_derep"] <- length(fwd_derep$uniques)
      saveRDS(fwd_derep, file = paste0(output_derep, "/", i, "_fwd_derep.rds"))
      rev_derep <- derepFastq(rev_filt, verbose = T)
      log_tab[i , "reads_rev_derep"] <- length(rev_derep$uniques)
      saveRDS(rev_derep, file = paste0(output_derep, "/", i, "_rev_derep.rds"))

      # infer sample composition
      fwd_dada <- dada(fwd_derep, err = fwd_err, multithread = F)
      log_tab[i , "fwd_asv"] <- length(fwd_dada$sequence)
      saveRDS(fwd_dada, file = paste0(output_dada, "/", i, "_fwd_dada.rds"))
      rev_dada <- dada(rev_derep, err = rev_err, multithread = F)
      log_tab[i , "rev_asv"] <- length(rev_dada$sequence)
      saveRDS(rev_dada, file = paste0(output_dada, "/", i, "_rev_dada.rds"))

      # remove objects and free-up memory for next lopp iteration
      rm(temp, x, fwd_filt, rev_filt, fwd_derep, rev_derep, fwd_err, rev_err,
        fwd_dada, rev_dada, merged_reads, merged_reads_nochim)
      gc()

    }, error = function(e) e)

    # record time for analysis end
    log_tab[i, "timestamp_end"] <- Sys.time()
  }


  ## create sequence table
  fwd_files <- list.files(output_dada, pattern = "_fwd_dada.rds",
    full.names = T)
  names(fwd_files) <- gsub('^.*ASVs/\\s*|\\s*_fwd_dada.*$', '', fwd_files)
  fwd_dada <- vector("list", length(fwd_files))
  names(fwd_dada) <- names(fwd_files)
  for(i in names(fwd_dada)) {
    fwd_dada[[i]] <- readRDS(fwd_files[i])
  }
  seqtab_fwd <- makeSequenceTable(fwd_dada)
  saveRDS(seqtab_fwd, file = paste0(output_main, "/", "seqtab_fwd.rds"))

  rev_files <- list.files(output_dada, pattern = "_rev_dada.rds",
    full.names = T)
  names(rev_files) <- gsub('^.*ASVs/\\s*|\\s*_rev_dada.*$', '', rev_files)
  rev_dada <- vector("list", length(rev_files))
  names(rev_dada) <- names(rev_files)
  for(i in names(rev_dada)) {
    rev_dada[[i]] <- readRDS(rev_files[i])
  }
  seqtab_rev <- makeSequenceTable(rev_dada)
  saveRDS(seqtab_rev, file = paste0(output_main, "/", "seqtab_rev.rds"))


  ## remove chimeras
  seqtab_fwd_nochim <- removeBimeraDenovo(seqtab_fwd, method = "consensus",
    multithread = F, verbose = T)
  saveRDS(seqtab_fwd_nochim, file = paste0(output_main, "/",
    "seqtab_fwd_nochim.rds"))

  seqtab_rev_nochim <- removeBimeraDenovo(seqtab_rev, method = "consensus",
    multithread = F, verbose = T)
  saveRDS(seqtab_rev_nochim, file = paste0(output_main, "/",
    "seqtab_rev_nochim.rds"))


  ## modify seqtable to add failed samples and extract ASV sequences
  failed <- names(fwd_reads)[!(names(fwd_reads) %in% names(fwd_dada))]
  failed_reads <- matrix(0, nrow = length(failed), ncol = ncol(seqtab_fwd),
    dimnames = list(failed, colnames(seqtab_fwd)))
  cdm_fwd <- rbind(seqtab_fwd, failed_reads)
  ASV_fwd_seqs <- colnames(cdm_fwd)
  names(ASV_fwd_seqs) <- paste0("ASV", sprintf("%05d", 1:length(ASV_fwd_seqs)))
  colnames(cdm_fwd) <- names(ASV_fwd_seqs)
  write(paste0(">", names(ASV_fwd_seqs), "\n", ASV_fwd_seqs),
    file = paste0(output_main, "/", "ASV_fwd_seqs.fasta"))
  write.table(cdm_fwd, file = paste0(output_main, "/", "cdm_fwd.csv"),
    col.names = NA)

  failed <- names(rev_reads)[!(names(rev_reads) %in% names(rev_dada))]
  failed_reads <- matrix(0, nrow = length(failed), ncol = ncol(seqtab_rev),
    dimnames = list(failed, colnames(seqtab_rev)))
  cdm_rev <- rbind(seqtab_rev, failed_reads)
  ASV_rev_seqs <- colnames(cdm_rev)
  names(ASV_rev_seqs) <- paste0("ASV", sprintf("%05d", 1:length(ASV_rev_seqs)))
  colnames(cdm_rev) <- names(ASV_rev_seqs)
  write(paste0(">", names(ASV_rev_seqs), "\n", ASV_rev_seqs),
    file = paste0(output_main, "/", "ASV_rev_seqs.fasta"))
  write.table(cdm_rev, file = paste0(output_main, "/", "cdm_rev.csv"),
    col.names = NA)


  ## save summary table
  write.table(log_tab, file = paste0(output_main, "/", "log.csv"),
    col.names = NA, sep = "\t")

}
