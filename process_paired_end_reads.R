process_paired_end_reads <- function(data_folder, output_string = "output",
  trim_right = 0, fwd_primer, rev_primer, unite_ref = NULL) {

  require(dada2)
  require(ShortRead)
  require(Biostrings)


  ## load cutadapt
  cutadapt <- "/usr/local/bin/cutadapt"
  system2(cutadapt, args = "--version")

  ## create output folders
  output_main <- paste(output_string, Sys.Date(), sep = "_")
  if (!dir.exists(output_main)) dir.create(output_main)

  output_filtN <- paste(output_main, "filtN", sep = "/")
  if (!dir.exists(output_filtN)) dir.create(output_filtN)

  output_cut <- paste(output_main, "primers_cut", sep = "/")
  if (!dir.exists(output_cut)) dir.create(output_cut)

  output_QC <- paste(output_main, "QC", sep = "/")
  if (!dir.exists(output_QC)) dir.create(output_QC)

  output_derep <- paste(output_main, "dereplicated_reads", sep = "/")
  if (!dir.exists(output_derep)) dir.create(output_derep)

  output_err <- paste(output_main, "learned_errors", sep = "/")
  if (!dir.exists(output_err)) dir.create(output_err)

  output_dada <- paste(output_main, "ASVs", sep = "/")
  if (!dir.exists(output_dada)) dir.create(output_dada)

  output_merged <- paste(output_main, "merged_reads", sep = "/")
  if (!dir.exists(output_merged)) dir.create(output_merged)


  ## input files
  fwd_reads <- sort(list.files(data_folder, pattern = "_1.fastq.gz",
    full.names = T))
  names(fwd_reads) <- gsub('^.*19239_\\s*|\\s*_lib.*$', '', fwd_reads)

  rev_reads <- sort(list.files(data_folder, pattern = "_2.fastq.gz",
    full.names = T))
  names(rev_reads) <- gsub('^.*19239_\\s*|\\s*_lib.*$', '', rev_reads)


  ## functions to remove primers
  # re-orient primers
  fwd_primer_rc <- dada2:::rc(fwd_primer)
  rev_primer_rc <- dada2:::rc(rev_primer)
  fwd_flags <- paste("-g", fwd_primer, "-a", rev_primer_rc) 
  rev_flags <- paste("-g", rev_primer, "-a", fwd_primer_rc) 

  # path to sequence files with Ns filtered out
  fwd_reads_filtN <- file.path(output_filtN, basename(fwd_reads))
  names(fwd_reads_filtN) <- names(fwd_reads)
  rev_reads_filtN <- file.path(output_filtN, basename(rev_reads))
  names(rev_reads_filtN) <- names(rev_reads)

  # path to sequence files with primers cut
  fwd_reads_cut <- file.path(output_cut, basename(fwd_reads))
  names(fwd_reads_cut) <- names(fwd_reads)
  rev_reads_cut <- file.path(output_cut, basename(rev_reads))
  names(rev_reads_cut) <- names(rev_reads)

  # create all orientations for sequencing primers
  allOrients <- function(primer) {
      # Create all orientations of the input sequence
      require(Biostrings)
      dna <- DNAString(primer)
      orients <- c(Forward = dna, Complement = complement(dna),
        Reverse = reverse(dna), RevComp = reverseComplement(dna))
      return(sapply(orients, toString))
  }
  fwd_orients <- allOrients(fwd_primer)
  rev_orients <- allOrients(rev_primer)

  # map primers to sequences
  primerHits <- function(primer, fn) {
      # Counts number of reads in which the primer is found
      nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = F)
      return(sum(nhits > 0))
  }


  ## create log table with summary data
  log_tab <- data.frame(row.names = names(fwd_reads))
  log_tab$filename <- rep(NA, nrow(log_tab))
  log_tab$timestamp_start <- as.POSIXlt(rep(NA, nrow(log_tab)))
  log_tab$timestamp_end <- as.POSIXlt(rep(NA, nrow(log_tab)))
  log_tab$fwd_primer_match <- rep(NA, nrow(log_tab))
  log_tab$rev_primer_match <- rep(NA, nrow(log_tab))
  log_tab$fwd_primer_match_cut <- rep(NA, nrow(log_tab))
  log_tab$rev_primer_match_cut <- rep(NA, nrow(log_tab))
  log_tab$reads_filt_in <- rep(NA, nrow(log_tab))
  log_tab$reads_filt_out <- rep(NA, nrow(log_tab))
  log_tab$reads_fwd_derep <- rep(NA, nrow(log_tab))
  log_tab$reads_rev_derep <- rep(NA, nrow(log_tab))
  log_tab$fwd_asv <- rep(NA, nrow(log_tab))
  log_tab$rev_asv <- rep(NA, nrow(log_tab))
  log_tab$merged_reads <- rep(NA, nrow(log_tab))
  log_tab$merged_asv <- rep(NA, nrow(log_tab))


  ## Individual sample processing in loop
  for(i in names(fwd_reads)) {
    cat("Processing:", i, "\n")

    # record time for analysis end
    log_tab[i, "timestamp_start"] <- Sys.time()

    tryCatch({

      # remove primers
      filterAndTrim(fwd_reads[i], fwd_reads_filtN[i], rev_reads[i],
        rev_reads_filtN[i], maxN = 0, multithread = F)
      log_tab[i, "fwd_primer_match"] <- sum(sapply(fwd_orients, primerHits,
        fn = fwd_reads_filtN[[i]]), sapply(fwd_orients, primerHits,
        fn = rev_reads_filtN[[i]]))
      log_tab[i, "rev_primer_match"] <- sum(sapply(rev_orients, primerHits,
        fn = fwd_reads_filtN[[i]]), sapply(rev_orients, primerHits,
        fn = rev_reads_filtN[[i]]))
      system2(cutadapt, args = c(fwd_flags, "-n", 1, "-o",
        fwd_reads_cut[i], fwd_reads_filtN[i]))
      system2(cutadapt, args = c(rev_flags, "-n", 1, "-o",
        rev_reads_cut[i], rev_reads_filtN[i]))
      log_tab[i, "fwd_primer_match_cut"] <- sum(sapply(fwd_orients, primerHits,
        fn = fwd_reads_cut[[i]]), sapply(fwd_orients, primerHits,
        fn = rev_reads_cut[[i]]))
      log_tab[i, "rev_primer_match_cut"] <- sum(sapply(rev_orients, primerHits,
        fn = fwd_reads_cut[[i]]), sapply(rev_orients, primerHits,
        fn = rev_reads_cut[[i]]))

      # quality-filter sequences
      fwd_filt <- tempfile(fileext = ".fastq.gz")  
      rev_filt <- tempfile(fileext = ".fastq.gz")
      temp <- filterAndTrim(fwd = fwd_reads_cut[i], filt = fwd_filt,
        rev = rev_reads_cut[i], filt.rev = rev_filt, maxN = 0, maxEE = c(2, 2),
        truncQ = 2, minLen = 50, compress = T, verbose = T)
      log_tab[i , "filename"] <- rownames(temp)
      log_tab[i , "reads_filt_in"] <- temp[1]
      log_tab[i , "reads_filt_out"] <- temp[2]

      # plot QC plots
      pdf(paste0(output_QC, "/", i, "_fwd.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(fwd_reads[i])
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_rev.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(rev_reads[i])
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_fwd_cut.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(fwd_reads_cut[i])
      print(x)
      dev.off()

      pdf(paste0(output_QC, "/", i, "_rev_cut.pdf"), w = 4, h = 3)
      x <- plotQualityProfile(rev_reads_cut[i])
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

      # merge forward and reverse reads
      merged_reads <- mergePairs(fwd_dada, fwd_derep, rev_dada, rev_derep,
        verbose = T)
      log_tab[i , "merged_reads"] <- sum(merged_reads$abundance)
      log_tab[i , "merged_asv"] <- length(merged_reads$sequence)
      saveRDS(merged_reads, file = paste0(output_merged, "/", i, "_merged.rds"))
  
      # remove objects and free-up memory for next loop iteration
      rm(temp, x, fwd_filt, rev_filt, fwd_derep, rev_derep, fwd_err, rev_err,
        fwd_dada, rev_dada, merged_reads, merged_reads_nochim)
      gc()

    }, error = function(e) e)

    # record time for analysis end
    log_tab[i, "timestamp_end"] <- Sys.time()
  }


  ## create sequence table
  seqtab_process_time <- as.POSIXlt(rep(NA, 2))
  names(seqtab_process_time) <- c("timestamp_start", "timestamp_end")
  seqtab_process_time["timestamp_start"] <- Sys.time()
 
  merged_files <- list.files(output_merged, pattern = "_merged.rds",
    full.names = T)
  names(merged_files) <- gsub('^.*reads/\\s*|\\s*_merged.*$', '', merged_files)
  merged_reads <- vector("list", length(merged_files))
  names(merged_reads) <- names(merged_files)
  for(i in names(merged_reads)) {
    merged_reads[[i]] <- readRDS(merged_files[i])
  }
  seqtab <- makeSequenceTable(merged_reads)
  saveRDS(seqtab, file = paste0(output_main, "/", "seqtab.rds"))


  ## remove chimeras
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
    multithread = F, verbose = T)
  saveRDS(seqtab_nochim, file = paste0(output_main, "/", "seqtab_nochim.rds"))
  seqtab_process_time["timestamp_end"] <- Sys.time()
  write(seqtab_process_time, file = paste0(output_main, "/",
    "seqtab_process_time.txt"))


  ## modify seqtable to add failed samples and extract ASV sequences
  failed <- names(fwd_reads)[!(names(fwd_reads) %in% names(merged_reads))]
  failed_reads <- matrix(0, nrow = length(failed), ncol = ncol(seqtab),
    dimnames = list(failed, colnames(seqtab)))
  cdm <- rbind(seqtab, failed_reads)
  ASV_seqs <- colnames(cdm)
  names(ASV_seqs) <- paste0("ASV", sprintf("%05d", 1:length(ASV_seqs)))
  colnames(cdm) <- names(ASV_seqs)
  write(paste0(">", names(ASV_seqs), "\n", ASV_seqs),
    file = paste0(output_main, "/", "ASV_seqs.fasta"))
  write.table(cdm, file = paste0(output_main, "/", "cdm.csv"), col.names = NA)


  ## assign taxonomy (SKIPPED: too slow!!!)
  #taxonomy_process_time <- as.POSIXlt(rep(NA, 2))
  #names(taxonomy_process_time) <- c("timestamp_start", "timestamp_end")
  #taxonomy_process_time["timestamp_start"] <- Sys.time()
  #taxonomy <- assignTaxonomy(ASV_seqs, unite_ref, multithread = F,
  #  tryRC = T, verbose = T, minBoot = 80)
  #saveRDS(taxonomy, file = paste0(output_main, "/", "taxonomy.rds"))
  #taxonomy_process_time["timestamp_end"] <- Sys.time()
  #write(taxonomy_process_time, file = paste0(output_main, "/",
  #  "taxonomy_process_time.txt"))


  ## save summary table
  write.table(log_tab, file = paste0(output_main, "/", "log.csv"), col.names = NA,
    sep = "\t")

}