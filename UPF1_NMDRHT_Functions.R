# General Functions ---------------------------------------------------------------

# Define function for plotting "n" in some plots
n_fun <- function(x){
  return(data.frame(y = my_n,
                    label = paste0("n = ",
                                   length(x)
                    )))
}


# factR -------------------------------------------------------------------

.extractCDSchecks <- function(cds, fasta, argnames, ...) {
  # define global variables
  exonorder <- NULL
  
  if (!has_consistentSeqlevels(cds, fasta, verbose = FALSE)) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchChromosomes(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]
    ))
  }
  # catch wrong cds class
  if (is_gtf(cds)) {
    cds <- S4Vectors::split(cds[cds$type == "CDS"], ~transcript_id)
    if (length(cds) == 0) {
      rlang::abort(sprintf(
        "`%s` do not contain CDS information", argnames[1]
      ))
    }
  }
  
  if (!is(cds, "GRangesList")) {
    rlang::abort("cds class type is not GRanges GTF or GRangesList")
  }
  
  cds <- filtereach(cds, ...)
  if (length(cds) == 0) {
    rlang::abort("No CDS to display")
  }
  return(sorteach(cds, exonorder))
}


.getSequence <- function(cds, fasta) {
  x <- y <- instop <- NULL
  
  rlang::inform("Checking CDSs and translating protein sequences")
  cdsSeq <- GenomicFeatures::extractTranscriptSeqs(fasta, cds)
  aaSeq <- suppressWarnings(
    Biostrings::translate(cdsSeq, if.fuzzy.codon = "solve")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y = strsplit(x, split = "")) %>%
    dplyr::mutate(noATG = ifelse(y[[1]] != "M", TRUE, FALSE)) %>%
    dplyr::mutate(instop = ifelse("*" %in% y, TRUE, FALSE)) %>%
    dplyr::ungroup()
  
  # check for ATG and internal stop_codon, truncate proteins with internal 
  # stop codon
  ## and remove entries without proteins after truncation
  if (TRUE %in% aaSeq$noATG) {
    rlang::warn(sprintf("%s CDSs do not begin with ATG", sum(aaSeq$noATG)))
  }
  if (TRUE %in% aaSeq$instop) {
    aaSeq <- suppressWarnings(aaSeq %>%
                                dplyr::rowwise() %>%
                                dplyr::mutate(x = ifelse(instop == TRUE,
                                                         paste(y[seq_len(which(y == "*") - 1)], collapse = ""),
                                                         x
                                )) %>%
                                dplyr::mutate(y = strsplit(x, split = "")) %>%
                                dplyr::ungroup())
    
    rlang::warn(sprintf(paste0("%s CDSs contain internal stop codon. ",
                               "Truncating CDS sequence to retain ORF"), 
                        sum(aaSeq$instop)))
    if ("" %in% aaSeq$x) {
      rlang::warn(sprintf(paste0(
        "After truncation, %s cds have no ",
        "coding sequences. These CDSs were not analyzed"), 
        sum(aaSeq$x == "")))
      aaSeq <- aaSeq[aaSeq$x != "", ]
    }
  }
  rlang::inform(sprintf(
    "Predicting domain families for %s proteins", nrow(aaSeq)))
  return(aaSeq)
}

is_gtf <- function(...) {
  type <- gene_id <- transcript_id <- NULL
  return(unlist(lapply(list(...), function(x) {
    if (is(x, "GRanges")) {
      x <- as.data.frame(x)
      if (all(c("type", "gene_id", "transcript_id") %in% names(x))) {
        x <- x %>%
          dplyr::select(type, gene_id, transcript_id) %>%
          dplyr::distinct()
        if (nrow(x) > 1) {
          return(TRUE)
        }
      }
    }
    return(FALSE)
  })))
}


# Proteomics --------------------------------------------------------------

# Comment: Many of these functions were adapted from Oliver Popp (MDC)

create_peptidome_VB <- function(proteome_AAString, missed_cleavages, aa_range){
  # adapted from PeptideRanger 
  # PMID: 36701129 | https://github.com/rr-2/PeptideRanger/
  # proteome_namedVector: vector of protein sequence with transcript_id as names
  # missed_cleavages: range of missed cleavages considered in in-silico digestion
  # aa_range: range of number of amino acids in peptides that are filtered for
  
  size <- symbol <- uniprot <- start <- end <- NULL
  
  if(missing(missed_cleavages)){
    
    # default missed cleavages
    missed_cleavages <- c(0,2)
    
    print('default: 0-2 missed cleavages included in digestion')
  }
  
  
  if(missing(aa_range)){
    
    # default peptide size range
    aa_range <- c(7,52)
    
    print('default: 7-52 amino acid range')
  }
  
  # directly accesses fasta file using directory and loads in proteome as large AAStringSet
  proteome <- proteome_AAString
  
  # creates empty peptidome df
  peptidome <- data.frame()
  
  # loops through missed cleavage range
  for(i in missed_cleavages[1]:missed_cleavages[2]){
    
    # performs in silico trypsin digestion of proteome yielding peptides with i missed cleavages
    curr_cleavage <- BiocGenerics::as.data.frame(cleaver::cleave(proteome, enzym="trypsin", missedCleavages = i, unique = FALSE)) 
    
    # provides locations of peptides within parent proteins
    curr_positions <- BiocGenerics::as.data.frame(cleaver::cleavageRanges(proteome, enzym="trypsin", missedCleavages = i))
    
    # adds peptide start and end locations
    curr_cleavage$start <- curr_positions$start
    curr_cleavage$end <- curr_positions$end
    # adds the number of amino acids in peptides
    curr_cleavage$size <- curr_positions$width
    
    # adds the number of missed cleavages that resulted in the peptide
    curr_cleavage$missed_cleavages <- i
    # adds peptides from i missed cleavages and info to peptidome
    peptidome <- rbind(peptidome, curr_cleavage)
    
  }
  
  
  colnames(peptidome)[1] <- 'tx_row_id'
  colnames(peptidome)[2] <- 'transcript_id'
  colnames(peptidome)[3] <- 'sequence'
  
  
  # selects a subset of columns
  peptidome <- dplyr::select(peptidome, tx_row_id, transcript_id, missed_cleavages, sequence, start, end, size) %>% 
    # Flag stop codons (marked by *), classify as terminal or internal and remove those stop codons
    mutate(terminal_peptide = case_when(str_detect(sequence, "\\*$") ~ TRUE,
                                        TRUE ~ FALSE)) %>% 
    mutate(internal_stop = case_when(str_detect(sequence, "\\*") & terminal_peptide == FALSE ~ TRUE,
                                     TRUE ~ FALSE)) %>% 
    mutate(sequence = stringr::str_replace(sequence, '\\*', '')) %>% 
    # adjust size and end depending on stop codon presence
    mutate(end = case_when(terminal_peptide == TRUE ~ end-1,
                           terminal_peptide == FALSE ~ end),
           size = case_when(terminal_peptide == TRUE ~ size-1,
                            terminal_peptide == FALSE ~ size))
  
  # filters for peptides in the amino acid range specified
  peptidome = dplyr::filter(peptidome, size >= aa_range[1], size <= aa_range[2])
  
  # removes peptides with sequences containing "U"
  peptidome <- dplyr::filter(peptidome, grepl("U", peptidome$sequence) != TRUE)
  
  return(peptidome)
  
}

create_peptidome_VB_simple <- function(proteome_AAString, missed_cleavages, aa_range){
  # adapted from PeptideRanger 
  # PMID: 36701129 | https://github.com/rr-2/PeptideRanger/
  # proteome_namedVector: vector of protein sequence with transcript_id as names
  # missed_cleavages: range of missed cleavages considered in in-silico digestion
  # aa_range: range of number of amino acids in peptides that are filtered for
  
  
  size <- symbol <- uniprot <- start <- end <- NULL
  
  if(missing(missed_cleavages)){
    
    # default missed cleavages
    missed_cleavages <- c(0,2)
    
    print('default: 0-2 missed cleavages included in digestion')
  }
  
  
  if(missing(aa_range)){
    
    # default peptide size range
    aa_range <- c(7,52)
    
    print('default: 7-52 amino acid range')
  }
  
  
  # directly accesses fasta file using directory and loads in proteome as large AAStringSet
  proteome <- proteome_AAString
  
  # creates empty peptidome df
  peptidome <- data.frame()
  
  cut_sites <- "^M|K|R"
  
  # loops through missed cleavage range
  for(i in missed_cleavages[1]:missed_cleavages[2]){
    
    # performs in silico trypsin digestion of proteome yielding peptides with i missed cleavages
    curr_cleavage <- BiocGenerics::as.data.frame(cleaver::cleave(proteome, custom = cut_sites, missedCleavages = i, unique = FALSE)) 
    
    # provides locations of peptides within parent proteins
    curr_positions <- BiocGenerics::as.data.frame(cleaver::cleavageRanges(proteome, custom = cut_sites, missedCleavages = i))
    
    # adds peptide start and end locations
    curr_cleavage$start <- curr_positions$start
    curr_cleavage$end <- curr_positions$end
    # adds the number of amino acids in peptides
    curr_cleavage$size <- curr_positions$width
    
    # adds the number of missed cleavages that resulted in the peptide
    curr_cleavage$missed_cleavages <- i
    # adds peptides from i missed cleavages and info to peptidome
    peptidome <- rbind(peptidome, curr_cleavage)
    
  }
  
  
  colnames(peptidome)[1] <- 'tx_row_id'
  colnames(peptidome)[2] <- 'transcript_id'
  colnames(peptidome)[3] <- 'sequence'
  
  
  # selects a subset of columns
  peptidome <- dplyr::select(peptidome, tx_row_id, transcript_id, missed_cleavages, sequence, start, end, size) %>% 
    # Flag stop codons (marked by *), classify as terminal or internal and remove those stop codons
    mutate(terminal_peptide = case_when(str_detect(sequence, "\\*$") ~ TRUE,
                                        TRUE ~ FALSE)) %>% 
    mutate(internal_stop = case_when(str_detect(sequence, "\\*") & terminal_peptide == FALSE ~ TRUE,
                                     TRUE ~ FALSE)) %>% 
    mutate(sequence = stringr::str_replace(sequence, '\\*', '')) %>% 
    # adjust size and end depending on stop codon presence
    mutate(end = case_when(terminal_peptide == TRUE ~ end-1,
                           terminal_peptide == FALSE ~ end),
           size = case_when(terminal_peptide == TRUE ~ size-1,
                            terminal_peptide == FALSE ~ size))
  
  # filters for peptides in the amino acid range specified
  peptidome = dplyr::filter(peptidome, size >= aa_range[1], size <= aa_range[2])
  
  # removes peptides with sequences containing "U"
  peptidome <- dplyr::filter(peptidome, grepl("U", peptidome$sequence) != TRUE)
  
  return(peptidome)
  
}

mylog <- function(x, base=2) {
  tmp <- log(x, base = base)
  tmp[is.infinite(tmp)] <- 0/0 # NA
  return(tmp)
} 

valid.n <- function(x) sum(!is.na(x))

myboxplot <- function( matrix, adj=0.965, srt=35, cex=0.9, main=NA, pdfout=F ) {
  
  if (pdfout) {
    nam <- ifelse(!is.na(main), paste0(main, ".pdf"), timepaste(".pdf"))
    try(dev.off(), silent = T)
    mywidth <- ncol(matrix)/exp(1)
    if(mywidth < 8) {
      mywidth <- 8
    }
    pdf(nam, width = mywidth)
  }
  ## Draw boxplot with no axes.
  boxplot(matrix, xaxt = "n", yaxt="n", frame.plot=FALSE, main=main, border="darkblue", col="#F2A93B", col.main="darkblue")
  
  ## Draw x-axis without labels.
  # axis(side = 1, labels = T)
  
  ## Draw y-axis.
  axis(side = 2, cex.axis = 1.2,
       ## Rotate labels perpendicular to y-axis.
       las = 2, family = "serif",
       ## Adjust y-axis label positions.
       mgp = c(3, 0.75, 0), col.axis="darkblue", col="darkblue")
  
  ## Draw the x-axis labels.
  text(x = 1:ncol(matrix),
       ## Move labels to just below bottom of chart.
       y = par("usr")[3] - 0.45,
       ## Use names from the data list.
       labels = colnames(matrix),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = srt,
       ## Adjust the labels to almost 100% right-justified.
       adj = adj,
       ## Increase label size.
       cex = cex,
       col = "darkblue",
       family = "serif")
  if (pdfout) {
    dev.off()
  }
}


mydensityplot <- function( matrix, xlab="log[2] intensity", width=NA, title="Density plot", show.legend=T ) {
  print( as.data.frame(matrix) %>% rownames_to_column("ID") %>% gather(var, val, -ID) %>% ggplot(aes(val, col=var)) +
           geom_density(show.legend = show.legend) +
           labs(y="", x=xlab, title = title))
  if ( title=="" ) {
    ggsave(timepaste(paste0(xlab, "_Density.pdf")), width = width)
  } else {
    ggsave(paste(title, paste0(xlab, "_Density.pdf"), sep = "_"), width = width)
  }
}


modF <- function(d, columnsForTest, class.vector, id.col="ID", FDR=0.05 ) {
  # ed <- rt("expDesign.txt")
  # columnsForTest <- names(df) %in% ed$Column.Name[ed$Experiment!=""]
  # class.vector   <- ed$Experiment[ed$Experiment != ""]
  id <- d[,id.col]
  data <-  d [,columnsForTest]
  if ( length(class.vector)!=ncol(data) ) { stop("Class vector and columns to be analyzed are not the same!")}
  
  # cat('\n-- modF.test --\n')
  
  f <- factor (class.vector)
  design <- model.matrix ( ~ 0 + f )
  
  # moderated F test
  fit <- lmFit (data, design)
  fit <- eBayes (fit)
  
  sig <- topTable (fit, number=nrow(data), sort.by='none')
  mod.sig <- sig [,'adj.P.Val'] < FDR
  non.na.n <- apply (data, 1, function (x) { sum (is.finite (x)) })
  
  
  mod.f <- data.frame ( cbind (id=id, sig, significant=mod.sig, total.n=non.na.n, stringsAsFactors=F) )
  # mod.f <- data.frame ( cbind (sig, significant=mod.sig, total.n=non.na.n, stringsAsFactors=F) )
  
  colnames(mod.f) <- sub("^f", "", colnames(mod.f))
  colnames(mod.f) <- paste0("modF_", colnames(mod.f))
  # volcano plot possible?
  
  # cat('\n-- modF.test exit --\n')
  # return( list( input=d, output=final.results, groups=class.vector)  )
  return(mod.f)
}


oneWayModAnova <- function (mat, group, repeated = FALSE, subject,
                            adjust.method = "BH", 
                            sort.by = "none", trend=TRUE ) {
  # function updated on 2024-09-03: uses now group as factors
  ## might render rename_anova_comparisons function obsolete
  if (missing(mat)) 
    stop("'mat' is missing")
  if (!is.matrix(mat)) 
    stop("'mat' must be a matrix")
  if (missing(group)) 
    stop("'group' is missing")
  if (!is.factor(group)) 
    group <- factor(group)
  if (ncol(mat) != length(group)) 
    stop("length of group must be equal to ncol(mat)")
  
  nlev <- nlevels(group)
  if (nlev < 3) 
    stop("'group' has less than three levels, use 'mod.t.test' instead")
  
  # Retain original factor levels instead of using LETTERS
  group.tmp <- group
  
  if (repeated) {
    if (missing(subject)) 
      stop("'subject' is missing")
    if (!is.factor(subject)) 
      subject <- factor(subject)
    
    design <- model.matrix(~group.tmp + subject)
    fit1 <- lmFit(mat, design)
    fit2 <- eBayes(fit1)
    res <- topTable(fit2, coef = 2:(1 + nlevels(group.tmp)), 
                    adjust.method = adjust.method, number = Inf, sort.by = sort.by)
    res <- res[, c("AveExpr", "F", "P.Value", "adj.P.Val")]
    names(res) <- c("grand mean", "F", "P.value", "adj.P.value")
  } else {
    design <- model.matrix(~0 + group.tmp)
    colnames(design) <- levels(group.tmp)
    fit1 <- lmFit(mat, design)
    
    combs <- apply(combn(levels(group.tmp), 2), 2, paste0, collapse = "-")
    cont.matrix <- makeContrasts(contrasts = combs, levels = design)
    fit2 <- contrasts.fit(fit1, cont.matrix)
    fit3 <- eBayes(fit2, trend = trend, robust = TRUE)
    
    res <- topTable(fit3, adjust.method = adjust.method, 
                    number = Inf, sort.by = sort.by)
    
    levs <- levels(group)
    combs <- apply(combn(levs, 2), 2, paste0, collapse = " vs ")
    names(res) <- c(combs, "grand mean", "F", "P.value", "adj.P.value")
  }
  
  return(res)
}



# Simplified Adaption of replaceMissingFromGaussian for Gaussian imputation by column -- original function from Matthias Ziehm
Impute <- function( data, width=0.3, shift=1.8, setSeed=T, plotHist=T, export_summary=T ) {
  
  if (export_summary) {
    # Generate summary table
    summary_table_rows <- data.frame(
      Rows          = rownames(data),
      Imputed_values = rowSums(is.na(data))
    )
    summary_table_columns <- data.frame(
      Columns       = colnames(data),
      Imputed_values = colSums(is.na(data))
    )
  }
  
  # automatically generates histograms in the working directory
  if ( plotHist ) { pdf(paste0("histograms_Imputation_", length(list.files(pattern = "histograms_Imputation")), ".pdf")) }
  
  data=apply(data, 2, function(x) {
    x[is.infinite(x)]=0/0
    return(x)
  }) # substitute +-Inf with NA, which will then be substituted by values from normal distribution
  
  
  lengths <- sapply(colnames(data), function(x) {
    x <- data[, which(colnames(data) == x)]
    x[is.na(x)] <- rnorm(length(x[is.na(x)]), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
    return(length(x))
  })
  print(lengths)
  
  
  data=sapply(colnames(data), function(x) {
    nam <- x
    x <- data[,which(colnames(data)==x)]
    tempSD=sd(x, na.rm=T)
    if ( plotHist ) { try(hist(x, breaks = floor(length(x)/100), main=nam, xlab = "", col=adjustcolor("darkblue", alpha.f = 0.3)), silent = T) }
    where <- which(is.na(x))
    if( setSeed ) { set.seed(28091980) }
    x[is.na(x)]=rnorm(length(x[is.na(x)]), mean(x, na.rm=T)-shift*tempSD, width*tempSD)
    if ( plotHist ) { try(hist(x[where], floor(length(where)/10), main=nam, xlab = "", col=adjustcolor("darkred", alpha.f = 0.5), add=T),
                          silent = T) }
    return(x)
  })
  
  if ( plotHist ) { dev.off() }
  
  return(data)
} # end Impute
