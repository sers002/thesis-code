#' Read in a genotype matrix stored in HAPGEN2 format
#'
#' @param filename Filename to read in, e.g. "test.haps.gz". This function also
#' expects a file named the same, but with an "info.gz" extension with metadata
#' (first five columns of a .gen file)
#'
#' @return The genotype matrix as a data.table
#' @export
read_genotype_hapgen <- function(filename) {
  hapgen_haps <- fread(filename)

  hapgen_meta <- fread(
    gsub("\\.haps\\.", ".info.", filename),
    col.names = c("i", "snp", "pos", "a1", "a2")
  )

  ## Add consecutive pairs of columns (i.e. col1+col2, col3+col4, col5+col6, ...)
  hapgen_dip <- hapgen_haps[, .SD, .SDcols=seq(1, ncol(hapgen_haps) - 1, 2)] +
    hapgen_haps[, .SD, .SDcols=seq(2, ncol(hapgen_haps), 2)]

  hapgen_dip <- 2 - hapgen_dip

  hapgen_dip <- cbind(hapgen_meta, hapgen_dip)

  hapgen_dip
}

#' Read in a genotype matrix stored in VCF format
#'
#' @param filename Filename to read in, e.g. "test.vcf.gz".
#' @param skip Number of rows of metadata to skip at beginning of file
#'
#' @return The genotype matrix as a data.table
#' @export
read_genotype_vcf <- function(filename) {
  test_genotype = fread(filename, skip = 100)

  ## VCFs have "0|0" for the genotype for each participant, so turn these into 0,1,2
  test_genotype[, 10:ncol(test_genotype)] <- lapply(
    test_genotype[, 10:ncol(test_genotype)],
    function(col) {
      as.numeric(substr(col, 1, 1)) + as.numeric(substr(col, 3, 3))
    }
  )

  setnames(
    test_genotype,
    c("ID",  "ALT", "REF"),
    c("snp", "a1", "a2")
  )

  test_genotype
}

#' Simulate gene expression data
#'
#' @param genotype_mat_input A genotype matrix in data.table format with rows
#' being variants and columns being participants. Metadata is stored in the
#' first columns (denoted by the first \code{genotype_start_col} columns).
#' @param betas A data.table with beta estimates.
#' @param genotype_start_col How many columns in the genotype matrix are
#' non-genotype data (i.e. rsid, ref/alt allele, etc)
#' @param seed A seed set just before errors are generated for expression data
#'
#' @return A sparse matrix containing gene expression for each gene and
#' participant.
#' @export
generate_expression <- function(
  genotype_mat_input,
  betas,
  genotype_start_col = 10,
  seed = 12345
) {
  message("  Starting expression function")
  ## Subset to only SNPs common to eQTL dataset and genotype
  in_both <- intersect(betas$snp, genotype_mat_input$snp)
  genotype_mat_input <- genotype_mat_input[snp %in% in_both]
  betas <- betas[snp %in% in_both]

  message("  Checking flipped SNPs")
  ## See if any SNPs need to be flipped
  genotype_flip.df <- merge(
    genotype_mat_input[, c("snp", "a1", "a2")],
    unique(betas[, c("snp", "a1", "a2")]),
    by = "snp"
  )

  to_remove = genotype_flip.df[a1.x != a1.y]$snp

  genotype_mat_input = genotype_mat_input[!snp %in% to_remove]
  betas = betas[!snp %in% to_remove]

  #check for duplicates and remove
  to_remove = genotype_mat_input[duplicated(snp)]$snp
  genotype_mat_input = genotype_mat_input[!snp %in% to_remove]
  betas = betas[!snp %in% to_remove]

  rm(to_remove, in_both, genotype_flip.df)

  message("  Creating coefficient matrix")
  ## Change betas data frame to be a matrix with rows=genes and cols=snps
  ## where entries are the betas
  all_assoc_genes <- dMcast(
    betas[, c("snp", "pid", "beta_shrink")],
    pid ~ snp,
    value.var = "beta_shrink",
    fun.aggregate = "sum"
  )

  ## dMcast appends the variable name (ie "snp") to column names, so remove this
  colnames(all_assoc_genes) <- sub("^snp", "", colnames(all_assoc_genes))

  ## Set order of coefficient matrix so it lines up with genotype matrix
  all_assoc_genes <- all_assoc_genes[, order(colnames(all_assoc_genes)), drop = FALSE]

  ## Reorder SNPs so they match up with gene-snp beta matrix
  genotype_mat_input <- genotype_mat_input[order(snp)]

  ## All SNPs are in same order between the two datasets
  all.equal(genotype_mat_input$snp, colnames(all_assoc_genes))
  #setdiff(colnames(all_assoc_genes), genotype_mat$snp )

  dim(genotype_mat_input)

  message("  Creating sparse genotype matrix")
  ## Rearrange datasets to get required matrices
  #previously started at col 6 - chsnge to 10 for vcf
  genotype_mat <- Matrix(
    as.matrix(genotype_mat_input[, genotype_start_col:ncol(genotype_mat_input)]),
    sparse = TRUE
  )
  genotype_mat <- t(genotype_mat)
  colnames(genotype_mat) <- genotype_mat_input$snp
  all_assoc_genes <- t(all_assoc_genes)

  rm(genotype_mat_input)
  gc()

  ## Check dimensions are as required for matrix multiplication
  dim(genotype_mat)
  dim(all_assoc_genes)

  message("  Generating gene expression")
  gene_expression <- genotype_mat %*% all_assoc_genes

  dim(gene_expression)

  message("  Adding random error")
  set.seed(seed)
  gene_expression = gene_expression +
    rnorm(n = ncol(gene_expression) * nrow(gene_expression), 0, 1)

  gene_expression
}

#' Export data for use by MatrixEQTL
#'
#' @param genotype_mat Genotype matrix with SNPs as rows and participants as columns
#' @param gene_expression Gene expression with genes as rows and participants as columns
#' @param combined
#' @param directory Which directory to save in
#'
#' @return A list containing the filenames where the data is stored
#' @export
write_matrixeqtl_files <- function(
  genotype_mat = NULL,
  gene_expression = NULL,
  combined = NULL,
  genotype_start_col = 6,
  directory = "Output/matrixEQTL"
) {
  filenames <- list()

  if (!is.null(genotype_mat)) {
    g_out = data.table(
      id = genotype_mat$snp,
      as.matrix(genotype_mat[, genotype_start_col:ncol(genotype_mat)])
    )
    fwrite(g_out, file.path(directory, "SNP.txt"), sep = "\t")

    filenames[["snp"]] <- file.path(directory, "SNP.txt")
  }

  if (!is.null(gene_expression)) {
    g_out = t(gene_expression)
    g_out = data.table(id = row.names(g_out), as.matrix(g_out))
    fwrite(g_out, file.path(directory, "GE.txt"), sep = "\t")

    filenames[["expression"]] <-  file.path(directory, "GE.txt")
  }

  if (!is.null(combined)) {
    g_out = unique(
      data.table(snp = combined$snp, chr = "chr19", pos = combined$sid_pos)
    )
    fwrite(g_out, file.path(directory, "snpsloc.txt"), sep = "\t")

    filenames[["snp_loc"]] <- file.path(directory, "snpsloc.txt")

    g_out = unique(
      data.table(
        geneid = combined$pid,
        chr = "chr19",
        s1 = combined$gene_start,
        s2 = combined$gene_end
      )
    )
    fwrite(g_out, file.path(directory, "geneloc.txt"), sep = "\t")

    filenames[["gene_loc"]] <- file.path(directory, "geneloc.txt")
  }

  filenames
}


run_matrixeqtl <- function(
  SNP_file_name,
  snps_location_file_name,
  expression_file_name,
  gene_location_file_name,
  output_file_name_cis,
  output_file_name_tra,
  pvOutputThreshold_cis = 1,
  pvOutputThreshold_tra = 1e-3,
  cisDist = 1e6
) {
  library(MatrixEQTL)
  library(data.table)

  ## Settings
  base.dir = "Output"
  useModel = modelLINEAR;
  covariates_file_name = character();
  errorCovariance = numeric();

  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);

  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);

  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name) > 0) {
    cvrt$LoadFile(covariates_file_name);
  }

  ## Run the analysis
  snpspos = fread(snps_location_file_name);
  genepos = fread(gene_location_file_name);

  ## If output cis file ends in "gz", then open a gzip connection
  if (grepl(".gz$", output_file_name_cis)) {
    output_file_name_cis <- gzcon(file(output_file_name_cis, "wb"), level = 9)
    cis_con <- TRUE
  } else {
    cis_con <- FALSE
  }

  if (grepl(".gz$", output_file_name_tra)) {
    output_file_name_tra <- gzcon(file(output_file_name_tra, "wb"), level = 9)
    tra_con <- TRUE
  } else {
    tra_con <- FALSE
  }

  me = Matrix_eQTL_main(
    snps                  = snps,
    gene                  = gene,
    cvrt                  = cvrt,
    output_file_name      = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel              = useModel,
    errorCovariance       = errorCovariance,
    verbose               = FALSE,
    output_file_name.cis  = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos               = snpspos,
    genepos               = genepos,
    cisDist               = cisDist,
    pvalue.hist           = "qqplot",
    min.pv.by.genesnp     = FALSE,
    noFDRsaveMemory       = TRUE
  );

  ## gzip connection needs to be closed otherwise some data is not written to
  ## the file
  if (cis_con) {
    close(output_file_name_cis)
  }

  if (tra_con) {
    close(output_file_name_tra)
  }

  me
}

calculate_power <- function(
  results_cis,
  combined,
  n_filtered = 382146,
  n_all = 11102563
) {
  library(data.table)

  ## Merge the results back onto the original SNP dataset
  results_merge = merge(
    combined,
    results_cis,
    by.x = c("snp", "pid"),
    by.y = c("snps", "gene"),
    all.x = TRUE
  )

  results_merge = results_merge[!is.na(statistic)]

  ## cis SNPs are within 1000000bp of gene start/end
  results_merge[, cis_flag := FALSE]
  results_merge[
    sid_pos >= (gene_start - 1000000 - 1) & sid_pos <= (gene_end + 1000000),
    cis_flag := TRUE
  ]

  results_merge[, causal_T := qtl_tag]

  ## Calculate FDR values and Bonferroni-adjusted p-values (on all variants,
  ## not just causal variants)
  results_merge[, fdr_filtered := p.adjust(pvalue, method = "BH", n = n_filtered)]
  results_merge[, fdr_all      := p.adjust(pvalue, method = "BH", n = n_all)]

  p_all <- c(
    results_merge$pvalue,
    sample(results_merge[qtl_tag == FALSE]$pvalue, size = n_all - nrow(results_merge), replace = TRUE)
  )
  p_all <- p.adjust(p_all, method = "BH")
  message("p_all: ", length(p_all), "; based on ", length(results_merge[qtl_tag == FALSE]$pvalue))
  results_merge[, fdr_all2      := p_all[1:nrow(results_merge)]]

  results_merge[, by_filtered  := p.adjust(pvalue, method = "BY", n = n_filtered)]
  results_merge[, by_all       := p.adjust(pvalue, method = "BY", n = n_all)]

  results_merge[, bon_filtered := p.adjust(pvalue, method = "bonferroni", n = n_filtered)]
  results_merge[, bon_all      := p.adjust(pvalue, method = "bonferroni", n = n_all)]

  ## Subset to only causal variants
  # results_merge <- results_merge[qtl_tag == TRUE]

  ## Determine which variants are discoveries
  results_merge[, causal_M_filt := fdr_filtered < 0.05]
  results_merge[, causal_M_all  := fdr_all < 0.05]
  results_merge[, causal_M_all2  := fdr_all2 < 0.05]

  results_merge[, causal_M_filt_by := by_filtered < 0.05]
  results_merge[, causal_M_all_by  := by_all < 0.05]

  results_merge[, causal_M_filt_bon := bon_filtered < 0.05]
  results_merge[, causal_M_all_bon  := bon_all < 0.05]

  # results_merge[, got_right_filtered := (causal_M_filt == 1 & causal_T == 1)]
  # results_merge[, got_right_all      := (causal_M_all  == 1 & causal_T == 1)]
  #
  # results_merge[, got_right_filtered_by := (causal_M_filt_by == 1 & causal_T == 1)]
  #
  # results_merge[, got_right_filtered_bon := (causal_M_filt_bon == 1 & causal_T == 1)]
  # results_merge[, got_right_all_bon      := (causal_M_all_bon  == 1 & causal_T == 1)]

  ## Create bins for effect size
  es_bin_cutoffs <- c(-Inf, 0.025, 0.057, 0.13, 0.265, Inf)
  results_merge[, es_bin := cut(
    abs(alpha),
    es_bin_cutoffs,
    right = FALSE,
    labels = c("< 0.025", "0.025-0.057", "0.057-0.13", "0.13-0.265", "0.265+")
  )]

  ## Create bins for MAF
  results_merge[, maf_bin := "10-20"]
  results_merge[maf > 0.20, maf_bin := "20-30"]
  results_merge[maf > 0.30, maf_bin := "30-40"]
  results_merge[maf > 0.40, maf_bin := "40-50"]

  results_merge[, num_correct_bin_filtered := sum(causal_M_filt), by = c("es_bin", "maf_bin", "qtl_tag")]
  results_merge[, num_correct_bin_all      := sum(causal_M_all), by = c("es_bin", "maf_bin", "qtl_tag")]
  results_merge[, num_correct_bin_all2      := sum(causal_M_all2), by = c("es_bin", "maf_bin", "qtl_tag")]

  results_merge[, num_correct_bin_filtered_by := sum(causal_M_filt_by), by = c("es_bin", "maf_bin", "qtl_tag")]
  results_merge[, num_correct_bin_all_by      := sum(causal_M_all_by), by = c("es_bin", "maf_bin", "qtl_tag")]

  results_merge[, num_correct_bin_filtered_bon := sum(causal_M_filt_bon), by = c("es_bin", "maf_bin", "qtl_tag")]
  results_merge[, num_correct_bin_all_bon      := sum(causal_M_all_bon), by = c("es_bin", "maf_bin", "qtl_tag")]

  results_merge[, num_bin := .N, by = c("es_bin", "maf_bin", "qtl_tag") ]

  sub = unique(
    results_merge[, c(
      "num_correct_bin_filtered",
      "num_correct_bin_all",
      "num_correct_bin_all2",
      "num_correct_bin_filtered_by",
      "num_correct_bin_all_by",
      "num_correct_bin_filtered_bon",
      "num_correct_bin_all_bon",
      "num_bin",
      "es_bin",
      "maf_bin",
      "qtl_tag"
    )]
  )

  sub[, power_filtered_fdr := num_correct_bin_filtered/num_bin]
  sub[, power_all_fdr      := num_correct_bin_all/num_bin]
  sub[, power_all_fdr2      := num_correct_bin_all2/num_bin]

  sub[, power_filtered_by := num_correct_bin_filtered_by/num_bin]
  sub[, power_all_by      := num_correct_bin_all_by/num_bin]

  sub[, power_filtered_bon := num_correct_bin_filtered_bon/num_bin]
  sub[, power_all_bon      := num_correct_bin_all_bon/num_bin]

  sub
}

calculate_power_egenes <- function(
  results_cis,
  n_filtered = 382146,
  n_all = 11102563
) {
  library(data.table)

  ## Merge the results back onto the original SNP dataset
  results_merge = merge(
    combined,
    results_cis,
    by.x = c("snp", "pid"),
    by.y = c("snps", "gene"),
    all.x = TRUE
  )

  results_merge = results_merge[!is.na(statistic)]

  ## cis SNPs are within 1000000bp of gene start/end
  results_merge[, cis_flag := FALSE]
  results_merge[
    sid_pos >= (gene_start - 1000000 - 1) & sid_pos <= (gene_end + 1000000),
    cis_flag := TRUE
  ]

  results_merge[, causal_T := qtl_tag]

  message(nrow(results_merge), ", ", ncol(results_merge))

  results_merge[, fdr_filtered := p.adjust(pvalue, method = "BH", n = n_filtered)]
  results_merge[, by_filtered := p.adjust(pvalue, method = "BY", n = n_filtered)]

  results_merge[, fdr_all      := p.adjust(pvalue, method = "BH", n = n_all)]

  results_merge[, bon_filtered := p.adjust(pvalue, method = "bonferroni", n = n_filtered)]
  results_merge[, bon_all      := p.adjust(pvalue, method = "bonferroni", n = n_all)]

  results_merge[, egene := any(beta_shrink != 0), by = pid]

  results_merge[, egene_sim_filt := any(fdr_filtered < 0.05 & qtl_tag == TRUE), by = pid]
  results_merge[, egene_sim_all  := any(fdr_all < 0.05 & qtl_tag == TRUE), by = pid]

  results_merge[, egene_sim_filt_bon := any(bon_filtered < 0.05 & qtl_tag == TRUE), by = pid]
  results_merge[, egene_sim_all_bon  := any(bon_all < 0.05 & qtl_tag == TRUE), by = pid]

  genes = unique(
    results_merge[, c("pid", "egene", "egene_sim_filt", "egene_sim_all", "egene_sim_filt_bon", "egene_sim_all_bon")]
  )

  # return(genes)

  message(" egenes: ", sum(genes$egene))

  genes[, got_right_filtered := (egene_sim_filt == 1 & egene == 1)]
  genes[, got_right_all      := (egene_sim_all  == 1 & egene == 1)]

  genes[, got_right_filtered_bon := (egene_sim_filt_bon == 1 & egene == 1)]
  genes[, got_right_all_bon      := (egene_sim_all_bon  == 1 & egene == 1)]

  genes[, num_correct_bin_filtered := sum(got_right_filtered), by = c("egene")]
  genes[, num_correct_bin_all      := sum(got_right_all), by = c("egene")]

  # genes[, num_correct_bin_filtered_by := sum(got_right_filtered_by), by = c("egene")]

  genes[, num_correct_bin_filtered_bon := sum(got_right_filtered_bon), by = c("egene")]
  genes[, num_correct_bin_all_bon      := sum(got_right_all_bon), by = c("egene")]

  genes[, num_bin := .N, by = c("egene") ]

  # return(genes)

  sub = unique(
    genes[, c(
      "egene",
      "num_bin",
      "num_correct_bin_filtered",
      "num_correct_bin_all",
      "num_correct_bin_filtered_bon",
      "num_correct_bin_all_bon"
    )]
  )

  sub[, power_filtered_fdr := num_correct_bin_filtered/num_bin]
  sub[, power_all_fdr      := num_correct_bin_all/num_bin]

  # sub[, power_filtered_by := num_correct_bin_filtered_by/num_bin]

  sub[, power_filtered_bon := num_correct_bin_filtered_bon/num_bin]
  sub[, power_all_bon      := num_correct_bin_all_bon/num_bin]

  sub[egene == TRUE]
}
