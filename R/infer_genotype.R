#' Infer Parental Genotypes from Hybrid Seed Endosperm
#'
#' This function reads a VCF file generated from hybrid seed endosperm,
#' uses a Beta-binomial model to test the allele depth ratio (2:1),
#' and infers the maternal and paternal genotypes for each hybrid sample.
#'
#' @param vcf_file A character string specifying the path to the input VCF file.
#'   Note: This single VCF file should/can contain genotype data for MULTIPLE F1 hybrid endosperm samples.
#' @param output_file A character string specifying the path to the output CSV file. Default is "parental_genotype.csv".
#' @param min_depth An integer specifying the minimum total read depth (REF + ALT) required to consider a site. Default is 150.
#' @param rho A numeric value specifying the overdispersion parameter for the Beta-binomial model. Default is 0.001.
#' @param alpha A numeric value specifying the significance level. Default is 0.05.
#'
#' @return A data frame containing the inferred maternal and paternal genotypes for each hybrid.
#' @export
#' @importFrom vcfR read.vcfR getCHROM getPOS
#' @importFrom VGAM pbetabinom
#' @importFrom utils write.table
infer_parental_genotype <- function(vcf_file,
                                    output_file = "parental_genotype.csv",
                                    min_depth = 150,
                                    rho = 0.001,
                                    alpha = 0.05) {

  message("Step 1: Reading VCF file...")
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  message("Step 2: Extracting Allele Depth (AD) from FORMAT column...")
  dat <- vcf@gt
  #
  extAD <- function(X){
    x <- strsplit(X, ":")
    return(unlist(lapply(x, function(x) x[2])))
  }

  dat <- apply(dat, 2, extAD)
  #
  if (is.null(dim(dat))) {
    stop("Error: Failed to extract AD. Please ensure the VCF FORMAT column contains 'AD'.")
  }

  dat <- data.frame(Chr = vcfR::getCHROM(vcf), Pos = vcfR::getPOS(vcf), dat)
  dat$Chr <- paste0(dat$Chr, "_", dat$Pos)
  dat <- data.frame(dat)

  #
  dat[] <- lapply(dat, function(x) gsub("\\.", ".,.", x))
  dat <- dat[, -c(2, 3)] # 删除原有的 Pos 列

  data <- data.frame(markerID = dat$Chr)
  hybrids <- colnames(dat)[2:ncol(dat)]

  message("Step 3: Inferring genotypes for ", length(hybrids), " samples using Beta-binomial model...")

  for (hybrid in hybrids) {
    # 提取 REF 和 ALT reads
    f1 <- unlist(strsplit(dat[, hybrid], ","))
    f1 <- data.frame(Pos = dat[, 1],
                     REF = as.numeric(f1[2*c(1:(length(f1)/2))-1]),
                     ALT = as.numeric(f1[2*c(1:(length(f1)/2))]))

    # 过滤深度不够的位点
    f1 <- f1[f1$REF + f1$ALT > min_depth, ]
    dats <- f1[!is.na(f1$REF), ]

    if (nrow(dats) == 0) {
      warning("Sample ", hybrid, " has no sites passing the depth filter (min_depth = ", min_depth, ").")
      next
    }

    # 定义 Beta-binomial 检验函数
    bbinom_test <- function(REF, ALT, p0, rho) {
      n <- REF + ALT
      k <- ALT
      if (n == 0) return(NA_real_)

      p_lower <- VGAM::pbetabinom(q = k, size = n, prob = p0, rho = rho)
      p_upper <- if (k == 0) 1 else 1 - VGAM::pbetabinom(q = k - 1, size = n, prob = p0, rho = rho)
      p_value <- min(1, 2 * min(p_lower, p_upper))
      return(p_value)
    }

    # 检验 1:3 和 2:3 比例
    dats$P_1_3 <- mapply(bbinom_test, REF = dats$REF, ALT = dats$ALT, MoreArgs = list(p0 = 1/3, rho = rho))
    dats$P_2_3 <- mapply(bbinom_test, REF = dats$REF, ALT = dats$ALT, MoreArgs = list(p0 = 2/3, rho = rho))

    dats$Fit_1_3 <- ifelse(is.na(dats$P_1_3), NA, dats$P_1_3 >= alpha)
    dats$Fit_2_3 <- ifelse(is.na(dats$P_2_3), NA, dats$P_2_3 >= alpha)

    # 判定亲本基因型
    dats$paternal <- ifelse(dats$ALT == 0, 0,
                            ifelse(dats$REF == 0, 2,
                                   ifelse(dats$Fit_1_3 == TRUE, 2,
                                          ifelse(dats$Fit_2_3 == TRUE, 0, NA))))

    dats$maternal <- ifelse(dats$ALT == 0, 0,
                            ifelse(dats$REF == 0, 2,
                                   ifelse(dats$Fit_1_3 == TRUE, 0,
                                          ifelse(dats$Fit_2_3 == TRUE, 2, NA))))

    # 匹配回主表
    data[, paste0(hybrid, "_maternal")] <- dats$maternal[match(data$markerID, dats$Pos)]
    data[, paste0(hybrid, "_paternal")] <- dats$paternal[match(data$markerID, dats$Pos)]
  }

  sample_cols <- 2:ncol(data)
  valid_rows <- rowSums(!is.na(data[, sample_cols, drop = FALSE])) > 0
  data <- data[valid_rows, ]

  message("Filtering complete. Remaining valid markers: ", nrow(data))
  
  message("Step 4: Saving results to ", output_file, "...")
  utils::write.table(data, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)

  return(data)
}
