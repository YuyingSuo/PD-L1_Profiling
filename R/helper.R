library(GSVA)
library(impute)
library(qs)

pdl1types=c('GBM','HNSCC','LSCC','LUAD','OV','PDAC','SCLC')
err_txt ='Error: Quantitative PD-L1 protein data are insufficient and failed to meet the 50% missing value threshold, so the result cannot be provided.'

HCC_err='Error: Transcriptomic data of HCC have not been obtained yet.'
rna_noNAT=c('AML','COAD','GBM','OV')
war_txt ='Warning: Transcriptomic data of NAT samples are unavailable for this cancer type.'


data_dir <- file.path("Data")



#### read data for exp, DEG, PPI  ####

dataprocess1 <- function(cn, type1) {
  if (type1 == "Proteomic") {
    if (!(cn %in% pdl1types)) {
      err <- err_txt
      return(list(expr_data = NULL, err = err))
    } else {
      file_name <- paste0("total_", cn, "_forprocess.qs")
      file_path <- file.path(data_dir, file_name)
      expr_data <- qread(file_path)
      gc(FALSE, TRUE)
      return(list(expr_data = expr_data, err = NULL))
    }
    
  } else {
    t_name <- paste0("RNAseq_", cn, "_log2_T.qs")
    n_name <- paste0("RNAseq_", cn, "_log2_N.qs")
    if (cn %in% rna_noNAT) {
      war <- war_txt
      rna_t <- qread(file.path(data_dir, t_name))
      gc(FALSE, TRUE)
      return(list(expr_data = rna_t, err = war))
      
    } else {
      t_path <- file.path(data_dir, t_name)
      n_path <- file.path(data_dir, n_name)
      rna_t <- qread(t_path)
      rna_n <- qread(n_path)
      colnames(rna_n) <- paste0(colnames(rna_n), "_N")
      colnames(rna_t) <- paste0(colnames(rna_t), "_T")
      
      expr_data <- do.call(cbind, list(rna_n, rna_t))
      rm(rna_n, rna_t)
      gc(FALSE, TRUE)
      
      return(list(expr_data = expr_data, err = NULL))
    }
  }
}




#### read data for pathway  ####
dataprocess_path <- function(cn, type) {
  if (type == "Proteomic") {
    file_name <- paste0("total_", cn, "_forprocess.qs")
    file_path <- file.path(data_dir, file_name)
    expr_data <- qread(file_path)
    gc(FALSE, TRUE)
    
    n <- expr_data[, grep("_N", colnames(expr_data)), drop = FALSE]
    t <- expr_data[, -grep("_N", colnames(expr_data)), drop = FALSE]
    rm(expr_data); gc(FALSE, TRUE)
    
    return(list(n = n, t = t, err = NULL))
    
  } else {
    t_name <- paste0("RNAseq_", cn, "_log2_T.qs")
    n_name <- paste0("RNAseq_", cn, "_log2_N.qs")
    
    if (cn %in% rna_noNAT) {
      war <- war_txt
      rna_t <- qread(file.path(data_dir, t_name))
      gc(FALSE, TRUE)
      return(list(n = NULL, t = rna_t, err = war))
      
    } else {
      t_path <- file.path(data_dir, t_name)
      n_path <- file.path(data_dir, n_name)
      rna_t <- qread(t_path)
      rna_n <- qread(n_path)
      
      colnames(rna_n) <- paste0(colnames(rna_n), "_N")
      colnames(rna_t) <- paste0(colnames(rna_t), "_T")
      gc(FALSE, TRUE)
      
      return(list(n = rna_n, t = rna_t, err = NULL))
    }
  }
}











#### DEG process  ####
deg_dataprocess <- function(expr_t) {
  gene_id <- "CD274"
  gi <- match(gene_id, rownames(expr_t))
  
  cd274_vec <- as.numeric(expr_t[gi, ])
  qs <- quantile(cd274_vec, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE)
  cut_1 <- qs[1]; cut_2 <- qs[2]
  print(cut_1); print(cut_2)
  
  low_idx  <- which(cd274_vec <= cut_1)
  high_idx <- which(cd274_vec >  cut_2)

  expr_mat <- as.matrix(expr_t)
  storage.mode(expr_mat) <- "double"

  dt_low  <- expr_mat[, low_idx,  drop = FALSE]
  dt_high <- expr_mat[, high_idx, drop = FALSE]
  rm(expr_mat); gc(FALSE, TRUE) 
  
  bad_low  <- rowMeans(is.na(dt_low))  > 0.5
  bad_high <- rowMeans(is.na(dt_high)) > 0.5
  keep <- !(bad_low & bad_high)
  if (!any(keep)) stop("All rows filtered by NA>50% in both groups.")
  
  dt_low  <- dt_low[keep, ,  drop = FALSE]
  dt_high <- dt_high[keep, , drop = FALSE]
  kept_genes <- rownames(dt_low)
  

  r <- nrow(dt_low)
  med_low <- numeric(r)
  med_high <- numeric(r)
  pvalue <- numeric(r)
  
  block_size <- 1000L
  for (start in seq(1, r, by = block_size)) {
    end <- min(start + block_size - 1L, r)
    rng <- start:end
    med_low[rng]  <- matrixStats::rowMedians(dt_low[rng, , drop = FALSE],  na.rm = TRUE)
    med_high[rng] <- matrixStats::rowMedians(dt_high[rng, , drop = FALSE], na.rm = TRUE)
    pvalue[rng]   <- matrixTests::row_wilcoxon_twosample(
      dt_high[rng, , drop = FALSE],
      dt_low[rng, , drop = FALSE],
      alternative = "two.sided", exact = FALSE
    )$pvalue
    gc(FALSE, TRUE)
  }
  
  log2fc <- med_high - med_low
  padj <- p.adjust(pvalue, method = "BH")
  
  thr <- log2(1.5)
  change_padj <- factor(
    ifelse(padj < 0.05 & abs(log2fc) >= thr,
           ifelse(log2fc >= thr, "Up", "Down"),
           "NoSignifi"),
    levels = c("Up","Down","NoSignifi")
  )
  colnames(dt_low)  <- paste0(colnames(expr_t)[low_idx],  "_low")
  colnames(dt_high) <- paste0(colnames(expr_t)[high_idx], "_high")
  
  dt_bind <- cbind(dt_low, dt_high)
  rm(dt_low, dt_high); gc(FALSE, TRUE)
  
  res <- data.frame(
    Log2Fc_high_low = log2fc,
    p = pvalue,
    p.adj = padj,
    change_padj = change_padj,
    row.names = kept_genes,
    check.names = FALSE
  )
  
  res <- cbind(dt_bind, res)
  res
}




#### cor process  ####
cor_dataprocess <- function(item, expr_t) {
  cd274 <- as.numeric(expr_t[which(rownames(expr_t) == item), ])
  l <- ncol(expr_t)
  n <- nrow(expr_t)
  M <- as.matrix(expr_t)
  storage.mode(M) <- "double"
  
  Rho <- numeric(n)
  pval <- numeric(n)
  
  block_size <- 1000L
  for (start in seq(1, n, by = block_size)) {
    end <- min(start + block_size - 1L, n)
    for (i in start:end) {
      tmp <- M[i, 1:l]
      ct <- suppressWarnings(cor.test(cd274, tmp, method = "spearman"))
      Rho[i] <- ct$estimate
      pval[i] <- ct$p.value
    }
    gc(FALSE, TRUE)
  }
  
  expr_t$Rho <- Rho
  expr_t$p <- pval
  expr_t$p.adj <- p.adjust(pval, method = "BH")
  
  expr_t$cortype <- factor(
    ifelse(expr_t$p.adj < 0.05 & abs(expr_t$Rho) >= 0.4,
           ifelse(expr_t$Rho >= 0.4, "Pos", "Neg"),
           "NoSig"),
    levels = c("Pos", "Neg", "NoSig")
  )
  
  return(expr_t)
}













#### ssgsea  ####
.load_kegg_once <- local({
  .gmt <- NULL; .union <- NULL
  function(p) {
    if (is.null(.gmt)) {
      load(p); g <- dt; rm(dt)
      if (is.null(names(g))) names(g) <- paste0("gs", seq_along(g))
      .gmt <<- g; .union <<- unique(unlist(g, use.names = FALSE))
    }
    list(gmt = .gmt, union = .union)
  }
})

ssgsea_prepare <- function(expr, union, na_row_thresh = 0.5) {
  expr <- as.matrix(expr[intersect(rownames(expr), union), , drop = FALSE])
  if (anyNA(expr)) {
    bad <- rowMeans(is.na(expr)) > na_row_thresh
    if (any(bad)) expr <- expr[!bad, , drop = FALSE]
    if (anyNA(expr)) {
      expr <- impute::impute.knn(
        expr,
        k = min(10L, max(3L, ncol(expr) - 1L)),
        rowmax = 0.5,  colmax = 0.8, maxp = Inf,rng.seed = 1
      )$data
    }
  }
  expr
}

.ssgsea_chunk <- function(expr, gmt, minSize = 5) {
  if (length(gmt) == 0L) return(NULL)
  keep <- vapply(gmt, \(x) {
    k <- sum(rownames(expr) %in% x); k >= minSize
  }, logical(1))
  gmt <- gmt[keep]; if (!length(gmt)) return(NULL)
  
  res <- GSVA::gsva(expr, gmt, method = "ssgsea", abs.ranking = TRUE,
                    min.sz = minSize, parallel.sz = 1, verbose = FALSE)
  if (is.null(dim(res)))
    res <- matrix(res, 1, dimnames = list(names(gmt)[1], colnames(expr)))
  res
}

ssgsea_run_chunked <- function(expr, gmt, gsize = 6, ssize = 30, minSize = 5) {
  gsize <- max(1L, as.integer(gsize)); ssize <- max(1L, as.integer(ssize))
  ns <- ncol(expr); if (ns < 1L) return(matrix(0,0,0))
  gnames <- if (is.null(names(gmt))) paste0("gs", seq_along(gmt)) else names(gmt)
  gchunks <- split(seq_along(gnames), ceiling(seq_along(gnames)/gsize))
  schunks <- split(seq_len(ns), ceiling(seq_len(ns)/ssize))
  res_list <- vector("list", length(schunks))
  
  for (i in seq_along(schunks)) {
    e <- expr[, schunks[[i]], drop = FALSE]
    tmp <- NULL
    for (j in seq_along(gchunks)) {
      g <- gmt[gnames[gchunks[[j]]]]
      r <- .ssgsea_chunk(e, g, minSize)
      if (!is.null(r)) tmp <- if (is.null(tmp)) r else rbind(tmp, r)
      rm(g, r); gc()
    }
    res_list[[i]] <- tmp; rm(e, tmp); gc()
  }
  
  valid <- Filter(\(m) !is.null(m) && nrow(m)>0, res_list)
  if (!length(valid)) return(matrix(0,0,ns))
  paths <- Reduce(intersect, lapply(valid, rownames))
  out <- NULL
  for (m in valid) {
    mb <- m[paths,,drop=FALSE]
    out <- if (is.null(out)) mb else cbind(out, mb)
    rm(mb); gc()
  }
  out
}

ssgsea_process <- function(expr, path = file.path(data_dir,"KEGG_1.42.0_genes.RData"),
                           gsize = 6, ssize = 30, minSize = 5) {
  kg <- .load_kegg_once(path)
  expr2 <- ssgsea_prepare(expr, kg$union)
  res <- ssgsea_run_chunked(expr2, kg$gmt, gsize, ssize, minSize)
  rm(expr2); gc(); res
}

ssgsea_testprocess <- function(cn, type, path = file.path(data_dir,"KEGG_1.42.0_genes.RData"),
                               gsize = 6, ssize = 30, minSize = 5) {
  d <- dataprocess_path(cn, type)
  if (is.null(d$n)) {
    t <- ssgsea_process(d$t, path, gsize, ssize, minSize)
    list(cn=cn,type=type,ssgsea_data=t,ssgsea_t=t,err=war_txt)
  } else {
    n <- ssgsea_process(d$n, path, gsize, ssize, minSize)
    t <- ssgsea_process(d$t, path, gsize, ssize, minSize)
    ov <- intersect(rownames(n), rownames(t))
    n2 <- n[ov,,drop=FALSE]; t2 <- t[ov,,drop=FALSE]
    list(cn=cn,type=type,ssgsea_data=cbind(n2,t2),ssgsea_t=t,err=NULL)
  }
}





#### pho ####
pho_dataprocess <- function(cn, data) {
  message("carried pho_test")
  df <- data.frame(data, sample = rownames(data))
  df <- df[!grepl("_N", rownames(df)), , drop = FALSE]
  df <- df[order(df$ssgsea_PD_L1_path), , drop = FALSE]
  n <- nrow(df); k <- floor(n / 3)
  low  <- rownames(df)[1:k]
  high <- rownames(df)[(n - k + 1):n]

  pho_path <- file.path(data_dir, paste0("pho_", cn, "_forprocess.qs"))
  pho <- qs::qread(pho_path)
  pho$pos <- sub(".*_", "", pho$Index)
  rn <- paste(pho$Gene, pho$pos, sep = "_"); rownames(pho) <- rn

  low_exist  <- low[low %in% colnames(pho)]
  high_exist <- high[high %in% colnames(pho)]
  
  if (length(low_exist) == 0L && length(high_exist) == 0L) {
    warning("No overlapping samples found in phosphorylation data.")
    return(data.frame())
  }
  
  message("block calculating")
  r <- nrow(pho)
  log2fc <- numeric(r); p <- numeric(r)
  block_size <- 400L
  
  for (s in seq(1L, r, by = block_size)) {
    e <- min(s + block_size - 1L, r)
    low_m  <- as.matrix(pho[s:e, low_exist,  drop = FALSE])
    high_m <- as.matrix(pho[s:e, high_exist, drop = FALSE])
    log2fc[s:e] <- matrixStats::rowMedians(high_m, na.rm = TRUE) - 
      matrixStats::rowMedians(low_m, na.rm = TRUE)
    for (i in seq_len(nrow(low_m))) {
      p[s + i - 1L] <- suppressWarnings(
        wilcox.test(high_m[i, ], low_m[i, ])$p.value
      )
    }
    rm(low_m, high_m); gc(FALSE, TRUE)
  }

  padj <- p.adjust(p, "BH")
  change <- ifelse(padj < 0.05 & abs(log2fc) >= log2(1.5),
                   ifelse(log2fc >= log2(1.5), "Up", "Down"), "NoSignifi")
  

  pho_low  <- pho[, low_exist,  drop = FALSE]
  pho_high <- pho[, high_exist, drop = FALSE]
  colnames(pho_low)  <- paste0(colnames(pho_low),  "_PDL1_ssgsea_Inactive")
  colnames(pho_high) <- paste0(colnames(pho_high), "_PDL1_ssgsea_Active")
  
  out <- data.frame(
    pho_low,
    pho_high,
    Log2Fc_High_Low = log2fc,
    p = p,
    padj = padj,
    change = factor(change, levels = c("Up", "Down", "NoSignifi")),
    row.names = rn, check.names = FALSE
  )
  
  rm(pho, pho_low, pho_high); gc(FALSE, TRUE)
  return(out)
}