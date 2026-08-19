library(Matrix)
library(anndata)

aggregate_group <- function(x, grouping, func = colMeans){
  grouping <- factor(grouping)
  t(sapply(levels(grouping), function(g) func(x[which(grouping == g),])))
}

rename_gene_symbol <- function(mat, rename){
  mat <- mat[intersect(rownames(mat), names(rename)),]
  rownames(mat) <- rename[rownames(mat)]
  return(mat)
}

rowNorm <- function(mat){
  diag_mat <- Diagonal(x = 1/rowSums(mat))
  res <- diag_mat %*% mat
  dimnames(res) <- dimnames(mat)
  return(res)
}

colNorm <- function(mat){
  diag_mat <- Diagonal(x = 1/colSums(mat))
  res <- mat %*% diag_mat
  dimnames(res) <- dimnames(mat)
  return(res)
}

expCounts <- function(mat){
  if (is(mat, "Matrix")){
    df_mat <- summary(mat)
    df_mat$x <- exp(df_mat$x)-1
    mat <- sparseMatrix(i = df_mat$i, j = df_mat$j, x = df_mat$x, dims = dim(mat), dimnames = dimnames(mat))
  } else{
    mat <- exp(mat)-1
  }
  return(mat)
}

colMeans_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), mean)
}

rowMeans_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  res <- rep(0, nrow(dgCMat))
  res[sort(unique(RowInd))] <- sapply(split(dgCMat@x, RowInd), mean)
  return(setNames(res, rownames(dgCMat)))
}

merge_count_matrices <- function(counts_1,
                                 counts_2,
                                 rename_features = NULL,
                                 only_intersect = FALSE,
                                 suffix_1 = NULL,
                                 suffix_2 = NULL)
{
  if (! is.null(suffix_1))
    colnames(counts_1) <- paste0(colnames(counts_1), suffix_1)
  if (! is.null(suffix_2))
    colnames(counts_2) <- paste0(colnames(counts_2), suffix_2)
  
  if (length(intersect(colnames(counts_1), colnames(counts_2))) > 0){
    new_colnames <- make.unique(c(colnames(counts_1), colnames(counts_2)))
    colnames(counts_1) <- new_colnames[1:ncol(counts_1)]
    colnames(counts_2) <- new_colnames[(ncol(counts_1)+1):(ncol(counts_1)+ncol(counts_2))]
  }
  
  if (! is.null(rename_features)){
    if (length(intersect(names(rename_features), rownames(counts_1)))>0)
      rownames(counts_1)[rownames(counts_1) %in% names(rename_features)] <- rename_features[rownames(counts_1)[rownames(counts_1) %in% names(rename_features)]]
    
    if (length(intersect(names(rename_features), rownames(counts_2)))>0)
      rownames(counts_2)[rownames(counts_2) %in% names(rename_features)] <- rename_features[rownames(counts_2)[rownames(counts_2) %in% names(rename_features)]]
  }
  
  if (only_intersect){
    features <- intersect(rownames(counts_1), rownames(counts_2))
    counts_1 <- counts_1[features,]
    counts_2 <- counts_2[features,]
  } else{
    features <- union(rownames(counts_1), rownames(counts_2))
    
    newidx_1 <- setNames(1:length(features), features)[rownames(counts_1)]
    df_counts_1 <- summary(counts_1)
    df_counts_1$i <- newidx_1[df_counts_1$i]
    counts_1 <- sparseMatrix(i = df_counts_1$i, j = df_counts_1$j, x = df_counts_1$x, dims = c(length(features), ncol(counts_1)), dimnames = list(features, colnames(counts_1)))
    
    newidx_2 <- setNames(1:length(features), features)[rownames(counts_2)]
    df_counts_2 <- summary(counts_2)
    df_counts_2$i <- newidx_2[df_counts_2$i]
    counts_2 <- sparseMatrix(i = df_counts_2$i, j = df_counts_2$j, x = df_counts_2$x, dims = c(length(features), ncol(counts_2)), dimnames = list(features, colnames(counts_2)))
  }
  
  counts <- cbind(counts_1, counts_2)
  return(counts)
}

summary_guides_to_targets <- function(guide_mat, guide_target){
  mat <- t(sapply(unique(guide_target), function(target){
    guides <- names(guide_target)[which(guide_target == target)]
    if (length(guides) == 1)
      return(setNames(guide_mat[guides,], colnames(guide_mat)))
    colSums(guide_mat[guides,])
  }))
  mat <- Matrix(mat, sparse = T, dimnames = list(unique(guide_target), colnames(guide_mat)))
  return(mat)
}


summarize_data_to_groups.Matrix <- function(object, groups, bycol = T)
{
  groups <- factor(groups)
  summfunc <- Matrix::rowMeans
  if (bycol) summfunc <- Matrix::colMeans
  
  summ_mat <- sapply(levels(groups), function(x){
    if (bycol) obj <- object[which(groups==x),]
    if (!bycol) obj <- object[,which(groups==x)]
    summfunc(obj)
  })
  if (bycol) summ_mat <- t(summ_mat)
  summ_mat <- Matrix(summ_mat, sparse = T)
  return(summ_mat)
}

summarize_data_to_groups.data.frame <- function(object, groups)
{
  groups <- factor(groups)
  ident_mat_groups <- sparseMatrix(i = which(!is.na(groups)),
                                   j = as.numeric(groups)[!is.na(groups)],
                                   x = 1,
                                   dims = c(length(groups), length(levels(groups))), dimnames = list(names(groups), levels(groups)))
  ident_mat_groups_norm <- colNorm(ident_mat_groups)
  
  df_summ <- do.call(cbind.data.frame, lapply(1:ncol(object), function(i){
    dat <- object[,i]
    if (is.numeric(dat)){
      if (sum(! is.na(dat)) == 0)
        return(setNames(rep(NA, length(levels(groups))), levels(groups)))
      dat <- matrix(dat, nrow = 1, dimnames = list(colnames(object)[i], rownames(object)))
      dat[is.na(dat)] <- 0
      summ_dat <- dat %*% ident_mat_groups_norm
      return(setNames(as.numeric(summ_dat), colnames(ident_mat_groups_norm)))
    } else{
      if (sum(! is.na(dat))==0){
        if (is.null(levels(dat)))
          return(setNames(rep(NA, length(levels(groups))), levels(groups)))
        return(setNames(factor(rep(NA, length(levels(groups))), levels = levels(dat)), levels(groups)))
      }
      dat <- factor(dat)
      lab_mat_groups <- sparseMatrix(i = which(!is.na(dat)),
                                     j = as.numeric(dat)[which(!is.na(dat))],
                                     x = 1,
                                     dims = c(length(dat), length(levels(dat))), dimnames = list(rownames(object), levels(dat)))
      summ_mat <- t(lab_mat_groups) %*% ident_mat_groups
      summ_dat <- factor(setNames(rownames(summ_mat)[apply(summ_mat, 2, which.max)], levels(groups)), levels = levels(dat))
      return(summ_dat)
    }
  }))
  colnames(df_summ) <- colnames(object)
  
  return(df_summ)
}

summarize_data_to_groups <- function(object, ...) {
  UseMethod(generic = 'summarize_data_to_groups', object = object)
}





# convert data to h5ad format
data_to_h5ad.Seurat <- function(object,
                                assay = DefaultAssay(object),
                                savefile = NULL,
                                verbose = F)
{
  if (verbose)
    cat(">> start to create the anndata object...\n")
  adata <- anndata::AnnData(X = t(object[[assay]]@data),
                            obs = object@meta.data,
                            var = object[[assay]]@meta.features,
                            layers = list(count = t(object[[assay]]@counts)),
                            obsm = setNames(lapply(names(object@reductions), function(x) Embeddings(object, x)), names(object@reductions))
                            )
  if (verbose)
    cat(">> done.\n")
  
  if (!is.null(savefile)){
    if (verbose)
      cat(paste0(">> saving the anndata object to file: ",savefile,"\n"))
    adata$write_h5ad(savefile)
    if (verbose)
      cat(">> done.\n")
  }
  
  return(adata)
}

data_to_h5ad.default <- function(object,
                                 vars = NULL,
                                 obs = NULL,
                                 obsm = list(),
                                 layers = list(),
                                 savefile = NULL,
                                 verbose=F)
{
  if (verbose)
    cat(">> start to create the anndata object...\n")
  if (is.null(vars))
    vars <- data.frame(id = colnames(object), row.names = colnames(object))
  if (is.null(obs))
    obs <- data.frame(id = rownames(object), row.names = rownames(object))
  adata <- anndata::AnnData(X = object,
                            obs = obs,
                            var = vars,
                            obsm = obsm,
                            layers = layers
                            )
  if (verbose)
    cat(">> done.\n")
  
  if (!is.null(savefile)){
    if (verbose)
      cat(paste0(">> saving the anndata object to file: ",savefile,"\n"))
    adata$write_h5ad(savefile)
    if (verbose)
      cat(">> done.\n")
  }
  
  return(adata)
}

data_to_h5ad <- function(object, ...) {
  UseMethod(generic = 'data_to_h5ad', object = object)
}


# impute the node values by means of their knn, except for nodes marked to be remained as provided
# use to assign values to nodes with NA values at the beginning
# can also be used to smooth the data
# multiple values per node (object is a dense/sparse matrix)
diffusion_means.Matrix <- function(object,
                                   diffusion_mat = NULL,
                                   idx_init = NULL,
                                   knn_mat = NULL,
                                   mutual = FALSE,
                                   num_iter = 20,
                                   verbose = T){
  if (is.null(diffusion_mat)){
    if (verbose)
      message("Diffusion matrix is missing. Start to generate from the kNN and the initial nodes...")
    if (is.null(idx_init) | is.null(knn_mat))
      stop("Either diffusion_mat or both idx_init and knn_mat should not be NULL")
    if (mutual){
      diffusion_mat <- ((knn_mat != 0) & (t(knn_mat != 0))) * 1.0
    } else{
      diffusion_mat <- ((knn_mat != 0) | (t(knn_mat != 0))) * 1.0
    }
    diffusion_mat <- as(diffusion_mat, "Matrix")
    
    if (verbose == 2)
      message("Incorporating the initial nodes into the kNN matrix...")
    diffusion_mat[idx_init,] <- 0
    diag(diffusion_mat) <- 1
    diffusion_mat <- drop0(diffusion_mat)
    diffusion_mat <- rowNorm(diffusion_mat)
    
    if (verbose)
      message("Diffusion matrix generated.")
  }
  
  diffused <- object
  diffused[is.na(diffused)] <- 0
  diffused <- drop0(diffused)
  
  if (verbose)
    message("Start diffusion...")
  for(i in 1:num_iter){
    diffused <- diffusion_mat %*% diffused
    if (verbose == 2)
      message(paste0("Done diffusion for ",i," times..."))
  }
  
  if (verbose)
    message("Diffused mean calculated.")
  return(diffused)
}

diffusion_means.matrix <- function(object, ...){
  diffusion_means.Matrix(object, ...)
}

diffusion_means <- function(object, ...) {
  UseMethod(generic = 'diffusion_means', object = object)
}

# random walk with restart to smooth node values in a network
# it is a generalized version of the diffusion_means function (alpha = 0)
# but only applicable when all the initial values are non-NAs
# only works for one value per node (i.e. init is a vector not matrix)
random_walk_with_restart <- function(init, # initial vortex weight
                                     knn_idx = NULL, # knn index matrix (e.g. RANN::nn2$nn.idx)
                                     knn_mat = NULL, # alternative, the adj matrix
                                     symm = TRUE, # symmetric the matrix
                                     normalize = TRUE, # normalize rows (i.e. diffusion probability)
                                     alpha = 0.2, # prob of restart
                                     num_rounds = 100) # num of random walks
{
  if (is.matrix(init)) dim(init) <- c(1, nrow(init)*ncol(init))
  if (is.vector(init)) dim(init) <- c(1, length(init))
  
  if (is.null(knn_mat)){
    knn_df <- data.frame(i = rep(1:nrow(knn_idx), ncol(knn_idx)), j = as.numeric(knn_idx))
    knn_mat <- sparseMatrix(i = knn_df$i, j = knn_df$j, dims = c(nrow(knn_idx), nrow(knn_idx)))
  }
  
  if (symm)
    knn_mat <- (knn_mat + t(knn_mat)) > 0
  
  if (normalize){
    rowsums_mat <- rowSums(abs(knn_mat))
    knn_df <- summary(knn_mat)
    knn_df$x <- knn_df$x / rowsums_mat[knn_df$i]
    knn_df <- knn_df[!is.na(knn_df$x),]
    knn_mat <- sparseMatrix(i = knn_df$i, j = knn_df$j, x = knn_df$x, dims = dim(knn_mat), dimnames = dimnames(knn_mat))
  }
  
  walked <- init
  for (i in 1:num_rounds){
    walked <- alpha * init + (1 - alpha) * (walked %*% knn_mat)
  }
  return(walked[1,])
}


# group genes into categories based on their expression detection across cells
# five categories: undetected, low/medium/high/ultrahigh detection
grp_det_rates <- function(expr = NULL,
                          pseudocount = TRUE,
                          det_rates = NULL,
                          num_breaks = 40,
                          df_sm = 10,
                          thres_det_rates = c(0.05, 0.3, 0.7, 0.95),
                          verbose = FALSE)
{
  if (length(thres_det_rates) != 4){
    stop("The length of thres_det_rates should be 4.")
  }
  
  if (is.null(det_rates)){
    det_rates <- rowMeans(expr > 0)
    if (pseudocount)
      det_rates <- det_rates + 1 / ncol(expr)
    x <- log(det_rates)
  } else{
    x <- log(det_rates)
  }
  
  h <- hist(x, breaks = num_breaks, plot = FALSE)
  m <- smooth.spline(x = (h$breaks[-1]+h$breaks[-length(h$breaks)])/2, y = h$counts, df = df_sm)
  pred <- predict(m, seq(min(x[!is.infinite(x)]), max(x[!is.infinite(x)]), length.out=50))
  
  mean_det <- mean(pred$x[max(which(pred$y[-length(pred$y)]-pred$y[-1] < 0)):(max(which(pred$y[-length(pred$y)]-pred$y[-1] < 0))+1)])
  idx <- which(x>=mean_det)
  sd_det <- sqrt(sum((x[idx]-mean_det)^2) / (length(idx)-1))
  
  p_det <- pnorm(x, mean = mean_det, sd = sd_det)
  det_grp <- ifelse(p_det < thres_det_rates[1],
                    "ultralow",
                    ifelse(p_det < thres_det_rates[2],
                           "low",
                           ifelse(p_det < thres_det_rates[3],
                                  "medium",
                                  ifelse(p_det < thres_det_rates[4],
                                         "high",
                                         "ultrahigh"))
                    )
  )
  names(det_grp) <- names(det_rates)
  det_grp <- factor(det_grp, levels=c("ultralow", "low", "medium", "high", "ultrahigh"))
  
  if (verbose){
    max_det_rates <- tapply(det_rates, det_grp, max)
    
    message("estimation done. the estimated ranges:")
    message(paste0(" * ultralow: [0, ", max_det_rates["ultralow"],"]"))
    message(paste0(" * low: (", max_det_rates["ultralow"],", ", max_det_rates["low"],"]"))
    message(paste0(" * medium: (", max_det_rates["low"],", ", max_det_rates["medium"],"]"))
    message(paste0(" * high: (", max_det_rates["medium"],", ", max_det_rates["high"],"]"))
    message(paste0(" * ultrahigh: (", max_det_rates["high"],", ", max_det_rates["ultrahigh"],"]"))
  }
  
  return(det_grp)
}

