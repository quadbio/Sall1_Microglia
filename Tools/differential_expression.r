# natural spline based pseudotime test, to identify genes with expression changes along a continuous trajectory (e.g. pseudotime)
pt_dep_test_ns <- function(expr, # expression matrix, rows as genes
                           pt, # pseudotime or other continuous variables
                           covar = NULL, # covariates
                           df = 5, # natural spline degree of freedom
                           num_threads = 20) # number of threads
{
  library(doParallel)
  library(splines)
  
  registerDoParallel(num_threads)
  p <- foreach(i = 1:nrow(expr), .combine = rbind) %dopar%{
    e <- as.numeric(expr[i,])
    idx <- which(!is.na(e) & !is.na(pt))
    
    if (var(e[idx]) == 0)
      return(rep(NaN, 3))
    
    if (is.null(covar)){
      m0 <- lm(e[idx] ~ 1)
      m1 <- lm(e[idx] ~ ns(pt[idx], df = df))
    } else{
      data <- data.frame(y = e[idx], covar[idx,])
      m0 <- lm(y ~ ., data = data)
      data <- data.frame(y = e[idx], covar[idx,])
      m1 <- lm(y ~ ns(pt[idx], df = df) + ., data = data)
    }
    a0 <- anova(m0)
    a1 <- anova(m1)
    a01 <- anova(m0, m1)
    
    p1 <- a1$Pr[1]
    p2 <- a01$Pr[2]
    p3 <- pf(a0["Residuals","Mean Sq"] / a1["Residuals","Mean Sq"], df1 = a0["Residuals","Df"], df2 = a1["Residuals","Df"], lower.tail = F)
    return(c(p1, p2, p3))
  }
  p <- data.frame(p)
  rownames(p) <- rownames(expr)
  colnames(p) <- c("ANOVA", "ANCOVA", "Resi")
  return(p)
}

# ANOVA-based test for expression changed linearly related to a given factor/numeric variable
ancova_group_test <- function(expr, # expression matrix, rows as genes
                              group, # group values, factor or numeric
                              covar = NULL, # covariate data frame
                              num_threads = 20, # number of threads to use
                              return_coef_group = NULL) # the level of factor in the group to return the coefficient
{
  require(doParallel, quietly = T)
  idx <- which(!is.na(group))
  expr <- expr[,idx]
  
  group <- group[idx]
  if (is.factor(group) | is.character(group)){
    group <- factor(group)
    g_lev <- levels(group)
    if (!is.null(return_coef_group) & sum(g_lev == return_coef_group)>0){
      g_lev <- c(setdiff(g_lev, return_coef_group), return_coef_group)
    } else{
      return_coef_group <- NULL
    }
    group <- factor(group, levels = g_lev)
  }
  
  if(!is.null(covar)) covar <- setNames(data.frame(covar[idx,]), colnames(covar))
  
  registerDoParallel(num_threads)
  res <- foreach(i = 1:nrow(expr), .combine = rbind) %dopar%{
    e <- as.numeric(expr[i,])
    
    if (is.null(covar)){
      m0 <- lm(e ~ 1)
      m1 <- lm(e ~ group)
    } else{
      data <- data.frame(y = e, covar, group = group)
      m0 <- lm(y ~ . - group, data = data)
      m1 <- lm(y ~ ., data = data)
    }
    
    a0 <- anova(m0)
    a1 <- anova(m1)
    a01 <- anova(m0, m1)
    
    p1 <- a1$Pr[nrow(a1)-1]
    p2 <- a01$Pr[2]
    p3 <- pf(a0["Residuals","Mean Sq"] / a1["Residuals","Mean Sq"], df1 = a0["Residuals","Df"], df2 = a1["Residuals","Df"], lower.tail = F)
    
    var_tot <- sum(a0[,2])
    var_covar <- sum(a0[-nrow(a0),2])
    var_g <- a1[nrow(a1)-1,2]
    
    if (is.factor(group)){
      coef_g <- ifelse(is.null(return_coef_group), NA, coef(m1)[paste0("group",return_coef_group)])
    } else{
      coef_g <- coef(m1)["group"]
    }
    
    return(c(var_g, var_covar, var_tot, p1, p2, p3, coef_g))
  }
  stopImplicitCluster()
  res <- data.frame(res)
  rownames(res) <- rownames(expr)
  colnames(res) <- c("var_group","var_covar","var_total","p_ANOVA", "p_ANCOVA", "p_Resi", "coef")
  return(res)
}

if(! require(glmGamPoi, quietly = T)){
  BiocManager::install("glmGamPoi")
  library(glmGamPoi)
}
# Gamma poisson GLM based DE test
glm_gp_group_test <- function(expr, # expression matrix, rows as genes
                              group, # group values, factor or numeric
                              covar = NULL, # covariate data frame
                              is_log = T, # whether the expression matrix is log-transformed
                              force_dense = T, # always convert the expression matrix to dense matrix
                              tmpfile = tempfile(), # temp h5 file to store the HDF5 array when in need
                              subsample = T, # doing subsampling for the model to estimate dispersion
                              verbose = F, # verbose
                              ..., # parameters passed to glm_gp
                              return_coef_group = NULL) # the level of factor in the group to return the coefficient
{
  idx <- which(!is.na(group))
  expr <- expr[,idx]
  
  group <- group[idx]
  if (is.factor(group) | is.character(group)){
    group <- factor(group)
    g_lev <- levels(group)
    if (!is.null(return_coef_group) & sum(g_lev == return_coef_group)>0){
      g_lev <- c(setdiff(g_lev, return_coef_group), return_coef_group)
    } else{
      return_coef_group <- NULL
    }
    group <- factor(group, levels = g_lev)
  }
  
  if(!is.null(covar)) covar <- setNames(data.frame(covar[idx,]), colnames(covar))
  
  if (is(expr, "Matrix")){
    require(Matrix, quietly = T)
    if (is_log){
      df_expr <- summary(expr)
      df_expr$x <- exp(df_expr$x) - 1
      df_expr$x[df_expr$x < 0] <- 0
      expr <- drop0(sparseMatrix(i = df_expr$i, j = df_expr$j, x = df_expr$x, dims = dim(expr), dimnames = dimnames(expr)))
    }
    
    if ((! force_dense) & require(HDF5Array, quietly = T)){
      if (file.exists(tmpfile)) file.remove(tmpfile)
      expr <- HDF5Array::writeHDF5Array(expr, filepath = tmpfile, name = "expr", with.dimnames = T, as.sparse = F)
    } else{
      if (! force_dense)
        cat("HDF5Array package is required to use sparse matrix as input.\n")
      expr <- as.matrix(expr)
      cat("Converted the matrix to dense matrix.\n")
    }
  } else if (is_log){
    expr <- exp(expr) - 1
  }
  
  if (is.null(covar)){
    meta <- data.frame(group = group, covar = 1, row.names = colnames(expr))
  } else{
    meta <- data.frame(group = group, covar, row.names = colnames(expr))
  }
  
  fit <- try(glm_gp(expr, design = ~ ., col_data = meta, subsample = subsample, ..., verbose = verbose), silent = T)
  if (class(fit) == "try-error"){
    warnings(paste("Model fitting error: ", fit))
    res <- data.frame(name = rownames(expr), pval = NA, adj_pval = 1, f_statistic = NA, df1 = NA, df2 = NA, lfc = NA, coef = NA)
  } else{
    if (is.null(return_coef_group) & is.factor(group)){
      coef_g <- NA
    } else if (is.factor(group)){
      coef_g <- fit$Beta[,paste0("group",return_coef_group)]
    } else{
      coef_g <- fit$Beta[,"group"]
    }
    res <- test_de(fit, reduced_design = ~ . - group, verbose = verbose)
    res <- data.frame(res, coef = coef_g)
  }
  return(res)
}

# Binomial GLM ANOVA to test for detection rate differences
binom_group_test <- function(expr,
                             group,
                             coverage,
                             covar = NULL,
                             test = c("Chisq","F","Rao","Cp"),
                             num_threads = 50,
                             return_coef_group = NULL)
{
  test <- match.arg(test)
  require(doParallel, quietly = T)
  
  if (is.factor(group) | is.character(group)){
    group <- factor(group)
    g_lev <- levels(group)
    if (!is.null(return_coef_group) & sum(g_lev %in% return_coef_group)>0){
      g_lev <- c(setdiff(g_lev, return_coef_group), return_coef_group)
    } else{
      return_coef_group <- NULL
    }
    group <- factor(group, levels = g_lev)
  }
  
  meta <- data.frame(group = group, coverage = coverage, row.names = colnames(expr))
  if (!is.null(covar)){
    meta <- data.frame(meta, covar)
  } else{
    meta <- data.frame(meta, covar = 1)
  }
  
  registerDoParallel(num_threads)
  res <- foreach(i = 1:nrow(expr), .combine = rbind) %dopar%{
    e <- as.numeric(expr[i,]) > 0
    data <- data.frame(y = e, meta)
    m1 <- glm(y ~ ., data = data, family = binomial())
    m0 <- glm(y ~ . - group, data = data, family = binomial())
    p_anova <- anova(m1,m0, test = test)$Pr[2]
    
    a0 <- anova(m0)
    a1 <- anova(m1)
    p_resi <- pf((a0[nrow(a0),4]/a0[nrow(a0),3]) / (a1[nrow(a1),4]/a1[nrow(a1),3]), df1 = a0[nrow(a0),3], df2 = a1[nrow(a1),3], lower.tail = F)
    
    if (is.factor(group)){
      if (is.null(return_coef_group)){
        coef_names <- paste0("group",levels(group))
      } else{
        coef_names <- paste0("group",return_coef_group)
      }
      coef_g <- ifelse(coef_names %in% names(coef(m1)), coef(m1)[coef_names], 0)
    } else{
      coef_g <- coef(m1)["group"]
    }
    
    return(c(p_anova, p_resi, coef_g))
  }
  stopImplicitCluster()
  res <- data.frame(res)
  rownames(res) <- rownames(expr)
  if (is.numeric(group)){
    colnames(res) <- c("p_ANOVA", "p_Resi", "coef")
  } else{
    if (is.null(return_coef_group)){
      coef_names <- paste0("coef_",levels(group))
    } else{
      coef_names <- paste0("coef_",return_coef_group)
    }
    colnames(res) <- c("p_ANOVA", "p_Resi", coef_names)
  }
  return(res)
}