options(drop=FALSE)

run_consensus_clust2=function (norm.dat, select.cells = colnames(norm.dat), niter = 100, 
                               sample.frac = 0.8, output_dir = "subsample_result", mc.cores = 1, 
                               de.param = de_param(), merge.type = c("undirectional", "directional"), 
                               override = FALSE, init.result = NULL, cut.method = "auto", 
                               ...) 
{
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  all.cells = select.cells
  if (!is.null(init.result)) {
    all.cells = intersect(all.cells, names(init.result$cl))
  }
  run <- function(i, ...) {
    prefix = paste("iter", i, sep = ".")
    print(prefix)
    outfile = file.path(output_dir, paste0("result.", i, 
                                           ".rda"))
    if (file.exists(outfile) & !override) {
      return(NULL)
    }
    #result=NULL
    #while (is.null(result)) {
    select.cells = sample(all.cells, round(length(all.cells) * 
                                             sample.frac))
    save(select.cells, file = file.path(output_dir, paste0("cells.", 
                                                           i, ".rda")))
    result <- iter_clust(norm.dat = norm.dat, select.cells = select.cells, 
                         prefix = prefix, de.param = de.param, merge.type = merge.type, 
                         result = init.result, ...)
    #}
    if (is.null(result)) {
      result=list()
      result$cl=rep(1,length(select.cells))
      names(result$cl)=select.cells
    }
    save(result, file = outfile)
  }
  if (mc.cores == 1) {
    sapply(1:niter, function(i) {
      run(i, ...)
    })
  }
  else {
    require(foreach)
    require(doParallel)
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    foreach(i = 1:niter, .combine = "c") %dopar% run(i)
    stopCluster(cl)
  }
  result.files = file.path(output_dir, dir(output_dir, "result.*.rda"))
  load(result.files[[1]])
  cl.size = table(result$cl)
  graph.size = sum(cl.size^2)
  if (graph.size < 10^8) {
    co.result <- collect_co_matrix_sparseM(norm.dat, result.files, 
                                           all.cells)
    co.ratio = co.result$co.ratio
    consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, 
                                            norm.dat, select.cells = all.cells, de.param = de.param, 
                                            merge.type = merge.type, method = cut.method, result = init.result)
    refine.result = refine_cl(consensus.result$cl, co.ratio = co.ratio, 
                              tol.th = 0.01, confusion.th = 0.6, min.cells = de.param$min.cells)
    markers = consensus.result$markers
  }
  else {
    result <- iter_clust(norm.dat = norm.dat, select.cells = all.cells, 
                         de.param = de.param, merge.type = merge.type, result = init.result, 
                         ...)
    co.result <- collect_subsample_cl_matrix(norm.dat, result.files, 
                                             all.cells)
    cl = merge_cl_by_co(result$cl, co.ratio = NULL, cl.mat = co.result$cl.mat, 
                        diff.th = 0.25)
    refine.result = refine_cl(cl, cl.mat = co.result$cl.mat, 
                              tol.th = 0.01, confusion.th = 0.6, min.cells = de.param$min.cells)
    markers = result$markers
  }
  cl = refine.result$cl
  merge.result = merge_cl(norm.dat = norm.dat, cl = cl, rd.dat = Matrix::t(norm.dat[markers, 
                                                                                    ,drop=F]), de.param = de.param, merge.type = merge.type, return.markers = FALSE)
  return(list(co.result = co.result, cl.result = merge.result))
}

collect_co_matrix_sparseM=function (norm.dat, result.files, all.cells, max.cl.size = 100000) 
{
  tmp = collect_subsample_cl_matrix(norm.dat, result.files, 
                                    all.cells, max.cl.size = max.cl.size)
  cl.list = tmp$cl.list
  cl.mat = tmp$cl.mat
  co.ratio = Matrix::crossprod(Matrix::t(cl.mat))
  co.ratio@x = co.ratio@x/length(result.files)
  return(list(co.ratio = co.ratio, cl.mat = cl.mat, cl.list = cl.list))
}

collect_subsample_cl_matrix=function (norm.dat, result.files, all.cells, max.cl.size = NULL) 
{
  select.cells = c()
  cl.list = list()
  for (f in result.files) {
    print(f)
    tmp = load(f)
    cl = result$cl
    test.cells = setdiff(all.cells, names(cl))
    markers = unique(result$markers)
    map.df = map_by_cor(norm.dat[markers, names(cl),drop=F], cl, 
                        norm.dat[markers, test.cells,drop=F], method = "means")$pred.df
    test.cl = setNames(map.df$pred.cl, row.names(map.df))
    all.cl = c(setNames(as.character(cl), names(cl)), setNames(as.character(test.cl), 
                                                               names(test.cl)))
    cl.list[[f]] = all.cl
  }
  if (!is.null(max.cl.size)) {
    select.cells = sample_cl_list(cl.list, max.cl.size = max.cl.size)
  }
  else {
    select.cells = all.cells
  }
  cl.mat = do.call("cbind", sapply(cl.list, function(cl) {
    get_cl_mat(cl[select.cells])
  }, simplify = F))
  return(list(cl.list = cl.list, cl.mat = cl.mat))
}

DE_genes_pairs=function (norm.dat, cl, pairs, method = "limma", low.th = 1, 
                         cl.present = NULL, use.voom = FALSE, counts = NULL) 
{
  require(limma)
  select.cl = unique(c(pairs[, 1], pairs[, 2]))
  cl = cl[cl %in% select.cl]
  cl.means = as.data.frame(get_cl_means(norm.dat, cl))
  norm.dat = as.matrix(norm.dat[, names(cl),drop=F])
  if (length(low.th) == 1) {
    low.th = setNames(rep(low.th, nrow(norm.dat)), row.names(norm.dat))
  }
  if (is.null(cl.present)) {
    cl.present = as.data.frame(get_cl_means(norm.dat >= 
                                              low.th[row.names(norm.dat)], cl))
  }
  cl.size = table(cl)
  de.df = list()
  fit = NULL
  if (method == "limma") {
    cl = setNames(as.factor(paste0("cl", cl)), names(cl))
    design = model.matrix(~0 + cl)
    colnames(design) = levels(as.factor(cl))
    if (use.voom & !is.null(counts)) {
      v = voom(as.matrix(counts[row.names(norm.dat), names(cl),drop=F]), 
               design)
      fit = lmFit(v, design)
    }
    else {
      fit = lmFit(norm.dat[, names(cl),drop=F], design = design)
    }
  }
  for (i in 1:nrow(pairs)) {
    x = as.character(pairs[i, 1])
    y = as.character(pairs[i, 2])
    pair = paste(x, y, sep = "_")
    if (method == "limma") {
      ctr <<- paste(paste0("cl", x), "- ", paste0("cl", 
                                                  y))
      contrasts.matrix <- makeContrasts(ctr, levels = design)
      fit2 = contrasts.fit(fit, contrasts.matrix)
      fit2 = eBayes(fit2)
      pval = fit2$p.value[, 1]
      padj = p.adjust(pval)
      lfc = coef(fit2)[, 1]
    }
    else if (method == "chisq") {
      lfc = cl.means[, x] - cl.means[, y]
      df = vec_chisq_test(cl.present[, x] * cl.size[[x]], 
                          cl.size[[x]], cl.present[, y] * cl.size[[y]], 
                          cl.size[[y]])
      pval = df[, "pval"]
      padj = p.adjust(pval)
    }
    df = data.frame(padj = padj, pval = pval, lfc = lfc, 
                    meanA = cl.means[[x]], meanB = cl.means[[y]], q1 = cl.present[[x]], 
                    q2 = cl.present[[y]])
    row.names(df) = row.names(norm.dat)
    de.df[[pair]] = df
  }
  return(de.df)
}

DE_genes_cat_by_cl=function (norm.dat, cl, binary.cat, ...) 
{
  cl = droplevels(as.factor(cl))
  cl.cat = setNames(paste0(cl, binary.cat[names(cl)]), names(cl))
  tmp = levels(cl)
  cl.cat.pairs = data.frame(cat1 = paste0(tmp, binary.cat[1]), 
                            cat2 = paste(tmp, binary.cat[2]), stringsAsFactors = FALSE)
  cl.cat.de.genes = deScore.pairs(norm.dat[, names(cl.cat),drop=F], 
                                  cl = cl.cat, pairs = cl.cat.pairs, ...)
  cat1.de.num = sapply(cl.cat.de.genes, function(x) {
    if (length(x) == 0) 
      return(0)
    length(x$up.genes)
  })
  cat2.de.num = sapply(cl.cat.de.genes, function(x) {
    if (length(x) == 0) 
      return(0)
    length(x$down.genes)
  })
  cat1.de.genes = sapply(cl.cat.de.genes, function(x) {
    if (is.null(x)) 
      return("")
    return(paste(head(x$up.genes, 8), collapse = " "))
  })
  cat2.de.genes = sapply(cl.cat.de.genes, function(x) {
    if (is.null(x)) 
      return("")
    return(paste(head(x$down.genes, 8), collapse = " "))
  })
  cl.cat.de.df = data.frame(cat1.de.num, cat2.de.num, cat1.de.genes, 
                            cat2.de.genes)
  return(list(cl.cat.de.df = cl.cat.de.df, cl.cat.de.genes = cl.cat.de.genes))
}

DE_genes_pw=function (norm.dat, cl, ...) 
{
  cn = as.character(sort(unique(cl)))
  cl.n = length(cn)
  pairs = cbind(rep(cn, rep(cl.n, cl.n)), rep(cn, cl.n))
  pairs = pairs[pairs[, 1] < pairs[, 2], , drop = F]
  de.df = DE_genes_pairs(norm.dat = norm.dat, cl = cl, pairs = pairs, 
                         ...)
  return(de.df)
}

de_pair=function (df, de.param = de_param(), cl.size1 = NULL, cl.size2 = NULL) 
{
  df = df[order(df$pval, -abs(df$lfc)), ]
  select = with(df, which(padj < de.param$padj.th & abs(lfc) > 
                            de.param$lfc.th))
  select = row.names(df)[select]
  if (is.null(select) | length(select) == 0) {
    return(list())
  }
  up = select[df[select, "lfc"] > 0]
  down = select[df[select, "lfc"] < 0]
  df = df[select, ]
  if (!is.null(de.param$q.diff.th) & is.null(df$q.diff)) {
    df$q.diff = with(df, abs(q1 - q2)/pmax(q1, q2))
    df$q.diff[is.na(df$q.diff)] = 0
  }
  if (!is.null(de.param$q1.th)) {
    up = with(df[up, , drop = F], up[q1 > de.param$q1.th])
    if (!is.null(cl.size1)) {
      up = with(df[up, , drop = F], up[q1 * cl.size1 >= 
                                         de.param$min.cells])
    }
    down = with(df[down, , drop = F], down[q2 > de.param$q1.th])
    if (!is.null(cl.size2)) {
      down = with(df[down, , drop = F], down[q2 * cl.size2 >= 
                                               de.param$min.cells])
    }
  }
  if (!is.null(de.param$q2.th)) {
    up = with(df[up, , drop = F], up[q2 < de.param$q2.th])
    down = with(df[down, , drop = F], down[q1 < de.param$q2.th])
  }
  if (!is.null(de.param$q.diff.th)) {
    up = with(df[up, , drop = F], up[abs(q.diff) > de.param$q.diff.th])
    down = with(df[down, , drop = F], down[abs(q.diff) > 
                                             de.param$q.diff.th])
  }
  select = c(up, down)
  if (length(select) == 0) {
    return(list())
  }
  else {
    df$padj[df$padj < 10^-20] = 10^-20
    up.score = sum(-log10(df[up, "padj"]))
    down.score = sum(-log10(df[down, "padj"]))
    if (length(up) == 0) {
      up.score = 0
    }
    if (length(down) == 0) {
      down.score = 0
    }
    tmp = list(score = up.score + down.score, up.score = up.score, 
               down.score = down.score, num = length(select), up.num = length(up), 
               down.num = length(down), genes = select, up.genes = up, 
               down.genes = down, de.df = df[df$padj < de.param$padj.th, 
                                             ])
  }
}

de_param=function (low.th = 1, padj.th = 0.01, lfc.th = 1, q1.th = 0.5, 
          q2.th = NULL, q.diff.th = 0.7, de.score.th = 100, min.cells = 4) 
{
  list(low.th = low.th, padj.th = padj.th, lfc.th = lfc.th, 
       q1.th = q1.th, q2.th = q2.th, q.diff.th = q.diff.th, 
       de.score.th = de.score.th, min.cells = min.cells)
}



de_score=function (norm.dat, cl, de.param = de_param(), method = "limma", 
                   de.genes = NULL, ...) 
{
  if (is.factor(cl)) {
    cl = droplevels(cl)
  }
  cn = as.character(sort(unique(cl)))
  cl.n = length(cn)
  pairs = cbind(rep(cn, rep(cl.n, cl.n)), rep(cn, cl.n))
  pairs = pairs[pairs[, 1] < pairs[, 2], , drop = F]
  row.names(pairs) = paste(pairs[, 1], pairs[, 2], sep = "_")
  if (!is.null(de.genes)) {
    pairs = pairs[!row.names(pairs) %in% names(de.genes), 
                  , drop = F]
  }
  de.result = de_score_pairs(norm.dat, cl = cl, pairs = pairs, 
                             de.param = de.param, method = method, ...)
  de.genes = c(de.genes, de.result$de.genes)
  return(de.genes)
}

de_score_pairs=function (norm.dat, cl, pairs, de.df = NULL, de.param = de_param(), 
                         method = "limma") 
{
  row.names(pairs) = paste(pairs[, 1], pairs[, 2], sep = "_")
  select.cl = unique(c(pairs[, 1], pairs[, 2]))
  cl = cl[cl %in% select.cl]
  norm.dat = as.matrix(norm.dat[, names(cl),drop=F])
  if (is.factor(cl)) {
    cl = droplevels(cl)
  }
  cl.size = table(cl)
  cl.n = names(cl.size)
  cl.small = cl.n[cl.size < de.param$min.cells]
  cl.big = setdiff(cl.n, cl.small)
  select.pair = pairs[, 1] %in% cl.big & pairs[, 2] %in% cl.big
  de.genes = list()
  if (sum(select.pair) > 0) {
    cl = cl[cl %in% c(pairs[select.pair, 1], pairs[select.pair, 
                                                   2])]
    select.cells = names(cl)
    low.th = de.param$low.th
    if (length(low.th) == 1) {
      low.th = setNames(rep(low.th, nrow(norm.dat)), row.names(norm.dat))
    }
    select.genes = row.names(norm.dat)[rowSums(norm.dat[, 
                                                        select.cells,drop=F] >= low.th[row.names(norm.dat)]) >= 
                                         de.param$min.cells]
    if (is.null(de.df)) {
      de.df = DE_genes_pairs(norm.dat[select.genes, select.cells], 
                             cl[select.cells], pairs[select.pair, , drop = F], 
                             low.th = low.th, method = method)
    }
    de.genes = sapply(names(de.df), function(x) {
      if (is.null(de.df[[x]])) {
        return(list())
      }
      df = de.df[[x]]
      if (!is.null(de.param$min.cells)) {
        cl.size1 = cl.size[as.character(pairs[x, 1])]
        cl.size2 = cl.size[as.character(pairs[x, 2])]
      }
      else {
        cl.size1 = cl.size2 = NULL
      }
      de_pair(df, de.param = de.param, cl.size1, cl.size2)
    }, simplify = F)
  }
  else {
    de.df = list()
  }
  for (i in which(!select.pair)) {
    pair = paste(pairs[i, 1], pairs[i, 2], sep = "_")
    de.genes[[pair]] = list()
    de.df[[pair]] = list()
  }
  list(de.df = de.df, de.genes = de.genes)
}

dend_lca=function (dend, l1, l2, l = rep(attr(dend, "label"), length(l1))) 
{
  node.height = setNames(get_nodes_attr(dend, "height"), get_nodes_attr(dend, 
                                                                        "label"))
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      tmp.l = attr(dend[[i]], "label")
      cat(tmp.l, i, "\n")
      labels = get_subtree_label(dend[[i]])
      select = l1 %in% labels & l2 %in% labels
      if (sum(select) > 0) {
        select = which(select)[node.height[l[select]] > 
                                 node.height[tmp.l]]
        l[select] = tmp.l
        l = dend_lca(dend[[i]], l1, l2, l)
      }
    }
  }
  return(l)
}

dend_list=function (dend) 
{
  l = list()
  l[[attr(dend, "label")]] = dend
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      l = c(l, dend_list(dend[[i]]))
    }
  }
  return(l)
}

findVG=function (dat, plot_file = NULL, return_type = "data") 
{
  require(Matrix)
  require(ggplot2)
  if (!is.matrix(dat)) {
    sample_totals <- Matrix::colSums(dat)
  }
  else {
    sample_totals <- colSums(dat)
  }
  scaling_factor <- sample_totals/median(sample_totals)
  scaled_data <- dat/scaling_factor[col(dat)]
  if (is.matrix(dat)) {
    gene_means <- rowMeans(scaled_data)
    gene_vars <- rowMeans(scaled_data^2) - gene_means^2
  }
  else {
    gene_means <- Matrix::rowMeans(scaled_data)
    gene_vars <- Matrix::rowMeans(scaled_data^2) - scaled_data^2
  }
  dispersion <- log10(gene_vars/gene_means)
  IQR <- quantile(dispersion, probs = c(0.25, 0.75), na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1])/(qnorm(0.75) - qnorm(0.25))
  z <- (dispersion - m)/delta
  gene_var_stats <- data.frame(gene = rownames(dat), g.means = gene_means, 
                               g.vars = gene_vars, dispersion = dispersion, z = z, 
                               pval = 1 - pnorm(z), padj = p.adjust(1 - pnorm(z), method = "fdr"))
  select <- !is.na(gene_var_stats$dispersion) & gene_var_stats$dispersion > 
    0
  fit <- with(gene_var_stats, loess(dispersion ~ log10(g.means), 
                                    subset = select))
  residual <- resid(fit)
  base <- min(predict(fit))
  diff <- gene_var_stats$dispersion - base
  diff[select] <- residual
  IQR <- quantile(diff, probs = c(0.25, 0.75), na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1])/(qnorm(0.75) - qnorm(0.25))
  gene_var_stats$loess.z <- (diff - m)/delta
  gene_var_stats$loess.pval = 1 - pnorm(gene_var_stats$loess.z)
  gene_var_stats$loess.padj = p.adjust(gene_var_stats$loess.pval, 
                                       method = "fdr")
  if (!is.null(plot_file) | return_type %in% c("plots", "both")) {
    qq_z_plot <- ggplot(gene_var_stats, aes(sample = z)) + 
      stat_qq() + stat_qq_line()
    qq_loess.z_plot <- ggplot(gene_var_stats, aes(sample = loess.z)) + 
      stat_qq() + stat_qq_line()
    z_density_plot <- ggplot(gene_var_stats, aes(z)) + geom_density()
    loess.z_density_plot <- ggplot(gene_var_stats, aes(loess.z)) + 
      geom_density()
    fit.df <- data.frame(x = fit$x[, 1], y = fit$fitted)
    fit.df <- fit.df[order(fit.df$x), ]
    dispersion_fit_plot <- ggplot() + geom_point(data = gene_var_stats, 
                                                 aes(x = log10(g.means), y = dispersion)) + geom_line(data = fit.df, 
                                                                                                      aes(x = x, y = y), color = "blue")
  }
  if (!is.null(plot_file)) {
    pdf(plot_file)
    qq_z_plot
    qq_loess.z_plot
    z_density_plot
    loess.z_density_plot
    dispersion_fit_plot
    dev.off()
  }
  if (return_type == "data") {
    return(gene_var_stats)
  }
  else if (return_type == "plots") {
    return(list(qq_z_plot = qq_z_plot, qq_loess.z_plot = qq_loess.z_plot, 
                z_density_plot = z_density_plot, loess.z_density_plot = loess.z_density_plot, 
                dispersion_fit_plot = dispersion_fit_plot))
  }
  else if (return_type == "both") {
    return(list(gene_var_stats = gene_var_stats, plots = list(qq_z_plot = qq_z_plot, 
                                                              qq_loess.z_plot = qq_loess.z_plot, z_density_plot = z_density_plot, 
                                                              loess.z_density_plot = loess.z_density_plot, dispersion_fit_plot = dispersion_fit_plot)))
  }
}

get_cl_mat=function (cl) 
{
  library(Matrix)
  if (!is.factor(cl)) {
    cl <- as.factor(cl)
  }
  cl = droplevels(cl)
  cl.mat <- sparseMatrix(i = 1:length(cl), j = as.integer(cl), 
                         x = 1)
  rownames(cl.mat) <- names(cl)
  colnames(cl.mat) <- levels(cl)
  return(cl.mat)
}

get_cl_means=function (mat, cl) 
{
  library(Matrix)
  cl.sums <- get_cl_sums(mat, cl)
  cl.size <- table(cl)
  cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
  return(cl.means)
}

get_cl_sums=function (mat, cl) 
{
  library(Matrix)
  cl.mat <- get_cl_mat(cl)
  cl.sums <- Matrix::tcrossprod(mat[, rownames(cl.mat),drop=F], Matrix::t(cl.mat))
  cl.sums <- as.matrix(cl.sums)
  return(cl.sums)
}

get_cl_medians=function (mat, cl) 
{
  library(Matrix)
  library(matrixStats)
  cl.med <- do.call("cbind", tapply(names(cl), cl, function(x) {
    rowMedians(as.matrix(mat[, x]))
  }))
  rownames(cl.med) <- rownames(mat)
  return(cl.med)
}

get_cell.cl.co.ratio=function (cl, co.ratio = NULL, cl.mat = NULL) 
{
  require(Matrix)
  if (!is.null(co.ratio)) {
    cell.cl.co.ratio = get_cl_means(co.ratio, cl)
    return(cell.cl.co.ratio)
  }
  if (!is.null(cl.mat)) {
    tmp = cl.mat %*% get_cl_sums(Matrix::t(cl.mat), cl)
    tmp = tmp/Matrix::rowSums(cl.mat)
    cl.size = table(cl)
    cell.cl.co.ratio = as.matrix(t(t(tmp)/as.vector(cl.size[colnames(tmp)])))
    return(cell.cl.co.ratio)
  }
  stop("Either co.ratio or cl.mat should not be NULL")
}

filter_gene_mod=function (norm.dat, select.cells, gene.mod, minModuleSize = 10, 
                          min.deScore = 40, de.param = de_param(), max.cl.size = NULL, 
                          rm.eigen = NULL, rm.th = 0.7, maxSize = 200, prefix = "cl", 
                          max.mod = NULL) 
{
  require(matrixStats)
  eigen = get_eigen(gene.mod, norm.dat, select.cells)[[1]]
  if (!is.null(rm.eigen)) {
    rm.cor = cor(eigen, rm.eigen[select.cells, ])
    rm.cor[is.na(rm.cor)] = 0
    rm.score = setNames(rowMaxs(abs(rm.cor)), colnames(eigen))
    print(rm.score)
    select = rm.score < rm.th
    if (sum(!select)) {
      print("Remove module")
      print(rm.score[!select, drop = F])
    }
    eigen = eigen[, select, drop = F]
    gene.mod = gene.mod[select]
  }
  if (length(select.cells) > 4000) {
    method = "ward.D"
  }
  else {
    method = c("ward.D", "kmeans")
  }
  nmod = min(20, length(gene.mod))
  if (nmod == 0) {
    return(NULL)
  }
  gene.mod = head(gene.mod, nmod)
  eigen = eigen[, 1:nmod, drop = F]
  mod.score = setNames(rep(0, length(gene.mod)), names(gene.mod))
  not.selected = 1:length(gene.mod)
  for (m in method) {
    tmp = score_gene_mod(norm.dat, select.cells, gene.mod = gene.mod[not.selected], 
                         eigen = eigen[select.cells, not.selected, drop = F], 
                         method = m, de.param = de.param, max.cl.size = max.cl.size)
    x = do.call("cbind", sapply(tmp, function(x) x[[1]], 
                                simplify = F))
    tmp = x["sc", ] > mod.score[not.selected]
    mod.score[not.selected[tmp]] = x["sc", tmp]
    tmp = x["sc", ] > min.deScore
    not.selected = not.selected[!tmp]
  }
  ord = order(mod.score, decreasing = T)
  gene.mod = gene.mod[ord]
  mod.score = mod.score[ord]
  eigen = eigen[, ord, drop = F]
  select.mod = head(which(mod.score > min.deScore), max.mod)
  if (length(select.mod) == 0) {
    select.mod = head(which(mod.score > min.deScore/2), 
                      2)
  }
  if (length(select.mod) > 0) {
    gene.mod = gene.mod[select.mod]
    mod.score = mod.score[select.mod]
    eigen = eigen[, select.mod, drop = F]
    return(list(gene.mod = gene.mod, eigen = eigen, gene.mod.val = mod.score))
  }
  return(NULL)
}


get_cl_co_stats=function (cl, co.ratio = NULL, cl.mat = NULL) 
{
  require(matrixStats)
  cell.cl.co.ratio = get_cell.cl.co.ratio(cl, co.ratio = co.ratio, 
                                          cl.mat = cl.mat)
  cl.co.ratio <- get_cl_means(t(cell.cl.co.ratio), cl)
  cell.co.stats <- sapply(1:ncol(cell.cl.co.ratio), function(i) {
    select.cells = names(cl)[cl == colnames(cell.cl.co.ratio)[i]]
    cohesion = setNames(cell.cl.co.ratio[select.cells, i, 
                                         drop = F], select.cells)
    best.between = rowMaxs(cell.cl.co.ratio[select.cells, 
                                            -i, drop = F])
    confusion = best.between/cohesion
    separability = cohesion - best.between
    df = data.frame(cohesion, separability, confusion)
    colnames(df) = c("cohesion", "separability", "confusion")
    df
  }, simplify = F)
  cell.co.stats = do.call("rbind", cell.co.stats)
  cl.co.stats = as.data.frame(do.call("rbind", tapply(1:nrow(cell.co.stats), 
                                                      cl[row.names(cell.co.stats)], function(x) {
                                                        sapply(cell.co.stats[x, ], median)
                                                      })))
  return(list(cell.cl.co.ratio = cell.cl.co.ratio, cl.co.ratio = cl.co.ratio, 
              cell.co.stats = cell.co.stats, cl.co.stats = cl.co.stats))
}


get_color=function (color.mapping) 
{
  while (1) {
    tmp.color = which(duplicated(color.mapping))
    if (length(tmp.color) == 0) {
      break
    }
    rgb = col2rgb(color.mapping)
    for (x in tmp.color) {
      if (x < length(tmp.color)) {
        tmp = round(rgb[, x - 1] * 0.8 + sample(20, 
                                                3) + rgb[, x + 1] * 0.2)
      }
      else {
        tmp = round(rgb[, x - 1] * 0.8 + sample(40, 
                                                3))
      }
      rgb[, x] = tmp
    }
    rgb[rgb > 255] = 255
    color.mapping = rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
  }
  return(color.mapping)
}

get_cols=function (big.dat, cols) 
{
  p = big.dat@p
  if (is.character(cols)) {
    cols = match(cols, big.dat@col_id)
  }
  select = lapply(cols, function(col) {
    if (p[col] + 1 <= p[col + 1]) {
      (p[col] + 1):(p[col + 1])
    }
    else {
      NULL
    }
  })
  select.index = do.call("c", select)
  l = sapply(select, length)
  p = c(0, cumsum(l))
  i = (big.dat@i)[select.index]
  x = (big.dat@x)[select.index]
  mat = sparseMatrix(i = as.integer(i + 1), x = as.double(x), 
                     p = p, dims = c(big.dat@dim[1], length(l)))
  colnames(mat) = big.dat@col_id[cols]
  row.names(mat) = big.dat@row_id
  return(mat)
}

get_core_intermediate=function (norm.dat, cl, select.markers, n.bin = 5, n.iter = 100, 
                                mc.cores = 10, method = "median") 
{
  require(parallel)
  cl.cv <- mclapply(1:n.iter, function(i) {
    tmp = map_cv(norm.dat, cl, select.markers, n.bin = n.bin, 
                 method = method)
  }, mc.cores = mc.cores)
  cell.cl.cor.map = do.call("rbind", sapply(cl.cv, function(x) {
    df = data.frame(cell = names(x), cl = x)
  }, simplify = F))
  cell.cl.cor.map = table(cell.cl.cor.map[, 1], cell.cl.cor.map[, 
                                                                2])
  cell.cl.cor.map = cell.cl.cor.map/rowSums(cell.cl.cor.map)
  cell.cl.map.df = data.frame(org.cl = as.character(cl[row.names(cell.cl.cor.map)]), 
                              best.score = rowMaxs(cell.cl.cor.map), best.cl = colnames(cell.cl.cor.map)[apply(cell.cl.cor.map, 
                                                                                                               1, which.max)], stringsAsFactors = FALSE)
  row.names(cell.cl.map.df) = row.names(cell.cl.cor.map)
  tmp = get_pair_matrix_coor(cell.cl.cor.map, row.names(cell.cl.map.df), 
                             as.character(cell.cl.map.df$best.cl))
  tmp1 = cell.cl.cor.map
  tmp1[tmp] = 0
  cell.cl.map.df$second.score = rowMaxs(tmp1)
  cell.cl.map.df$second.cl = colnames(tmp1)[apply(tmp1, 1, 
                                                  which.max)]
  cell.cl.map.df$second.cl[cell.cl.map.df$second.score == 
                             0] = NA
  cell.cl.map.df$transition.cl = NA
  tmp = with(cell.cl.map.df, org.cl != best.cl | best.score < 
               0.9)
  cell.cl.map.df[tmp, "transition.cl"] = as.character(cell.cl.map.df[tmp, 
                                                                     "best.cl"])
  tmp = with(cell.cl.map.df, which(org.cl == transition.cl))
  cell.cl.map.df$transition.cl[tmp] = as.character(cell.cl.map.df[tmp, 
                                                                  "second.cl"])
  cl.med <- do.call("cbind", tapply(names(cl), cl, function(x) {
    rowMedians(as.matrix(norm.dat[select.markers, x]))
  }))
  row.names(cl.med) = select.markers
  cell.cl.cor = cor(as.matrix(norm.dat[select.markers, row.names(cell.cl.map.df),drop=F]), 
                    cl.med[select.markers, ,drop=F])
  cell.cl.map.df$cor = with(cell.cl.map.df, get_pair_matrix(cell.cl.cor, 
                                                            row.names(cell.cl.map.df), as.character(org.cl)))
  cell.cl.map.df$core = is.na(cell.cl.map.df$transition.cl)
  cl.transition.df = with(cell.cl.map.df, as.data.frame(table(org.cl, 
                                                              transition.cl)))
  colnames(cl.transition.df)[1:2] = c("cl1", "cl2")
  cl.transition.df = cl.transition.df[cl.transition.df$Freq > 
                                        0, ]
  cl.transition.df$org.cl = as.character(cl.transition.df$org.cl)
  cl.transition.df$transition.cl = as.character(cl.transition.df$transition.cl)
  cl.transition.df$cl.min = pmin(cl.transition.df$org.cl, 
                                 cl.transition.df$transition.cl)
  cl.transition.df$cl.max = pmax(cl.transition.df$org.cl, 
                                 cl.transition.df$transition.cl)
  cl.transition.df$cl.pair = paste(cl.transition.df$cl.min, 
                                   cl.transition.df$cl.max)
  trainsition.df.bi = do.call("rbind", tapply(1:nrow(cl.transition.df), 
                                              cl.transition.df$cl.pair, function(x) {
                                                tmp = cl.transition.df[x, ][1, ]
                                                tmp$Freq = sum(cl.transition.df[x, "Freq"])
                                                tmp[, c(4, 5, 3)]
                                              }))
  cl.size = table(cl)
  cl.transition.df.bi$cl.min.size = cl.size[cl.transition.df.bi$cl.min]
  cl.transition.df.bi$cl.max.size = cl.size[cl.transition.df.bi$cl.max]
  cl.transition.df.bi$ratio = with(cl.transition.df.bi, Freq/pmin(cl.min.size, 
                                                                  cl.max.size))
  cl.transition.df = cl.transition.df.bi
  return(list(cell.cl.map.df = cell.cl.map.df, cl.transition.df = cl.transition.df))
}

get_de_genes_sym=function (de.genes) 
{
  de.genes.sym = de.genes
  pairs = names(de.genes)
  pairs = as.data.frame(do.call("rbind", strsplit(pairs, "_")))
  row.names(pairs) = names(de.genes)
  pairs$tmp = paste0(pairs[, 2], "_", pairs[, 1])
  de.genes.sym[pairs$tmp] = de.genes[row.names(pairs)]
  for (i in 1:nrow(pairs)) {
    de.genes.sym[[pairs[i, "tmp"]]]$up.genes = de.genes[[row.names(pairs)[i]]]$down.genes
    de.genes.sym[[pairs[i, "tmp"]]]$down.genes = de.genes[[row.names(pairs)[i]]]$up.genes
  }
  return(de.genes.sym)
}

get_de_score=function (de.df, top.genes, upper = 20) 
{
  gene.score = -log10(de.df[top.genes, "padj"])
  gene.score[gene.score > upper] = upper
  gene.score = sum(gene.score)
  return(gene.score)
}

get_dend_markers=function (dend, n = 20) 
{
  if (length(dend) > 1) {
    m = head(names(sort(-attr(dend, "markers"))), n)
    markers = list()
    markers[[attr(dend, "label")]] = m
    for (i in 1:length(dend)) {
      markers = c(markers, get_dend_markers(dend[[i]]))
    }
    return(markers)
  }
  return(NULL)
}

get_dend_markers_direction=function (dend, n = 20) 
{
  if (length(dend) > 1) {
    m = unique(unlist(lapply(attr(dend, "markers.byCl"), 
                             head, n)))
    markers = list()
    markers[[attr(dend, "label")]] = m
    for (i in 1:length(dend)) {
      markers = c(markers, get_dend_markers(dend[[i]]))
    }
    return(markers)
  }
  return(NULL)
}

get_dend_parent=function (dend) 
{
  if (length(dend) > 1) {
    p = attr(dend, "label")
    c = sapply(1:length(dend), function(i) attr(dend[[i]], 
                                                "label"))
    edges = data.frame(parent = rep(p, length(c)), child = c)
    tmp = do.call("rbind", sapply(1:length(dend), function(i) get_dend_parent(dend[[i]]), 
                                  simplify = F))
    edges = rbind(edges, tmp)
    return(edges)
  }
  return(NULL)
}

get_eigen=function (gene.mod, norm.dat, select.cells = colnames(norm.dat), 
                    prefix = NULL, method = "ward.D", hc = NULL, ...) 
{
  tmp.dat = as.matrix(norm.dat[unlist(gene.mod), select.cells, 
                               drop = F])
  eigen = sapply(gene.mod, function(x) {
    tmp.num = sum(rowSums(tmp.dat[x, select.cells] > 0) > 
                    0)
    if (tmp.num < 3) {
      return(rep(0, length(select.cells)))
    }
    pr.result = prcomp(t(tmp.dat[x, select.cells,drop=F]), tol = 0.8)
    pc1 = pr.result$x[, 1]
    rot = pr.result$rotatio[, 1, drop = F]
    if (sum(rot > 0) < length(x)/2) {
      pc1 = -pc1
    }
    pc1
  })
  colnames(eigen) = paste0("ME", names(gene.mod))
  row.names(eigen) = select.cells
  kME = cor(t(as.matrix(norm.dat[unlist(gene.mod), select.cells,drop=F])), 
            eigen)
  hub = sapply(names(gene.mod), function(i) {
    x = gene.mod[[i]]
    x[which.max(kME[x, paste0("ME", i)])]
  })
  hub = setNames(hub, paste0("ME", names(gene.mod)))
  colnames(eigen) = paste(colnames(eigen), hub[colnames(eigen)])
  row.names(eigen) = select.cells
  if (!is.null(prefix) & ncol(eigen) > 1) {
    if (is.null(hc)) {
      hc = hclust(dist(eigen), method = method)
    }
    colv = as.dendrogram(hc)
    pdf(paste0(prefix, ".eigen.pdf"), height = 6, width = 7)
    heatmap.3(t(eigen), Colv = colv, col = blue.red(100), 
              trace = "none", dendrogram = "column", ...)
    dev.off()
  }
  return(list(eigen, hc))
}

get_gene_score=function (de.genes, top.n = 50, max.num = 1000, bin.th = 4) 
{
  select.genes <- sapply(de.genes, function(x) {
    up.genes = x$up.genes
    up.binary.genes = up.genes[x$q.stats[up.genes, 3] < 
                                 1 & x$q.stats[up.genes, 2] > bin.th]
    up.genes = up.genes[up.genes %in% c(head(up.genes, top.n), 
                                        up.binary.genes)]
    down.genes = x$down.genes
    down.binary.genes = down.genes[x$q.stats[down.genes, 
                                             1] < 1 & x$q.stats[down.genes, 4] > bin.th]
    down.genes = down.genes[down.genes %in% c(head(down.genes, 
                                                   top.n), down.binary.genes)]
    list(up = up.genes, down = down.genes)
  }, simplify = F)
  all.genes = unique(unlist(select.genes))
  up.gene.score = sapply(select.genes, function(x) {
    tmp = match(all.genes, x$up)
    tmp[is.na(tmp)] = max.num
    tmp
  })
  down.gene.score = sapply(select.genes, function(x) {
    tmp = match(all.genes, x$down)
    tmp[is.na(tmp)] = max.num
    tmp
  })
  row.names(up.gene.score) = row.names(down.gene.score) = all.genes
  return(list(up.gene.score = up.gene.score, down.gene.score = down.gene.score))
}

get_pair_matrix=function (m, rows, cols) 
{
  v <- as.vector(m)
  coor <- get_pair_matrix_coor(m, rows, cols)
  return(v[coor])
}

get_pair_matrix_coor=function (m, rows, cols) 
{
  if (!is.numeric(rows)) {
    rows <- match(rows, row.names(m))
  }
  if (!is.numeric(cols)) {
    cols <- match(cols, colnames(m))
  }
  coor <- (cols - 1) * nrow(m) + rows
  return(coor)
}

getGroupSpecificMarkers=function (cl.g, norm.dat, cl, cl.present.counts = NULL, low.th = 1, 
                                  q1.th = 0.5, q.diff.th = 0.7, n.markers = 5) 
{
  cl = droplevels(cl)
  select.cells = names(cl)[cl %in% cl.g]
  not.select.cells = setdiff(names(cl), select.cells)
  if (is.null(cl.present.counts)) {
    fg = rowSums(norm.dat[, select.cells,drop=F] > low.th)/length(select.cells)
    bg = rowSums(norm.dat[, not.select.cells,drop=F] > low.th)/length(not.select.cells)
  }
  else {
    fg = rowSums(cl.present.counts[, cl.g, drop = F])
    bg = rowSums(cl.present.counts[, levels(cl), drop = F]) - 
      fg
    bg.freq = bg/length(not.select.cells)
    fg.freq = fg/length(select.cells)
  }
  tau = (fg.freq - bg.freq)/pmax(bg.freq, fg.freq)
  g = names(fg.freq)[fg.freq > q1.th & tau > q.diff.th]
  g = g[order(tau[g] + fg.freq[g], decreasing = T)]
  g = head(g, n.markers)
  if (length(g > 0)) {
    df = data.frame(g = g, specifity = tau[g], fg.freq = fg.freq[g], 
                    bg.freq = bg.freq[g], fg.counts = fg[g], bg.counts = bg[g])
    return(df)
  }
  return(NULL)
}

getNodeSpecificMarkers=function (dend.list, norm.dat, cl, cl.df, ...) 
{
  do.call("rbind", sapply(names(dend.list), function(x) {
    print(x)
    cl.g = row.names(cl.df)[cl.df$cluster_label %in% labels(dend.list[[x]])]
    df = getGroupSpecificMarkers(cl.g, norm.dat, cl, ...)
    if (!is.null(df)) {
      df$cl = x
    }
    df
  }, simplify = F))
}

getNodeVsSiblingMarkers=function (dend.list, norm.dat, cl, cl.df, ...) 
{
  do.call("rbind", sapply(names(dend.list), function(x) {
    dend = dend.list[[x]]
    all.cl = droplevels(cl[cl %in% row.names(cl.df)[cl.df$cluster_label %in% 
                                                      labels(dend)]])
    if (length(dend) > 1) {
      do.call("rbind", sapply(1:length(dend), function(i) {
        cl.g = row.names(cl.df)[cl.df$cluster_label %in% 
                                  labels(dend[[i]])]
        df = getGroupSpecificMarkers(cl.g, norm.dat, 
                                     all.cl, ...)
        if (!is.null(df)) {
          df$cl = attr(dend[[i]], "label")
        }
        df
      }, simplify = F))
    }
    else {
      NULL
    }
  }, simplify = F))
}

getTopMarkersByPropNew=function (propExpr, medianExpr, propDiff = 0, propMin = 0.5, 
                                 medianFC = 1, excludeGenes = NULL, sortByMedian = TRUE) 
{
  specGenes <- rep("none", dim(propExpr)[2])
  names(specGenes) <- colnames(propExpr)
  propSort <- t(apply(propExpr, 1, function(x) return(-sort(-x))))
  propWhich <- t(apply(propExpr, 1, function(x, y) return(y[order(-x)]), 
                       colnames(propExpr)))
  medianDif <- apply(data.frame(propWhich[, 1], medianExpr), 
                     1, function(x, y) {
                       wIn <- y == as.character(x[1])
                       mIn <- as.numeric(x[2:length(x)])[wIn]
                       mOut <- max(as.numeric(x[2:length(x)])[!wIn])
                       return(mIn - mOut)
                     }, colnames(propExpr))
  keepProp <- (propSort[, 1] >= propMin) & ((propSort[, 1] - 
                                               propSort[, 2]) > propDiff) & (medianDif >= medianFC) & 
    (!is.element(rownames(propExpr), excludeGenes))
  propSort <- propSort[keepProp, ]
  propWhich <- propWhich[keepProp, ]
  ord <- order(-medianDif[keepProp] * ifelse(sortByMedian, 
                                             1, 0), propSort[, 2] - propSort[, 1])
  propSort <- propSort[ord, ]
  propWhich <- propWhich[ord, ]
  while (sum(keepProp) > 1) {
    keepProp <- !is.element(propWhich[, 1], names(specGenes)[specGenes != 
                                                               "none"])
    if (sum(keepProp) <= 1) 
      break
    tmp <- propWhich[keepProp, ]
    specGenes[tmp[1, 1]] <- rownames(tmp)[1]
  }
  specGenes
}

group_cl=function (anno, cluster, norm.dat, select.cells = colnames(norm.dat), 
                   neun.thresh = 0.5, neun.colname = "facs_population_plan", 
                   outlier.cl = NULL, donor.cl = NULL) 
{
  select.cells <- intersect(names(cluster), select.cells)
  select.id <- match(select.cells, colnames(norm.dat))
  norm.dat <- norm.dat[, select.id,drop=F]
  anno <- droplevels(anno[select.id, ])
  anno$cluster <- cluster[select.cells]
  neuronal.cl <- check_neun(anno, cluster, neun.thresh = neun.thresh, 
                            neun.colname = neun.colname)
  test.genes <- c("GAD1", "GAD2", "SLC17A7", "SLC17A6")
  expr.cnt <- apply(norm.dat[test.genes, ,drop=F], 1, function(x) {
    expr.thresh <- 3
    tapply(x, cluster, function(x) sum(x > expr.thresh)/length(x))
  })
  cl.class <- list()
  cl.class[["inh"]] <- setdiff(names(neuronal.cl)[neuronal.cl & 
                                                    (expr.cnt[, "GAD1"] > 0.5 | expr.cnt[, "GAD2"] > 0.5)], 
                               c(outlier.cl, donor.cl))
  cl.class[["exc"]] <- setdiff(names(neuronal.cl)[neuronal.cl & 
                                                    (expr.cnt[, "SLC17A7"] > 0.2 | expr.cnt[, "SLC17A6"] > 
                                                       0.2)], c(outlier.cl, donor.cl, cl.class[["inh"]]))
  cl.class[["glia"]] <- setdiff(names(neuronal.cl)[!neuronal.cl], 
                                c(outlier.cl, donor.cl))
  cl.class[["donor"]] <- setdiff(donor.cl, outlier.cl)
  cl.class[["outlier"]] <- setdiff(unique(cluster), unlist(cl.class[c("inh", 
                                                                      "exc", "glia", "donor")]))
  return(cl.class)
}

group_specific_markers=function (cl.g, norm.dat, cl, de.param, n.markers = 5, cl.present.counts = NULL) 
{
  cl = droplevels(cl)
  select.cells = names(cl)[cl %in% cl.g]
  not.select.cells = setdiff(names(cl), select.cells)
  if (is.null(cl.present.counts)) {
    fg = rowSums(norm.dat[, select.cells,drop=F] > de.param$low.th)
    bg = rowSums(norm.dat[, not.select.cells,drop=F] > de.param$low.th)
  }
  else {
    fg = rowSums(cl.present.counts[, cl.g, drop = F])
    bg = rowSums(cl.present.counts[, levels(cl), drop = F]) - 
      fg
  }
  bg.freq = bg/length(not.select.cells)
  fg.freq = fg/length(select.cells)
  tau = (fg.freq - bg.freq)/pmax(bg.freq, fg.freq)
  ratio = fg/(fg + bg)
  stats <- vec_chisq_test(fg, rep(length(select.cells), length(fg)), 
                          bg, rep(length(not.select.cells), length(bg)))
  g = names(fg.freq)[fg.freq > de.param$q1.th & tau > de.param$q.diff.th]
  g = g[order(tau[g] + ratio[g]/4 + fg.freq[g]/5, decreasing = T)]
  select.g = c(g[tau[g] > 0.95], head(g, n.markers))
  g = g[g %in% select.g]
  if (length(g > 0)) {
    df = data.frame(g = g, specificity = round(tau[g], digits = 2), 
                    fg.freq = round(fg.freq[g], digits = 2), bg.freq = round(bg.freq[g], 
                                                                             digits = 2), fg.counts = fg[g], bg.counts = bg[g], 
                    pval = stats[g, "pval"])
    return(df)
  }
  return(NULL)
}

heatmap.3=function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                       "row", "column", "none"), symm = FALSE, scale = c("none", 
                                                                                                                         "row", "column"), na.rm = TRUE, revC = identical(Colv, 
                                                                                                                                                                          "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
                      scale != "none", col = "heat.colors", colsep, rowsep, 
                    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
                    notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("column", 
                                                                                   "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
                    vline = median(breaks), linecol = tracecol, margins = c(5, 
                                                                            5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
                    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
                    key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                                                "density", "none"), denscol = tracecol, symkey = min(x < 
                                                                                                                       0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
                    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
                    ...) 
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) 
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", 
            "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv)) 
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) 
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv)) 
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) 
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) 
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) 
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  if (exists("hcc")) {
    retval$hcc <- hcc
  }
  if (exists("hcr")) {
    retval$hcr <- hcr
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
      1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) 
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid)) 
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors)) 
        stop("'ColSideColors' must be a character ")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                      1)
      if (is.vector(ColSideColors)) 
        nnn = 1
      else nnn = nrow(ColSideColors)
      lhei <- c(lhei[1], nnn * 0.1, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors)) 
        stop("'RowSideColors' must be a character ")
      if (is.vector(RowSideColors)) 
        nnn = 1
      else nnn = ncol(RowSideColors)
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                           1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], nnn * 0.1, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat)) 
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) 
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    if (is.vector(RowSideColors)) {
      image(rbind(1:nr), col = RowSideColors[rowInd], 
            axes = FALSE)
    }
    if (is.matrix(RowSideColors)) {
      jk.row = RowSideColors
      jk.xy = matrix(which(jk.row != "0"), dim(jk.row))
      colnames(jk.xy) <- colnames(jk.row)
      image(t(jk.xy), col = jk.row[rowInd, ], xaxt = "n", 
            yaxt = "n")
      axis(1, at = seq(0, 1, 1/(ncol(jk.xy) - 1)), labels = colnames(jk.xy), 
           las = 2, cex.axis = cexCol, tick = 0)
    }
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    if (is.vector(ColSideColors)) {
      image(cbind(1:nc), col = ColSideColors[colInd], 
            axes = FALSE)
    }
    if (is.matrix(ColSideColors)) {
      jk.col = ColSideColors
      jk.xy = matrix(which(jk.col != "0"), dim(jk.col))
      image(t(jk.xy), col = jk.col[, colInd], axes = FALSE)
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
    cellnote <- t(cellnote)
  }
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) 
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) 
    retval$rowDendrogram <- ddr
  if (exists("ddc")) 
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(xlab)) 
    mtext(xlab, side = 3, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  if (!missing(colsep)) 
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1], 
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep)) 
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, 
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) 
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x))
      tmpbreaks[length(tmpbreaks)] <- max(abs(x))
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row") 
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column") 
      mtext(side = 1, "Column Z-Score", line = 2)
    else {
      mtext(side = 1, "", line = 2)
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, 
           pretty(dens$y))
      title("")
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("", cex = 0.25)
      par(cex = 0.25)
      mtext(side = 2, "", line = 2)
    }
    else {
      title("", cex = 0.25)
    }
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

incrementHex=function (col, r = 0, g = 0, b = 0) 
{
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(pmin(1, 
                                                          pmax(0, x[1] + r/255)), pmin(1, pmax(0, x[2] + g/255)), 
                                                     pmin(1, pmax(0, x[3] + b/255))))
}

init_cut=function (co.ratio, select.cells, cl.list, min.cells = 4, th = 0.3, 
                   method = "ward.D", verbose = FALSE) 
{
  avg.cl.num = mean(sapply(cl.list, function(cl) {
    sum(table(cl[select.cells]) >= min.cells)
  }))
  tmp.dat = co.ratio[select.cells, select.cells]
  hc = hclust(as.dist(1 - as.matrix(crossprod(tmp.dat))), 
              method = "ward.D")
  tmp.cl = cutree(hc, ceiling(avg.cl.num) + 2)
  tmp.cl = refine_cl(tmp.cl, co.ratio = co.ratio, min.cells = min.cells, 
                     niter = 1, confusion.th = 1)$cl
  if (length(unique(tmp.cl)) == 1) {
    return(NULL)
  }
  tmp.cl = merge_cl_by_co(tmp.cl, tmp.dat, diff.th = th)
  if (length(unique(tmp.cl)) == 1) {
    return(NULL)
  }
  return(cl = tmp.cl)
}

iter_clust=function (norm.dat, select.cells = colnames(norm.dat), prefix = NULL, 
                     split.size = 10, result = NULL, method = "auto", ...) 
{
  if (!is.null(prefix)) {
    print(prefix)
  }
  if (method == "auto") {
    if (length(select.cells) > 3000) {
      select.method = "louvain"
    }
    else {
      select.method = "ward.D"
    }
  }
  else {
    select.method = method
  }
  if (length(select.cells) <= 3000) {
    if (!is.matrix(norm.dat)) {
      norm.dat = as.matrix(norm.dat[, select.cells,drop=F])
    }
  }
  if (is.null(result)) {
    result = onestep_clust(norm.dat, select.cells = select.cells, 
                           prefix = prefix, method = select.method, ...)
    gc()
  }
  if (!is.null(result)) {
    select.cells = intersect(select.cells, names(result$cl))
    cl = result$cl[select.cells]
    gene.mod = result$gene.mod
    markers = result$markers
    cl = setNames(as.integer(cl), names(cl))
    new.cl = cl
    cl.size = table(cl)
    to.split = names(cl.size)[cl.size >= split.size]
    if (length(to.split) > 0) {
      n.cl = 1
      for (x in sort(unique(cl))) {
        tmp.cells = names(cl)[cl == x]
        if (!x %in% to.split) {
          new.cl[tmp.cells] = n.cl
        }
        else {
          tmp.prefix = paste(prefix, x, sep = ".")
          tmp.result = iter_clust(norm.dat = norm.dat, 
                                  select.cells = tmp.cells, prefix = tmp.prefix, 
                                  split.size = split.size, method = method, 
                                  ...)
          gc()
          if (is.null(tmp.result)) {
            new.cl[tmp.cells] = n.cl
          }
          else {
            tmp.cl = tmp.result$cl
            if (length(unique(tmp.cl) > 1)) {
              new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
              markers = union(markers, tmp.result$markers)
            }
          }
        }
        n.cl = max(new.cl) + 1
      }
      cl = new.cl
    }
    result = list(cl = cl, markers = markers)
    return(result)
  }
  return(NULL)
}

iter_consensus_clust=function (co.ratio, cl.list, norm.dat, select.cells = colnames(co.ratio), 
                               de.param = de_param(), merge.type = c("undirectional", "directional"), 
                               all.col = NULL, diff.th = 0.25, prefix = NULL, method = c("auto", 
                                                                                         "louvain", "ward.D"), verbose = FALSE, result = NULL, 
                               rd.dat = NULL) 
{
  method = method[1]
  require(igraph)
  if (verbose) {
    print(prefix)
  }
  if (!is.null(result)) {
    markers = result$markers
    cl = setNames(as.integer(as.character(result$cl)), names(result$cl))
  }
  else {
    markers = NULL
    if (length(select.cells) < 2 * de.param$min.cells) {
      return(NULL)
    }
    if (method == "auto") {
      if (length(select.cells) > 3000) {
        select.method = "louvain"
      }
      else {
        if (!is.matrix(co.ratio)) {
          co.ratio = as.matrix(co.ratio[select.cells, 
                                        select.cells])
        }
        select.method = "ward.D"
      }
    }
    else {
      select.method = method
    }
    if (select.method == "ward.D") {
      tmp.cl = init_cut(co.ratio, select.cells, cl.list, 
                        min.cells = de.param$min.cells, th = diff.th, 
                        method = select.method)
      if (is.null(tmp.cl)) {
        return(NULL)
      }
    }
    else {
      adj.mat = co.ratio[select.cells, select.cells]
      gr = graph.adjacency(adj.mat, mode = "undirected", 
                           weighted = TRUE)
      comm = cluster_louvain(gr)
      rm(gr)
      gc()
      if (pass_louvain(modularity(comm), adj.mat)) {
        tmp.cl = setNames(comm$membership, select.cells)
      }
      else {
        return(NULL)
      }
    }
    if (length(unique(tmp.cl)) == 1) {
      return(NULL)
    }
    if (verbose) {
      print(table(tmp.cl))
    }
    if (is.null(rd.dat)) {
      rd.dat = get_cell.cl.co.ratio(tmp.cl, co.ratio)
    }
    tmp = merge_cl(norm.dat = norm.dat, cl = tmp.cl, rd.dat = rd.dat, 
                   verbose = verbose, de.param = de.param, merge.type = merge.type)
    if (is.null(tmp) | !is.list(tmp)) 
      return(NULL)
    if (length(unique(tmp$cl)) == 1) 
      return(NULL)
    tmp.cl = tmp$cl
    tmp.cl = merge_cl_by_co(tmp.cl, co.ratio, diff.th)
    tmp.cl = setNames(as.integer(tmp.cl), names(tmp.cl))
    if (length(unique(tmp.cl)) == 1) {
      return(NULL)
    }
    cl = tmp.cl
    markers = tmp$markers
    if (verbose) {
      print(table(cl))
    }
  }
  cell.cl.co.ratio = get_cl_means(co.ratio, cl)
  n.cl = max(cl)
  new.cl = cl
  for (i in sort(unique(cl))) {
    tmp.prefix = paste0(prefix, ".", i)
    tmp.cells = names(cl)[cl == i]
    uncertain.cells = sum(cell.cl.co.ratio[tmp.cells, as.character(i)] < 
                            1 - diff.th)
    if (uncertain.cells < de.param$min.cells) {
      next
    }
    result = iter_consensus_clust(co.ratio = co.ratio, cl.list = cl.list, 
                                  norm.dat = norm.dat, select.cells = tmp.cells, prefix = tmp.prefix, 
                                  all.col = all.col, diff.th = diff.th, method = method, 
                                  de.param = de.param, merge.type = merge.type, rd.dat = rd.dat, 
                                  verbose = verbose)
    if (is.null(result)) {
      next
    }
    tmp.cl = result$cl
    new.cl[names(tmp.cl)] = tmp.cl + n.cl
    n.cl = max(new.cl)
    markers = union(markers, result$markers)
  }
  cl = new.cl
  cl = setNames(as.integer(as.factor(cl)), names(cl))
  return(list(cl = cl, markers = markers))
}

jaccard=function (m) 
{
  library(Matrix)
  A <- tcrossprod(m)
  A <- as(A, "dgTMatrix")
  b <- Matrix::rowSums(m)
  x = A@x/(b[A@i + 1] + b[A@j + 1] - A@x)
  A@x = x
  return(A)
}

knn_jaccard=function (knn) 
{
  knn.df = data.frame(i = rep(1:nrow(knn), ncol(knn)), j = as.vector(knn))
  knn.mat = sparseMatrix(i = knn.df[[1]], j = knn.df[[2]], 
                         x = 1)
  jaccard(knn.mat)
}

label_dend=function (dend, n = 1) 
{
  if (is.null(attr(dend, "label"))) {
    attr(dend, "label") = paste0("n", n)
    n = n + 1
  }
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      tmp = label_dend(dend[[i]], n)
      dend[[i]] = tmp[[1]]
      n = tmp[[2]]
    }
  }
  return(list(dend = dend, n))
}

lm_matrix=function (dat, x) 
{
  coef <- colSums(t(dat) * x)/sum(x^2)
  fit <- matrix(coef, ncol = 1) %*% matrix(x, nrow = 1)
  resid <- dat - fit
  se_r <- rowVars(dat)
  se_f <- rowVars(resid)
  R_2 <- 1 - se_f/se_r
  R_2[is.na(R_2)] <- 0
  R_2 <- setNames(R_2, row.names(dat))
  n <- ncol(dat)
  F_stats <- (se_r - se_f)/(se_f/(n - 2))
  F_stats <- setNames(F_stats, row.names(dat))
  p.value <- df(F_stats, 1, n - 2)
  p.value[is.na(p.value)] <- 1
  padj <- p.adjust(p.value)
  fit <- data.frame(coef, R_2, F_stats, p.value, padj)
  return(list(fit = fit, resid = resid))
}

lm_normalize=function (dat, x, R_2.th = 0.2, padj.th = 0.01, min.genes = 5) 
{
  m <- rowMedians(dat)
  m <- setNames(m, row.names(dat))
  q.max <- rowMaxs(dat)
  q.max <- setNames(q.max, row.names(dat))
  dat <- dat - m
  x <- x - median(x)
  norm.result <- lm_matrix(dat, x)
  fit <- norm.result$fit
  resid <- norm.result$resid
  select.genes <- row.names(fit)[fit$padj < padj.th & fit$R_2 > 
                                   R_2.th]
  if (length(select.genes) >= min.genes) {
    dat[select.genes, ] <- resid[select.genes, ]
    dat <- dat + m
    dat[dat < 0] <- 0
    dat[select.genes, ] <- apply(dat[select.genes, ], 2, 
                                 function(x) {
                                   select <- x > q.max[select.genes]
                                   x[select] <- q.max[select]
                                   return(x)
                                 })
  }
  else {
    dat <- dat + m
  }
  return(list(dat, fit))
}


logcPM=function (counts) 
{
  norm.dat = cpm(counts)
  norm.dat@x = log2(norm.dat@x + 1)
  norm.dat
}

makeColorsUnique=function (colorVector, seed = 1) 
{
  set.seed(seed)
  while (max(table(colorVector)) > 1) {
    cl <- names(table(colorVector))[table(colorVector) > 
                                      1]
    for (i in which(is.element(colorVector, cl))) colorVector[i] <- as.character(incrementHex(colorVector[i], 
                                                                                              sample(-2:2, 1), sample(-2:2, 1), sample(-2:2, 1)))
  }
  colorVector
}

makeLayerLabel=function (Samp.dat, column = "roi_label") 
{
  val <- as.character(as.matrix(Samp.dat[, column]))
  val <- gsub("2_3", "3", val)
  val <- as.character(sapply(val, function(x) return(strsplit(x, 
                                                              "_")[[1]][2])))
  val <- substr(val, 2, nchar(val))
  val[val == "ZZ_Missing"] <- NA
  val <- gsub("1-6", "0", val)
  val <- as.numeric(substr(val, 1, 1))
  val[val == 0] <- NA
  val
}

makeRegionLabel=function (Samp.dat, column = "roi_label") 
{
  val <- as.character(as.matrix(Samp.dat[, column]))
  val[nchar(val) <= 4] <- paste0("LGN_", val[nchar(val) <= 
                                               4])
  val <- gsub("-", "_", val)
  val <- as.character(sapply(val, function(x) return(strsplit(x, 
                                                              "_")[[1]][1])))
  val
}

map_by_cor=function (train.dat, train.cl, test.dat, method = "median") 
{
  if (method == "median") {
    library(matrixStats)
    cl.meds <- tapply(names(train.cl), train.cl, function(x) {
      train.mat <- train.dat[, x, drop = F]
      train.mat <- as.matrix(train.mat)
      rowMedians(train.mat)
    })
    cl.dat <- do.call("cbind", cl.meds)
  }
  else {
    cl.dat <- get_cl_means(train.dat, train.cl)
  }
  row.names(cl.dat) <- row.names(train.dat)
  test.cl.cor <- cor(as.matrix(test.dat), cl.dat)
  test.cl.cor[is.na(test.cl.cor)] <- 0
  max.cl.cor <- apply(test.cl.cor, 1, which.max)
  pred.cl <- colnames(test.cl.cor)[max.cl.cor]
  pred.cl <- setNames(pred.cl, row.names(test.cl.cor))
  pred.score <- apply(test.cl.cor, 1, max)
  if (is.factor(train.cl)) {
    pred.cl <- setNames(factor(pred.cl, levels = levels(train.cl)), 
                        names(pred.cl))
  }
  pred.df <- data.frame(pred.cl = pred.cl, pred.score = pred.score)
  out_list <- list(pred.df = pred.df, cor.matrix = test.cl.cor)
  return(out_list)
}

map_cl_summary=function (ref.dat, ref.cl, map.dat, map.cl) 
{
  map.result <- map_by_cor(ref.dat, ref.cl, map.dat)
  cor.matrix <- map.result$cor.matrix
  map.df <- map.result$pred.df
  colnames(map.df)[1] <- "map.cl"
  map.df$org.cl <- map.cl[row.names(map.df)]
  cl.size <- table(map.cl)
  cl.map.df <- as.data.frame(with(map.df, table(org.cl, map.cl)))
  cl.map.df$Prob <- round(cl.map.df$Freq/cl.size[as.character(cl.map.df$org.cl)], 
                          digits = 2)
  cl.map.df$pred.score <- 0
  for (i in 1:nrow(cl.map.df)) {
    select <- names(map.cl)[map.cl == as.character(cl.map.df$org.cl[i])]
    cl.map.df$pred.score[i] <- mean(cor.matrix[select, as.character(cl.map.df$map.cl[i])])
  }
  cl.map.df$pred.score <- round(cl.map.df$pred.score, digits = 2)
  cl.map.df <- cl.map.df[cl.map.df$Freq > 0, ]
  out_list <- list(map.df = map.df, cl.map.df = cl.map.df)
  return(out_list)
}

map_cv=function (norm.dat, cl, markers, n.bin = 5, g.perc = 1, method = "median") 
{
  bins <- tapply(names(cl), cl, function(x) {
    if (length(x) > n.bin) {
      tmp <- rep_len(1:n.bin, length(x))
    }
    else {
      tmp <- sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  })
  bins <- unlist(bins)
  names(bins) <- gsub(".*\\.", "", names(bins))
  bins <- bins[names(cl)]
  pred.cl <- setNames(rep(NA, length(cl)), names(cl))
  for (i in 1:n.bin) {
    print(i)
    train.cells <- names(cl)[bins != i]
    test.cells <- names(cl)[bins == i]
    select.markers <- sample(markers, round(length(markers) * 
                                              g.perc))
    map.result <- map_by_cor(norm.dat[select.markers, ,drop=F], 
                             cl[train.cells], norm.dat[select.markers, test.cells,drop=F], 
                             method = method)$pred.df
    pred.cl[test.cells] <- as.character(map.result[test.cells, 
                                                   "pred.cl"])
  }
  return(pred.cl)
}

map_sampling=function (train.dat, train.cl, test.dat, markers, markers.perc = 0.8, 
                       iter = 100, method = "median", verbose = TRUE) 
{
  map.result <- sapply(1:iter, function(i) {
    if (verbose) {
      cat("\r", paste0("Running iteration ", i, " of ", 
                       iter, ".        "))
      flush.console()
    }
    tmp.markers <- sample(markers, round(length(markers) * 
                                           markers.perc))
    map_by_cor(train.dat[tmp.markers, ], train.cl, test.dat[tmp.markers, 
                                                            ], method = method)
  }, simplify = F)
  map.cl <- sapply(map.result, function(x) {
    x$pred.df$pred.cl
  })
  row.names(map.cl) <- colnames(test.dat)
  map <- as.data.frame(as.table(as.matrix(map.cl)))
  map.freq <- table(map$Var1, map$Freq)
  max.freq <- apply(map.freq, 1, which.max)
  pred.cl <- colnames(map.freq)[max.freq]
  pred.cl <- setNames(pred.cl, row.names(map.freq))
  map.df <- data.frame(pred.cl = pred.cl, prob = rowMaxs(map.freq)/iter)
  out_list <- list(map.df = map.df, map.freq = map.freq)
  return(out_list)
}

mapDendMarkers=function (dend.list, map.dat, select.cells, th = 0.5) 
{
  map.gene.num <- matrix(0, nrow = ncol(map.dat), ncol = length(dend.list))
  row.names(map.gene.num) = colnames(map.dat)
  colnames(map.gene.num) = names(dend.list)
  for (x in names(dend.list)) {
    node = dend.list[[x]]
    if (length(node) > 1) {
      for (i in 1:length(node)) {
        l = attr(node[[i]], "label")
        print(l)
        tmp.genes = intersect(attr(node, "markers.byCl")[[i]], 
                              row.names(map.dat))
        map.gene.num[select.cells, l] = colSums(map.dat[tmp.genes, 
                                                        select.cells] > th)
      }
    }
  }
  return(map.gene.num)
}

markers_max_tau=function (cl.dat, th = 0.5, tau.th = 0.7, n.markers = 3) 
{
  tau = calc_tau(cl.dat)
  tmp = split(row.names(cl.dat), colnames(cl.dat)[apply(cl.dat, 
                                                        1, which.max)])
  tau.df = do.call("rbind", sapply(names(tmp), function(x) {
    g = tmp[[x]]
    g = g[cl.dat[g, x] > th]
    g = g[order(tau[g], decreasing = T)]
    if (length(g) > 0) {
      df = data.frame(g = g, cl = x, cl.dat = cl.dat[g, 
                                                     x], tau = tau[g])
      df = df[df$tau > tau.th, ]
      head(df, n = n.markers)
    }
    else {
      NULL
    }
  }, simplify = F))
}

markers_tau_one_vs_other=function (cl, cl.present.counts, present.th = 0.4, tau.th = 0.8, 
                                   top.n = 10) 
{
  all.g = row.names(cl.present.counts)
  all.g = setdiff(all.g, c(grep("LOC", all.g, value = T), 
                           grep("Rik$", all.g, value = T)))
  cl.size = table(cl)[colnames(cl.present.counts)]
  cl.present.prob = t(t(cl.present.counts)/as.vector(cl.size))
  tau.genes = sapply(colnames(cl.present.counts), function(x) {
    fg = cl.present.counts[, x]/sum(cl == x)
    bg = (rowSums(cl.present.counts) - cl.present.counts[, 
                                                         x])/(length(cl) - sum(cl == x))
    tau = (fg - bg)/pmax(bg, fg)
    g = all.g[cl.present.prob[all.g, x] > 0.5]
    g = g[order(tau[g], decreasing = T)]
    g = g[tau[g] > 0.97 | (tau[g] > 0.7 & g %in% head(g, 
                                                      top.n))]
  }, simplify = F)
}

merge_cl=function (norm.dat, cl, rd.dat, de.param = de_param(), merge.type = c("undirectional", 
                                                                               "directional"), max.cl.size = 300, de.method = "limma", 
                   de.genes = NULL, return.markers = TRUE, verbose = 0) 
{
  merge.type = merge.type[1]
  print(merge.type)
  cl = setNames(as.integer(as.character(cl)), names(cl))
  de.df = list()
  select.cells = names(cl)
  pairs = NULL
  if (!is.null(de.genes)) {
    pairs = do.call("rbind", strsplit(names(de.genes), "_"))
    row.names(pairs) = names(de.genes)
  }
  if (is.matrix(norm.dat)) {
    cell.gene.counts = colSums(norm.dat[, select.cells,drop=F] > 
                                 0)
  }
  else {
    cell.gene.counts = Matrix::colSums(norm.dat[, select.cells,drop=F] > 
                                         0)
  }
  cell.weights = cell.gene.counts - min(cell.gene.counts) + 
    200
  while (length(unique(cl)) > 1) {
    cl.size = table(cl)
    tmp.cells = sample_cells(cl, max.cl.size, weights = cell.weights)
    tmp.dat = as.matrix(norm.dat[, tmp.cells,drop=F])
    cl.rd = Matrix::t(get_cl_means(Matrix::t(rd.dat), cl))
    while (TRUE) {
      cl.size = table(cl)
      if (ncol(cl.rd) > 2) {
        cl.sim = cor(t(cl.rd))
      }
      else {
        cl.diff = as.matrix(dist(cl.rd/as.vector(cl.size[row.names(cl.rd)])))
        cl.sim = 1 - cl.diff/max(cl.diff)
      }
      cl.small = names(cl.size)[cl.size < de.param$min.cells]
      if (length(cl.small) > 0) {
        tmp = as.data.frame(as.table(cl.sim[cl.small, 
                                            , drop = F]))
        tmp[, 1] = as.integer(as.character(tmp[, 1]))
        tmp[, 2] = as.integer(as.character(tmp[, 2]))
        tmp = tmp[tmp[, 1] != tmp[, 2], , drop = F]
        closest.pair = which.max(tmp$Freq)
        x = tmp[closest.pair, 1]
        y = tmp[closest.pair, 2]
        if (verbose > 0) {
          cat("Merge: ", x, y, "sim:", tmp[closest.pair, 
                                           3], "\n")
        }
        cl[cl == x] = y
        cl.rd[as.character(y), ] = cl.rd[as.character(y), 
                                         ] + cl.rd[as.character(x), ]
        cl.rd = cl.rd[row.names(cl.rd) != x, , drop = F]
      }
      else {
        break
      }
    }
    if (length(cl.size) == 1) {
      return(NULL)
    }
    if (length(cl.size) < 10 & length(de.genes) == 0) {
      de.genes = de_score(tmp.dat, cl = cl[tmp.cells], 
                          de.param = de.param, method = de.method, de.genes = de.genes)
      sc = sapply(de.genes, function(x) {
        if (length(x) > 0) {
          x$score
        }
        else {
          0
        }
      })
      merge.pairs = do.call("rbind", strsplit(names(de.genes), 
                                              "_"))
      row.names(merge.pairs) = names(de.genes)
      pairs = merge.pairs
    }
    else {
      if (nrow(cl.rd) > 2) {
        k = 3
        knn.matrix = t(sapply(1:nrow(cl.sim), function(i) {
          colnames(cl.sim)[-i][head(order(cl.sim[i, 
                                                 -i], decreasing = T), k)]
        }))
        row.names(knn.matrix) = row.names(cl.sim)
        merge.pairs = do.call("rbind", apply(knn.matrix, 
                                             2, function(x) data.frame(c1 = row.names(knn.matrix), 
                                                                       c2 = x)))
        merge.pairs[, 1] = as.integer(as.character(merge.pairs[, 
                                                               1]))
        merge.pairs[, 2] = as.integer(as.character(merge.pairs[, 
                                                               2]))
        tmp1 = pmin(merge.pairs[, 1], merge.pairs[, 
                                                  2])
        tmp2 = pmax(merge.pairs[, 1], merge.pairs[, 
                                                  2])
        merge.pairs = cbind(tmp1, tmp2)
      }
      else {
        merge.pairs = matrix(as.integer(row.names(cl.rd)), 
                             nrow = 1)
      }
      row.names(merge.pairs) = paste(merge.pairs[, 1], 
                                     merge.pairs[, 2], sep = "_")
      merge.pairs = merge.pairs[!duplicated(row.names(merge.pairs)), 
                                , drop = F]
      if (nrow(merge.pairs) == 0) {
        break
      }
      new.pairs = setdiff(row.names(merge.pairs), names(de.genes))
      pairs = rbind(pairs, merge.pairs[new.pairs, , drop = F])
      tmp.de.genes = de_score_pairs(tmp.dat, cl = cl[tmp.cells], 
                                    pairs = merge.pairs[new.pairs, , drop = F], 
                                    de.param = de.param, method = de.method)$de.genes
      de.genes[names(tmp.de.genes)] = tmp.de.genes
      sc = sapply(de.genes[row.names(merge.pairs)], function(x) {
        if (length(x) > 0) {
          x$score
        }
        else {
          0
        }
      })
    }
    sc = sort(sc)
    to.merge = sapply(names(sc), function(p) {
      x = de.genes[[p]]
      if (length(x) == 0) {
        to.merge = TRUE
      }
      else if (merge.type == "undirectional") {
        to.merge = x$score < de.param$de.score.th
      }
      else {
        to.merge = x$up.score < de.param$de.score.th | 
          x$down.score < de.param$de.score.th
      }
      to.merge
    })
    if (sum(to.merge) == 0) {
      break
    }
    sc = sc[to.merge]
    to.merge = merge.pairs[names(sc), , drop = FALSE]
    merged = c()
    for (i in 1:nrow(to.merge)) {
      p = to.merge[i, ]
      if (i == 1 | sc[i] < de.param$de.score.th/2 & length(intersect(p, 
                                                                     merged)) == 0) {
        if (verbose > 0) {
          cat("Merge ", p[1], p[2], sc[i], "\n")
        }
        cl[cl == p[2]] = p[1]
        rm.pairs = row.names(pairs)[pairs[, 1] %in% 
                                      p | pairs[, 2] %in% p]
        de.genes = de.genes[setdiff(names(de.genes), 
                                    rm.pairs)]
      }
      merged = c(merged, p)
    }
    pairs = pairs[names(de.genes), , drop = F]
  }
  if (length(unique(cl)) < 2) {
    return(NULL)
  }
  if (verbose > 0) {
    print(table(cl))
  }
  markers = NULL
  if (return.markers) {
    de.genes = de_score(as.matrix(norm.dat[, tmp.cells,drop=F]), 
                        cl[tmp.cells], de.genes = de.genes, de.param = de.param)
    markers = select_markers(norm.dat, cl, de.genes = de.genes, 
                             n.markers = 50)$markers
  }
  sc = sapply(de.genes, function(x) {
    if (length(x) > 0) {
      x$score
    }
    else {
      0
    }
  })
  return(list(cl = cl, de.genes = de.genes, sc = sc, markers = markers))
}

merge_cl_by_co=function (cl, co.ratio = NULL, cl.mat = NULL, diff.th = 0.25, 
                         verbose = 0) 
{
  cell.cl.co.ratio = get_cell.cl.co.ratio(cl, co.ratio = co.ratio, 
                                          cl.mat = cl.mat)
  cl.co.ratio <- do.call("rbind", tapply(names(cl), cl, function(x) colMeans(cell.cl.co.ratio[x, 
                                                                                              , drop = F])))
  co.within = diag(cl.co.ratio)
  co.df <- as.data.frame(as.table(cl.co.ratio), stringsAsFactors = FALSE)
  co.df = co.df[co.df[, 1] < co.df[, 2] & co.df[, 3] > 0.1, 
                ]
  co.df$within1 = co.within[co.df[, 1]]
  co.df$within2 = co.within[co.df[, 2]]
  co.df$diff = pmax(co.df$within1, co.df$within2) - co.df[, 
                                                          3]
  co.df = co.df[co.df$diff < diff.th, ]
  co.df = co.df[order(co.df[, 1], decreasing = T), ]
  if (verbose > 0) {
    print(co.df)
  }
  for (i in 1:nrow(co.df)) {
    cl[cl == co.df[i, 2]] = co.df[i, 1]
  }
  cl = setNames(as.integer(as.character(cl)), names(cl))
  return(cl)
}

merge_cl_fast=function (norm.dat, cl, rd.dat.t, de.param = de_param(), merge.type = c("undirectional", 
                                                                                      "directional"), max.cl.size = 300, de.method = "limma", 
                        de.genes = NULL, return.markers = FALSE, pairBatch = 40, 
                        sampled = FALSE, verbose = 0) 
{
  if (!is.integer(cl)) {
    cl = setNames(as.integer(as.character(cl)), names(cl))
  }
  merge.type = merge.type[1]
  de.df = list()
  pairs = NULL
  if (!is.null(de.genes)) {
    pairs = do.call("rbind", strsplit(names(de.genes), "_"))
    row.names(pairs) = names(de.genes)
  }
  cl.rd = get_cl_means(rd.dat.t, cl[names(cl) %in% colnames(rd.dat.t)])
  while (TRUE) {
    cl.size = table(cl)
    if (length(cl.size) == 1) {
      break
    }
    cl.small = names(cl.size)[cl.size < de.param$min.cells]
    if (length(cl.small) == 0) {
      break
    }
    if (ncol(cl.rd) > 2 & nrow(cl.rd) > 2) {
      cl.sim = cor(cl.rd)
    }
    else {
      cl.diff = as.matrix(dist(t(cl.rd)))
      cl.sim = 1 - cl.diff/max(cl.diff)
    }
    tmp = as.data.frame(as.table(cl.sim[cl.small, , drop = F]))
    tmp[, 1] = as.integer(as.character(tmp[, 1]))
    tmp[, 2] = as.integer(as.character(tmp[, 2]))
    tmp = tmp[tmp[, 1] != tmp[, 2], , drop = F]
    closest.pair = which.max(tmp$Freq)
    x = tmp[closest.pair, 1]
    y = tmp[closest.pair, 2]
    if (verbose > 0) {
      cat("Merge: ", x, y, "sim:", tmp[closest.pair, 3], 
          "\n")
    }
    cl[cl == x] = y
    tmp = rowMeans(rd.dat.t[, names(cl)[cl == y], drop = FALSE])
    cl.rd[, as.character(y)] = tmp
    cl.rd = cl.rd[, colnames(cl.rd) != x, drop = F]
  }
  while (length(unique(cl)) > 1) {
    if (!is.null(max.cl.size) & !sampled) {
      sampled.cells = sample_cells(cl, max.cl.size)
      tmp.dat = norm.dat[, sampled.cells,drop=F]
    }
    else {
      tmp.dat = norm.dat
    }
    if (length(unique(cl)) == 2) {
      merge.pairs = as.data.frame(matrix(as.integer(colnames(cl.rd)), 
                                         nrow = 1))
      merge.pairs$sim = cor(cl.rd[, 1], cl.rd[, 2])
      row.names(merge.pairs) = paste(merge.pairs[, 1], 
                                     merge.pairs[, 2], sep = "_")
    }
    else {
      if (ncol(cl.rd) > 2 & nrow(cl.rd) > 2) {
        cl.sim = cor(cl.rd)
      }
      else {
        cl.diff = as.matrix(dist(t(cl.rd)))
        cl.sim = 1 - cl.diff/max(cl.diff)
      }
      knn.matrix = t(sapply(1:nrow(cl.sim), function(i) {
        colnames(cl.sim)[-i][order(cl.sim[i, -i], decreasing = T)]
      }))
      row.names(knn.matrix) = row.names(cl.sim)
      knn.matrix = knn.matrix[, 1:min(3, ncol(knn.matrix))]
      merge.pairs = do.call("rbind", apply(knn.matrix, 
                                           2, function(x) data.frame(c1 = row.names(knn.matrix), 
                                                                     c2 = x, stringsAsFactors = FALSE)))
      merge.pairs$sim = get_pair_matrix(cl.sim, merge.pairs[, 
                                                            1], merge.pairs[, 2])
      merge.pairs[, 1] = as.integer(merge.pairs[, 1])
      merge.pairs[, 2] = as.integer(merge.pairs[, 2])
      tmp1 = pmin(merge.pairs[, 1], merge.pairs[, 2])
      tmp2 = pmax(merge.pairs[, 1], merge.pairs[, 2])
      merge.pairs[, 1:2] = cbind(tmp1, tmp2)
      p = paste(merge.pairs[, 1], merge.pairs[, 2], sep = "_")
      merge.pairs = merge.pairs[!duplicated(p), , drop = F]
      row.names(merge.pairs) = p[!duplicated(p)]
      merge.pairs = merge.pairs[order(merge.pairs$sim, 
                                      decreasing = T), ]
    }
    if (nrow(merge.pairs) == 0) {
      break
    }
    new.pairs = setdiff(row.names(merge.pairs), names(de.genes))
    while (length(new.pairs) > 0) {
      new.pairs = new.pairs[head(order(merge.pairs[new.pairs, 
                                                   "sim"], decreasing = T), pairBatch)]
      pairs = rbind(pairs, merge.pairs[new.pairs, , drop = F])
      tmp.de.genes = de_score_pairs(tmp.dat, cl = cl[colnames(tmp.dat)], 
                                    pairs = merge.pairs[new.pairs, , drop = F], 
                                    de.param = de.param, method = de.method)$de.genes
      de.genes[names(tmp.de.genes)] = tmp.de.genes
      tmp.pairs = intersect(names(de.genes), row.names(merge.pairs))
      sc = sapply(de.genes[tmp.pairs], function(x) {
        if (length(x) > 0) {
          x$score
        }
        else {
          0
        }
      })
      sc = sort(sc)
      to.merge = sapply(names(sc), function(p) {
        x = de.genes[[p]]
        if (length(x) == 0) {
          to.merge = TRUE
        }
        else if (merge.type == "undirectional") {
          to.merge = x$score < de.param$de.score.th
        }
        else {
          to.merge = x$up.score < de.param$de.score.th | 
            x$down.score < de.param$de.score.th
        }
        to.merge
      })
      if (sum(to.merge) > 0) {
        sc = sc[to.merge]
        to.merge = merge.pairs[names(sc), , drop = FALSE]
        to.merge$sc = sc
        break
      }
      new.pairs = setdiff(row.names(merge.pairs), names(de.genes))
    }
    if (length(new.pairs) == 0) {
      break
    }
    merged = c()
    for (i in 1:nrow(to.merge)) {
      p = c(to.merge[i, 1], to.merge[i, 2])
      if (i == 1 | sc[i] < de.param$de.score.th/2 & length(intersect(p, 
                                                                     merged)) == 0) {
        if (verbose > 0) {
          cat("Merge ", p[1], p[2], to.merge[i, "sc"], 
              to.merge[i, "sim"], "\n")
        }
        cl[cl == p[2]] = p[1]
        rm.pairs = row.names(pairs)[pairs[, 1] %in% 
                                      p | pairs[, 2] %in% p]
        de.genes = de.genes[setdiff(names(de.genes), 
                                    rm.pairs)]
        tmp.cells = intersect(names(cl)[cl == p[1]], 
                              colnames(rd.dat.t))
        tmp = Matrix::rowMeans(rd.dat.t[, tmp.cells, 
                                        drop = FALSE])
        cl.rd[, as.character(p[1])] = tmp
        cl.rd = cl.rd[, colnames(cl.rd) != p[2], drop = F]
      }
      merged = c(merged, p)
    }
    pairs = pairs[names(de.genes), , drop = F]
  }
  if (length(unique(cl)) < 2) {
    return(NULL)
  }
  if (verbose > 0) {
    print(table(cl))
  }
  markers = NULL
  if (return.markers) {
    de.genes = de_score(tmp.dat, cl[colnames(tmp.dat)], 
                        de.genes = de.genes, de.param = de.param)
  }
  markers = select_markers(norm.dat, cl, de.genes = de.genes, 
                           n.markers = 50)$markers
  sc = sapply(de.genes, function(x) {
    if (length(x) > 0) {
      x$score
    }
    else {
      0
    }
  })
  return(list(cl = cl, de.genes = de.genes, sc = sc, markers = markers))
}

merge_co_matrix=function (co.ratio1, co.ratio2) 
{
  all.cells = c(row.names(co.ratio1), row.names(co.ratio2))
  tmp.co = Matrix(0, nrow = ncol(co.ratio1), ncol = ncol(co.ratio2))
  colnames(tmp.co) = colnames(co.ratio2)
  row.names(tmp.co) = row.names(co.ratio1)
  tmp.co1 = rbind(co.ratio1, t(tmp.co))
  tmp.co2 = rbind(tmp.co, co.ratio2)
  co.ratio = cbind(tmp.co1, tmp.co2)
  return(co.ratio)
}

multiplot=function (co.ratio1, co.ratio2) 
{
  all.cells = c(row.names(co.ratio1), row.names(co.ratio2))
  tmp.co = Matrix(0, nrow = ncol(co.ratio1), ncol = ncol(co.ratio2))
  colnames(tmp.co) = colnames(co.ratio2)
  row.names(tmp.co) = row.names(co.ratio1)
  tmp.co1 = rbind(co.ratio1, t(tmp.co))
  tmp.co2 = rbind(tmp.co, co.ratio2)
  co.ratio = cbind(tmp.co1, tmp.co2)
  return(co.ratio)
}

node_specific_markers=function (dend.list, norm.dat, cl, ...) 
{
  do.call("rbind", sapply(names(dend.list), function(x) {
    print(x)
    df = group_specific_markers(labels(dend.list[[x]]), 
                                norm.dat, cl, ...)
    if (!is.null(df)) {
      df$cl = x
    }
    df
  }, simplify = F))
}

node_vs_sibling_markers=function (dend.list, norm.dat, cl, cl.df, ...) 
{
  do.call("rbind", sapply(names(dend.list), function(x) {
    dend = dend.list[[x]]
    print(labels(dend))
    all.cl = droplevels(cl[cl %in% labels(dend)])
    if (length(dend) > 1) {
      do.call("rbind", sapply(1:length(dend), function(i) {
        cl.g = labels(dend[[i]])
        df = group_specific_markers(cl.g, norm.dat, 
                                    all.cl, ...)
        if (!is.null(df)) {
          df$node = attr(dend[[i]], "label")
          df$parent = attr(dend, "label")
        }
        df
      }, simplify = F))
    }
    else {
      NULL
    }
  }, simplify = F))
}

onestep_clust=function (norm.dat, select.cells = colnames(norm.dat), counts = NULL, 
                        method = c("louvain", "ward.D", "kmeans"), vg.padj.th = 0.5, 
                        dim.method = c("pca", "WGCNA"), max.dim = 20, rm.eigen = NULL, 
                        rm.th = 0.7, de.param = de_param(), min.genes = 5, merge.type = c("undirectional", 
                                                                                          "directional"), maxGenes = 3000, sampleSize = 4000, 
                        max.cl.size = 300, prefix = NULL, verbose = FALSE) 
{
  library(matrixStats)
  method = method[1]
  dim.method = dim.method[1]
  merge.type = merge.type[1]
  if (length(select.cells) > sampleSize) {
    sampled.cells = sample(select.cells, pmin(length(select.cells), 
                                              sampleSize))
  }
  else {
    sampled.cells = select.cells
  }
  if (is.matrix(norm.dat)) {
    select.genes = row.names(norm.dat)[which(rowSums(norm.dat[, 
                                                              select.cells,drop=F] > de.param$low.th) >= de.param$min.cells)]
  }
  else {
    select.genes = row.names(norm.dat)[which(Matrix::rowSums(norm.dat[, 
                                                                      select.cells,drop=F] > de.param$low.th) >= de.param$min.cells)]
  }
  if (is.null(counts)) {
    counts = 2^(norm.dat[select.genes, sampled.cells,drop=F]) - 
      1
  }
  plot_file = NULL
  if (verbose & !is.null(prefix)) {
    plot_file = paste0(prefix, ".vg.pdf")
  }
  vg = findVG(as.matrix(counts[select.genes, sampled.cells]), 
              plot_file = plot_file)
  if (dim.method == "auto") {
    if (length(select.cells) > 1000) {
      dim.method = "pca"
    }
    else {
      dim.method = "WGCNA"
    }
  }
  if (dim.method == "WGCNA") {
    select.genes = row.names(vg)[which(vg$loess.padj < 1)]
    select.genes = head(select.genes[order(vg[select.genes, 
                                              "padj"], -vg[select.genes, "z"])], maxGenes)
    rd.dat = rd_WGCNA(norm.dat, select.genes = select.genes, 
                      select.cells = select.cells, sampled.cells = sampled.cells, 
                      de.param = de.param, max.mod = max.dim, max.cl.size = max.cl.size)
  }
  else {
    select.genes = row.names(vg)[which(vg$loess.padj < vg.padj.th | 
                                         vg$dispersion > 3)]
    select.genes = head(select.genes[order(vg[select.genes, 
                                              "padj"], -vg[select.genes, "z"])], maxGenes)
    if (verbose) {
      cat("Num high variance genes:", length(select.genes), 
          "\n")
    }
    if (length(select.genes) < min.genes) {
      return(NULL)
    }
    rd.dat = rd_PCA(norm.dat, select.genes, select.cells, 
                    sampled.cells = sampled.cells, max.pca = max.dim)
  }
  if (is.null(rd.dat) || ncol(rd.dat) == 0) {
    return(NULL)
  }
  if (!is.null(rm.eigen)) {
    rm.cor = cor(rd.dat, rm.eigen[row.names(rd.dat), ])
    rm.cor[is.na(rm.cor)] = 0
    rm.score = rowMaxs(abs(rm.cor))
    select = rm.score < rm.th
    if (sum(!select) > 0 & verbose) {
      print("Remove dimension:")
      print(rm.score[!select])
    }
    if (sum(select) == 0) {
      return(NULL)
    }
    rd.dat = rd.dat[, select, drop = F]
  }
  if (verbose) {
    print(method)
  }
  max.cl = ncol(rd.dat) * 2 + 1
  if (method == "louvain") {
    k = pmin(15, round(nrow(rd.dat)/2))
    tmp = type(rd.dat, k)
    if (is.null(tmp)) {
      return(NULL)
    }
    cl = tmp$cl
    if (length(unique(cl)) > max.cl) {
      tmp.means = do.call("cbind", tapply(names(cl), cl, 
                                          function(x) {
                                            colMeans(rd.dat[x, , drop = F])
                                          }, simplify = F))
      tmp.hc = hclust(dist(t(tmp.means)), method = "average")
      tmp.cl = cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
      cl = setNames(tmp.cl[as.character(cl)], names(cl))
    }
  }
  else if (method == "ward.D") {
    hc = hclust(dist(rd.dat), method = "ward.D")
    cl = cutree(hc, max.cl)
  }
  else if (method == "kmeans") {
    cl = kmeans(rd.dat, max.cl)$cluster
  }
  else {
    stop(paste("Unknown clustering method", method))
  }
  merge.result = merge_cl(norm.dat, cl = cl, rd.dat = rd.dat, 
                          merge.type = merge.type, de.param = de.param, max.cl.size = max.cl.size)
  gc()
  if (is.null(merge.result)) 
    return(NULL)
  sc = merge.result$sc
  cl = merge.result$cl
  if (length(unique(cl)) > 1) {
    if (verbose) {
      cat("Expand", prefix, "\n")
      cl.size = table(cl)
      print(cl.size)
      save(cl, file = paste0(prefix, ".cl.rda"))
    }
    de.genes = merge.result$de.genes
    markers = merge.result$markers
    cl.dat = get_cl_means(norm.dat[markers, ,drop=F], cl[sample_cells(cl, 
                                                               max.cl.size)])
    cl.hc = hclust(dist(t(cl.dat)), method = "average")
    cl = setNames(factor(as.character(cl), levels = colnames(cl.dat)[cl.hc$order]), 
                  names(cl))
    if (verbose & !is.null(prefix)) {
      display_cl(norm.dat, cl, prefix = prefix, markers = markers, 
                 max.cl.size = max.cl.size)
    }
    levels(cl) = 1:length(levels(cl))
    result = list(cl = cl, markers = markers)
    return(result)
  }
  return(NULL)
}

pass_louvain=function (mod.sc, adj.mat) 
{
  library(Matrix)
  p <- mean(Matrix::colSums(adj.mat > 0) - 1)/nrow(adj.mat)
  n <- ncol(adj.mat)
  rand.mod1 <- 0.97 * sqrt((1 - p)/(p * n))
  rand.mod2 <- (1 - 2/sqrt(n)) * (2/(p * n))^(2/3)
  rand.mod.max <- max(rand.mod1, rand.mod2, na.rm = TRUE)
  return(mod.sc > rand.mod.max)
}

plot_cell_cl_co_matrix=function (co.ratio, cl, max.cl.size = 100, col = NULL) 
{
  blue.red <- colorRampPalette(c("blue", "white", "red"))
  select.cells = sample_cells(cl, max.cl.size)
  co.stats = get_cl_co_stats(cl, co.ratio)
  mat = co.stats$cell.cl.co.ratio
  tom = Matrix::tcrossprod(mat[select.cells, ])
  row.names(tom) = colnames(tom) = select.cells
  all.hc = hclust(as.dist(1 - tom), method = "average")
  ord1 = all.hc$labels[all.hc$order]
  ord = ord1[order(cl[ord1])]
  sep = cl[ord]
  sep = which(sep[-1] != sep[-length(sep)])
  if (is.null(col)) {
    heatmap.3(mat[ord, ], col = blue.red(150)[50:150], trace = "none", 
              Rowv = NULL, Colv = NULL, rowsep = sep, sepcolor = "black", 
              dendrogram = "none", labRow = "")
  }
  else {
    heatmap.3(mat[ord, ], col = blue.red(150)[50:150], trace = "none", 
              Rowv = NULL, Colv = NULL, rowsep = sep, sepcolor = "black", 
              ColSideColors = col[, ord], dendogram = "none", 
              labRow = "")
  }
}

plot_cl_heatmap=function (norm.dat, cl, markers, prefix = NULL, hc = NULL, gene.hc = NULL, 
                          centered = FALSE, labels = names(cl), sorted = FALSE, by.cl = TRUE, 
                          ColSideColors = NULL, maxValue = 5, min.sep = 4, main = "", 
                          height = 13, width = 9) 
{
  library(matrixStats)
  blue.red <- colorRampPalette(c("blue", "white", "red"))
  select.cells = names(cl)
  tmp.dat = as.matrix(norm.dat[markers, select.cells, drop = F])
  if (!is.null(ColSideColors)) {
    ColSideColors = ColSideColors[, select.cells, drop = F]
  }
  if (centered) {
    tmp.dat = tmp.dat - rowMeans(tmp.dat)
    breaks = c(min(min(tmp.dat) - 0.1, -maxValue), seq(-maxValue, 
                                                       maxValue, length.out = 99), max(max(tmp.dat) + 1))
  }
  else {
    tmp.dat = tmp.dat/pmax(rowMaxs(tmp.dat), 1)
    breaks = c(0, seq(0.05, 1, length.out = 100))
  }
  colnames(tmp.dat) = labels
  cexCol = min(70/ncol(tmp.dat), 1)
  cexRow = min(60/nrow(tmp.dat), 1)
  if (is.null(gene.hc)) {
    gene.hc = hclust(dist(tmp.dat), method = "ward.D")
  }
  if (is.null(hc) & !sorted & length(select.cells) < 2000) {
    hc = hclust(dist(t(tmp.dat)), method = "ward.D")
  }
  col = blue.red(150)[51:150]
  if (!is.null(prefix)) {
    pdf(paste(prefix, "pdf", sep = "."), height = height, 
        width = width)
  }
  if (by.cl) {
    if (sorted) {
      ord = 1:length(cl)
    }
    else {
      if (!is.null(hc)) {
        ord = order(cl, order(hc$order))
      }
      else {
        ord = order(cl)
      }
    }
    sep = cl[ord]
    sep = which(sep[-1] != sep[-length(sep)])
    sep = c(sep[1], sep[which(sep[-1] - sep[-length(sep)] >= 
                                min.sep) + 1])
    heatmap.3(tmp.dat[, ord], Rowv = as.dendrogram(gene.hc), 
              Colv = NULL, col = col, trace = "none", dendrogram = "none", 
              cexCol = cexCol, cexRow = cexRow, ColSideColors = ColSideColors[, 
                                                                              ord], breaks = breaks, colsep = sep, sepcolor = "black", 
              main = main)
    cells.order = colnames(tmp.dat)[ord]
  }
  else {
    heatmap.3(tmp.dat, Rowv = as.dendrogram(gene.hc), Colv = as.dendrogram(hc), 
              col = col, trace = "none", dendrogram = "none", 
              cexCol = cexCol, cexRow = cexRow, ColSideColors = ColSideColors, 
              breaks = breaks, main = main)
    cells.order = colnames(tmp.dat)[hc$order]
  }
  if (!is.null(prefix)) {
    dev.off()
  }
  return(cells.order)
}

plot_cl_meta_barplot=function (cluster, meta, col = NULL) 
{
  library(ggplot2)
  meta = as.factor(meta)
  final.tbl <- table(cluster, meta)
  final.tbl = final.tbl/rowSums(final.tbl)
  tb.df = droplevels(as.data.frame(final.tbl))
  g = ggplot(data = tb.df, aes(x = cluster, y = Freq, fill = meta)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, 
                                                                   hjust = 1, vjust = 0.5), panel.grid.major = element_blank(), 
                                        panel.background = element_blank())
  if (!is.null(col)) {
    g = g + scale_fill_manual(values = col)
  }
  return(g)
}

plot_co_matrix=function (co.ratio, cl, max.cl.size = 100, col = NULL) 
{
  blue.red <- colorRampPalette(c("blue", "white", "red"))
  select.cells = names(cl)
  select.cells = sample_cells(cl, max.cl.size)
  tom = Matrix::crossprod(co.ratio[select.cells, select.cells])
  row.names(tom) = colnames(tom) = select.cells
  all.hc = hclust(as.dist(1 - tom), method = "average")
  ord1 = all.hc$labels[all.hc$order]
  ord1 = ord1[ord1 %in% select.cells]
  ord = ord1[order(cl[ord1])]
  sep = cl[ord]
  sep = which(sep[-1] != sep[-length(sep)])
  if (is.null(col)) {
    heatmap.3(as.matrix(co.ratio[ord, ord]), col = blue.red(150)[50:150], 
              trace = "none", Rowv = NULL, Colv = NULL, colsep = sep, 
              sepcolor = "black", labRow = "")
  }
  else {
    heatmap.3(as.matrix(co.ratio[ord, ord]), col = blue.red(150)[50:150], 
              trace = "none", Rowv = NULL, Colv = NULL, colsep = sep, 
              sepcolor = "black", ColSideColors = col[, ord], 
              labRow = "")
  }
}

plot_de_lfc_num=function (de.genes, top.n = 100, select.pair = NULL, cl.label = NULL) 
{
  de.score <- sapply(de.genes, function(x) {
    x$score
  })
  de.up.score <- sapply(de.genes, function(x) {
    x$up.score
  })
  de.down.score <- sapply(de.genes, function(x) {
    x$down.score
  })
  de.num <- sapply(de.genes, function(x) {
    x$num
  })
  de.lfc = sapply(de.genes, function(x) {
    top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])], 
                     top.n)
    mean(abs(x$de.df[top.genes, "lfc"]))
  })
  de.q.diff = sapply(de.genes, function(x) {
    top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])], 
                     top.n)
    q.diff = with(x$de.df[top.genes, ], abs(q1 - q2)/pmax(q1, 
                                                          q2))
    mean(q.diff)
  })
  de.summary = data.frame(de.num, de.lfc, de.q.diff, de.score, 
                          de.up.score, de.down.score)
  row.names(de.summary) = names(de.genes)
  tmp = do.call("rbind", strsplit(row.names(de.summary), "_"))
  de.summary$cl1 = tmp[, 1]
  de.summary$cl2 = tmp[, 2]
  g = ggplot(de.summary, aes(de.num, de.lfc, color = de.q.diff)) + 
    geom_point() + scale_color_gradient2(midpoint = 0.85) + 
    scale_x_log10()
  if (!is.null(select.pair)) {
    select.df = de.summary[select.pair, ]
    if (!is.null(cl.label)) {
      select.df$pair.label = with(select.df, paste(cl.label[as.character(cl1)], 
                                                   cl.label[as.character(cl2)], sep = ":"))
    }
    g = g + geom_text(data = select.df, aes(de.num - 0.02, 
                                            de.lfc, label = pair.label), size = 2, color = "black")
    g = g + geom_point(data = select.df, aes(de.num, de.lfc), 
                       color = "red", pch = 1)
  }
  g = g + xlab("Number of DE genes")
  g = g + ylab("Mean log2(FC) of top 100 DE.genes")
  return(list(g = g, de.summary = de.summary))
}

plot_de_num=function (de.genes, dend, cl.label = NULL, directed = FALSE, 
                      file = "log10.de.num.pdf", ...) 
{
  up.de.num = sapply(de.genes, function(x) length(x$up.genes))
  down.de.num = sapply(de.genes, function(x) length(x$down.genes))
  head(sort(pmin(up.de.num, down.de.num)))
  names(down.de.num) = with(pairs[names(down.de.num), ], paste0(P2, 
                                                                "_", P1))
  label = as.hclust(dend)$label
  de.num = c(up.de.num, down.de.num)
  de.num.matrix <- convert_pair_matrix(de.num, directed = directed)
  de.num.matrix <- de.num.matrix[label, label]
  save(de.num.matrix, file = "de.num.matrix.rda")
  breaks = c(-1, seq(0.2, 4, length.out = 100))
  if (!is.null(cl.label)) {
    colnames(de.num.matrix) = row.names(de.num.matrix) = cl.label[row.names(de.num.matrix)]
  }
  tmp.dat = log10(de.num.matrix + 1)
  pdf(file, ...)
  heatmap.3(tmp.dat, col = jet.colors(100), breaks = breaks, 
            trace = "none", Colv = dend, Rowv = dend, dendrogram = "row", 
            cexRow = 0.3, cexCol = 0.3)
  dev.off()
}

plot_dend=function (dend, dendro_data = NULL, node_size = 1, r = c(-0.1, 
                                                                   1)) 
{
  require(dendextend)
  require(ggplot2)
  if (is.null(dendro_data)) {
    dendro_data = as.ggdend(dend)
    dendro_data$nodes$label = get_nodes_attr(dend, "label")
    dendro_data$nodes = dendro_data$nodes[is.na(dendro_data$nodes$leaf), 
                                          ]
  }
  node_data = dendro_data$nodes
  label_data <- dendro_data$labels
  segment_data <- dendro_data$segments
  if (is.null(node_data$node_color)) {
    node_data$node_color = "black"
  }
  ggplot() + geom_text(data = node_data, aes(x = x, y = y, 
                                             label = label, color = node_color), size = node_size, 
                       vjust = 1) + geom_segment(data = segment_data, aes(x = x, 
                                                                          xend = xend, y = y, yend = yend), color = "gray50") + 
    geom_text(data = label_data, aes(x = x, y = -0.01, label = label, 
                                     color = col), size = node_size, angle = 90, hjust = 1) + 
    scale_color_identity() + theme_dendro() + scale_y_continuous(limits = r)
}

plot_river=function (anno.df1, anno.df2, common.cells, ...) 
{
  source("~/zizhen/My_R/map_river_plot.R")
  source("~/zizhen/My_R/sankey_functions.R")
  cols = c("sample_id", "cluster_id", "cluster_label", "cluster_color")
  df1 = anno.df1 %>% filter(sample_id %in% common.cells) %>% 
    select(cols)
  df2 = anno.df2 %>% filter(sample_id %in% common.cells) %>% 
    select(cols)
  colnames(df2)[-1] = paste0("map_", colnames(df2)[-1])
  df = left_join(df1, df2)
  g = river_plot(df, ...)
}

plot_tsne_cl=function (norm.dat, select.genes, cl, cl.df, tsne.df = NULL, 
                       show.legend = FALSE, cex = 0.15, fn.size = 2, alpha.val = 1, 
                       ...) 
{
  library(ggplot2)
  require(Rtsne)
  if (is.null(tsne.df)) {
    tsne.result = Rtsne(t(as.matrix(norm.dat[select.genes, 
                                             names(cl),drop=F])), ...)$Y
    row.names(tsne.result) = names(cl)
    tsne.df = as.data.frame(tsne.result[names(cl), ])
    colnames(tsne.df) = c("Lim1", "Lim2")
  }
  tsne.df$cl = cl[row.names(tsne.df)]
  tsne.df$cl_label = factor(cl.df[as.character(tsne.df$cl), 
                                  "cluster_label"], levels = as.character(cl.df$cluster_label))
  tsne.df$cl_label = droplevels(tsne.df$cl_label)
  cl.center = do.call("rbind", tapply(1:nrow(tsne.df), tsne.df$cl, 
                                      function(x) {
                                        x = sample(x, pmin(length(x), 500))
                                        center = c(median(tsne.df[x, 1]), median(tsne.df[x, 
                                                                                         2]))
                                        dist = as.matrix(dist(tsne.df[x, 1:2]))
                                        tmp = x[which.min(rowSums(dist))]
                                        c(x = tsne.df[tmp, 1], y = tsne.df[tmp, 2])
                                      }))
  row.names(cl.center) = cl.df[row.names(cl.center), "cluster_label"]
  cl.col = setNames(as.character(cl.df$cluster_color), cl.df$cluster_label)
  shape = setNames(1:length(levels(tsne.df$cl_label))%%20 + 
                     1, levels(tsne.df$cl_label))
  g = ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color = cl_label, 
                                                        shape = cl_label), size = cex)
  g = g + scale_color_manual(values = alpha(as.vector(cl.col[levels(tsne.df$cl_label)]), 
                                            alpha.val)) + scale_shape_manual(values = as.vector(shape[levels(tsne.df$cl_label)]))
  for (i in 1:nrow(cl.center)) {
    g = g + annotate("text", label = row.names(cl.center)[i], 
                     x = cl.center[i, 1], y = cl.center[i, 2], size = fn.size, 
                     color = "black")
  }
  g = g + geom_point(data = as.data.frame(cl.center), aes(x = x, 
                                                          y = y), size = cex * 1.5)
  g = g + theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
                axis.line.y = element_line(colour = "black"))
  if (show.legend) {
    g = g + guides(colour = guide_legend(override.aes = list(shape = shape[levels(tsne.df$cl_label)])), 
                   ncol = 5)
    g = g + theme(legend.position = "bottom")
  }
  else {
    g = g + theme(legend.position = "none")
  }
  return(list(tsne.df = tsne.df, g = g))
}

plot_tsne_gene=function (tsne.df, norm.dat, genes, cex = 0.15) 
{
  library(ggplot2)
  plots = list()
  for (g in genes) {
    tsne.df$expr = norm.dat[g, row.names(tsne.df),drop=F]
    p = ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color = expr), 
                                                      size = cex)
    p = p + scale_color_gradient(low = "gray", high = "red") + 
      xlab(g)
    p = p + theme(legend.position = "none")
    plots[[g]] = p
  }
  return(plots)
}

plot_tsne_meta=function (tsne.df, meta, meta.col = NULL, show.legend = TRUE, 
                         cex = 0.15, legend.size = 5) 
{
  library(ggplot2)
  tsne.df$meta = meta
  p = ggplot(tsne.df, aes(Lim1, Lim2)) + geom_point(aes(color = meta), 
                                                    size = cex)
  if (is.factor(meta)) {
    if (is.null(meta.col)) {
      meta.col = setNames(jet.colors(length(levels(meta))), 
                          levels(meta))
    }
    p = p + scale_color_manual(values = as.vector(meta.col[levels(tsne.df$meta)]))
    p = p + theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
                  axis.line.y = element_line(colour = "black"))
  }
  else {
    p = p + scale_color_gradient(low = "blue", high = "red")
  }
  if (!show.legend) {
    p = p + theme(legend.position = "none")
  }
  else {
    p = p + guides(colour = guide_legend(override.aes = list(size = legend.size)))
  }
  return(p)
}


predict_annotate_cor=function (cl, norm.dat, ref.markers, ref.cl, ref.cl.df, ref.norm.dat, 
                               method = "median", reorder = FALSE) 
{
  map_results <- map_by_cor(ref.norm.dat[ref.markers, ,drop=F], ref.cl, 
                            norm.dat[ref.markers, names(cl),drop=F], method = method)
  pred.cl <- setNames(factor(as.character(map_result$pred.df$pred.cl), 
                             levels = row.names(ref.cl.df)), row.names(map_results$pred.df))
  results <- compare_annotate(cl, pred.cl, ref.cl.df, reorder = reorder)
  return(results)
}

prune_dend=function (dend, rm.labels, top.level = TRUE) 
{
  if (length(dend) > 1) {
    new_dend = list()
    for (i in 1:length(dend)) {
      new_dend[[i]] = prune_dend(dend[[i]], rm.labels, 
                                 top.level = FALSE)
    }
    new_dend = new_dend[!sapply(new_dend, is.null)]
    if (length(new_dend) > 1) {
      member = sum(sapply(new_dend, function(x) attr(x, 
                                                     "member")))
      class(new_dend) = "dendrogram"
      attr(new_dend, "height") = attr(dend, "height")
      attr(new_dend, "members") = member
      attr(new_dend, "edgePar") = attr(dend, "edgePar")
      attr(new_dend, "label") = attr(dend, "label")
      dend = new_dend
    }
    else if (length(new_dend) == 0) {
      dend = NULL
    }
    else if (length(new_dend) == 1) {
      dend = new_dend[[1]]
    }
  }
  else {
    if (labels(dend) %in% rm.labels) {
      cat("Remove nodes", labels(dend), "\n")
      dend = NULL
    }
  }
  if (top.level & !is.null(dend)) {
    dend = collapse_branch(dend)
  }
  return(dend)
}

pvclust_show_signif_gradient=function (dend, pvclust_obj, signif_type = c("bp", "au"), signif_col_fun = colorRampPalette(c("black", 
                                                                                                                           "darkred", "red")), ...) 
{
  require(dendextend)
  require(dplyr)
  signif_type <- match.arg(signif_type)
  pvalue_per_node <- pvclust_obj$edges[[signif_type]]
  ord <- rank(get_branches_heights(dend, sort = FALSE))
  pvalue_per_node <- pvalue_per_node[ord]
  signif_col <- signif_col_fun(100)
  pvalue_by_all_nodes <- rep(NA, nnodes(dend))
  ss_leaf <- which_leaf(dend)
  pvalue_by_all_nodes[!ss_leaf] <- pvalue_per_node
  pvalue_by_all_nodes <- na_locf(pvalue_by_all_nodes)
  the_cols <- signif_col[round(pvalue_by_all_nodes * 100)]
  signif_lwd = seq(0.5, 2, length.out = 100)
  the_lwds = signif_lwd[round(pvalue_by_all_nodes * 100)]
  dend = dend %>% assign_values_to_branches_edgePar(the_cols, 
                                                    "col") %>% assign_values_to_branches_edgePar(the_lwds, 
                                                                                                 "lwd") %>% assign_values_to_branches_edgePar(pvalue_by_all_nodes, 
                                                                                                                                              "conf")
}

rd_PCA=function (norm.dat, select.genes, select.cells, sampled.cells = select.cells, 
                 max.pca = 10, th = 2, w = NULL) 
{
  require(Matrix)
  if (is.null(w)) {
    pca = prcomp(t(as.matrix(norm.dat[select.genes, sampled.cells,drop=F])), 
                 tol = 0.01)
    pca.importance = summary(pca)$importance
    v = pca.importance[2, ]
  }
  else {
    dat <- scale(norm.dat[select.genes, select.cells,drop=F], center = TRUE, 
                 scale = TRUE)
    w = w[select.genes, select.cells]
    for (x in 1:nsamp) {
      for (y in 1:nsamp) {
        wt1 <- w[, x] * w[, y]
        cov1 <- cov.wt(e[, c(x, y)], wt = wt1, center = FALSE)$cov
        e.cov[x, y] <- cov1[1, 2]
      }
    }
    nsamp <- ncol(e)
    wt <- crossprod(w)
    cov = cov.wt(dat, wt, center = FALSE)$cov
    eig1 <- eigen(cov, symmetric = TRUE)
    eig.val <- eig1$values
    eig.val[is.na(eig.val) | eig.val < 0] <- 0
    eig.vec <- eig1$vectors
    dimnames(eig.vec) <- list(colnames(e), paste0("PC", 
                                                  1:ncol(eig.vec)))
    pca1 <- list()
    pca1$sdev <- sqrt(eig.val)
    pca1$rotation <- eig.vec
    pca1$x <- e %*% eig.vec
  }
  select = which((v - mean(v))/sd(v) > th)
  tmp = head(select, max.pca)
  if (length(sampled.cells) < length(select.cells)) {
    rot = pca$rotatio[, tmp]
    tmp.dat = norm.dat[row.names(rot), select.cells,drop=F]
    rd.dat = as.matrix(Matrix::t(tmp.dat) %*% rot)
  }
  else {
    rd.dat = pca$x[, tmp, drop = F]
  }
  return(rd.dat)
}

rd_WGCNA=function (norm.dat, select.genes, select.cells, sampled.cells = select.cells, 
                   minModuleSize = 10, cutHeight = 0.99, type = "unsigned", 
                   softPower = 4, rm.gene.mod = NULL, rm.eigen = NULL, ...) 
{
  suppressPackageStartupMessages(require(WGCNA))
  dat = as.matrix(norm.dat[select.genes, sampled.cells,drop=F])
  adj = adjacency(t(dat), power = softPower, type = type)
  adj[is.na(adj)] = 0
  TOM = TOMsimilarity(adj, TOMType = type, verbose = 0)
  dissTOM = as.matrix(1 - TOM)
  row.names(dissTOM) = colnames(dissTOM) = row.names(dat)
  rm(dat)
  gc()
  geneTree = hclust(as.dist(dissTOM), method = "average")
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, 
                                              distM = dissTOM, cutHeight = cutHeight, deepSplit = 2, 
                                              pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
  gene.mod = split(row.names(dissTOM), dynamicMods)
  gene.mod = gene.mod[setdiff(names(gene.mod), "0")]
  if (!is.null(rm.gene.mod)) {
    rm.eigen = get_eigen(rm.gene.mod, norm.dat, select.cells)[[1]]
  }
  else {
    rm.eigen = NULL
  }
  if (is.null(gene.mod) | length(gene.mod) == 0) {
    return(NULL)
  }
  gm = filter_gene_mod(norm.dat, select.cells, gene.mod, minModuleSize = minModuleSize, 
                       rm.eigen = rm.eigen, ...)
  if (is.null(gm)) {
    return(NULL)
  }
  rd.dat = gm$eigen
  return(rd.dat)
}

refine_cl=function (cl, co.ratio = NULL, cl.mat = NULL, confusion.th = 0.6, 
                    min.cells = 4, niter = 50, tol.th = 0.02, verbose = 0) 
{
  cl = setNames(as.integer(as.character(cl)), names(cl))
  while (TRUE) {
    correct = 0
    iter.num = 0
    while (iter.num < niter) {
      co.stats <- get_cl_co_stats(cl, co.ratio = co.ratio, 
                                  cl.mat = cl.mat)
      cell.cl.co.ratio <- co.stats$cell.cl.co.ratio
      tmp.dat = cell.cl.co.ratio[names(cl), as.character(sort(unique(cl)))]
      pred.cl <- setNames(colnames(tmp.dat)[apply(tmp.dat, 
                                                  1, which.max)], row.names(tmp.dat))
      if (sum(cl == pred.cl) <= correct) {
        break
      }
      correct = sum(cl == pred.cl)
      correct.frac = correct/length(cl)
      if (verbose) {
        print(correct.frac)
      }
      if (1 - correct.frac < tol.th) {
        break
      }
      tmp.cells = names(pred.cl)[pred.cl != cl[names(pred.cl)]]
      cl[tmp.cells] = pred.cl[tmp.cells]
      iter.num = iter.num + 1
    }
    cl.size = table(cl)
    cl.confusion = setNames(co.stats$cl.co.stats$confusion, 
                            row.names(co.stats$cl.co.stats))
    cl.small = names(cl.size)[cl.size < min.cells]
    rm.cl = union(names(cl.confusion)[cl.confusion > confusion.th], 
                  cl.small)
    if (length(rm.cl) == 0) {
      break
    }
    if (length(rm.cl) == length(cl.size)) {
      cl[names(cl)] = min(cl)
      return(list(cl = cl, co.stats = co.stats))
    }
    tmp.cells = names(cl)[cl %in% rm.cl]
    tmp.dat = cell.cl.co.ratio[tmp.cells, as.character(setdiff(unique(cl), 
                                                               rm.cl)), drop = F]
    pred.cl = setNames(colnames(tmp.dat)[apply(tmp.dat, 
                                               1, which.max)], row.names(tmp.dat))
    cl[tmp.cells] = pred.cl[tmp.cells]
    if (length(unique(cl)) == 1) {
      break
    }
  }
  return(list(cl = cl, co.stats = co.stats))
}

renameAndOrderClusters=function (sampleInfo, classNameColumn = "cluster_type_label", 
                                 classGenes = c("GAD1", "SLC17A7", "SLC1A3"), classLevels = c("inh", 
                                                                                              "exc", "glia"), layerNameColumn = "layer_label", regionNameColumn = "Region_label", 
                                 matchNameColumn = "cellmap_label", newColorNameColumn = "cellmap_color", 
                                 otherColumns = NULL, propLayer = 0.3, dend = NULL, orderbyColumns = c("layer", 
                                                                                                       "region", "topMatch"), includeClusterCounts = FALSE, 
                                 includeBroadGenes = FALSE, broadGenes = NULL, includeSpecificGenes = FALSE, 
                                 propExpr = NULL, medianExpr = NULL, propDiff = 0, propMin = 0.5, 
                                 medianFC = 1, excludeGenes = NULL, sortByMedian = TRUE, 
                                 sep = "_") 
{
  sampleInfo <- as.data.frame(sampleInfo)
  kpColumns <- unique(c("cluster_id", "cluster_label", "cluster_color", 
                        classNameColumn, matchNameColumn, regionNameColumn, 
                        otherColumns))
  kpColumns <- intersect(kpColumns, colnames(sampleInfo))
  clusterInfo <- t(sampleInfo[, kpColumns])
  colnames(clusterInfo) <- sampleInfo[, "cluster_label"]
  clusterInfo <- clusterInfo[, unique(colnames(clusterInfo))]
  clusterInfo <- t(clusterInfo)
  clusterInfo <- as.data.frame(clusterInfo)
  clusterInfo$old_cluster_label <- clusterInfo$cluster_label
  clusterInfo <- clusterInfo[order(clusterInfo[, "cluster_id"]), 
                             ]
  rownames(clusterInfo) <- 1:dim(clusterInfo)[1]
  if (!is.null(classNameColumn)) {
    if (!is.na(classLevels[1])) {
      clusterInfo[, classNameColumn] <- factor(clusterInfo[, 
                                                           classNameColumn], levels = classLevels)
    }
    classLabel <- clusterInfo[, classNameColumn]
    if (is.factor(classLabel)) 
      classLabel <- droplevels(classLabel)
    clusterInfo$class <- classLabel
    if (is.null(classGenes)) 
      classGenes <- NA
    if (is.na(classGenes[1])) {
      classLab <- rep("none", length(classLabel))
    }
    else {
      classLab <- classGenes[apply(propExpr[classGenes, 
                                            ], 2, which.max)]
    }
  }
  else {
    classLab <- classGenes[apply(propExpr[classGenes, ], 
                                 2, which.max)]
    names(classLevels) <- classGenes
    clusterInfo$class <- factor(classLevels[classLab], levels = classLevels)
    classNameColumn <- "class"
  }
  if (includeBroadGenes) {
    broadLab <- broadGenes[apply(propExpr[broadGenes, ], 
                                 2, which.max)]
    broadProp <- apply(propExpr[broadGenes, ], 2, max)
    broadLab[broadProp < propMin] <- classLab[broadProp < 
                                                propMin]
    names(broadLab) <- colnames(propExpr)
    clusterInfo$broadGene <- broadLab[clusterInfo$old_cluster_label]
  }
  if (includeSpecificGenes) {
    kpGn <- rep(TRUE, dim(propExpr)[1])
    specGenes <- getTopMarkersByPropNew(propExpr = propExpr[kpGn, 
                                                            ], medianExpr = medianExpr[kpGn, ], propDiff = propDiff, 
                                        propMin = propMin, medianFC = medianFC, excludeGenes = excludeGenes)
    specGenes0 <- getTopMarkersByPropNew(propExpr = propExpr[kpGn, 
                                                             ], medianExpr = medianExpr[kpGn, ], propDiff = 0, 
                                         propMin = propMin, medianFC = 0, excludeGenes = excludeGenes)
    for (s in colnames(propExpr)[(specGenes == "none")]) {
      if ((specGenes0[s] != "none") & (specGenes[s] == 
                                       "none")) {
        specGenes[s] <- specGenes0[s]
      }
    }
    clusterInfo$specificGene <- specGenes[clusterInfo$old_cluster_label]
  }
  cl3 <- sampleInfo[, "cluster_label"]
  names(cl3) <- sampleInfo$sample_id
  cl3 <- factor(cl3, levels = clusterInfo$cluster_label)
  if (is.null(newColorNameColumn)) 
    newColorNameColumn <- NA
  if (is.na(newColorNameColumn)) 
    newColorNameColumn <- "cluster_color"
  colorVec <- as.character(tapply(names(cl3), cl3, function(x) {
    col <- as.factor(sampleInfo[, newColorNameColumn])
    names(col) <- sampleInfo$sample_id
    return(names(sort(-table(col[x])))[1])
  }))
  clusterInfo$cluster_color <- makeColorsUnique(colorVec)
  if (is.null(regionNameColumn)) 
    regionNameColumn <- NA
  if (!is.na(regionNameColumn)) {
    regionVec <- as.character(tapply(names(cl3), cl3, function(x) {
      rg <- as.factor(sampleInfo[, regionNameColumn])
      names(rg) <- sampleInfo$sample_id
      rg <- rg[x]
      rg <- table(rg)/table(sampleInfo[, regionNameColumn])
      rg <- -sort(-round(100 * rg/sum(rg)))[1]
      return(paste(names(rg), rg, sep = "~"))
    }))
    clusterInfo$region <- regionVec
  }
  if (!is.na(layerNameColumn)) {
    clLayer <- sampleInfo[, layerNameColumn]
    names(clLayer) <- names(cl3)
    layerVec <- (tapply(names(cl3), cl3, function(x) {
      lyy <- factor(clLayer)[x]
      if (mean(is.na(lyy)) >= 0.5) 
        return(c(0, 0, 0, 0, 0, 0))
      lyy <- lyy[!is.na(lyy)]
      layTab <- cbind(as.numeric(names(table(lyy))), table(lyy), 
                      table(clLayer))
      return(((layTab[, 2]/layTab[, 3])/max(layTab[, 2]/layTab[, 
                                                               3])))
    }))
    rn <- names(layerVec)
    layerVec <- matrix(unlist(layerVec), ncol = 6, byrow = TRUE)
    rownames(layerVec) <- rn
    colnames(layerVec) <- 1:6
    layLab <- apply(layerVec, 1, function(x, y) {
      z <- as.numeric(colnames(layerVec)[x >= y])
      if (length(z) == 0) 
        return("x")
      if (length(z) == 1) 
        return(z)
      return(paste(range(z), collapse = "-"))
    }, propLayer)
    layLab <- paste0("L", layLab)
    clusterInfo$layer <- layLab
  }
  if (is.null(matchNameColumn)) 
    matchNameColumn <- NA
  if (!is.na(matchNameColumn)) {
    matchVec <- as.character(tapply(names(cl3), cl3, function(x) {
      y <- is.element(sampleInfo$sample_id, x)
      nm <- -sort(-table(sampleInfo[y, matchNameColumn]))
      return(names(nm)[1])
    }))
    clusterInfo$topMatch <- matchVec
  }
  if (includeClusterCounts) {
    clusterInfo$cellCount <- table(factor(sampleInfo$cluster_id, 
                                          levels = as.numeric(clusterInfo$cluster_id)))
  }
  clusterInfo$cluster_label <- as.character(clusterInfo$cluster_label)
  for (i in 1:dim(clusterInfo)[1]) {
    lab <- NULL
    if (!is.null(clusterInfo$class)) 
      lab <- paste(lab, clusterInfo$class[i], sep = sep)
    if (!is.null(clusterInfo$layer)) 
      lab <- paste(lab, clusterInfo$layer[i], sep = sep)
    if (!is.null(clusterInfo$broadGene)) 
      lab <- paste(lab, clusterInfo$broadGene[i], sep = sep)
    if (!is.null(clusterInfo$specificGene)) 
      lab <- paste(lab, clusterInfo$specificGene[i], sep = sep)
    if (!is.null(clusterInfo$region)) 
      lab <- paste(lab, clusterInfo$region[i], sep = sep)
    if (!is.null(clusterInfo$topMatch)) 
      lab <- paste(lab, clusterInfo$topMatch[i], sep = sep)
    if (!is.null(clusterInfo$cellCount)) 
      lab <- paste0(lab, sep, "N=", clusterInfo$cellCount[i])
    lab <- substr(lab, 2, nchar(lab))
    clusterInfo[i, "cluster_label"] <- lab
  }
  clNames <- clusterInfo[, "cluster_label"]
  clNames <- gsub("Exc L1", "Exc L2", clNames)
  clNames <- gsub("exc L1", "exc L2", clNames)
  clNames <- gsub("L2-2", "L2", clNames)
  clusterInfo[, "cluster_label"] <- clNames
  ordNew <- 1:dim(clusterInfo)[1]
  ordCols <- intersect(intersect(orderbyColumns, colnames(clusterInfo)), 
                       c("topMatch", "layer", "region"))
  ordVal <- "ordNew = order(clusterInfo[,classNameColumn],"
  if (length(ordCols) >= 1) 
    for (i in 1:length(ordCols)) {
      ordVal <- paste0(ordVal, "clusterInfo[,\"", ordCols[i], 
                       "\"],")
    }
  ordVal <- paste0(ordVal, "clusterInfo[,\"cluster_id\"])")
  eval(parse(text = ordVal))
  clusterInfo <- clusterInfo[ordNew, ]
  if (!is.null(dend)) {
    lab <- c(labels(dend), setdiff(labels(dend), clusterInfo$old_cluster_label))
    clusterInfo <- clusterInfo[match(lab, clusterInfo$old_cluster_label), 
                               ]
    print("NOTE: You will need to regenerate the dendrogram with these new cluster labels now.")
  }
  rownames(clusterInfo) <- clusterInfo$lrank <- clusterInfo$cluster_id <- 1:dim(clusterInfo)[1]
  clusterInfo
}

reorder_cl=function (cl, dat) 
{
  cl.means = get_cl_means(dat, cl)
  cl.hc = hclust(as.dist(1 - cor(cl.means)), method = "average")
  cl = setNames(factor(as.character(cl), levels = cl.hc$labels[cl.hc$order]), 
                names(cl))
  cl = setNames(as.integer(cl), names(cl))
}

reorder_dend=function (dend, l.rank, top.level = TRUE) 
{
  tmp.dend = dend
  sc = sapply(1:length(dend), function(i) {
    l = dend[[i]] %>% labels
    mean(l.rank[dend[[i]] %>% labels])
  })
  ord = order(sc)
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      if (ord[i] != i) {
        dend[[i]] = tmp.dend[[ord[i]]]
      }
      if (length(dend[[i]]) > 1) {
        dend[[i]] = reorder_dend(dend[[i]], l.rank, 
                                 top.level = FALSE)
      }
    }
  }
  if (top.level) {
    dend = collapse_branch(dend, 10^-10)
  }
  return(dend)
}

reorder_factor=function (clMatch, clRef) 
{
  out <- table(clMatch, clRef)
  out <- t(out)/colSums(out)
  fac <- apply(out, 2, function(x) which.max(x) + 0.01 * sum(cumsum(x)))
  fac <- names(sort(fac))
  factor(clMatch, levels = fac)
}

run_consensus_clust=function (norm.dat, select.cells = colnames(norm.dat), niter = 100, 
                              sample.frac = 0.8, output_dir = "subsample_result", mc.cores = 1, 
                              de.param = de_param(), merge.type = c("undirectional", "directional"), 
                              override = FALSE, init.result = NULL, cut.method = "auto", 
                              ...) 
{
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  all.cells = select.cells
  if (!is.null(init.result)) {
    all.cells = intersect(all.cells, names(init.result$cl))
  }
  run <- function(i, ...) {
    prefix = paste("iter", i, sep = ".")
    print(prefix)
    outfile = file.path(output_dir, paste0("result.", i, 
                                           ".rda"))
    if (file.exists(outfile) & !override) {
      return(NULL)
    }
    select.cells = sample(all.cells, round(length(all.cells) * 
                                             sample.frac))
    save(select.cells, file = file.path(output_dir, paste0("cells.", 
                                                           i, ".rda")))
    result <- iter_clust(norm.dat = norm.dat, select.cells = select.cells, 
                         prefix = prefix, de.param = de.param, merge.type = merge.type, 
                         result = init.result, ...)
    save(result, file = outfile)
  }
  if (mc.cores == 1) {
    sapply(1:niter, function(i) {
      run(i, ...)
    })
  }
  else {
    require(foreach)
    require(doParallel)
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    foreach(i = 1:niter, .combine = "c") %dopar% run(i)
    stopCluster(cl)
  }
  result.files = file.path(output_dir, dir(output_dir, "result.*.rda"))
  load(result.files[[1]])
  cl.size = table(result$cl)
  graph.size = sum(cl.size^2)
  if (graph.size < 10^8) {
    co.result <- collect_co_matrix_sparseM(norm.dat, result.files, 
                                           all.cells)
    co.ratio = co.result$co.ratio
    consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, 
                                            norm.dat, select.cells = all.cells, de.param = de.param, 
                                            merge.type = merge.type, method = cut.method, result = init.result)
    refine.result = refine_cl(consensus.result$cl, co.ratio = co.ratio, 
                              tol.th = 0.01, confusion.th = 0.6, min.cells = de.param$min.cells)
    markers = consensus.result$markers
  }
  else {
    result <- iter_clust(norm.dat = norm.dat, select.cells = all.cells, 
                         de.param = de.param, merge.type = merge.type, result = init.result, 
                         ...)
    co.result <- collect_subsample_cl_matrix(norm.dat, result.files, 
                                             all.cells)
    cl = merge_cl_by_co(result$cl, co.ratio = NULL, cl.mat = co.result$cl.mat, 
                        diff.th = 0.25)
    refine.result = refine_cl(cl, cl.mat = co.result$cl.mat, 
                              tol.th = 0.01, confusion.th = 0.6, min.cells = de.param$min.cells)
    markers = result$markers
  }
  cl = refine.result$cl
  merge.result = merge_cl(norm.dat = norm.dat, cl = cl, rd.dat = Matrix::t(norm.dat[markers, ,drop=F
                                                                                    ]), de.param = de.param, merge.type = merge.type, return.markers = FALSE)
  return(list(co.result = co.result, cl.result = merge.result))
}

sample_cells=function (cl, sample.size, weights = NULL) 
{
  n_cl <- unique(cl)
  if (length(sample.size) == 1) {
    sample.size <- setNames(rep(sample.size, length(n_cl)), 
                            n_cl)
  }
  cl.cells <- split(names(cl), cl)
  sampled.cells <- unlist(sapply(names(cl.cells), function(x) {
    cells <- cl.cells[[x]]
    if (sample.size[[x]] == length(cells)) {
      return(cells)
    }
    to.sample <- pmin(sample.size[[x]], length(cells))
    if (!is.null(weights)) {
      sampled <- sample(cells, to.sample, prob = weights[cells])
    }
    else {
      sampled <- sample(cells, to.sample)
    }
    sampled
  }, simplify = FALSE))
  return(sampled.cells)
}

sample_cl_list=function (cl.list, max.cl.size = 500) 
{
  select.cells = c()
  for (cl in cl.list) {
    cl.size = table(cl)
    more.cl = cl[setdiff(names(cl), select.cells)]
    more.size = table(more.cl)
    add.cells = unlist(lapply(names(more.size), function(x) {
      sample(names(more.cl)[more.cl == x], min(more.size[[x]], 
                                               max.cl.size))
    }))
    select.cells = c(select.cells, add.cells)
  }
  return(select.cells)
}

score_gene_mod=function (norm.dat, select.cells, gene.mod, eigen = NULL, method = "average", 
                         de.param = de_param(), max.cl.size = NULL) 
{
  if (length(gene.mod) == 0) {
    return(NULL)
  }
  if (is.null(eigen)) {
    eigen = get_eigen(gene.mod, norm.dat[, names(select.cells),drop=F])
  }
  colnames(eigen) = names(gene.mod)
  gene.mod.val = sapply(names(gene.mod), function(y) {
    x = gene.mod[[y]]
    if (is.null(max.cl.size)) {
      tmp.dat = norm.dat[x, select.cells,drop=F]
    }
    else {
      v = eigen[select.cells, y]
      ord = order(v)
      tmp.cells = unique(select.cells[c(head(ord, max.cl.size), 
                                        tail(ord, max.cl.size))])
      tmp.dat = norm.dat[x, tmp.cells,drop=F]
    }
    tmp.dat = as.matrix(tmp.dat)
    tmp.dat = tmp.dat - rowMeans(tmp.dat)
    if (method == "average") {
      tmp.cl = cutree(hclust(dist(t(tmp.dat)), method = "average"), 
                      2)
    }
    else if (method == "ward.D") {
      tmp.cl = cutree(hclust(dist(t(tmp.dat)), method = "ward.D"), 
                      2)
    }
    else if (method == "kmeans") {
      tmp.cl = kmeans(t(tmp.dat), 2)$cluster
    }
    else if (method == "louvain") {
      tmp = jaccard_louvain(t(tmp.dat), 15)
      if (is.null(tmp)) {
        return(list(c(0, 0), NULL))
      }
      tmp.cl = tmp$cl
    }
    else {
      stop(paste("Unknown method", method))
    }
    de.genes = de_score(as.matrix(norm.dat[, names(tmp.cl),drop=F]), 
                        cl = tmp.cl, de.param = de.param)
    de.genes = de.genes[sapply(de.genes, length) > 1]
    if (length(de.genes) > 0) {
      sc = max(sapply(de.genes, function(x) x$score))
      gene.num = max(sapply(de.genes, function(x) x$num))
    }
    else {
      sc = 0
      gene.num = 0
    }
    return(list(c(sc = sc, gene.num = gene.num), de.genes = de.genes))
  }, simplify = F)
  return(gene.mod.val)
}

select_markers=function (norm.dat, cl, n.markers = 20, de.genes = NULL, ...) 
{
  if (is.null(de.genes)) {
    de.genes = de_score(norm.dat, cl, ...)
  }
  pairs = names(de.genes)
  pairs.df = gsub("cl", "", do.call("rbind", strsplit(pairs, 
                                                      "_")))
  row.names(pairs.df) = pairs
  select.pairs = pairs[pairs.df[, 1] %in% cl & pairs.df[, 
                                                        2] %in% cl]
  de.markers = sapply(select.pairs, function(s) {
    tmp = de.genes[[s]]
    c(head(tmp$up.genes, n.markers), head(tmp$down.genes, 
                                          n.markers))
  }, simplify = F)
  markers = intersect(unlist(de.markers), row.names(norm.dat))
  return(list(markers = markers, de.genes = de.genes[select.pairs]))
}

select_markers_pair=function (norm.dat, de.genes, add.genes, gene.score = NULL, 
                              rm.genes = NULL, top.n = 50, max.num = 2000) 
{
  pairs = do.call("rbind", strsplit(gsub("cl", "", names(add.genes)), 
                                    "_"))
  row.names(pairs) = names(add.genes)
  de.genes.list = list()
  if (is.null(gene.score)) {
    de.genes = de.genes[names(add.genes)]
    tmp = get_gene_score(de.genes, top.n = top.n, max.num = max.num, 
                         bin.th = 4)
    up.gene.score = tmp$up.gene.score
    down.gene.score = tmp$down.gene.score
    gene.score = pmin(up.gene.score, down.gene.score)
    row.names(gene.score) = row.names(up.gene.score)
  }
  select.genes = setdiff(row.names(gene.score), rm.genes)
  gene.score = gene.score[select.genes, names(add.genes), 
                          drop = F]
  final.genes = list()
  while (sum(add.genes) > 0 & nrow(gene.score) > 0) {
    g = order(rowSums(gene.score))[1]
    g.n = row.names(gene.score)[g]
    select.pair = colnames(gene.score)[gene.score[g, ] < 
                                         top.n]
    add.genes[select.pair] = add.genes[select.pair] - 1
    for (p in select.pair) {
      de.genes.list[[p]] = union(de.genes.list[[p]], g.n)
    }
    add.genes = add.genes[add.genes > 0]
    gene.score = gene.score[-g, names(add.genes), drop = F]
  }
  return(de.genes.list)
}

select_markers_pair_direction=function (de.genes, add.up, add.down, up.gene.score = NULL, 
                                        down.gene.score = NULL, rm.genes = NULL, top.n = 50, max.num = 2000) 
{
  up.genes = down.genes = list()
  final.genes = c()
  if (is.null(up.gene.score)) {
    pairs.n = union(names(add.up), names(add.down))
    pairs = do.call("rbind", strsplit(gsub("cl", "", pairs.n), 
                                      "_"))
    row.names(pairs) = pairs.n
    de.genes = de.genes[pairs.n]
    tmp = get_gene_score(de.genes, top.n = top.n, max.num = max.num, 
                         bin.th = 4)
    up.gene.score = tmp$up.gene.score
    down.gene.score = tmp$down.gene.score
  }
  select.genes = setdiff(row.names(up.gene.score), rm.genes)
  while ((sum(add.up) + sum(add.down) > 0) & length(select.genes) > 
         0) {
    sc = 0
    if (length(add.up) > 0) {
      sc = rowSums(up.gene.score[select.genes, names(add.up), 
                                 drop = F])
    }
    if (length(add.down) > 0) {
      sc = sc + rowSums(down.gene.score[select.genes, 
                                        names(add.down), drop = F])
    }
    g = select.genes[which.min(sc)]
    select.up.pair = names(add.up)[up.gene.score[g, names(add.up)] < 
                                     top.n]
    select.down.pair = names(add.down)[down.gene.score[g, 
                                                       names(add.down)] < top.n]
    if (length(select.up.pair) + length(select.down.pair) == 
        0) {
      break
    }
    add.up[select.up.pair] = add.up[select.up.pair] - 1
    add.down[select.down.pair] = add.down[select.down.pair] - 
      1
    add.up = add.up[add.up > 0, drop = F]
    add.down = add.down[add.down > 0, drop = F]
    for (p in select.up.pair) {
      up.genes[[p]] = c(up.genes[[p]], g)
    }
    for (p in select.down.pair) {
      down.genes[[p]] = c(down.genes[[p]], g)
    }
    final.genes = c(final.genes, g)
    select.genes = setdiff(select.genes, final.genes)
    if (length(select.genes) > 0) {
      select = rep(TRUE, length(select.genes))
      if (length(add.up) > 0) {
        select = rowSums(up.gene.score[select.genes, 
                                       names(add.up), drop = F] < top.n) > 0
      }
      if (length(add.down) > 0) {
        select = select | rowSums(down.gene.score[select.genes, 
                                                  names(add.down), drop = F] < top.n) > 0
      }
      select.genes = select.genes[select]
    }
    cat(g, length(add.up), length(add.down), "\n")
  }
  markers = final.genes
  return(list(markers = markers, up.genes = up.genes, down.genes = down.genes))
}

select_N_markers=function (de.genes, up.gene.score = NULL, down.gene.score = NULL, 
                           default.markers = NULL, pair.num = 1, rm.genes = NULL) 
{
  add.up = add.down = setNames(rep(pair.num, length(de.genes)), 
                               names(de.genes))
  if (!is.null(default.markers)) {
    up.default = sapply(de.genes, function(x) {
      intersect(x$up.genes, default.markers)
    }, simplify = F)
    down.default = sapply(de.genes, function(x) {
      intersect(x$down.genes, default.markers)
    }, simplify = F)
    add.up = pmax(add.up - sapply(up.default, length), 0)
    add.down = pmax(add.down - sapply(down.default, length), 
                    0)
  }
  add.up = add.up[add.up > 0, drop = F]
  add.down = add.down[add.down > 0, drop = F]
  result = select_markers_pair_direction(add.up, add.down, 
                                         de.genes = de.genes, up.gene.score = up.gene.score, 
                                         down.gene.score = down.gene.score, rm.genes = c(rm.genes, 
                                                                                         default.markers), top.n = 50, max.num = 2000)
  up.genes = up.default
  down.genes = down.default
  for (x in names(result$up.genes)) {
    up.genes[[x]] = c(up.genes[[x]], result$up.genes[[x]])
  }
  for (x in names(result$down.genes)) {
    down.genes[[x]] = c(down.genes[[x]], result$down.genes[[x]])
  }
  return(list(up.genes = up.genes, down.genes = down.genes, 
              markers = c(default.markers, result$markers)))
}

selectMarkersPair=function (norm.dat, add.genes, de.genes = NULL, gene.score = NULL, 
                            rm.genes = NULL, top.n = 50, max.num = 2000) 
{
  pairs = do.call("rbind", strsplit(gsub("cl", "", names(add.genes)), 
                                    "_"))
  row.names(pairs) = names(add.genes)
  de.genes.list = list()
  if (is.null(gene.score)) {
    de.genes = de.genes[names(add.genes)]
    tmp = get_gene_score(de.genes, top.n = top.n, max.num = max.num, 
                         bin.th = 4)
    up.gene.score = tmp$up.gene.score
    down.gene.score = tmp$down.gene.score
    gene.score = pmin(up.gene.score, down.gene.score)
    row.names(gene.score) = all.genes
  }
  select.genes = setdiff(row.names(gene.score), rm.genes)
  gene.score = gene.score[select.genes, names(add.genes), 
                          drop = F]
  final.genes = list()
  while (sum(add.genes) > 0 & nrow(gene.score) > 0) {
    g = order(rowSums(gene.score))[1]
    g.n = row.names(gene.score)[g]
    select.pair = colnames(gene.score)[gene.score[g, ] < 
                                         top.n]
    add.genes[select.pair] = add.genes[select.pair] - 1
    for (p in select.pair) {
      de.genes.list[[p]] = union(de.genes.list[[p]], g.n)
    }
    add.genes = add.genes[add.genes > 0]
    gene.score = gene.score[-g, names(add.genes), drop = F]
  }
  return(de.genes.list)
}

selectMarkersPairDirection=function (add.up, add.down, de.genes = NULL, up.gene.score = NULL, 
                                     down.gene.score = NULL, rm.genes = NULL, top.n = 50, max.num = 2000) 
{
  up.genes = down.genes = list()
  final.genes = c()
  if (is.null(up.gene.score)) {
    pairs.n = union(names(add.up), names(add.down))
    pairs = do.call("rbind", strsplit(gsub("cl", "", pairs.n), 
                                      "_"))
    row.names(pairs) = pairs.n
    de.genes = de.genes[pairs.n]
    tmp = get_gene_score(de.genes, top.n = top.n, max.num = max.num, 
                         bin.th = 4)
    up.gene.score = tmp$up.gene.score
    down.gene.score = tmp$down.gene.score
  }
  select.genes = setdiff(row.names(up.gene.score), rm.genes)
  while ((sum(add.up) + sum(add.down) > 0) & length(select.genes) > 
         0) {
    sc = 0
    if (length(add.up) > 0) {
      sc = rowSums(up.gene.score[select.genes, names(add.up), 
                                 drop = F])
    }
    if (length(add.down) > 0) {
      sc = sc + rowSums(down.gene.score[select.genes, 
                                        names(add.down), drop = F])
    }
    g = select.genes[which.min(sc)]
    select.up.pair = names(add.up)[up.gene.score[g, names(add.up)] < 
                                     top.n]
    select.down.pair = names(add.down)[down.gene.score[g, 
                                                       names(add.down)] < top.n]
    if (length(select.up.pair) + length(select.down.pair) == 
        0) {
      break
    }
    add.up[select.up.pair] = add.up[select.up.pair] - 1
    add.down[select.down.pair] = add.down[select.down.pair] - 
      1
    add.up = add.up[add.up > 0, drop = F]
    add.down = add.down[add.down > 0, drop = F]
    for (p in select.up.pair) {
      up.genes[[p]] = c(up.genes[[p]], g)
    }
    for (p in select.down.pair) {
      down.genes[[p]] = c(down.genes[[p]], g)
    }
    final.genes = c(final.genes, g)
    select.genes = setdiff(select.genes, final.genes)
    if (length(select.genes) > 0) {
      select = rep(TRUE, length(select.genes))
      if (length(add.up) > 0) {
        select = rowSums(up.gene.score[select.genes, 
                                       names(add.up), drop = F] < top.n) > 0
      }
      if (length(add.down) > 0) {
        select = select | rowSums(down.gene.score[select.genes, 
                                                  names(add.down), drop = F] < top.n) > 0
      }
      select.genes = select.genes[select]
    }
    cat(g, length(add.up), length(add.down), "\n")
  }
  markers = final.genes
  return(list(markers = markers, up.genes = up.genes, down.genes = down.genes))
}

selectMarkersPairGroup=function (norm.dat, cl, g1, g2, de.genes, top.n = 50, max.num = 1000, 
                                 n.markers = 20, up.gene.score = NULL, down.gene.score = NULL) 
{
  pairs = do.call("rbind", strsplit(names(de.genes), "_"))
  pairs = gsub("cl", "", pairs)
  row.names(pairs) = names(de.genes)
  up.pairs = row.names(pairs)[pairs[, 1] %in% g1 & pairs[, 
                                                         2] %in% g2]
  down.pairs = row.names(pairs)[pairs[, 1] %in% g2 & pairs[, 
                                                           2] %in% g1]
  select.pairs = c(up.pairs, down.pairs)
  if (is.null(up.gene.score)) {
    tmp = get_gene_score(de.genes, top.n = top.n, max.num = max.num, 
                         bin.th = 4)
    up.gene.score = tmp$up.gene.score
    down.gene.score = tmp$down.gene.score
    row.names(down.gene.score) = all.genes
  }
  all.genes = row.names(up.gene.score)
  tmp.up.gene.score = cbind(up.gene.score[, up.pairs, drop = F], 
                            down.gene.score[, down.pairs, drop = F])
  tmp.down.gene.score = cbind(down.gene.score[, up.pairs, 
                                              drop = F], up.gene.score[, down.pairs, drop = F])
  up.genes = row.names(tmp.up.gene.score)[head(order(rowSums(tmp.up.gene.score)), 
                                               n.markers)]
  down.genes = row.names(tmp.down.gene.score)[head(order(rowSums(tmp.down.gene.score)), 
                                                   n.markers)]
  up.num = colSums(tmp.up.gene.score[up.genes, , drop = F] < 
                     max.num)
  down.num = colSums(tmp.down.gene.score[down.genes, , drop = F] < 
                       max.num)
  total.num = up.num + down.num
  add.genes = setNames(rep(n.markers, ncol(tmp.up.gene.score)), 
                       colnames(tmp.up.gene.score)) - total.num
  add.genes = add.genes[add.genes > 0]
  up.genes = up.genes[rowMins(tmp.up.gene.score[up.genes, 
                                                , drop = F]) < max.num]
  down.genes = down.genes[rowMins(tmp.down.gene.score[down.genes, 
                                                      , drop = F]) < max.num]
  genes = union(up.genes, down.genes)
  if (length(add.genes) > 0) {
    tmp = selectMarkersPair(norm.dat, add.genes = add.genes, 
                            de.genes = de.genes, gene.score = pmin(tmp.up.gene.score, 
                                                                   tmp.down.gene.score), rm.genes = c(up.genes, 
                                                                                                      down.genes), top.n = top.n)
    genes = union(genes, unlist(tmp))
  }
  return(genes)
}

selectNMarkers=function (de.genes, up.gene.score = NULL, down.gene.score = NULL, 
                         default.markers = NULL, pair.num = 1, rm.genes = NULL) 
{
  add.up = add.down = setNames(rep(pair.num, length(de.genes)), 
                               names(de.genes))
  if (!is.null(default.markers)) {
    up.default = sapply(de.genes, function(x) {
      intersect(x$up.genes, default.markers)
    }, simplify = F)
    down.default = sapply(de.genes, function(x) {
      intersect(x$down.genes, default.markers)
    }, simplify = F)
    add.up = pmax(add.up - sapply(up.default, length), 0)
    add.down = pmax(add.down - sapply(down.default, length), 
                    0)
  }
  add.up = add.up[add.up > 0, drop = F]
  add.down = add.down[add.down > 0, drop = F]
  result = selectMarkersPairDirection(add.up, add.down, de.genes = de.genes, 
                                      up.gene.score = up.gene.score, down.gene.score = down.gene.score, 
                                      rm.genes = c(rm.genes, default.markers), top.n = 50, 
                                      max.num = 2000)
  up.genes = up.default
  down.genes = down.default
  for (x in names(result$up.genes)) {
    up.genes[[x]] = c(up.genes[[x]], result$up.genes[[x]])
  }
  for (x in names(result$down.genes)) {
    down.genes[[x]] = c(down.genes[[x]], result$down.genes[[x]])
  }
  return(list(up.genes = up.genes, down.genes = down.genes, 
              markers = c(default.markers, result$markers)))
}

set_pair_matrix=function (m, rows, cols, vals) 
{
  coor <- get_pair_matrix_coor(m, rows, cols)
  m[coor] <- vals
  return(m)
}

sparse_cor=function (m) 
{
  library(Matrix)
  n_rows <- nrow(m)
  n_cols <- ncol(m)
  ii <- unique(m@i) + 1
  Ex <- colMeans(m)
  nozero <- as.vector(m[ii, ]) - rep(Ex, each = length(ii))
  covmat <- (crossprod(matrix(nozero, ncol = n_cols)) + crossprod(t(Ex)) * 
               (n_rows - length(ii)))/(n_rows - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/crossprod(t(sdvec))
  return(cormat)
}

tau_one_vs_other_genes=function (cl, cl.present.counts, present.th = 0.4, tau.th = 0.8, 
                                 top.n = 10) 
{
  all.g = row.names(cl.present.counts)
  all.g = setdiff(all.g, c(grep("LOC", all.g, value = T), 
                           grep("Rik$", all.g, value = T)))
  cl.size = table(cl)[colnames(cl.present.counts)]
  cl.present.prob = t(t(cl.present.counts)/as.vector(cl.size))
  tau.genes = sapply(colnames(cl.present.counts), function(x) {
    fg = cl.present.counts[, x]/sum(cl == x)
    bg = (rowSums(cl.present.counts) - cl.present.counts[, 
                                                         x])/(length(cl) - sum(cl == x))
    tau = (fg - bg)/pmax(bg, fg)
    g = all.g[cl.present.prob[all.g, x] > 0.5]
    g = g[order(tau[g], decreasing = T)]
    g = g[tau[g] > 0.97 | (tau[g] > 0.7 & g %in% head(g, 
                                                      top.n))]
  }, simplify = F)
}

type=function (dat, k = 10) 
{
  suppressPackageStartupMessages(library(Rphenograph))
  rpheno <- Rphenograph(dat, k = k)
  cl <- setNames(rpheno[[2]]$membership, row.names(dat)[as.integer(rpheno[[2]]$names)])
  return(list(cl = cl, result = rpheno))
}

type.RANN=function (dat, k = 10) 
{
  library(igraph)
  library(matrixStats)
  library(RANN)
  knn.matrix = nn2(dat, k = k)[[1]]
  jaccard.adj <- knn_jaccard(knn.matrix)
  jaccard.gr <- igraph::graph.adjacency(jaccard.adj, mode = "undirected", 
                                        weighted = TRUE)
  louvain.result <- igraph::cluster_louvain(jaccard.gr)
  mod.sc <- igraph::modularity(louvain.result)
  if (pass_louvain(mod.sc, jaccard.adj)) {
    cl <- setNames(louvain.result$membership, row.names(dat))
    return(list(cl = cl, result = louvain.result))
  }
  else {
    return(NULL)
  }
}

unbranch_by_conf=function (dend, conf.th) 
{
  if (length(dend) > 1) {
    conf = c()
    for (i in 1:length(dend)) {
      if (is.null(attr(dend[[i]], "edgePar"))) {
        conf[i] = 1
      }
      else {
        conf[i] = attr(dend[[i]], "edgePar")$conf
      }
      dend[[i]] = unbranch_by_conf(dend[[i]], conf.th)
    }
    select = conf < conf.th
    select.children = which(select)
    if (length(select.children) > 0) {
      unchanged = which(!select)
      new_dend = dend[unchanged]
      names(new_dend) = unchanged
      for (i in select.children) {
        if (length(dend[[i]]) > 1) {
          for (j in 1:length(dend[[i]])) {
            ind = sprintf("%02d", j)
            attr(dend[[i]][[j]], "edgePar") = attr(dend[[i]], 
                                                   "edgePar")
            new_dend[[paste(i, ind, sep = ".")]] = dend[[i]][[j]]
          }
        }
        else {
          new_dend[[as.character(i)]] = dend[[i]]
        }
      }
      new_dend = new_dend[order(as.integer(names(new_dend)))]
      class(new_dend) = "dendrogram"
      attr(new_dend, "height") = attr(dend, "height")
      attr(new_dend, "members") = attr(dend, "members")
      attr(new_dend, "midpoint") = attr(dend, "midpoint")
      attr(new_dend, "edgePar") = attr(dend, "edgePar")
      attr(new_dend, "label") = attr(dend, "label")
      dend = new_dend
    }
  }
  return(dend)
}

unbranch_by_length=function (dend, length.th) 
{
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      dend[[i]] = unbranch_by_length(dend[[i]], length.th)
    }
    child.h = get_childrens_heights(dend)
    h = attr(dend, "height")
    select = h - child.h < length.th
    select.children = which(select)
    if (length(select.children) > 0) {
      unchanged = which(!select)
      new_dend = dend[unchanged]
      names(new_dend) = unchanged
      for (i in select.children) {
        if (length(dend[[i]]) > 1) {
          for (j in 1:length(dend[[i]])) {
            ind = sprintf("%02d", j)
            idx = paste(i, ind, sep = ".")
            new_dend[[idx]] = dend[[i]][[j]]
            attr(new_dend[[idx]], "edgePar") = attr(dend[[i]], 
                                                    "edgePar")
          }
        }
        else {
          new_dend[[as.character(i)]] = dend[[i]]
        }
      }
      new_dend = new_dend[order(as.integer(names(new_dend)))]
      class(new_dend) = "dendrogram"
      attr(new_dend, "height") = attr(dend, "height")
      attr(new_dend, "members") = attr(dend, "members")
      attr(new_dend, "midpoint") = attr(dend, "midpoint")
      attr(new_dend, "edgePar") = attr(end, "edgePar")
      dend = new_dend
    }
  }
  return(dend)
}

updateSampDat=function (Samp.dat, clusterInfo) 
{
  lab <- as.character(Samp.dat$cluster_label)
  for (i in 1:dim(clusterInfo)[1]) {
    kp <- lab == clusterInfo$old_cluster_label[i]
    Samp.dat$cluster_label[kp] <- clusterInfo$cluster_label[i]
    Samp.dat$cluster_id[kp] <- clusterInfo$cluster_id[i]
    Samp.dat$cluster_color[kp] <- clusterInfo$cluster_color[i]
  }
  Samp.dat
}

vec_chisq_test=function (x, x.total, y, y.total) 
{
  total <- x.total + y.total
  present = x + y
  absent = total - x - y
  o <- cbind(x, x.total - x, y, y.total - y)
  e <- cbind(present * x.total, absent * x.total, present * 
               y.total, absent * y.total)
  e <- e/as.vector(total)
  stat <- rowSums(pmax(0, abs(o - e) - 0.5)^2/e)
  tmp <- cbind(stats = stat, pval = pchisq(stat, 1, lower.tail = FALSE), 
               logFC = log2((x * y.total)/(y * x.total)), diff = x/x.total - 
                 y/y.total)
  tmp
}

within_group_specific_markers=function (cl.g, norm.dat, cl, ...) 
{
  cl = droplevels(cl[cl %in% cl.g])
  df = do.call("rbind", sapply(cl.g, function(x) {
    df = group_specific_markers(x, norm.dat, cl, ...)
    if (!is.null(df)) {
      df$cl = x
      df
    }
    else {
      NULL
    }
  }, simplify = F))
  return(df)
}



plot_RD_meta <- function(rd.dat, meta, meta.col=NULL,show.legend=TRUE, cex=0.15, legend.size=5,alpha.val=1)
  {
    rd.dat = as.data.frame(rd.dat)
    colnames(rd.dat)[1:2] = c("Dim1","Dim2")
    library(ggplot2)
    rd.dat$meta = meta
    p=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=meta),size=cex)
    if(is.factor(meta)){
      rd.dat = droplevels(rd.dat)
      if(is.null(meta.col)){
        if(length(levels(meta)) > 2){
          meta.col = setNames(jet.colors(length(levels(meta))), levels(meta))
        }
        else{
          meta.col = setNames(c("blue", "orange"), levels(meta))
        }
      }      
      p = p+ scale_color_manual(values=alpha(as.vector(meta.col[levels(rd.dat$meta)]),alpha.val))
      p = p+ theme(panel.background=element_blank(),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"))
    }
    else{
      p = p+ scale_color_gradient(low="blue",high="red")
    }
    if(!show.legend){
      p = p + theme(legend.position="none") 
    }
    else{
      if(is.factor(meta)){
        p = p + guides(colour = guide_legend(override.aes = list(size=legend.size)))
      }
    }
    #p = p + coord_fixed(ratio=1)
    return(p)
  }



plot_RD_gene <- function(rd.dat, norm.dat, genes, cex=0.15)
  {
    library(ggplot2)
    plots=list()
    rd.dat = as.data.frame(rd.dat)
    colnames(rd.dat)[1:2] = c("Dim1","Dim2")
    for(g in genes){
      rd.dat$expr = norm.dat[g,row.names(rd.dat)]
      p=ggplot(rd.dat, aes(Dim1, Dim2)) + geom_point(aes(color=expr),size=cex)
      p = p+ scale_color_gradient(low="gray80",high="red") 
      p = p + theme_void() + theme(legend.position="none")
      p = p + coord_fixed(ratio=1)
      p = p + ggtitle(g)
      plots[[g]]= p
    }
    return(plots)
  }




 varibow <- function(n_colors) {
  sats <- rep_len(c(0.55,0.7,0.85,1),length.out = n_colors)
  vals <- rep_len(c(1,0.8,0.6),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}