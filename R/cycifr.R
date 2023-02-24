#' Fit Gaussian Mixture Model
#'
#' Fits a Gaussian Mixture Model to a given input data array, and returns a list 
#' with a ggplot object and the cutoff points. This internally uses the python 
#' sklearn.mixture package
#'
#' @param array Input data array to fit the GMM to
#' @param n_components Number of Gaussian components to use
#' 
#' @return A list containing the following components:
#' \item{plot}{A ggplot2 object showing the GMM fit and histogram of input data}
#' \item{cutoffs}{A numeric vector of the cutoff points of the GMM components}
#'
#' @importFrom reticulate import np_array
#' @import ggplot2
#' @import data.table
#' 
#' @examples 
#' fit_gmm(iris$Petal.Length, n_components=2)
#' 
#' @export
fit_gmm <- function(array, n_components = 2, n_steps=500) {
  
  # setup python
  sklearn <- reticulate::import("sklearn.mixture")
  #np <- import("numpy", convert=FALSE)
  
  # fit the model
  pygmm = sklearn$GaussianMixture(n_components=2L, covariance_type="spherical",random_state=0L)
  pygmm$fit(reticulate::np_array(array)$reshape(-1L,1L))
  
  # Compute a set of cutoff points that can be used to threshold the data
  low  <- min(pygmm$means_) - 2*sqrt(pygmm$covariances_[which.min(pygmm$means_)])
  high <- max(pygmm$means_) + 2*sqrt(pygmm$covariances_[which.max(pygmm$means_)])
  ref_space = seq(low, high, length.out = 5000)
  result = pygmm$predict(reticulate::np_array(ref_space)$reshape(-1L, 1L))
  idx = which(diff(result) != 0)
  cutoffs = ref_space[idx]
  
  # Evaluate the PDF of the GMM at a set of points
  x = seq(min(array), max(array), length.out=n_steps)
  pyx = reticulate::np_array(x)$reshape(-1L, 1L)
  log_prob = pygmm$score_samples(pyx)
  responsibilities = pygmm$predict_proba(pyx)
  pdf = exp(log_prob)
  pdf_individual = data.table(sweep(responsibilities, MARGIN=1, STATS=pdf, FUN="*"))
  pdf_individual$x = x
  pdf_individual = melt(pdf_individual, id.vars = "x", measure.vars=c("V1","V2"))
  pdf = data=data.table(x=x,pdf=pdf)
  
  # plot
  #g <- ggplot() + geom_histogram(data=data.table(x=array), binwidth=0.01, alpha=0.3, aes(x=x, y=after_stat(density))) +
  #  geom_line(data=pdf_individual, aes(x=x, y=value, color=variable)) + 
  #  geom_line(data=pdf, aes(x=x,y=pdf), color="black", linetype="dashed") + 
  #  geom_vline(xintercept=cutoffs, color="red") + 
  #  theme_bw() + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y = element_blank()) 
  
  pdf_individual$pdf = rep(pdf$pdf, n_components)
  return(list(cutoffs=cutoffs, pdf=pdf_individual))
  
}

#' Fit GMMs for each marker in a cell table
#'
#' This function fits a Gaussian mixture model (GMM) to the expression values of each marker in a cell table. The input
#' cell table is a data.table object with required columns: CellID, marker, and value (expression level). The function
#' first checks that the cell table has a "CellID" column and a manageable size (< 3 million cells). It then extracts
#' the unique markers and fits a GMM for each marker using the \code{fit_gmm} function. The output is a list with 
#' elements for the GMM fit (pdf), cutoff points, histogram of original data, and the ggplot of the GMM fit for each marker.
#'
#' @param dt A data.table with required columns: CellID, marker, and value (expression level).
#' @param markers Optional character vector of markers to use, if only a subset of markers is to be analyzed.
#' @return A list containing the ggplot, cutoff points, histogram of original data, and the GMM fit for each marker.
#' @export
#' @import ggplot2 aes color coord_cartesian element_blank facet_wrap geom_area geom_line geom_vline ggplot theme_bw
#' @import data.table rbindlist setnames setkey
#' @importFrom utils stdout
fit_all_gmms <- function(dt, markers=NA,
                         n_components = 2, n_steps=500,
                         n_sample=NA) {
  
  # sanity check
  if (!all("CellID" %in% colnames(dt))) {
    stop("Should be a cell table (e.g. with column CellID")
  }
  
  if (nrow(dt) > 3e6) {
    stop("Cell table too big. Please subset to < 3 million cells. Did you accidentally provide all of the samples?")
  }
  
  # isolate and melt just this one sample
  useless <- c("Eccentricity","Solidity","Extent","Orientation","Area","MajorAxisLength","MinorAxisLength")
  id.vars <- intersect(colnames(dt),
                       c("CellID","X_centroid","Y_centroid","sample",
                         useless))
  cells.m <- melt(dt, id.vars=id.vars)[,.(variable,value)]
  setnames(cells.m, "variable","marker")
  setkey(cells.m, marker)
  
  # isolate the markers
  markers <- as.character(cells.m[,unique(marker)])
  
  # isolate the markers
  markers <- as.character(cells.m[,unique(marker)])
  
  # fit the GMM for each marker
  fit_list <- lapply(markers, function(v) {
    write(paste("...",v), stdout())
    
    # fit the GMM using sklearn.mixture
    set.seed(1337)
    val <- cells.m[v][value > 0, log10(value)]
    if (!is.na(n_sample))
      val <- sample(val, n_sample, replace=FALSE)
    
    fit <- fit_gmm(val,
                   n_components=n_components,
                   n_steps=n_steps)
    fit$pdf = fit$pdf[, marker := v]
    fit$cutoffs = data.table(marker=v, cutoff=fit$cutoffs)
    
    # get the histogram of the original data for ease of plotting
    h <- hist(val, breaks=fit$pdf[variable=="V1"]$x, plot=FALSE)
    b <- h$breaks
    fit$hist = data.table(breaks=(b[-1] + b[-length(b)])/2, density=h$density, marker=v)
    
    return(fit)
  })
  
  # aggregate the data for plotting
  dt.plot <- rbindlist(lapply(fit_list, function(x) x$pdf))
  dt.hist <- rbindlist(lapply(fit_list, function(x) x$hist))
  dt.cut  <- rbindlist(lapply(fit_list, function(x) x$cutoffs))
  
  # create plot object
  g <- ggplot() + 
    geom_area(data=dt.hist, alpha=0.3, aes(x=breaks, y=density)) +
    geom_line(data=dt.plot, aes(x=x, y=value, color=variable)) + 
    geom_line(data=dt.plot, aes(x=x,y=pdf), color="black", linetype="dashed") + 
    facet_wrap(~marker, scales="free_y") + 
    geom_vline(data=dt.cut, aes(xintercept=cutoff), color="red") + 
    theme_bw() + 
    theme(legend.position="none", axis.title.x=element_blank(), axis.title.y = element_blank()) +
    coord_cartesian(xlim=c(2,4))

  return(list(plot=g, cut=dt.cut, hist=dt.hist, pdf=dt.plot))
}
