#' Predict posterior probabilities for variants to be in 
#' each cluster in a fitted LDA model
#'
#' @param predict.data A data frame with character-valued columns 
#'        rs and cat and numeric-valued columns with annotations.
#'        Each row contains data for one SNP in one tissue.
#'        rs is an ID for the SNP, which need not be unique, and cat 
#'        is an ID for each tissue, which must match an ID
#'        for which training data was included when fitting the
#'        LDA model. Aannotation columns must have column names matching
#'        those supplied when fitting the LDA model.
#' @param summary.c A fitted LDA model created using 
#'        \code{\link{FitLDAModel}} or \code{\link{FitLDAModelNewTissues}}.
#' @return p.labeled, a data frame with one row per training variant
#'         with the posterior probability of each variant to be in each
#'         each cluster (with column names CLUSTER1,...) and columns
#'         cat and rs.
#' @examples
#' \dontrun{
#'   data(training)
#'   summary.c <- FitLDAModel(training.data=training, nclust=3,
#'                          kde.nbins=100, iters=50, inner.iters=50)
#'   pred <- Predict(predict.data=training, summary.c=summary.c)
#' }
#' @export
Predict <- function(predict.data, summary.c){
  predict.data <- as.data.frame(predict.data)
  z <- predict.data[, colnames(summary.c$data)]
  cat.dict <- summary.c$categories$ncat
  names(cat.dict) <- summary.c$categories$cat
  cat <- cat.dict[predict.data$cat]
  pred <- ebmme.Predict(data=z, cat=cat, 
                  block_id=summary.c$block.id, 
                  nclust=summary.c$nclust, 
                  levels=summary.c$levels.data, 
                  fs_binned=summary.c$fs_binned, 
                  a=summary.c$a)
  p.labeled <- as.data.frame(t(pred$p))
  colnames(p.labeled) <- paste0("CLUSTER", 1:ncol(p.labeled))
  p.labeled <- cbind(p.labeled, predict.data[, c("cat", "rs")])
  return(p.labeled)
}

#' Update the tissues in a fitted LDA model
#'
#' @param newtissue.data A data frame with character-valued columns 
#'        rs and cat and numeric-valued columns with annotations.
#'        Each row contains data for one SNP in one tissue
#'        rs is an ID for the SNP, which need not be unique, and
#'        cat is an ID for each tissue, which must match an ID
#'        for which training data was included when fitting the
#'        LDA model.  Annotation columns must have column names matching
#'        those supplied when fitting the LDA model.
#' @inheritParams Predict
#' @inheritParams FitLDAModel
#' @return a fitted LDA model, as returned by \code{\link{FitLDAModel}},
#'         that can be used with \code{\link{Predict}} for
#'         the tissues in newtissue.data
#' @examples 
#' \dontrun{
#'   data(training)
#'   summary.c <- FitLDAModel(training.data=training, nclust=3,
#'                          kde.nbins=100, iters=50, inner.iters=50)
#'   data(tissuenew)
#'   newtissues <- FitLDAModelNewTissues(newtissue.data=tissuenew,
#'                                      summary.c=summary.c,
#'                                      inner.iters=50)
#' }
#' @export
FitLDAModelNewTissues <- function(newtissue.data, summary.c, inner.iters=200){
  newtissue.data <- as.data.frame(newtissue.data)
  z <- newtissue.data[, colnames(summary.c$data)]
  cat <- newtissue.data$cat
  ncat <- as.numeric(as.factor(cat)) - 1
  df <- unique(data.frame(cat=cat, ncat=ncat))
  pred <- ebmme.FitNewTissue(data=z, cat=ncat,
                  block_id=summary.c$block.id,
                  nclust=summary.c$nclust,
                  levels=summary.c$levels.data,
                  fs_binned=summary.c$fs_binned, 
                  inner.iters=inner.iters)
  p.labeled <- as.data.frame(t(pred$p))
  colnames(p.labeled) <- paste0("CLUSTER", 1:ncol(p.labeled))
  p.labeled <- cbind(p.labeled, newtissue.data[, c("cat", "rs")])
  a.labeled <- as.data.frame(t(pred$a))
  colnames(a.labeled) <- paste0("CLUSTER", 1:ncol(a.labeled))
  a.labeled$ncat <- 0:max(df$ncat)
  a.labeled <- merge(a.labeled, df)
  a.labeled <- a.labeled[, !colnames(a.labeled)=="ncat"]
  summary.c$a <- pred$a
  summary.c$p <- pred$p
  summary.c$data <- z
  summary.c$cat <- ncat
  summary.c$all_as <- NULL
  summary.c$data.labeled <- NULL
  summary.c$categories <- df
  summary.c$a.labeled <- a.labeled
  summary.c$p.labeled <- p.labeled
  return(summary.c)
}

#' fit an LDA model
#' @param training.data A data frame with character-valued columns 
#'        rs and cat and numeric-valued columns with annotations.
#'        Each row is data for one SNP in one tissue.
#'        rs is an ID for the SNP, which need not be unique, and
#'        cat is an ID for each tissue.
#' @param nclust Integer specifying the number of clusters to estimate
#' @param kde.nbins Integer specifying how many bins to use for binning 
#'        each annotation
#' @param iters Integer specifying number of outer iterations
#' @param inner.iters Integer specifying number of inner iterations
#' @return A fitted LDA model, i.e., a list (apart from elements
#'         used internally) with elements
#'         \describe{
#'             \item{p.labeled}{a data frame with one row per 
#'                              training variant, with
#'                              the posterior probability for each 
#'                              variant to be in 
#'                              each cluster in columns CLUSTER1,... and
#'                              also with columns cat and rs}
#'             \item{a.labeled}{a data frame with one row per tissue
#'                              with membership vectors for each tissue
#'                              with columns cat and CLUSTER1,...}
#'         } 
#' @examples 
#' \dontrun{
#'   data(training)
#'   summary.c <- FitLDAModel(training.data=training, nclust=3,
#'                          kde.nbins=100, iters=50, inner.iters=50)
#' }
#' @export
FitLDAModel <- function(training.data, nclust=9, kde.nbins=1000, iters=250, 
                        inner.iters=200){
  training.data <- as.data.frame(training.data)
  cat <- training.data$cat
  ncat <- as.numeric(as.factor(cat)) - 1
  df <- unique(data.frame(cat=cat, ncat=ncat))
  training.data <- as.data.frame(training.data)
  z <- training.data[, !colnames(training.data) %in% c("rs", "cat")]
  summary.c <- ebmme.lda(data= z, cat=ncat, block.id = 1:ncol(z), 
                         iters = iters, 
                         inner.iters = inner.iters, nclust=nclust, 
                         kde.nbins=kde.nbins)
  # fix up p and a
  training.data$order <- 1:nrow(training.data)
  data.labeled <- merge(training.data, df)
  data.labeled <- data.labeled[order(data.labeled$order), ]
  data.labeled <- data.labeled[, !colnames(data.labeled) %in% "order"]
  summary.c$data.labeled <- data.labeled
  summary.c$categories <- df
  a.labeled <- as.data.frame(t(summary.c$a))
  colnames(a.labeled) <- paste0("CLUSTER", 1:ncol(a.labeled))
  a.labeled$ncat <- 0:max(df$ncat)
  a.labeled <- merge(a.labeled, df)
  a.labeled <- a.labeled[, !colnames(a.labeled)=="ncat"]
  summary.c$a.labeled <- a.labeled
  p.labeled <- as.data.frame(t(summary.c$p))
  colnames(p.labeled) <- paste0("CLUSTER", 1:ncol(p.labeled))
  p.labeled <- cbind(p.labeled, data.labeled[, c("cat", "rs")])
  summary.c$p.labeled <- p.labeled
  summary.c
}
