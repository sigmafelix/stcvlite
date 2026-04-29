

#' Generate manual rset object from spatiotemporal cross-validation indices
#' @param cvindex integer. Output of [`generate_cv_index`].
#' @param data data.frame from `stdt`. Should be the
#' same object as what was used for `covars` argument of [`generate_cv_index`]
#' @param cv_mode character(1). Spatiotemporal cross-validation indexing method.
#' See `cv_mode` description in [`generate_cv_index`].
#' @returns rset object of `rsample` package. A tibble with a list column of
#' training-test data.frames and a column of labels.
#' @author Insang Song
#' @importFrom rsample make_splits
#' @importFrom rsample manual_rset
#' @export
convert_cv_index_rset <- function(cvindex, data, cv_mode) {
  maxcvi <- max(cvindex)
  len_cvi <- seq_len(maxcvi)
  list_cvi <- split(len_cvi, len_cvi)
  list_cvi_rows <-
    lapply(
      list_cvi,
      function(x) {
        list(analysis = which(cvindex != x),
             assessment = which(cvindex == x))
      }
    )
  list_split_dfs <-
    lapply(
      list_cvi_rows,
      function(x) {
        rsample::make_splits(x = x, data = data)
      }
    )
  modename <- sprintf("cvfold_%s_%03d", cv_mode, len_cvi)
  rset_stcv <- rsample::manual_rset(list_split_dfs, modename)
  return(rset_stcv)

}