#' Reduce transcript models
#'
#' Combine identical transcript models and get non-redundent and reduced
#' transcript models
#'
#' @param transcript_models_ls A list of matrices where each row in the matrix
#' corresponds to a bin and each column is a transcript.
#' @param bin_operation Three different modes to deal with decimals in the
#' transript model (due to partial overlap of the first or last exon and bins).
#' Either "ceiling", "floor", or "round" (default: "round").
#'
#' @return  A list of matrices holding the reduced transcript models and a
#' dataframe holding each transcript belongs to which group and reduced model.
#' @export

reduce_transcript_models <-
    function(transcript_models_ls,
             bin_operation = c("round", "floor", "ceiling")) {
        # check transcript class
        tx_model_class <- unique(sapply(transcript_models_ls, class))
        if (length(tx_model_class) != 1 || tx_model_class != "matrix") {
            stop("invalid input for trascript models")
        }
        # check bin_operation
        if (length(bin_operation) > 1) bin_operation <- bin_operation[[1]]
        if (!bin_operation %in% c("ceiling", "floor", "round")) {
            stop("invalid arguement for bin_operation, it should be 'ceiling',
                 'floor' or 'round'.")
        }
        # deal with decimals in transcript models
        integer_models <- lapply(transcript_models_ls,
                                 function(x, fun) {
                                     do.call(fun, list(x))
                                 }, fun = bin_operation)
        # compute reduced transcript groups
        tx_groups <- lapply(integer_models, group_trancript_models)
        # get reduced tx models
        get_reduced_tx_models <- function(tx_models, tx_groups) {
            tx_models[, unlist(lapply(tx_groups, function(x) x[[1]])),
                      drop = FALSE]
        }
        reduced_tx_models <- mapply(get_reduced_tx_models,
                                    integer_models,
                                    tx_groups, SIMPLIFY = FALSE)
        # create group and model identifier for each transcript
        group_len <- lapply(tx_groups, lengths)
        tx_num <- sapply(group_len, sum)
        tx_group_model <- data.frame(tx_name = unlist(tx_groups),
                   group = rep(seq_along(tx_num), tx_num),
                   model = unlist(mapply(rep, lapply(group_len, seq_along),
                                         group_len, SIMPLIFY = FALSE)))
        rownames(tx_group_model) <- NULL
        return(list(reduced_tx_models, tx_group_model))
}

#' Group transcripts
#'
#' Group transripts if they have identical transcript model.
#'
#' @param tx_models A processed matrix without decimals in transcript model
#' @seealso \code{\link{reduce_transcript_models}},
#' @return  A list of grouped transcript names.

group_trancript_models <- function(tx_models) {
    pool <- colnames(tx_models)
    groups <- list()
    while (length(pool) > 1) {
        identical_models <- which(
            colSums(abs(tx_models[, pool[-1], drop = FALSE] -
                            tx_models[, pool[1]])) == 0
        )
        groups <- append(groups, list(c(pool[1], names(identical_models))))
        pool <- pool[-c(1, identical_models + 1)]
    }
    if (length(pool) == 1) groups <- append(groups, pool)
    return(groups)
}
