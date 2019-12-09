################################################################################
# build tx models for test
# test on indentical and different tx
tx_model_1 <- matrix(1, nrow = 5, ncol = 3,
                     dimnames = list(as.character(1:5), paste0("t1.", 1:3)))
tx_model_1[1, 2] <- 0

# test on partial overlap
tx_model_2 <- matrix(1, nrow = 5, ncol = 3,
                     dimnames = list(as.character(1:5), paste0("t2.", 1:3)))

tx_model_2[c(1, 5), 2] <- 0.3
tx_model_2[c(1, 5), 3] <- 0.7

# test on grouping multiple groups and name concatenate
tx_model_3 <- matrix(1, nrow = 5, ncol = 4,
                     dimnames = list(as.character(1:5)))

colnames(tx_model_3) <- as.vector(t(outer(c("t3.", "t4."), 1:2, paste0)))
tx_model_3[c(1, 5), 1:2] <- 0

tx_model_4 <- cbind(tx_model_1, tx_model_2)

# test inputs
test_that("Inputs are correct for reduce_transcript_models", {
    # not a list
    expect_error(reduce_transcript_models(tx_model_1, "invalid input"))
    # element in list is not matrix
    expect_error(reduce_transcript_models(list(tx_model_1, NA)),
                 "invalid input")
    expect_error(reduce_transcript_models(list(tx_model_1, 3)),
                 "invalid input")
    # invalid argument
    expect_error(reduce_transcript_models(list(tx_model_1), "Round"),
                 "invalid arguement")
})

test_that("Transcript models are reduced correctly", {
    true_reduced_modle_1 <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t1.1 = rep(1, 5),
            t1.2 = c(0, rep(1, 4))
        ))

    true_reduced_modle_2r <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t2.1 = rep(1, 5),
            t2.2 = c(0, rep(1, 3), 0)
        ))

    true_reduced_modle_2c <-
        as.matrix(data.frame(row.names = as.character(1:5),
                             t2.1 = rep(1, 5)))

    true_reduced_modle_2f <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t2.1 = rep(1, 5),
            t2.2 = c(0, rep(1, 3), 0)
        ))

    true_reduced_modle_3 <-
        as.matrix(data.frame(
            row.names = as.character(1:5),
            t3.1 = c(0, rep(1, 3), 0),
            t4.1 = rep(1, 5)
        ))

    expect_equal(
        reduce_transcript_models(list(tx_model_1, tx_model_2))[[1]],
        list(true_reduced_modle_1, true_reduced_modle_2r)
    )

    expect_equal(
        reduce_transcript_models(list(tx_model_1, tx_model_2), "ceiling")[[1]],
        list(true_reduced_modle_1, true_reduced_modle_2c)
    )

    expect_equal(
        reduce_transcript_models(list(tx_model_1, tx_model_3), "floor")[[1]],
        list(true_reduced_modle_1, true_reduced_modle_3)
    )

})

test_that("Transcript groupings are correct", {
    true_group_1r <-
        data.frame(
            TXNAME = c("t1.1", "t1.3", "t1.2", "t2.1", "t2.3", "t2.2"),
            group = c(rep(1, 3), rep(2, 3)),
            model = c(1, 1, 2, 1, 1, 2)
        )

    true_group_1f <-
        data.frame(
            TXNAME = c("t1.1", "t1.3", "t1.2", "t2.1", "t2.2", "t2.3"),
            group = c(rep(1, 3), rep(2, 3)),
            model = c(1, 1, 2, 1, 2, 2)
        )

    true_group_4 <-
        data.frame(
            TXNAME = c("t1.1", "t1.3", "t2.1", "t2.3", "t1.2", "t2.2"),
            group = c(rep(1, 6)),
            model = c(1, 1, 1, 1, 2, 3)
        )

    expect_equal(reduce_transcript_models(list(tx_model_1, tx_model_2))[[2]],
                 true_group_1r)

    expect_equal(
        reduce_transcript_models(list(tx_model_1, tx_model_2), "floor")[[2]],
        true_group_1f)

    expect_equal(reduce_transcript_models(list(tx_model_4))[[2]], true_group_4)
})
