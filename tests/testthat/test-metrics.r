test_that("Making triangles", {
    A <- c(0, 1)
    B <- c(1, 1)
    equi <- equilateral(A, B)
    ra <- right_angled(A,B)
    is <- isoceles(A, B)
    expect_equal(round(equi, 3), c(0.5, 1.866))
    expect_equal(ra, c(0, 2))
    expect_equal(round(is, 3), c(0.5, 1.866))
})
test_that("Metrics equilateral",{
    A <- c(0, 1)
    B <- c(1, 1)
    C <- equilateral(A, B)
    m <- do.call('rbind', metrics(A, B, C))
    expect_equal(as.numeric(m[13:18,1]), rep(1, 6))
})
test_that("mems", {
    if(requireNamespace("sf")){
    data(example_mesh, package = "mems")
    meshmetrics <- mems(example_mesh)
    expect_equal(round(meshmetrics$circumcircle_R, 2), c(1.25, 1.58, 1.12, 1.41))
    expect_equal(round(meshmetrics$incircle_r, 3), c(0.618, 0.744, 0.382, 0.586))
    expect_equal(round(meshmetrics$re, 3), c(0.625, 0.707, 1.118, 0.707))
    }
})
test_that("segments", {
    data(example_mesh, package = "mems")
    segs <- segments(example_mesh)
    expect_equal(round(segs$length[1], 3), 2.236)
})
test_that("mems plot", {
    if(requireNamespace("sf")){
    data(example_mesh, package = "mems")
    meshmetrics <- mems(example_mesh)
    expect_true(ggplot2::is.ggplot(gg_mems(meshmetrics)))
    }
})
