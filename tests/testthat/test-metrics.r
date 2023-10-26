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
    a <- dist(B,C); b <- dist(A, C); c <- dist(A, B)
    m <- do.call('rbind', metrics(a, b, c))
    expect_equal(as.numeric(m[,1]), rep(1, 6))
})
test_that("mems", {
    if(requireNamespace("sf")){
    data(example_mesh, package = "mems")
    meshmetrics <- mems(example_mesh)
    expect_equal(meshmetrics$circumcircle_R, rep(2, 4))
    expect_equal(round(meshmetrics$incircle_r, 3), rep(0.8280, 4))
    expect_equal(round(meshmetrics$radius_edge, 3), rep(0.7070, 4))
    }
})
test_that("segments", {
    data(example_mesh, package = "mems")
    segs <- segments(example_mesh)
    expect_equal(round(segs$length[1], 3), 2.828)
})
