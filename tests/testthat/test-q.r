test_that("Q blocks", {
    if(requireNamespace("fmesher")){
        data(example_mesh, package = "mems")
        fm <- fmesher::fm_fem(example_mesh)
        qb <- q_blocks(example_mesh)
        ## just checks @x slot entries
        expect_equal(as.matrix(fm$c0), as.matrix(qb$C0))
        expect_equal(as.matrix(fm$c1), as.matrix(qb$C1))
        expect_equal(as.matrix(fm$g1), as.matrix(qb$G1))
        expect_equal(as.matrix(fm$b1), as.matrix(qb$B1))
        expect_equal(fm$ta, as.matrix(qb$triangle_areas))
        }
})
