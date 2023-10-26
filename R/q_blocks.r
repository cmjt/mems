#' Calculates the precision matrix blocks contributions from each triangle/node
#' in the triangulation
#' @source \url{https://github.com/inlabru-org/fmesher/blob/d1dcbba1a42144a22b8e8056cb7236381e0b67ef/src/mesh.cc#L2278}
#' @details Simpler version of the \code{C++} versions implemented in
#' \code{fmesher::fm_fem}. Is only for demonstration.
#' @param mesh object of \code{fmesher::fm_mesh_2d} (i.e., 2D Delaunay triangulation, 'mesh').
q_blocks <- function(mesh)  {
    n_v <- mesh$n
    n_t <- nrow(mesh$graph$tv)
    ## Initialize matrices
    C1 <- C0 <- matrix(0, n_v, n_v)
    G1 <- matrix(0, n_v, n_v)
    B1 <- matrix(0, n_v, n_v)
    triangle_areas <- mems(mesh)$area
    ##
    TV <-  mesh$graph$tv ## which triangles share nodes/verticies
    TT <- mesh$graph$tt ## which triangles (t)ouch
    S <- mesh$loc
    for (t in 1:n_t) {
        tv <- TV[t, ]
        s0 <- S[tv[1], ]
        s1 <- S[tv[2], ]
        s2 <- S[tv[3], ]
        ## edges
        e <- list(
            s2 - s1,
            s0 - s2,
            s1 - s0
        )
        ## browser()
        eij <- matrix(0, 3, 3)
        for (i in 1:3) {
            eij[i, i] <- sum(e[[i]] * e[[i]])
            js <- 1:3
            js <- js[js > i & js <= 3]
            for(j in js){
                eij[i, j] <- sum(e[[i]] * e[[j]])
                eij[j, i] <- eij[i, j]
                
            }
        }
        ## C0, C1, G1
        a <- triangle_areas[t]
        fa <- abs(sum(crossprod(e[[1]], e[[2]]))) / 2
        for (i in 1:3) {
            C0[tv[i], tv[i]] <- C0[tv[i], tv[i]] + a / 3
            C1[tv[i], tv[i]] <- C1[tv[i], tv[i]] + a / 6
            G1[tv[i], tv[i]] <- G1[tv[i], tv[i]] + eij[i, i] / (4 * fa)
            js <- 1:3
            js <- js[js > i & js <= 3]
            for(j in js){
                C1[tv[i], tv[j]] <- C1[tv[i], tv[j]] + a / 12
                C1[tv[j], tv[i]] <- C1[tv[j], tv[i]] + a / 12
                vij <- eij[i, j] / (4 * fa)
                G1[tv[i], tv[j]] <- G1[tv[i], tv[j]] + vij
                G1[tv[j], tv[i]] <- G1[tv[j], tv[i]] + vij
            }
        }
        ## B1, boundary matrix
        b <- logical(3)
        for(boo in 1:3){
         if(is.na(TT[t, boo])) {b[boo] <- TRUE}
            }
        if (any(b)) {
            vij <- -1 / (4 * fa)
            for (i in 1:3) {
                for (j in 1:3) {
                    for (k in 1:3) {
                        if (b[k] && (i != k)) {
                            B1[tv[i], tv[j]] <- B1[tv[i], tv[j]] + eij[k, j] * vij
                        }
                    }
                }
            }
        }
    } ## matches t
    ## turn Q blocks into sparse Matrcies
    return(list(C0 = methods::as(methods::as(C0, "generalMatrix"), "TsparseMatrix"),
                C1 = methods::as(methods::as(C1, "generalMatrix"), "TsparseMatrix"),
                G1 = methods::as(methods::as(G1, "generalMatrix"), "TsparseMatrix"),
                B1 = methods::as(methods::as(B1, "generalMatrix"), "TsparseMatrix"),
                triangle_areas = triangle_areas))
}




    
