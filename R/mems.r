#' Calculate the distance between two points
#' @param A A numeric vector of length 2 specifying (first) location.
#' @param B A numeric vector of length 2 specifying (second) location.
#' @examples 
#' dist(c(1, 1), c(2, 1))
#' @export
dist <- function(A, B) {
  sqrt((A[1] - B[1])^2 + (A[2] - B[2])^2)
}
#' Find mid point on straight line between two points
#' @inheritParams dist
#' @examples
#' mid(c(1, 1), c(2, 1))
#' @export
mid <- function(A, B){
    xm <- (A[1] + B[1])/2
    ym <- (A[2] + B[2])/2
    return(c(xm, ym))
}
#' Cosine triangle function to calculate angle C (in degrees)
#' given edge lengths a, b, and c
#' 
#' @param a A numeric triangle edgelength (opposite angle A).
#' @param b A numeric triangle edgelength (opposite angle B).
#' @param c A numeric triangle edgelength (opposite angle C).
#' @examples
#' ang(1, 2, 1.5)
#' @export
ang <- function(a, b, c) {
  cosC <- (a^2 + b^2 - c^2) / (2 * a * b)
  angC <- acos(cosC)
  return(angC * (180 / pi))
}
#' Sin function to calculate the area of a triangle given
#' edge lengths a, b and angle C.
#' 
#' @param a A numeric triangle edgelength (opposite angle A).
#' @param b A numeric triangle edgelength (opposite angle B).
#' @param angC angle (in degrees)
#' @examples 
#' A <- c(1, 1)
#' B <- c(2, 1)
#' C <- c(2.5, 1.5)
#' a <- dist(B, C)
#' b <- dist(A, C)
#' c <- dist(A, B)
#' angC <- ang( a, b, c)
#' tri_area(a, b, angC)
#' @export
tri_area <- function(a, b, angC) {
    (1 / 2) * a * b * sin(angC / (180 / pi))
}
#' Calculate cicumcircle radius given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples
#' circum_R(c(1, 1), c(2, 1), c(3, 2.5))
#' @export
circum_R <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  abc <- a * b * c
  d1 <- a + b + c
  d2 <- b + c - a
  d3 <- c + a - b
  d4 <- a + b - c
  return((abc) / (sqrt(d1 * d2 * d3 * d4)))
}
#' Calculate cicumcircle centroid given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples
#' circum_O(c(1, 1), c(2, 1), c(3, 2.5))
#' @export
circum_O <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  angA <- ang(b, c, a) * (pi / 180)
  angB <- ang(a, c, b) * (pi / 180)
  angC <- ang(b, a, c) * (pi / 180)
  sumsins <- sin(2 * angA) + sin(2 * angB) + sin(2 * angC)
  xo <- (A[1] * sin(2 * angA) + B[1] * sin(2 * angB) + C[1] * sin(2 * angC)) / sumsins
  yo <- (A[2] * sin(2 * angA) + B[2] * sin(2 * angB) + C[2] * sin(2 * angC)) / sumsins
  return(c(xo, yo))
}
#' Calculate incumcircle radius given vertex A, B, C
#' locations of a triangle
#' 
#' @inheritParams incircle_r
#' @examples 
#' incircle_O(c(1, 1), c(2, 1), c(3, 2.5))
#' @export
incircle_O <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  abc <- a + b + c
  xc <- (a * A[1] + b * B[1] + c *  C[1]) / abc
  yc <- (a * A[2] + b * B[2] + c *  C[2]) / abc
  return(c(xc, yc))
}
#' Calculate incumcircle centroid given vertex A, B, C
#' locations of a triangle
#' 
#' @param A A numeric vector of length 2 specifying vertex location "A".
#' @param B A numeric vector of length 2 specifying vertex location "B".
#' @param C A numeric vector of length 2 specifying vertex location "C".
#' @examples 
#' incircle_r(c(1, 1), c(2, 1), c(3, 2.5))
#' @export
incircle_r <- function(A, B, C) {
  a <- dist(B, C)
  b <- dist(A, C)
  c <- dist(B, A)
  s <- (a + b + c) / 2
  sqrt(((s - a) * (s - b) * (s - c)) / s)
}
#' Extract a dataframe of mesh triangle segments (start and end locations)
#' 
#' @inheritParams mems
#' @examples 
#' data(example_mesh)
#' segments(example_mesh)
#' @export
segments <- function(mesh) {
  df <- rbind(data.frame(a = mesh$loc[mesh$graph$tv[, 1], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 2], c(1, 2)]),
              data.frame(a = mesh$loc[mesh$graph$tv[, 2], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 3], c(1, 2)]),
              data.frame(a = mesh$loc[mesh$graph$tv[, 3], c(1, 2)],
                         b = mesh$loc[mesh$graph$tv[, 1], c(1, 2)]))
  colnames(df) <- c("x", "y", "xend", "yend")
  df$length <- c(unlist(dist(df[, 1:2], df[, 3:4])))
  return(df)
}
#' Calculate all interior Delaunay triangulation triangle angles
#' 
#' @inheritParams mems
#' @examples 
#' data(example_mesh)
#' mesh_ang(example_mesh)
#' @export
mesh_ang <- function(mesh) {
  tv <- mesh$graph$tv
  angs <- matrix(numeric(3 * nrow(tv)), ncol = 3)
  for(i in 1:nrow(tv)) {
    ## the three verts of one triangle
    vs <- data.frame(x = rep(mesh$loc[tv[i,], 1], each = 2),
                     y = rep(mesh$loc[tv[i,], 2], each = 2))
    vs$xend <- vs$x[c(3, 5, 1, 5, 1, 3)];vs$yend <- vs$y[c(3, 5, 1, 5, 1, 3)]
    vs <- vs[c(1,2,4), ]
    dists <- c(unlist(dist(vs[, 1:2], vs[, 3:4])))
    angs[i, 1] <- ang(dists[1], dists[2], dists[3])
    angs[i, 2] <- ang(dists[2], dists[3], dists[1])
    angs[i, 3] <- ang(dists[3], dists[1], dists[2])
  }
  angs <- as.data.frame(angs)
  names(angs) <- c("angleA", "angleB", "angleC")
  angs$ID <- 1:nrow(angs)
  return(angs)
}
#' Calculate minimum edge length for each triangle in a given Delaunay triangulation
#' 
#' @inheritParams mems
#' @examples
#' data(example_mesh)
#' lmin(example_mesh)
#' @export
lmin <- function(mesh) {
  tv <- mesh$graph$tv
  lmin <- numeric(nrow(tv))
  for (i in 1:nrow(tv)) {
    A <- mesh$loc[tv[i, 1], 1:2]
    B <- mesh$loc[tv[i, 2], 1:2]
    C <- mesh$loc[tv[i, 3], 1:2]
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(A, B)
    lmin[i] <- min(a, b, c)
  }
  return(lmin)
}
#' Transform a \code{fmesher::fm_mesh_2d} into a \code{sf} object
#' 
#' @inheritParams mems
#' @source Modified from \code{sp} based function suggested by Finn in the
#' R-inla discussion Google Group
#' \url{https://groups.google.com/g/r-inla-discussion-group/c/z1n1exlZrKM}.
#' @return A simple features, \code{sf}, object.
#' @seealso \code{\link{mems}}
#' @examples
#' data(example_mesh)
#' sf <- mesh_2_sf(example_mesh)
#' @export
mesh_2_sf <- function(mesh) {
    st <- sf::st_sfc(lapply(
                  1:nrow(mesh$graph$tv),
                  function(x) {
                      tv <- mesh$graph$tv[x, , drop = TRUE]
                      sf::st_polygon(list(mesh$loc[tv[c(1, 3, 2, 1)],
                                                   1:2,
                                                   drop = FALSE]))
                  }
              )
              )
    dat <- as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE])
    dat$ID <- 1:length(st)
    res <- sf::st_sf(dat,
                     geometry = st)
    return(res)
}
#'
#' @source \url{https://gist.github.com/joeroe/43bb6efa233cb5c8994b6a6d41a4911a}
deldir_2_sf <- function(deldir) {
    deldir %>%
        purrr::map(~{cbind(x = .$x, y = .$y)} %>%
                       rbind(.[1,]) %>%
                       list() %>%
                       sf::st_polygon()) %>%
        sf::st_sfc() %>%
        sf::st_sf() %>%
        dplyr::mutate(id = dplyr::row_number()) %>%
        return()
}
#' Calculate a number of different geometric attributes of a Delaunay triangulation
#' 
#' Calculates a number of geometric attributes for a given
#' Delaunay triangulation based on the circumscribed and inscribed circle of each triangle.
#' 
#' @param mesh A \code{fmesher::fm_mesh_2d} object.
#' @param plot If \code{TRUE} then triangle quality metrics are plotted.
#' @param metric character specifying which metric of \code{x} to plot one of "metric_1",
#' "metric_2", "metric_3", "metric_4", "metric_5", "metric_6" (default "metric_1"), ignored if \code{plot = FALSE}.
#' @return An object of class \code{sf} with the following data for each triangle in the
#' triangulation
#' \itemize{
#' \item \code{V1}, \code{V2}, and \code{V3} corresponding vertices
#' of \code{mesh} matches \code{mesh$graph$tv};
#' \item \code{ID}, numeric triangle id;
#' \item \code{angleA}, \code{angleB}, and \code{angleC}, the
#' interior angles;
#' \item circumcircle radius, circumradius, \code{circumcircle_R} (\eqn{R});
#' \item incircle radius \code{incircle_r} (\eqn{r});
#' \item centroid locations of the circumcircle, circumcenter, (\code{c_Ox, c_Oy});
#' \item centroid locations of the incircle, incenter, (\code{i_Ox, i_Oy});
#' \item the radius-edge ratio \code{radius_edge} \eqn{\frac{R}{l_{min}}},
#' where \eqn{l_{min}} is the minimum edge length;
#' \item the radius ratio \code{radius_ratio} \eqn{\frac{r}{R}};
#' \item \code{area}, area (\eqn{A});
#' \item \code{metric_c1}:  sum of each triangle's nodes rowsums of the diagonal of the  C1 matrix 
#' \item \code{metric__g1}: sum of each triangle's nodes rowsums of the diagnonal of the G1 matrix 
#' \item \code{metric_b1}: summation of each triangle's nodes rowsums of the diagonal of the B1 matrix 
#' \item \code{metric_1}: \eqn{\frac{4\sqrt{3}|A|}{\Sigma_{i = 1}^3 L_i^2}},
#' where \eqn{L_i} is the length of edge \eqn{i}.
#' \item \code{metric_2}: \eqn{6 \sqrt{\frac{A}{\sqrt{3}}}}.
#' \item \code{metric_3}:\eqn{\frac{\text{min}(L_i)}{\text{max}(L_i)}}.
#' \item \code{metric_4}:  \eqn{3 \frac{\text{min}(\phi_i)}{\pi}}, where
#' \eqn{\phi_i} is the \eqn{i^{th}} internal angle.
#' \item \code{metric_5}: \eqn{\frac{4\sqrt{3} A}{\text{max}(L_i)} * \Sigma_{i = 1}^3 L_i^2}.
#' \item \code{metric_6}: \eqn{\frac{1}{4 q_b} + q_w \frac{\sqrt{3}}{16}}, where
#' \eqn{q_b = \frac{\Sigma_{i = 1}^3 L_i^2}{4 \sqrt{3} A}} and
#' \eqn{q_m = (\frac{1}{3A})(\Sigma_{i = 1}^3 L_i)^2}.
#' }
#' @details A triangle's circumcircle (circumscribed circle) is the unique circle that passes
#' through each of its three vertices. A triangle's incircle (inscribed circle) is the
#' largest circle that can be contained within it (i.e., touches it's three edges).
#' @examples
#' data(example_mesh, package = "mems")
#' metrics <- mems(example_mesh)
#' @export
mems <- function(mesh, plot = FALSE, metric = "metric_1") {
    tv <- mesh$graph$tv
    fm <- fmesher::fm_fem(mesh)
    c1 <- g1 <- b1 <- numeric(nrow(tv))
    c_R <- i_R <- area <- numeric(nrow(tv))
    c_O <- i_O <- matrix(rep(0, 2 * nrow(tv)), ncol = 2)
    quality_metrics  <- list()
    for (i in 1:nrow(tv)) {
        c1[i] <- sum(sum(fm$c1[tv[i, 1],]), sum(fm$c1[tv[i, 2],]), sum(fm$c1[tv[i, 3],]))
        g1[i] <- sum(fm$g1[tv[i, 1],tv[i, 1]], fm$g1[tv[i, 2],tv[i, 2]], fm$g1[tv[i, 3],tv[i, 3]])
        b1[i] <- sum(fm$b1[tv[i, 1], tv[i, 1]], fm$b1[tv[i, 2], tv[i, 2]], fm$b1[tv[i, 3],tv[i, 3]])
        A <- mesh$loc[tv[i, 1], 1:2]
        B <- mesh$loc[tv[i, 2], 1:2]
        C <- mesh$loc[tv[i, 3], 1:2]
        quality_metrics[[i]] <- metrics(A, B, C) |> as.data.frame()
    }
    mn <- lmin(mesh)
    sf <- mesh_2_sf(mesh)
    df <- data.frame(ID = 1:nrow(sf),
                     metric_c1 = c1,
                     metric_g1 = g1, metric_b1 = b1)
    df <- cbind(df, do.call('rbind', quality_metrics))
    sf <- dplyr::left_join(sf, df, by = "ID")
    if(plot) print(gg_mems(sf, metric = metric))
    return(sf)
}
#' A range of triangle quality metrics and measures
#'
#' @param A (x,y) coords of vertex A
#' @param B (x,y) coords of vertex B
#' @param C (x,y) coords of vertex C
#' @export
metrics <- function(A, B, C){
    a <-  dist(B, C)
    b <- dist(A, C)
    c <- dist(A, B)
    mn <- min(a, b, c)
    area <- abs(tri_area(a = a, b = b,
                         angC = ang(a, b, c)))
    c_R <- circum_R(A, B, C)
    i_R <- incircle_r(A, B, C)
    c_O <- circum_O(A, B, C)
    i_O <- incircle_O(A, B, C)
    m_1 <- (4 * sqrt(3) * abs(area)) / (sum(c(a, b, c)^2))
    m_2 <- (6*sqrt(area/sqrt(3))/sum(c(a, b, c)))^2
    m_3 <- min(a, b, c)/ max(a, b, c)
    m_4 <- 3 * (min(ang(a, b, c), ang(b, c, a), ang(c, a, b))/(180 / pi))/pi
    m_5 <- (4 * sqrt(3) * abs(area))/ ( max(a, b, c) * sum(c(a, b, c)))
    q_b <- (sum(c(a, b, c)^2))/(4 * sqrt(3) * area)
    q_w <- (1/(3*area))*(sum(c(a, b, c))^2)
    m_6 <- 1/(4*q_b) + q_w*(sqrt(3)/16)
    rr <- i_R/ c_R
    re <- c_R/mn
    return(list(area = area,
                incircle_r = i_R,
                circumcircle_R = c_R,
                c_Ox = c_O[1],  c_Oy = c_O[2],
                i_Ox = i_O[1],  i_Oy = i_O[2],
                angleA = ang(a, b, c),
                angleB = ang(b, c, a),
                angleC = ang(c, a, b),
                rr = rr, re = re,
                metric_1 = m_1, metric_2 = m_2,
                metric_3 = m_3, metric_4 = m_4,
                metric_5 = m_5, metric_6 = m_6))
}
#' Function to calculate third vertex of equilateral triangle
#' given two verticies
#' @inheritParams dist
#' @noRd
equilateral <- function(A, B){
    dX <- B[1] - A[1]
    dY <- B[2] - A[2]
    x3 <- (cos(60 / (180 / pi)) * dX - sin(60 / (180 / pi)) * dY) + A[1]
    y3 <- (sin(60 / (180 / pi)) * dX + cos(60 / (180 / pi)) * dY) + A[2]
    return(c(x3, y3))
}
#' Function to calculate third vertex of right angled triangle
#' given two verticies an height \code{h}
#' @inheritParams dist
#' @param h height
#' @noRd
right_angled <- function(A, B, h = 1){
    dX <- B[1] - A[1]
    dY <- B[2] - A[2]
    n <- dist(A, B)
    x3 <- (h*dY)/n + A[1]
    y3 <- (h*dX)/n + A[2]
    return(c(x3, y3))
}
#' Function to calculate third vertex of isoceles triangle
#' given two verticies an height \code{h}
#' @inheritParams dist
#' @param h height
#' @noRd
isoceles <- function(A, B, h = 1){
    n <- dist(A, B)
    phi <- ang(h, h, n)  / (180 / pi)
    x3 <- (n * cos(phi)) + A[1]
    y3 <- (n *sin(phi)) + A[2]
    return(c(x3, y3))
}
#' Function to get mesh triangle center points of a \code{mesh} object
#' @details Used for plotting label locations
#' @inheritParams mems
#' @export
cens <- function(mesh){
     mid <- t(apply(mesh$graph$tv, 1, function (x)
         incircle_O(mesh$loc[x[1],1:2],mesh$loc[x[2],1:2],mesh$loc[x[3],1:2])))
     df <- data.frame(triangle = paste("T",1:nrow(mid), sep = ""),
                       x = mid[,1], y = mid[,2])
     return(df)
}
#' Function to find the half way point of each
#' mesh edge and turn into an \code{sf} \code{LINESTRING}
#' object
#' @inheritParams mems
half_segments <- function(mesh){
    ## mesh$vv is the node (vertex) association
    ## mesh$tv is the triangle associations
    mm <- mesh$graph$vv
    dp <- diff(mm@p)
    idx <- cbind(mm@i+1,rep(seq_along(dp),dp))
    ## finding mid points from nodes
    nodes <- cbind(mesh$loc[,1], mesh$loc[,2])
    mp <- t(apply(idx, 1, function(x) mid(nodes[x[1], ], nodes[x[2],])))
    df <- data.frame(start_x = mp[,1], start_y = mp[,2],
                     end_x = nodes[idx[,1],1], end_y = nodes[idx[,1],2])
    df$ID <- seq(nrow(df))
    rows <- split(df, seq(nrow(df)))
    lines <- lapply(rows, function(row) {
        lmat <- matrix(unlist(row[1:4]), ncol = 2, byrow = TRUE)
        sf::st_linestring(lmat)
    })
    lines <- sf::st_sfc(lines)
    lines_sf <- sf::st_sf('ID' = df$ID, 'geometry' = lines)
    return(lines_sf)
}
#' Function to find which segments on boundary
#' @inheritParams mems
bound_half_segments <- function(mesh){
    fm <- fmesher::fm_fem(mesh)
    ## mesh verticies
    mm <- mesh$graph$vv
    ## Boundary segments
    tmp <- fm$b1*mm
    ## segments
    dp <- diff(tmp@p)
    idx <- cbind(tmp@i+1,rep(seq_along(dp),dp))
    ## finding mid points from nodes
    nodes <- cbind(mesh$loc[,1], mesh$loc[,2])
    mp <- t(apply(idx, 1, function(x) mid(nodes[x[1], ], nodes[x[2],])))
    df <- data.frame(start_x = mp[,1], start_y = mp[,2],
                     end_x = nodes[idx[,1],1], end_y = nodes[idx[,1],2])
    df$ID <- seq(nrow(df))
    rows <- split(df, seq(nrow(df)))
    lines <- lapply(rows, function(row) {
        lmat <- matrix(unlist(row[1:4]), ncol = 2, byrow = TRUE)
        sf::st_linestring(lmat)
    })
    lines <- sf::st_sfc(lines)
    lines_sf <- sf::st_sf('ID' = df$ID, 'geometry' = lines)
    return(lines_sf)
}
    
    
#' Internal function to construct the dual mesh
#'
#' @param mesh A spatial mesh of class \code{fmesher::fm_mesh_2d()}.
#' @return An\ simple features, code{sf}, object of the Voronoi tessellation
#' centered at each \code{mesh} node.
#' @source \url{http://www.r-inla.org/spde-book},  \url{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @noRd
dual_mesh <- function(mesh) {
    if (mesh$manifold == 'R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        pls <- lapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[, k] == i)
                if (length(j) > 0)
                    return(rbind(ce[j, , drop = FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1],
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]) / 2))
                else return(ce[j, , drop = FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[, 1] == i)
            j2 <- which(mesh$segm$bnd$idx[, 2] == i)
            if ((length(j1) > 0) | (length(j2) > 0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2,
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2))
                yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
                xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
            }
            else {
                yy <- p[, 2] - mesh$loc[i, 2]
                xx <- p[, 1] - mesh$loc[i, 1]
            }
            p <- p[order(atan2(yy, xx)), ]
            ## close polygons
            p <- rbind(p, p[1, ])
            sf::st_polygon(list(p))
        })
        geometry <- sf::st_sfc(pls)
        dat <- data.frame(ID = 1:length(geometry))
        return(sf::st_sf(dat,
                         geometry = geometry))
    }
    else stop("It only works for R2!")
}
