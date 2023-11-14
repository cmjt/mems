#' Plots the calculated metrics from a \code{mems} output
#' (i.e., each triangle metric).
#' @return A  \code{"patchwork"}, \code{"gg"}, \code{"ggplot"} object.
#' @export
gg_mems <- function(x){
    m1 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_1)) +
        ggplot2::labs(title = expression(paste("Metric = ",
                                               frac(4*sqrt(3)*abs(A),Sigma[i = 1]^3 *l[i]^2))))
    m2 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_2))+
        ggplot2::labs(title = expression(paste("Metric = ",
                                               6*sqrt(frac(A, sqrt(3))))))
    m3 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_3)) +
        ggplot2::labs(title = expression(paste("Metric = ",
                                               frac("min("~underline(bold(l))~")",
                                                    "max("~underline(bold(l))~")"))))
    m4 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_4)) +
        ggplot2::labs(title = expression(paste("Metric = ",
                                               3 * frac("min("~underline(bold(phi))~")", pi))))
    m5 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_5)) +
        ggplot2::labs(title = expression(paste("Metric = ",
                                               frac(4 *sqrt(3) * A,
                                                    "max("~underline(bold(l))~")") * Sigma[i = 1]^3 * l[i]^2)))
    m6 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_6))  +
        ggplot2::labs(title = expression(paste("Metric = ",
                                              frac(1, 4*q[b]) + q[w] * frac(sqrt(3),16))))
    patchwork::wrap_plots(m1, m2, m3, m4, m5, m6, nrow = 3) &
        ggplot2::theme_void()
}
#' Plots Q blocks contribution of nodes
gg_qblocks <- function(mesh, sf = TRUE){
    sf <- mesh_2_sf(mesh)
    qb <- q_blocks(mesh)
    df <- data.frame(ID = 1:nrow(sf), area = qb$triangle_area)
    sf <- dplyr::left_join(sf, df, by = "ID")
    mids <- cens(mesh) ## center of triangles
    plt_C0 <- ggplot2::ggplot(sf) + ggplot2::geom_sf(ggplot2::aes(fill = area)) + theme_void() +
        ggplot2::geom_text(data = mids, ggplot2::aes(x = x, y = y, label = triangle), inherit.aes = FALSE) +
        geom_point(data = data.frame(x = mesh$loc[ , 1], y = mesh$loc[ , 2]),
                   ggplot2::aes(x = x, y = y, size = diag(as.matrix(qb$C0)))) +
        labs(fill = "Area", size = "C0")
    plt_G1 <- ggplot2::ggplot(sf) + ggplot2::geom_sf(ggplot2::aes(fill = area)) + theme_void() +
        ggplot2::geom_text(data = mids, ggplot2::aes(x = x, y = y, label = triangle), inherit.aes = FALSE) +
        geom_point(data = data.frame(x = mesh$loc[ , 1], y = mesh$loc[ , 2]),
                   ggplot2::aes(x = x, y = y, size = rowSums(abs(as.matrix(qb$G1))))) +
        labs(fill = "Area", size = "G1")
    plt_B1 <- ggplot2::ggplot(sf) + ggplot2::geom_sf(ggplot2::aes(fill = area)) + theme_void() +
        ggplot2::geom_text(data = mids, ggplot2::aes(x = x, y = y, label = triangle), inherit.aes = FALSE) +
        geom_point(data = data.frame(x = mesh$loc[ , 1], y = mesh$loc[ , 2]),
                   ggplot2::aes(x = x, y = y, size = rowSums(abs(as.matrix(qb$B1))))) +
        labs(fill = "Area", size = "B1")
     patchwork::wrap_plots(plt_C0, plt_G1,plt_B1, nrow = 3)
}
