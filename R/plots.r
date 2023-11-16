#' Plots the calculated metrics from a \code{mems} output
#' (i.e., each triangle metric).
#' @param x output of \code{mems::mems}
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
#' only useful for small number of nodes
#' @inheritParams mems
#' @export
gg_qblocks <- function(mesh){
    fm <- fmesher::fm_fem(mesh) ## mems q_block to test
    nv <- mesh$n
    ## initialise matricies
    matc1 <- matg1 <- matb1 <- matrix(NA, ncol = nv, nrow = nv)
    ## C1
    matc1[cbind((fm$c1@i + 1), (fm$c1@j + 1))] <- fm$c1@x 
    c1.df <- reshape2::melt(matc1, c("x", "y"), value.name = "z")
    plt_C1 <- ggplot(data = c1.df, aes(x = y,y = -x,fill = z)) +
        geom_tile(size = 1.5, alpha = 0.8) + 
        scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                           expand = c(0, 0), position = "top") +
        scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
        labs(fill = "C1") + xlab("") + ylab("") +
        geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
        ggtitle("C1")
    ## G1
    matg1[cbind((fm$g1@i + 1), (fm$g1@j + 1))] <- fm$g1@x 
    g1.df <- reshape2::melt(matg1, c("x", "y"), value.name = "z")
    plt_G1 <- ggplot(data = g1.df, aes(x = y,y = -x,fill = z)) +
        geom_tile(size = 1.5, alpha = 0.8) + 
        scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                           expand = c(0, 0), position = "top") +
        scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
        labs(fill = "G1") + xlab("") + ylab("") +
        geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
        ggtitle("G1")

    ## B1
    matb1[cbind((fm$b1@i + 1), (fm$b1@j + 1))] <- fm$b1@x 
    b1.df <- reshape2::melt(matb1, c("x", "y"), value.name = "z")
    plt_B1 <- ggplot(data = b1.df, aes(x = y,y = -x,fill = z)) +
        geom_tile(size = 1.5, alpha = 0.8) + 
        scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                           expand = c(0, 0), position = "top") +
        scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
        labs(fill = "B1") + xlab("") + ylab("") +
        geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
        ggtitle("B1")

    ## plot all matricies
    patchwork::wrap_plots(plt_C1, plt_G1,plt_B1, nrow = 1) & 
        scale_fill_continuous(type = "viridis", na.value = NA) &
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA,color = "gray90", size = 0.5, linetype = "solid"),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill="gray90"),
              plot.background = element_rect(fill="gray90"),
              legend.position = "none", 
              axis.text = element_text(color="black", size=10),
              plot.title = element_text(size = 14))
}
#' Plot segment influence
#' only useful for small number of nodes
#' @inheritParams mems
#' @export
gg_seg_influence <- function(mesh){
    nv <- mesh$n
    nodes <- data.frame(x = mesh$loc[,1], y = mesh$loc[,2])
    fm <- fmesher::fm_fem(mesh)
    ## diag elements
    nodes$valc1 <- Matrix::diag(fm$c1)
    nodes$valg1 <- Matrix::diag(fm$g1)
    nodes$valb1 <- Matrix::diag(fm$b1)
    segs <- half_segments(mesh)
    bound_segs <- bound_half_segments(mesh)
    ## Voronoi
    tesselation <- deldir::deldir(nodes$x, nodes$y)
    tiles <- deldir::tile.list(tesselation) 
    vor <- deldir_2_sf(tiles)
    ## label points for line segments
    crd <- sf::st_coordinates(segs)
    lst <- split.data.frame(crd[, 1:2], crd[,3])
    cens <- lapply(lst, function(x) mid(x[1,], x[2,]))
    cents <- as.data.frame(do.call('rbind', cens))
    ## label points for boundary line segments
    crd <- sf::st_coordinates(bound_segs)
    lst <- split.data.frame(crd[, 1:2], crd[,3])
    cens <- lapply(lst, function(x) mid(x[1,], x[2,]))
    centsb1 <- as.data.frame(do.call('rbind', cens))
    ## dataframe C1
    dfc1 <- data.frame(row = (fm$c1@i + 1), col = (fm$c1@j + 1), val = fm$c1@x)
    ## remove off diag elements
    dfc1 <- dfc1[-which(dfc1$row - dfc1$col == 0),]
    plt_C1 <- ggplot() + geom_sf(data = vor,fill = NA,
                                 linetype = 2, alpha = 0.3) +
        geom_sf(data = segs, aes(col = dfc1$val), linewidth = 3, alpha = 0.7) +
        geom_text(data = cents, aes(x = X, y = Y, label = round(dfc1$val,3))) + 
        geom_point(data = nodes, aes(x = x, y = y, col = valc1), size = 10,) +
        geom_text(data = nodes, aes(x = x, y = y, label = round(valc1, 3))) +
        theme_void() + ggtitle("C1") + theme(legend.position = "none")
    ## dataframe G1
    dfg1 <- data.frame(row = (fm$g1@i + 1), col = (fm$g1@j + 1), val = fm$g1@x)
    ## remove off diag elements
    dfg1 <- dfg1[-which(dfg1$row - dfg1$col == 0),]
    plt_G1 <- ggplot() + geom_sf(data = vor,fill = NA,
                                 linetype = 2, alpha = 0.3) +
        geom_sf(data = segs, aes(col = dfg1$val), linewidth = 3, alpha = 0.7) +
        geom_text(data = cents, aes(x = X, y = Y, label = round(dfg1$val,3))) + 
        geom_point(data = nodes, aes(x = x, y = y, col = valg1), size = 10,) +
        geom_text(data = nodes, aes(x = x, y = y, label = round(valg1, 3))) +
        theme_void() + ggtitle("G1") + theme(legend.position = "none")
    ## dataframe B1
    dfb1 <- data.frame(row = (fm$b1@i + 1), col = (fm$b1@j + 1), val = fm$b1@x)
    ## remove diag elements
    dfb1 <- dfb1[-which(dfb1$row - dfb1$col == 0),]
    plt_B1 <- ggplot() + geom_sf(data = vor,fill = NA,
                                 linetype = 2, alpha = 0.3) +
        geom_sf(data = mesh_2_sf(mesh), fill = NA, col = "grey") +
        geom_point(data = nodes, aes(x = x, y = y), col = "grey") +
        geom_sf(data = bound_segs, aes(col = dfb1$val), linewidth = 3, alpha = 0.7) +
        geom_text(data = centsb1, aes(x = X, y = Y, label = round(dfb1$val,3))) + 
        geom_point(data = nodes[Matrix::diag(fm$b1) != 0, ], aes(x = x, y = y, col = valb1), size = 10) +
        geom_text(data = nodes[Matrix::diag(fm$b1) != 0, ], aes(x = x, y = y, label = round(valb1, 3))) +
        theme_void() + ggtitle("B1") + theme(legend.position = "none")
     patchwork::wrap_plots(plt_C1, plt_G1, plt_B1, nrow = 1) & 
        scale_fill_continuous(type = "viridis", na.value = NA)
}
#' Plot circle
circle <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
