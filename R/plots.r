#' Plots the calculated metrics from a \code{mems} output
#' (i.e., each triangle metric).
#' @param x output of \code{mems::mems}
#' @param metric character specifying which metric of \code{x} to plot one of "metric_1",
#' "metric_2", "metric_3", "metric_4", "metric_5", "metric_6"  (default "metric_1").
#' @return A  \code{"ggplot"} object.
#' @export
gg_mems <- function(x, metric = "metric_1"){
    opts <- c("metric_1", "metric_2", "metric_3", "metric_4", "metric_5", "metric_6" )
    if(!(metric %in% opts)) stop("Please provide a valide metric")
    if(metric == "metric_1"){
        m1 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_1)) +
            ggplot2::labs(title = expression(paste("Metric = ",
                                                   frac(4*sqrt(3)*abs(A),Sigma[i = 1]^3 *l[i]^2))))
        m1 + ggplot2::theme_void()
    }else{ if(metric == "metric_2") {
               m2 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_2))+
                   ggplot2::labs(title = expression(paste("Metric = ",
                                                          6*sqrt(frac(A, sqrt(3))))))
               m2 + ggplot2::theme_void()
           }else{if(metric == "metric_3") {
                     m3 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_3)) +
                         ggplot2::labs(title = expression(paste("Metric = ",
                                                                frac("min("~underline(bold(l))~")",
                                                                     "max("~underline(bold(l))~")"))))
                     m3 + ggplot2::theme_void()
                 }else{if(metric == "metric_4") {
                           m4 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_4)) +
                               ggplot2::labs(title = expression(paste("Metric = ",
                                                                      3 * frac("min("~underline(bold(phi))~")",
                                                                               pi))))
                           m4 + ggplot2::theme_void()
                       }else{if(metric == "metric_5") {
                                 m5 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_5)) +
                                     ggplot2::labs(title = expression(paste("Metric = ",
                                                                            frac(4 *sqrt(3) * A,
                                                                                 "max("~underline(bold(l))~")") * Sigma[i = 1]^3 * l[i]^2)))
                                 m5 + ggplot2::theme_void()
                             }else{if(metric == "metric_6") {
                                       m6 <- ggplot2::ggplot(x) + ggplot2::geom_sf(ggplot2::aes(fill = metric_6))  +
                                           ggplot2::labs(title = expression(paste("Metric = ",
                                                                                  frac(1, 4*q[b]) + q[w] * frac(sqrt(3),16))))
                                       m6 + ggplot2::theme_void()
                                   }
                             }
                       }
                 }
           }
    }
   
}
#' Plots Q blocks contribution of nodes
#' only useful for small number of nodes
#' @inheritParams mems
#' @param element One of "C1", "G1", "B1", corresponding matricies from \code{fmesher::fm_fem()},
#' default "C1".
#' @export
gg_qblocks <- function(mesh, element = "C1"){
    opts <- c("C1", "G1", "B1")
    if(!(element %in% opts)) stop("element should be one of C1, G1, or B1")
    fm <- fmesher::fm_fem(mesh) ## mems q_block to test
    nv <- mesh$n
    ## initialise matricies
    matc1 <- matg1 <- matb1 <- matrix(NA, ncol = nv, nrow = nv)
    if(element == "C1"){
        ## C1
        matc1[cbind((fm$c1@i + 1), (fm$c1@j + 1))] <- fm$c1@x 
        c1.df <- reshape2::melt(matc1, c("x", "y"), value.name = "z")
        plt <- ggplot(data = c1.df, aes(x = y,y = -x,fill = z)) +
            geom_tile(size = 1.5, alpha = 0.8) + 
            scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                               expand = c(0, 0), position = "top") +
            scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
            labs(fill = "C1") + xlab("") + ylab("") +
            geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
            ggtitle("C1")
    }else{if(element == "G1"){
              ## G1
              matg1[cbind((fm$g1@i + 1), (fm$g1@j + 1))] <- fm$g1@x 
              g1.df <- reshape2::melt(matg1, c("x", "y"), value.name = "z")
              plt <- ggplot(data = g1.df, aes(x = y,y = -x,fill = z)) +
                  geom_tile(size = 1.5, alpha = 0.8) + 
                  scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                                     expand = c(0, 0), position = "top") +
                  scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
                  labs(fill = "G1") + xlab("") + ylab("") +
                  geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
                  ggtitle("G1")
          }else{if(element == "B1"){
                    ## B1
                    matb1[cbind((fm$b1@i + 1), (fm$b1@j + 1))] <- fm$b1@x 
                    b1.df <- reshape2::melt(matb1, c("x", "y"), value.name = "z")
                    plt <- ggplot(data = b1.df, aes(x = y,y = -x,fill = z)) +
                        geom_tile(size = 1.5, alpha = 0.8) + 
                        scale_x_continuous(breaks = seq(1, nv, 1), labels = paste("n",1:nv,sep = ""),
                                           expand = c(0, 0), position = "top") +
                        scale_y_continuous(breaks = -nv:-1,labels = paste("n",nv:1,sep = ""), expand = c(0, 0)) +
                        labs(fill = "B1") + xlab("") + ylab("") +
                        geom_text( aes(label = round(z, 3)), color="black", size=rel(4)) +
                        ggtitle("B1")
                }
          }
    }
    ## plot 
    plt +
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
#' @inheritParams gg_qblocks
#' @export
gg_seg_influence <- function(mesh, element = "C1"){
    opts <- c("C1", "G1", "B1")
    if(!(element %in% opts)) stop("element should be one of C1, G1, or B1")
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
    vor <- dual_mesh(mesh)
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
    if(element == "C1") {
        ## dataframe C1
        dfc1 <- data.frame(row = (fm$c1@i + 1), col = (fm$c1@j + 1), val = fm$c1@x)
        ## remove off diag elements
        dfc1 <- dfc1[-which(dfc1$row - dfc1$col == 0),]
        plt <- ggplot() + geom_sf(data = vor,fill = NA,
                                  linetype = 2, alpha = 0.3) +
            geom_sf(data = segs, aes(col = dfc1$val), linewidth = 3, alpha = 0.7) +
            geom_text(data = cents, aes(x = X, y = Y, label = round(dfc1$val,3))) + 
            geom_point(data = nodes, aes(x = x, y = y, col = valc1), size = 10,) +
            geom_text(data = nodes, aes(x = x, y = y, label = round(valc1, 3))) +
            theme_void() + ggtitle("C1") + theme(legend.position = "none")
    }else{if(element == "G1"){
              ## dataframe G1
              dfg1 <- data.frame(row = (fm$g1@i + 1), col = (fm$g1@j + 1), val = fm$g1@x)
              ## remove off diag elements
              dfg1 <- dfg1[-which(dfg1$row - dfg1$col == 0),]
              plt <- ggplot() + geom_sf(data = vor,fill = NA,
                                        linetype = 2, alpha = 0.3) +
                  geom_sf(data = segs, aes(col = dfg1$val), linewidth = 3, alpha = 0.7) +
                  geom_text(data = cents, aes(x = X, y = Y, label = round(dfg1$val,3))) + 
                  geom_point(data = nodes, aes(x = x, y = y, col = valg1), size = 10,) +
                  geom_text(data = nodes, aes(x = x, y = y, label = round(valg1, 3))) +
                  theme_void() + ggtitle("G1") + theme(legend.position = "none")
          }else{if(element == "B1"){
                    ## dataframe B1
                    dfb1 <- data.frame(row = (fm$b1@i + 1), col = (fm$b1@j + 1), val = fm$b1@x)
                    ## remove diag elements
                    dfb1 <- dfb1[-which(dfb1$row - dfb1$col == 0),]
                    plt <- ggplot() + geom_sf(data = vor,fill = NA,
                                              linetype = 2, alpha = 0.3) +
                        geom_sf(data = mesh_2_sf(mesh), fill = NA, col = "grey") +
                        geom_point(data = nodes, aes(x = x, y = y), col = "grey") +
                        geom_sf(data = bound_segs, aes(col = dfb1$val), linewidth = 3, alpha = 0.7) +
                        geom_text(data = centsb1, aes(x = X, y = Y, label = round(dfb1$val,3))) + 
                        geom_point(data = nodes[Matrix::diag(fm$b1) != 0, ], aes(x = x, y = y, col = valb1), size = 10) +
                        geom_text(data = nodes[Matrix::diag(fm$b1) != 0, ], aes(x = x, y = y, label = round(valb1, 3))) +
                        theme_void() + ggtitle("B1") + theme(legend.position = "none")
                }
          }
    }
    
    plt +
        scale_fill_continuous(type = "viridis", na.value = NA)
}
#' Returns x and y coords of circle
circle <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
#' Plot triangle with all the bells
#' internal function only
gg_tri <- function(A, B, C, node_labels = c("A", "B", "C"),
                   edge_labels = c("a", "b", "c"), r = c(0.2, 0.2, 0.2),
                   ang_labels = c(expression(theta[a]), expression(theta[b]), expression(theta[c])),
                   circum  = FALSE, incircle = FALSE){
    tri <- data.frame(rbind(A, A, B, B, C, C), rbind(B, C, C, A, A, B))
    names(tri) <- c("x", "y","xend", "yend"); rownames(tri) <- NULL
    and <- with(tri, atan2(xend - x, yend - y))
    aA <- and[1:2]; aB <- and[3:4]; aC <- and[5:6]
    aA <-  ifelse(aA < 0, aA + 2 * pi, aA)
    aB <-  ifelse(aB < 0, aB + 2 * pi, aB)
    aC <-  ifelse(aC < 0, aC + 2 * pi, aC)
    ## direction arc centered at node A
    aA <- ifelse(abs(aA - aB) > pi & aA < aB, aA + 2 * pi, aA)
    ## ## direction arc centered at node B
    aB <- ifelse(abs(aB - aA) > pi & aB < aA, aB + 2 * pi, aB)
    ## ## direction arc centered at node C
    aC <- ifelse(abs(aC - aA) > pi & aC < aA, aC + 2 * pi, aC)
    edge  <- data.frame(edge = edge_labels)
    mids <- rbind(mems::mid(B, C), mems::mid(A, C), mems::mid(A, B))
    edge$x <- mids[,1]; edge$y <- mids[,2]
    plt <- ggplot(tri) +
        geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = 2) +
        theme_void() +
        ggforce::geom_arc(aes(x0 = A[1],y0 = A[2], start = aA[1], end = aA[2], r = r[1])) + 
        ggforce::geom_arc(aes(x0 = B[1],y0 = B[2], start = aB[1], end = aB[2], r = r[2])) + 
        ggforce::geom_arc(aes(x0 = C[1],y0 = C[2], start = aC[1], end = aC[2], r = r[3])) + 
        coord_fixed() + geom_point(data = unique(tri[,1:2]), aes(x = x, y = y), size = 7) +
        geom_text(data = tri[c(1, 3, 6),1:2], aes(x = x, y = y), label = node_labels, col = "white") +
        ggrepel::geom_text_repel(data = edge, aes(x = x , y = y ,
                                   label = edge), fontface = "bold")
    p <- ggplot_build(plt)
    m <- sapply(c(2, 3, 4), function(h) p$data[[h]] %>% dplyr::summarize(x = mean(x), y = mean(y))) |> t()
    m <- data.frame(x = unlist(m[,1]), y = unlist(m[,2]))
    mids_theta <- rbind(mems::mid(A, m[1,]), mems::mid(B, m[2,]), mems::mid(C, m[3,]))
    thetas <- data.frame(x = unlist(mids_theta[,1]), y = unlist(mids_theta[,2]))
    plt <- plt + geom_text(data = thetas, aes(x = x , y = y), 
                           label = ang_labels)
    if(circum | incircle){
        mets <- mems:::metrics(A, B, C) |> as.data.frame()
        if(circum){
            cc <- mems:::circle(c(mets$c_Ox, mets$c_Oy), mets$circumcircle_R*2, npoints = 100)
            plt <- plt + geom_path(data = cc, aes(x, y))  +
                geom_point(data = mets, aes(x = c_Ox, y = c_Oy), pch = 18, size = 3) +
                ggrepel::geom_text_repel(data = mets, aes(x = c_Ox, y = c_Oy, label = "C")) +
                geom_segment(data = mets, aes(x = c_Ox, y = c_Oy, xend = c_Ox,
                                              yend = c_Oy - circumcircle_R), linetype = 2) +
                ggrepel::geom_text_repel(data = mets, aes(x = c_Ox, y = c_Oy - circumcircle_R/2, label = "R"))
        }
        if(incircle){
            ic <- mems:::circle(c(mets$i_Ox, mets$i_Oy), mets$incircle_r*2, npoints = 100)
            plt <- plt +  geom_path(data = ic, aes(x, y), col = "grey") +
                geom_point(data = mets, aes(x = i_Ox, y = i_Oy), col = "grey") +
                ggrepel::geom_text_repel(data = mets, aes(x = i_Ox, y = i_Oy, label = "c"), col = "grey") +
                geom_segment(data = mets, aes(x = i_Ox, y = i_Oy, xend = i_Ox,
                                              yend = i_Oy - incircle_r), col = "grey", linetype = 2) +
                ggrepel::geom_text_repel(data = mets, aes(x = i_Ox, y = i_Oy - incircle_r/2, label = "r"),
                                         col = "grey")
        }
    }
   print(plt)     
}
