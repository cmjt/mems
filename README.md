## mems

Trying to make sense of mesh *influence*.

[![codecov](https://codecov.io/gh/cmjt/mems/graph/badge.svg?token=DRFASZHYAU)](https://codecov.io/gh/cmjt/mems)
[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/joethorley/stability-badges#experimental)

To install run

```r
devtools::install_github("mems")
```

## mesh triangles

Below is an illustration of mesh (Delaunay triangulation) elements. Each triangle (T1, ..., T4) is formed by the connection of three nodes (i.e., three of ni, i = 1,...,6) via an edge. Each node *influence* can be thought of as proportional to the area of the Voronoi cell (coloured polygon) centred at that node.  In fact, the diagonal elemens of **C0**, **C_{ii}** are equal to the sum of one third the area of each triangle that the node *i* is part of. The *influence* of each triangle (T1, ..., T4) on each node (n1, ..., n6) is illustrated by each coloured line segment (i.e., the *length* of the edged connected to each node within the surrounding Voronoi cell).


![](https://github.com/cmjt/mems/blob/main/docs/example_mesh_attributes.png?raw=true)

<details>
  <summary>Plot code</summary>
  
```r
## packages & data
require(sf)
require(ggplot2)
require(mems)
data(example_mesh, package = "mems")
mesh <- example_mesh
## manipulation
nodes <- data.frame(x = mesh$loc[,1], y = mesh$loc[,2])
vor <- mems:::dual_mesh(mesh)
sf <- mesh_2_sf(mesh)
lin_sf <- half_segments(mesh)
col_vor <- RColorBrewer::brewer.pal(6,"Dark2") 
col_lin_sf <- col_vor[(st_intersection(vor, lin_sf) %>%
    subset(., st_geometry_type(st_geometry(.)) == "LINESTRING"))$ID]
## plot
ggplot() +  geom_sf(data = vor,fill = col_vor,
                    linetype = 2, alpha = 0.3) +
    geom_sf(data = sf, fill = NA, linewidth = 1) +
    geom_sf(data = lin_sf, col = col_lin_sf, linewidth = 2, alpha = 0.5) + 
    theme_void() + geom_text(data = cens(mesh), aes(x = x, y = y,
    label = triangle),size = 7) +
    geom_point(data = nodes, aes(x = x, y = y), size = 10, col = col_vor) +
    geom_text(data = nodes, aes(x = x, y = y, label = paste("n", 1:nrow(nodes))))
```

</details>

## Triangle circumcircles and incircles

Below shows the circumcircles and incircles of each triangle in a mesh (Delaunay triangulation). Each triangle (T1, ..., T4) is formed by the connection of three nodes (i.e., three of ni, i = 1,...,6) via an edge. A circumcircle (ircumscribed circle) of a triangle is a circle that passes through all three of its nodes (vertices). An incircle (inscribed circle) is a circle that is tangent to each of the traingles's sides (edges).


![](https://github.com/cmjt/mems/blob/main/docs/circles.png?raw=true)

<details>
  <summary>Plot code</summary>

```r
mesh <- example_mesh
mem <- mems(mesh)
nodes <- data.frame(x = mesh$loc[,1], y = mesh$loc[,2])
npoints <- 100
## circumcircle
circum <- apply(cbind(mem$c_Ox, mem$c_Oy, mem$circumcircle_R), 1,
                function(x) circle(c(x[1], x[2]), x[3]*2, npoints = npoints))
circum <- do.call('rbind', circum)
circum$id <- rep(1:nrow(mem), each = npoints)
## incircle
incir <- apply(cbind(mem$i_Ox, mem$i_Oy, mem$incircle_r), 1,
                function(x) circle(c(x[1], x[2]), x[3]*2, npoints = npoints))
incir <- do.call('rbind', incir)
incir$id <- rep(1:nrow(mem), each = npoints)
## plot
ggplot(mem) + geom_sf(fill = NA, linewidth = 2, col = "grey", alpha = 0.7) +
    theme_void() +
    geom_segment(aes(x = c_Ox, y = c_Oy, xend = c_Ox + circumcircle_R,
                     yend = c_Oy, col = as.character(1:4))) +
    geom_segment(aes(x = i_Ox, y = i_Oy, xend = i_Ox + incircle_r,
                     yend = i_Oy,  col = as.character(1:4))) +
    geom_point(aes(x = c_Ox, y = c_Oy, col = as.character(1:4)), pch = 18, size = 3) +
    geom_point(aes(x = i_Ox, y = i_Oy, col = as.character(1:4))) +
    geom_path(data = circum, aes(x, y, group = id, col = as.character(id)))  +
    geom_path(data = incir, aes(x, y, group = id, col = as.character(id))) +
    theme(legend.position = "none") +
    geom_text(data = cens(mesh), aes(x = x, y = y,
    label = triangle, col = as.character(1:4)),size = 7) +
    geom_text(data = nodes, aes(x = x, y = y, label = paste("n", 1:nrow(nodes))), alpha = 0.8) +
    scale_color_manual(values =  RColorBrewer::brewer.pal(4, "Dark2") )
```

</details>