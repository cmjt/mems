## mems

[![codecov](https://codecov.io/gh/cmjt/mems/graph/badge.svg?token=DRFASZHYAU)](https://codecov.io/gh/cmjt/mems)

[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/joethorley/stability-badges#experimental)

To install run

```r
devtools::install_github("mems")
```

## mesh triangles

![](https://github.com/cmjt/mems/blob/main/docs/example_mesh_attributes.png?raw=true)

Illustration of mesh (Delaunay triangulation) elements. Each triangle (T1, ..., T4) is formed by the connection of three nodes (i.e., three of ni, i = 1,...,6) via an edge. Each node *influence* can be thought of as proportional to the area of the Voronoi cell (coloured polygon) centred at that node. The *influence* of each triangle (T1, ..., T4) on each node (n1, ..., n6) is illustrated by each coloured line segment (i.e., the *length* of the edged connected to each node within the surrounding Voronoi cell).

```r
require(mems)
data(example_mesh, package = "mems")
mesh <- example_mesh
require(deldir)
require(sf)
require(ggplot2)
nodes <- data.frame(x = mesh$loc[,1], y = mesh$loc[,2])
tesselation <- deldir(nodes$x, nodes$y)
tiles <- tile.list(tesselation) ## voronoi diagram (extends outside boundary)
vor <- deldir_2_sf(tiles)
sf <- mesh_2_sf(mesh)
lin_sf <- half_segments(mesh)
col_vor <- RColorBrewer::brewer.pal(6,"Dark2") 
col_lin_sf <- col_vor[(st_intersection(vor, lin_sf) %>%
    subset(., st_geometry_type(st_geometry(.)) == "LINESTRING"))$id]
## main plot
ggplot() +  geom_sf(data = deldir_2_sf(tiles),fill = col_vor,
                    linetype = 2, alpha = 0.3) +
    geom_sf(data = sf, fill = NA, linewidth = 1) +
    geom_sf(data = lin_sf, col = col_lin_sf, linewidth = 2, alpha = 0.5) + 
    theme_void() + geom_text(data = cens(mesh), aes(x = x, y = y, label = triangle),
                             size = 7) +
    geom_point(data = nodes, aes(x = x, y = y), size = 10, col = col_vor) +
    geom_text(data = nodes, aes(x = x, y = y, label = paste("n", 1:nrow(nodes))))
```		   