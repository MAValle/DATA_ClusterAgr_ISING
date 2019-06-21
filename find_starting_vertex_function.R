# # # # # # # # # # # #  # # # # # # FIND strongest vertex # # # # # # # # ## #  # # # # # # #
# June 11, 2019
# Input: igraph object the mst
# output: the name of teh starting vertex in numeric format
find_starting_vertex <- function(mst_g) {
  # strengh de los nodos en el MST
  #st <- strength(mst_g, weights = E(mst_g)$coupling)
  st <- strength(mst_g, weights = E(mst_g)$weight)
  id <- which(st == max(st))
  v <- as.numeric(as_ids(V(mst_g)[id]))
  v <- c(v)
  return(v)
}
# example:
#v <- find_starting_vertex(mst_g)
# Result: node 80
# # # # # # # # # # # #  # # # # # # FIND strongest vertex # # # # # # # # ## #  # # # # # # #
