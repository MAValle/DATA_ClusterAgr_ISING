# la funcion hierarchical_clustering_v4 ha tenido varios problemas de ejecucion
# al probarla con distintos J y D aleatorias generadas en simulation_hc_V2.R
# Se crea ahora una algoritmo mejorado de hierarchical_clustering_v4 
# que parte de la base de ir uniendo clusters con el menor numero de nodos
# asumiendo que mientras mas nodos tiene un cluster, la distancia de acople
# siempre es mayor.

# el resultado pretendido es la funcion: hierarchical_clustering_v5 


# El Algoritmo (el cerebro):
# Nota: el elgoritmo completo sigue el single linkage que se encuentra originalmente en hierarchical_clustering_v2
# Previamente se tiene que tener calculado el MST.
# 1. La matriz de distancia D nos indica que hay cluster A y B potencial para unir:
# NOta: A y B puede tener mas de 1 nodo.
      # se calcula la distancia de acople entre A y B =  d_c_original
      # para el cluster {A.B} tentativo, se calcula v (nodos de los clusters), nombres_num (nombre de los clusters)  y m (dist. ultrametrica de D)
# 2. Se calcula n_a numero de nodos de A / n_b = numero de nodos de B
# 3. Para cada cluster A y B:
#         se calcula n_a_cx que es el numero de nodos que tiene el cluster adyacentes a A o B en el MST.
# 4. Elegir el cluster C potencial de fusionarse a A o B segun el MINIMO de n_a_cx
# 5. calcular la distancia de acople entre C y A o B segun el caso = d_c_tentativo
# 6. Elegir fusionar {A,B} si d_c_original < d_c_tentativo
# # # a continuacion parte comun con hierarchical_clustering_v2:
# 7. obtener mn, v, nombres_num y m de los dos clusters que se van a fusionar
# 8. calcular d_mst y d_clp
# 9. agregar informacion a matriz merge 
# 10. seguir proceso de actualizacion de matriz de distancia D.

# notas:
# mn es la lista de los nodos que integran cada cluster fusionado en cada iteracion
# v son todos los nodos involucrados en los clusters que se van a fusionar 
# nombres_num son los dos nombres de los clusters que se van a fusionar
# m es la distancia ultrametrica (d_mst) entre los dos clusters utilizando single linkage.

# actual name: testing_hierarchical_clusteringV9.R
# 27.oct.19 : creation



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cargar funciones 
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R")
#source("network_functions_v1.R") # aqui nos interesa la funcion find_mst_barrier
source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function

source("create_coupling_function.R") # crea matriz J de acoples y D de distancias para simular
source("create_mst_function.R") # crea MST a partir de una matriz de acople.
source("acople_distance_sum_function.R") # suma las distantcias de acople dc
source("find_ady_dc_function.R") # encuentra los nodos adyacentes en un mst dado un nodo
source("get_num_nodes_function.R") # nos dice el numero de nodos que tiene el par de clusters a fusionar
source("get_num_nodes_of_ady_clusters_function.R") # nos dice informacion del numero de nodos que son adyacentes a un cluster
source("get_name_of_the_cluster_function.R") # dado el nombre de un nodo, nos dice a que cluster pertenece
source("get_nodes_of_the_cluster_function.R") # nod da los nodos involucrados en un cluster.
library(igraph)
library(Matrix)
load("test_borrar.RData") # aqui cargamos una matriz D y J que se genero de menra aletoria en 
# simulation_hc_V2.R y que sabemos que nos causa problemas y nor sirve por tanto para depurar
# hierarchical_clustering_v5.
mst_g <- create_mst(J=J, D=D)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


hierarchical_clustering_v5 <- function(D, J, mst_g) {
  # # # Inicializacion de la matriz merge que es el output
  N <- ncol(D)
  merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dmst", "dc")
  iteraciones <- N - 2
  mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 
  it = 1
  
  while (it <= iteraciones) {
    # # # BEGIN CEREBRO DEL ALGORITMO
    # Identificar cluster A y B
    rl <- cluster_find_name(D)
    id <- rl$id # 
    nombres <- rl$nombres  # 
    m <- rl$m  #
    # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
    nombres_num <- as.numeric(nombres) # nombres de los clusters
    
    # calculamos v, nombres_num, m y d_c para cluster original {A,B} (para comparar mas adelante)
    nombres_original <- nombres
    nombres_num_original <- nombres_num
    v_original <- c( get_nodes_of_the_cluster(name_cluster = nombres_num_original[1], lista=mn) , get_nodes_of_the_cluster(name_cluster = nombres_num_original[2], lista=mn)   )
    d_cpl_original <- acople_distance_sum2(J, v_original)
    m_original <- m
    
  
    # calculamos el numero de nodos que posee el cluster A y B (en nombres_num)
    # Nota: numbres_num puede ser: dos cluster con nombres positivos, dos cluster con nombres negativos, 
    #       dos clusters con un nombre positivo y otro negativo. 
    num_nodos <- get_num_nodes(name_cls = nombres_num, lista=mn)
    
    # identificamos los nodos adyacentes en el MST de los clusters A y B y su numero de nodos
    # creamos matriz temporal que tiene: 
    # cluster / node_ady / cluster_ady / num_nodes_clust_ady
    # donde cluster_ady es el nombre del cluster en donde se encuentra node_ady, y num_nodes_clust_ady
    # es el numero de nodos de cluster_ady .
    nodos_adyacentes <- get_num_nodes_of_ady_clusters(name_cls = nombres_num, lista=mn, mst_g = mst_g)
    # como cluster A y B se repetiran en la matriz de "nodos_adyacentes", debemos eliminar ambas filas:
    id1 =  which( ( nodos_adyacentes[,1] == nombres_num_original[1] ) &  ( nodos_adyacentes[,3] == nombres_num_original[2] )  )
    id2 =  which( ( nodos_adyacentes[,1] == nombres_num_original[2] ) & ( nodos_adyacentes[,3] == nombres_num_original[1] )   )
    nodos_adyacentes <- nodos_adyacentes[-c(id1, id2), , drop=FALSE]
    # ahora borramos las adyacencias internas de un mismo cluster
    id3 <- c( which(nodos_adyacentes[, 3] == nombres_num_original[1]),   which(nodos_adyacentes[, 3] == nombres_num_original[2]) )
    if ( ! is.integer0( id3) ) { nodos_adyacentes <- nodos_adyacentes[-c(id3), , drop=FALSE] }

    
    # vamos chequeando si la fusion del cluster original {A,B} es mejor que las tentativas en nodos_adyacentes
    # vamos a meter cluster original {A,B} en temporal:
    meter <- c(nombres_num_original[1], NA, nombres_num_original[2], sum(num_nodos) )
    nodos_adyacentes <- rbind(nodos_adyacentes, meter)
    id <- as.numeric(which( nodos_adyacentes[, 4] == min( nodos_adyacentes[, 4], na.rm = TRUE) ) )
    # borramos las filas de nodos_adyacentes que no son minimos en id
    temp <- (1:nrow(nodos_adyacentes) %in% id) # filas de nodos_adyacentes que no estan en id
    nodos_adyacentes <- nodos_adyacentes[temp, , drop=FALSE]
    # ahora tenemos que calcular distancia de acople dc a cada potencial cluster en "nodos_adyacentes"
    coupling_distances <- c()
    for (j in 1:length(id) ) {
      cluster_temporal <- as.numeric( nodos_adyacentes[j, c(1,3)] )
      v_temporal <- c( get_nodes_of_the_cluster(name_cluster = cluster_temporal[1], lista=mn) ,  
                       get_nodes_of_the_cluster(name_cluster = cluster_temporal[2], lista=mn)  )
      d_cpl_temporal <- acople_distance_sum2(J, v_temporal)
      coupling_distances <- c(coupling_distances, d_cpl_temporal)
    }
    # Ahora ponemos el vector "coupling_distances" como nueva columna en "nodos_adyacentes"
    nodos_adyacentes <- cbind(nodos_adyacentes, coupling_distances)
    id <- which.min(nodos_adyacentes[, 5])
    
    # final para seguir adelante:
    nombres_num <- c(nodos_adyacentes[id, 1], nodos_adyacentes[id, 3])
    nombres <- as.character(nombres_num)
    v <- c( get_nodes_of_the_cluster(name_cluster = nombres_num[1], lista=mn)  ,  
            get_nodes_of_the_cluster(name_cluster = nombres_num[2], lista=mn)  )
    mn <- c(mn, list( as.character(v) ) )
    m <- D[which(colnames(D)==nombres_num[1]), which(colnames(D)==nombres_num[2])]
    # actualizacion del id, necesario para borrar las filas y columnas en D
    id <- c( which(colnames(D) == nombres_num[1]) ,   which(colnames(D) == nombres_num[2])   )
    # # # END CEREBRO DEL ALGORITMO
    
    # # PUTTING INFO IN merge ANd UPGRADING DISTANCE MATRIX D
    # calculo de la distancia MST
    d_mst <- find_mst_barrier(mst_g, v)
    # calculo de las distancias de acople
    d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
    
    merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
    
    # # # # # # # # # Upgrading the distance matrix
    # fila id[1] y columna id[2] se deben borrar
    Dold <- D
    D <- D[-id, -id]
    
    # pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
    if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
      colno <- colnames(Dold) # nombres de columnas o fila de Dold
      #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
      #sel <- colno[!(colno %in% nombres)] # nodo identificado
      sel <- colno[!(colno %in% as.character(nombres_num) )] # nodo identificado 16-oct-19
      col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
      dit <- vector()
      dit <- c(dit,min(Dold[, col])) # la distancia minima
      nombres_ <- sel
    } else {
      #temp <- find_min_dist(D, Dold, nombres)
      temp <- find_min_dist(D = D, Dold = Dold, nombres = as.character(nombres_num)) # 16-oct-19
      dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
      nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
    }
    
    #ahora tengo que poner los valores en fila.columna de D.
    D <- put_dis(D, dit, nombres_, N, it)
    print(paste("Iteration: ", it))
    
    it <- it + 1
    
  } ### FIN DE LAS ITERACIONES
  
  
  # # # # LAST ITERATION # # # # 
  # Last iteration: when it = N - 1
  #it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres) # nuevo
  clusters_con_mas_de_un_nodo <- nombres_num[nombres_num > 0]
  if (length(clusters_con_mas_de_un_nodo ) > 1) { # En caso que se vayan a fusionar dos clusters que a su vez contienen varios nodos cada uno.
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(temp,  unlist(   mn[ clusters_con_mas_de_un_nodo[n] ] )  )
    }
    v <- temp
  } else {
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(nombres_num[nombres_num < 1],  unlist(   mn[ clusters_con_mas_de_un_nodo ] )  )
    }
    #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
    v <- temp
  }
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, v)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
  merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) 
  print(paste("Iteration: ", it))
  # retorno de los resultados matriz merge
  return(merge)

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 







# BORAR DESPUES
# 24-10-19   tratando de ver que pasa con V4 de hirarchical function.
#rm(list=setdiff(ls(), c("J", "D")))
rep = 100
Nn = 20
for (i in 1:rep) {
  print(paste("starting iteration", i))
  ot <- create_coupling(Nn=Nn, media=0, sj=1, lw=-3, up=3)
  J <- ot$J
  D <- ot$D
  mst_g <- create_mst(J=J, D=D)
  hierarchical_clustering_v5(D, J=J, mst_g)
  print(paste("ending iteration", i))
}




