### Project Cluster Aggregation - DATA_ClusterAgr

Utilizamos ISING inference para descubrir los acoples de una red de productos. Luego encontramos un MST. Luego podemos darnos un punto de partida (nodo) cualquier y comenzar a agregar productos de acuerdo a la energía de barrera dada por el MST y la energía que minimiza la energía de acoples.

Al respecto, el paper de conferencia para ENICC se tienen detalles.



#### Lista de códigos relevantes

**searching_v1.R**: Lo que hacemos aqui es ir buscando las energias de acople Ec y las barreras de energia Emst, de tal forma de minimizar alphaX*Ec + (1-alpha)*Emst.

Algoritmo:

\1. Se selecciona un vertice o nodo (o producto) de partida S = {nodo inicial} 
nota: el nodo inicial puede ser aquel con el mayor strength en el MST.

\2. subrutina energy_coupling_search:

​	 2.1 Para cada uno de los N-1 nodos restantes not in S, calcular la energia de acople y tener el vector de  	energias de acoples para cada agregacion de producto. (temp_Ec)

\3. subrutina de energy_mst_search:

​         3.1 Para cada uno de los N-1 nodos restantos not in S, calcular la energia de barreras o distancia a lo largo del MST que una S con el nodo not in S. Obtener el vector de energias mst para cada uno de los nodos. (temp_Emst)

 \4. sumar alpha*temp_Ec + (1-alpha)*temp_Ec para cada nodo y seleccionar nodo x que corresponda al minimo de alphaX*Ec + (1-alpha)*Emst.

\5. S = S Union {nodo x}

\6. repeat to 2 until S = {esten todos los nodos o hasta un numero determinado de nodos}



**mstdistances_and_energies_v2.R**: Con este archivo trabajamos actualmente. Aqui esta el código principal de agregacion de productos. Detalles:
Hacemos lo mismo que esta en mstdistances_and_energies.R, utilizando sus mismos codigos pero con datos de acoples inferidos con 25 productos: (ver new_inferring_parameters.R). Las inferencias se encuentran en: new_inferring_parameters_environment270219.RData. Dado que son 25 productos, el numero de estados es de app 33.5 millones, lo cual hace que generar la matriz de active_vertex sea impractico. Hay que hacer otro enfoque:

 \1. Se comienza por un vertice (producto) de posea alto strength en el MST (equivale a un producto de alta rotacion)

\2. A este vertice le encontramos su energia de acople Ec y distancia MST  (creamos matriz tableau con columnas: vertex, Ec , Emst, Et=Ec+Emst)

\3. inicio de loop: buscar otro vertex adyacente a v en el mst, que minimice DEc (ver paper congreso) agregar ese vertex a tableau

\4. volver a paso 3 hasta desired number of size of the cluster cluster seria el conjunto de todos los elementos de la columna vertgex en tableau.

Nota: La version 1 de  mst_distances_and_energies.R  viene de vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.



**new_inferring_parameters.R**: Este script es parte del proyecto que viene de mstdistances_and_energies.R y lo que hace es utilizar la maquina de boltzmann para inferir parametros,  pero esta vez con mayor numero de subcategorias de productos. En los proyectos anteriores hemos obtenidos 250 muestras de wide_basket para llevar a acabo 250 procesos de inferencia. Esta vez pretendemos utilizar solo una inferencia pero tomando en cuenta al menos el 50% de los productos mas vendidos en terminos de volumen.



Funciones:

**entropies_functions.R**: vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.

**network_functions_v1.R**: vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.

**ising_functions_v3.R**: vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.



Ambiente de datos generados:

**new_inferring_parameters_environment270219**: resultado de new_inerring_parameters.R





Para el congreso ENICC:

**paper_congreso_figure_replicav2.R**: replica de figura para el congreso ENICC.

**paper_congreso_figure_replica3.R**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica2v2**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica2.R**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica.R**: replica de figura para el paper del congreso a ENICC.