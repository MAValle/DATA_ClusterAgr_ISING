### Project Cluster Aggregation - DATA_ClusterAgr

Utilizamos ISING inference para descubrir los acoples de una red de productos. Luego encontramos un MST. Luego podemos darnos un punto de partida (nodo) cualquier y comenzar a agregar productos de acuerdo a la energía de barrera dada por el MST y la energía que minimiza la energía de acoples.

Al respecto, el paper de conferencia para ENICC se tienen detalles.



#### Lista de códigos relevantes

**clusteraggr_v1.R**: Lo que hacemos aqui es  utilizar el concepto de clustering aggregation basado en aglomeración. El procedimiento es así:

\1. Generamos matriz E = $\alpha*E_{c} + (1-\alpha)*E_{mst}$ para cada par de nodos.

\2. Search for a pair of nodes $i,j$ que produzca la menor E y lo agrupamos.

\3. Se repite el proceso de agregacion aglomerativo hasta que se terminen todos los nodos.

Si se hace el merge o aglomeracion con single linkage, obtendria el MST tree (con alpha=0).

**NOTA: actualmente aquí hacemos un HAC con solo una variable que es el E, pero podriamos hacerlo tomando Ec y Emst por separado**. **Para efectos de comparar, podemos comparar el HAC propuesto con energia con el de las correlaciones y tambien en terminos de complejidad de computo. **

**NOTA: Para efectos de comprar HAC tambien puedo utilizar sistema de particion, que segun el paper hierarchical clustering algorithms for document dataset, la particion es mejor que la aglomeracion.**

**clusteraggr_v2.R**: Actualmente, este script solo es una prueba artificial del 21.08.19 con n=4 nodos que hace clustering sin las funciones de *functions_hclust.R*. Hay que actualizar esto para comenzar a incorporar la energia de acople (de alguna manera) en el calculo de las distancias. (ver pag. 242 a 244 apuntes).



**monotonocity_checks.R**: Debido a que no es correcto mezclar distancias con energias en el calculo de $E$, vemos la posibilidad de utilizar distancias.



**testing_hierarchical_clustering.R**: Version primera de pruebas simples para hacer clustering jerarquico con mi propia programacion. Es la base para ``clusteraggr_v2.R``.

**testing_hierarchical_clusteringV2.R**: Version segunda de pruebas simples para hacer clustering jerarquico con mi propia programacion. Es la base para las funciones en ``functions_hclust.R`` y  ``clusteraggr_v2.R``.

**testing_hierarchical_clusteringV3.R**: En esta version 3, hacemos un test utilizando un ejemplo de matriz de distancias que se encuentra en: http://84.89.132.1/~michael/stanford/maeb7.pdf

**testing_hierarchical_clusteringV4.R**: En esta version 4, hacemos un test utilizando una matriz de distancia mas grande y con las funciones generales en functions_hclust.R desarrolladas en el 27.ago.19 basados en testing_hierarchical_clusteringV3.R

**testing_hierarchical_clusteringV5.R**: En esta version 5, probamos las funciones de *functions_hclust.R* para hacer HAC sobre la matriz de distancias del MST de los acoples para el paper y verificamos que el resultado del dendograma sea equivalente al del MST. Es decir, deseamos verificar lo indicado en pag. 237 (ver tambien pag. 241 de los apuntes.). Vemos que mi programacion de clustering jerarquico con single linkage está correcto y coincide con el MST. Implementado en la funcion hierarchical_clustering_v2.

**testing_hierarchical_clusteringV6.R**: En esta versión de single linkage, la matriz de distancia se va actualizando de acuerdo al factor gamma que declaramos en la pag. 242, en donde las distancias mst en la matriz se modifican de acuerdo a la energia de acople involucrada entre los nodos de los clusters que se fusionan. Los resultados están en pag.247. Se pierde ultrametricidad property. 

 **testing_hierarchical_clusteringV7b.R**: En esta version se implementa la idea de la pag.252, 243, 244, 250 en donde para cada fusion de cluster, se chequea primero la distancia de acople (que es la sumatoria de las energia de acople entre cada par de nodos convertidas a distancia), y se busca aquella fusion que tenga la menor energia de acople, pero bajo la restriccion de que la fusion se haga con algun nodo que esté en el MST. El proceso reconoce tres alternativas de fusion. La primera cuando se fusionan dos clusters en donde l numero de nodos de cada uno es 1. Aqui no e snecesario buscar distancia de acople porque es la misma que la distancoa mst. Segundo, cuando ambos cluster tienen más de un nodo. En este caso no se busca mejor distancia de acople y los dos clusters se unen. Tercero, cuando un cluster solo tiene 1 nodo, el el otro tiene mas de 1. En este caso se busca el mejor nodo de fusion al cluster grande que tenga la menor distancia de acople. Este algoritmo esta implementado en hierarchical_clustering_v3. Los resultados indican que no hay diferencias entre este forma de fusion y la de single linkage tradicional (implementado en hierarchical_clustering_v2).
Nota: la version testing_hierarchical_clusteringV7.R es una verison antigua de backup no vigente.

 **testing_hierarchical_clusteringV8.R**: Esta version es la misma que la V7b, pero en este caso, cuando hay dos cluster que se van a fusionar en donde los dos tienen más de 1 nodo, se buscar para cad acluster, otro nodo u otro cluster que logre la menor distancia de acople. En la version V7b, cuando teníamos este caso, los dos cluster se unían sin buscar minima distancia de acople.

**testing_hierarchical_clusteringV9.R**: Se crea ahora una algoritmo mejorado de hierarchical_clustering_v4 que parte de la base de ir uniendo clusters con el menor numero de nodos asumiendo que mientras mas nodos tiene un cluster, la distancia de acople siempre es mayor.

**testing_hierarchical_clusteringV10.R**: Este script consiste en ensayo y test para crear la funcion hierarchical_clustering_greedy que lo que hace es utilizar algoritmo greedy para ir fusionando pares de clusters basado en la distancia de acople, intentando minimizarla.

**testing_hierarchical_clusteringV11.R**: Este script consiste en ensayo y test para crear la funcion hierarchical_clustering_probabilistic_greedy que lo que hace es una aproximación greedy, pero probabilistica, es decir, para evitar que el greedy se quede atascado en un optimo local, con un cierto nivel de probabilidad vamos suitchiando entre una fusion greedy determinista y una fusion greedy probabilistica. Esta ultima consiste en determinar la fusion de dos cluster en forma aleatoria.

**testing_hierarchical_clusteringV12.R**: Este script lo que hace es buscar la manera de obtener la distancia de acople en cada iteracion del proceso clustering con modularidad utilizando la funcion cluster_fast_greedy de igraph. El proceso es basicamente ejecutar la función, y a partir del output de la funcion determinar los nodos que se van fusionando en  cada iteración, y luego sacarles la distancia de acople.

**testing_hierarchical_clusteringV13.R**: Basado en los resultados de simulation_hc_V6.R y simulation_hc_V7.R, vemos que al agregar mas randomness al algoritmo greedy, la media de distancia de acople $d_c$ va aumentando en relación al greedy deterministico. Asimismo, la varianza tambien aumenta, lo cual indica que a veces, el greedy probabilistico hace muy bien el trabajo buscando una distancia de acople menor al greedy deterministico, pero a veces tambien lo hace muy mal. Por lo tanto, es de suponer, que el greedy probabilistico necesita una ayuda. Lo que hacemos en este script es dar un punto de partida basado en el greedy deterministico, para mejor aun mas el desempeño del greedy probabilistico, disminuyendo su variabilidad en los resultado y tratando de estar siempre mejor que el deterministico.





**simulation_hc_V1.R**: En este script vamos a utilizar la funcion hierarchical_clustering_v2 y funcion hierarchical_clustering_v4 que estan en functions_hclust.R.

**simulation_hc_V2.R**: En este script vamos a utilizar la funcion hierarchical_clustering_v2 y funcion hierarchical_clustering_v4 que estan en functions_hclust.R

**simulation_hc_V3.R**: En este script vamos a simular clusterig jerarquicos de distintos tamanos de nodos de N= 20, 30, 50, 100, 250, y 500 con matrices de acoples con <J>=0. Luego haremos una grafica scatterplot en que graficamos en el eje X la distancia de acople del algoritmo normal, y en el eje Y la distancia de acople del algoritmo modificado para cada una de las iteraciones.

**simulation_hc_V4.R**: En este script vamos a simular clustering jerarquicos de distintos tamanos de nodos de N=  25, 50, 100 con matrices de acoples con <J>. La intencion es hacer un mapa de calor en que en el eje X va la iteracion de merge (para un determinado numero de nodos) y en el eje Y va <J>. Es decir, para un determinado numero de nodos, fabricamos una matriz con las iteraciones en un eje y <J> por otro, y la  en la superficie graficamos la distancia de acople. Podemos hacer el mapa de calor, uno para la distancia de acople de merge2, otra de merge5 y otra que sea la diferencia.

**simulation_hc_V5.R**: En este script simulamos una matriz J y comparamos las distancias de acople $d_c$ para distinto numero de clusters formados con greeddy algorithm con distintos niveles de aleatoriadad. Los resultados se puede ver en results_from_simulation_hc_V5_191119(b).pdf y también en results_from_simulation_hc_V5_191119.pdf.

**simulation_hc_V6.R**:  En este script simulamos varias matrices de acople N(0,1) y con N=20 para graficar las distancias de acople en cada iteraciones en cada uno de los algoritmos de prueba.

**simulation_hc_V7.R**: En este script simulamos 1 (un) ejemplo con N=25 nodos, para mostar como caso puntual y comprar entre todos los algoritmos. Los resultados se pueden ver en results_from_simulation_hc_V7_261219.pdf y en results_from_simulation_hc_V7_261219(2).pdf. Lo que vemos es que mientras le agregamos mas randomness al proceso greedy, la *media* de la distancia de acoples es cada vez mas mala en relacion al greedy deterministico, lo cual  tiene sentido puesto que al ser mas random, es como ir con  los ojos cerrados, y tendremos a veces mala suerte o buena suerte en mejorar la distancia de acople en relacion al greedy deterministico. También vemos que vemos que mientras mas randomnes le agregamos al proceso, habra mas varianza, lo que significa que podemos tener muy buenos resultados o muy malos resultados en referencia al grreedy deterministico. 

**simulation_hc_V8.R**: En este script analizamos el ejemplo real de la bases de datos transaccional que he utilizado en los papers de ICANN con 25 nodos, para calcular los dendogramas con MST y con greedy. Adicionalmente, vemos las simulaciones con greedy probabilistico.





**computing_silhouettes_V1.R**:  Tratamos de encontrar una distancia k de corte en el greedy deterministico (gamma=1) tal que el promedio de silhouette o el promedio de a(i) sea el mas grande posible. 

**computing_silhouettes_V2.R**: Dada una solucion a una distancia de acople h, encontramos consistencia interna utilizando la funcion my_own_silhouette en  my_own_silhouette_function.R. Ver pag. 275-277 apuntes.





**searching_v1.R**: Lo que hacemos aqui es ir buscando las energias de acople Ec y las barreras de energia Emst, de tal forma de minimizar $\alpha *E_c + (1-\alpha)*E_{mst}$. Este es para ir haciendo la agregación de productos desde un punto de partida.

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

\4. volver a paso 3 hasta desired number of size of the cluster cluster seria el conjunto de todos los elementos de la columna vertex en tableau.

Nota: La version 1 de  mst_distances_and_energies.R  viene de vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.



**new_inferring_parameters.R**: Este script es parte del proyecto que viene de mstdistances_and_energies.R y lo que hace es utilizar la maquina de boltzmann para inferir parametros,  pero esta vez con mayor numero de subcategorias de productos. En los proyectos anteriores hemos obtenidos 250 muestras de wide_basket para llevar a acabo 250 procesos de inferencia. Esta vez pretendemos utilizar solo una inferencia pero tomando en cuenta al menos el 50% de los productos mas vendidos en terminos de volumen.



####Figure replicas for paper

**density_couplings_forpaper.R**: Este script calcula el histigrama y densidad y sus respectiva grafica de la pdf de los acoples para efectos de publicacion en el paper.

**consistency_check_plot_forpaper.R**: Replica de la figura de recovered <si> , <sisj> y Cij para el paper. Seguimos el mismo procedimiento establecido en one_and_two_body_abalysisV2.R y exactamente el mismo para el paper de congreso en paper_congreso_figure_replica2-R

**networks_plot_forpaper.R**: Replica de las figuras de redes de acoples con distintos umbrales a la distirbución de los acoples. Esto se basa en paper_congreso_figure_replica3.R.

**mst_plot_forpaper.R**: Replcia de figuras de MST para paper.

**figure_netcoupling_for_IEEE**: 091220 -  En este script replicamos la figura "coupling_network_for_IEEE.eps" Aqui tomamos exactamente la figura que se hizo para la conferencia ICANN 2019 munich, hecha en paper_congreso_figure_replica3.R, pero con distinto orden y distinto colores para que no se vea  igual.

**heatmaps_figures.R**:  En este script nos enfocamos en hacer las graficas de  heatmap, que originalmente iniciamos en  simulation_hc_V8.R para el caso real de 25 nodos (datos de verdad) Este script reproduce los heatmaps que iran en el paper.

**figure_netcoupling_for_IEEE.R**:  En este script replicamos la figura "figure_coupling_network_for_IEEE.eps" Aqui tomamos exactamente la figura que se hizo para la conferencia ICANN 2019 munich, hecha en paper_congreso_figure_replica3.R, pero con distinto orden y distinto colores para que no se vea  igual.





####Funciones:

**entropies_functions.R**: vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.

**network_functions_v1.R**: vienen originalmente de PAPER_MBA by solving the inverse ISING problem.

**ising_functions_v3.R**: vienen originalmente de PAPER_MBA by solving teh inverse ISING problem.

**find_starting_vertex_function.R**: FIND strongest vertex

**find_mst_barrier_function.R**: Inspirado en la idea de la L164 en adelante en mstdistances_and_energies.R

**is_integer0_function.R**: función para testear si un vector esta vacio.

**functions_hclust.R**: Conjunto de functiones para hacer clustering jerarquico HAC con single linkage a partir de una matriz de distancia $D$. 

**acople_distance_sum_function.R**: esta funcion calcula la sumatoria de las acoples convertidas en distancia previamente.

**hierarchical_clustering_greedy_function.R**: Lleva a cabo un hiearchical clustering utilizando greedy approach. 

**hierarchical_clustering_probabilistic_greedy_function.R**: Lleva a cabo un hiearchical clustering utilizando greedy approach pero en forma probabilistica, es decir, en cada oportunidad de fusion de dos clusters, con cierta probabilidad se decide hacer la fusion basado en la minima distancia de acople o en forma aleatoria. Esta ultima modalidad nos permite salir eventualmente de un optimo local.

**pick_a_cluster_function.R**: esta funcion selecciona al azar un par de cluster de la matriz de distancia de acoples. Esta funcion se utiliza en la funcion hierarchical_clustering_greedy y hierarchical_clustering_probabilistic_greedy_function.

**find_min_distcp_function.R**: esta funcion lo que hace es actualizar la matriz de distancias de acople. Esta funcion se utiliza en las funciones hierarchical_clustering_greedy y hierarchical_clustering_probabilistic_greedy_function.

**create_coupling_function.R**: Funcion que crea matriz de acople y de distancia con media mu, desv estandar de acople sj y con Nn nodos.

**recovering_couplingdistances_from_fastgreedyigraph_function.R**: Set de funciones para recuperar la distancia de acople a partir  de community detection using modularity with cluster_fast_greedy of igraph. developed in testing_hierarchical_clusteringV12.R.

**post_hc_analysis_function.R**: Set de funciones para analizar los resultados de  varias simulaciones de los algoritmos de hiearchical clustering.

**recovering_couplingdistances_from_fastgreedyigraph_function.R**: Set de funciones para recuperar la distancia de acople a partir  de community detection using modularity with cluster_fast_greedy of igraph. Developed in testing_hierarchical_clusteringV12.R

**cophonetic_cor_function.R**: Funcion para calcular la correlacion cofonetica.

**silhouette_byhand_function.R**: Calculamos silhuette. He comprobado que esta funcion da lo mismo resultados que la funcion silhouette de la librerya cluster. Link:  https://en.wikipedia.org/wiki/Silhouette_(clustering) silhouette: a(i) = 1/(#Ci - 1) sum d(i,j)  suma de todas las distancias desde i a j, donde j son elementos intracluster. b(i) = min (1/#Ck) sum d(i,j) min de la suma de las distancias desce i a j donde j son elementos fuera del cluster donde esta i.

**my_own_silhouette_function.R**:  computing my own measure of internal consistency simmilar to the silhouette Given a solution with a coupling distance h, we have for each element of the system their corresponding cluster. From this solution we want compute a internal consistency to determine which  groups of clusters are more choesive and find whether teh solution has a good internal consistency.







####Ambiente de datos generados:

**new_inferring_parameters_environment270219.RData**: resultado de new_inerring_parameters.R

**plots_from_network_plot_forpaper310719.rds**: 

**fieldsandcouplings250_doing_parallel_for_physicaA_090718.rds**





####Figure replicas for conference ENICC:

**paper_congreso_figure_replicav2.R**: replica de figura para el congreso ENICC.

**paper_congreso_figure_replica3.R**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica2v2**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica2.R**: replica de figura para el paper del congreso a ENICC.

**paper_congreso_figure_replica.R**: replica de figura para el paper del congreso a ENICC.