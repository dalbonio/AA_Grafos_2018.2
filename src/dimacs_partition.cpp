#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/energybased/GEMLayout.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/LongestPathRanking.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/FastSimpleHierarchyLayout.h>
#include <ogdf/layered/OptimalHierarchyClusterLayout.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/Array.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf/basic/Queue.h>


#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterPlanarizationLayout.h>

#include <iostream>     // std::cout, std::ostream, std::ios
#include <fstream>      // std::filebuf
#include <string>

using namespace ogdf;
using namespace std;

#define EPSILON 0.0000000000000000001
#define QUANTIDADE_PARA_PERMUTAR 2
#define QUANTIDADE_CROMOSSOMOS 500
#define QUANT_FILHOS 16
#define FILHOS_TROCA 4
#define MAX_ITER 60000
#define MIN_CROSSING_OVER_RATE 20/(double)100
#define MAX_CROSSING_OVER_RATE 80/(double)100

#define RANDOM ((double)rand() / RAND_MAX)

static const int nmaxClusters = 8;
static const Color::Name clustCol[nmaxClusters] = {
		Color::Name::Green,
		Color::Name::Magenta,
		Color::Name::Red,
		Color::Name::Yellow,
		Color::Name::Brown,
		Color::Name::Violet,
		Color::Name::Lime,
		Color::Name::Turquoise
};

struct Cromossomo
{
	double fitness;

	node* genes;
};

static void getGraphSize(FILE* graphFile, size_t * n, size_t * m) {
	char type = ' ';
	char linestr[100];
	char* datastr;
	*n = 0;
	*m = 0;

	while (type != 'p') {
		type = fgetc(graphFile);
		if (type != EOF) {

			/* header */
			if (type == 'c') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr != NULL)
					printf("%s", linestr);
				else {
					*n = -1;
					return;
				}
			}

			/* Vertices */
			if (type == 'p') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr == NULL) {
					*n = -1;
					return;
				}
				datastr = strtok(linestr, " ");
				printf("\tdatastr:%s\n", datastr);
				if (!strcmp(datastr, "edge") || !strcmp(datastr, "col")) {
					datastr = strtok(NULL, " ");
					printf("\tdatastr:%s\n", datastr);
				}
				*n = atoi(datastr);

				datastr = strtok(NULL, " ");
				printf("\tdatastr:%s\n", datastr);
				*m = atoll(datastr);

				printf("Graph with (%ld) vertices and (%ld) edges.\n", *n, *m);
				printf("Density: %3.2f\n", 2.0 * ((double) *m) / (*n * (*n)));
			}
		}
	}
}


void readDimacsGraph(Graph& g, FILE* graphFile) {
	char type = ' ';
	char linestr[100];
	char* datastr;
	int i, j;
	size_t n;
	size_t m;

	printf("\t#Starts reading...\n");
	fflush(NULL);
	getGraphSize(graphFile, &n, &m);

	////
	// Graph variables
	////
	Array<node> V(n+1);
	for(int i = 1; i <= n; ++i)
		V[i] = g.newNode(i);

	type = fgetc(graphFile);
	while (type != EOF) {
		/* Edges */
		if (type == 'e') {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr == NULL)
				return;

			datastr = strtok(linestr, " ");
			i = atoi(datastr);

			datastr = strtok(NULL, " ");
			j = atoi(datastr);

			g.newEdge(V[i],V[j]);
			//g.newEdge(V[j],V[i]);
		} else {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr != NULL)
				printf(" %s\n", linestr);
			else
				return;
		}
		type = fgetc(graphFile);
	}
	printf("\t#Finishes reading graph\n");
	fflush(NULL);

	return;
}

static inline void ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) {
        return !isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](int ch) {
        return !isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(string &s) {
    ltrim(s);
    rtrim(s);
}

void readJostleGraph(Graph& g, FILE* graphFile)
{
	char type = ' ';
	char linestr[300];
	char* datastr;
	int i, j;
	size_t n;
	size_t m;

	printf("\t#Starts reading...\n");
	fflush(NULL);
	fgets(linestr, 300, graphFile);

	datastr = strtok(linestr, " ");
	printf("\tdatastr:%s\n", datastr);
	n = atoi(datastr);
	datastr = strtok(NULL, " ");
	printf("\tdatastr:%s\n", datastr);
	m = atoll(datastr);

	////
	// Graph variables
	////
	Array<node> V(n+1);
	for(int i = 1; i <= n; ++i)
		V[i] = g.newNode(i);

	i = 1;
	while (i <= n)
	{
		/* Edges */
		fgets(linestr, 300, graphFile);
		int length = strlen(linestr);
		if( linestr[length - 1] == '\n')
		{
			linestr[length - 1] = '\0';
		}

		datastr = strtok(linestr, " ");

		int j;
		while(datastr != NULL)
		{
			j = atoi(datastr);
			g.newEdge(V[i],V[j]);
			datastr = strtok(NULL, " ");
		}
		++i;
	}
	printf("\t#Finishes reading graph\n");
	fflush(NULL);

	return;
}



//					PARTE DO GENETICO			//


//genetic one
int partitionEdgesArray(Graph& G, node* sub_graph, int& array_size, int& bfs_times)
{
	int edge_count = 0;
	NodeArray<bool>visitados = NodeArray<bool>(G);
    visitados.fill(false);
    NodeArray<bool>visitados_total = NodeArray<bool>(G);
	visitados_total.fill(false);

    node iter_node = sub_graph[0];
	visitados[iter_node] = true;
    int visited_nodes = 1;
    int list_size = array_size;
    while(visited_nodes < list_size)
    {
        visitados.fill(false);
        Queue<node> q;
        q.append(iter_node);
        while(q.size()!=0)
        {
            node currentNode = q.pop();    //pop elements from the q

			if(visitados[currentNode] == true)
			{
				continue;
			}

			visitados[currentNode] = true;
			visitados_total[currentNode] = true;
			visited_nodes++;
            adjEntry adj = currentNode->firstAdj();

            for(; adj; adj = adj->succ())
    		{
    			node adjNode = adj->twinNode();
				int same_cluster = 0;
				for(int i = 0; i < list_size; ++i)
				{
					if(sub_graph[i]->index() == adjNode->index())
					{
						same_cluster = 1;
						break;
					}
				}

				//printf("Current: %d, Twin: %d\n", currentNode->index(), adjNode->index());
                if(same_cluster == 1)
                {
                    if(visitados[adjNode] == false)
                    {
						++edge_count;
                        q.append(adjNode);
                    }
                }
            }
        }

        iter_node = sub_graph[0];
		int i = 0;
		while(visited_nodes < list_size && i < list_size && visitados_total[iter_node] == true)
        {
        	iter_node = sub_graph[i];
			i++;
        }

		bfs_times += 1;
	}

	return edge_count;
}

void shuffle(void *array, size_t n, size_t size)
{
	char tmp[size];
	char *arr = (char*)array;
	size_t stride = size * sizeof(char);

	if (n > 1)
	{
		size_t i;
		for (i = 0; i < n - 1; ++i)
		{
			size_t rnd = (size_t) rand();
			size_t j = i + rnd / (RAND_MAX / (n - i) + 1);

			memcpy(tmp, arr + j * stride, size);
			memcpy(arr + j * stride, arr + i * stride, size);
			memcpy(arr + i * stride, tmp, size);
		}
	}
}

double sub_graph_density_array(Graph& g, node* subGraph, int& k, int& length, int& bfs_times)
{
	int pEdgeCount = partitionEdgesArray(g, subGraph, length, bfs_times);
	int nodeQuant = length;

	//printf("Edges: %d - Nodes: %d - ", pEdgeCount, nodeQuant);

	return ((double)(2 * pEdgeCount)) / ( (nodeQuant) * (nodeQuant - 1) );
}

double partition_density_array(Graph& g, node* partition, int& k, int& size1, int& len1, int& len2, int& bfs_times)
{
	double sum = 0;
	int i;
	for(i = 0; i < size1; ++i)
	{
		//printf("P[%d]: ", i);
		//printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[i * len1], k, len1));
		sum += sub_graph_density_array(g, &partition[i * len1], k, len1, bfs_times);
	}
	int starter = size1 * len1;
	int counter = 0;
	for(; i < k; ++i)
	{
		//printf("P[%d]: ", i);
		//printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[starter + counter * len2], k, len2));
		sum += sub_graph_density_array(g, &partition[starter + counter * len2], k, len2, bfs_times);
		counter++;
	}

	return sum;
}

node* initialize_genes(node* initial_array, int& size)
{
	node* genes = new node[size];

	shuffle(initial_array, size, sizeof(node));

	memcpy(genes, initial_array, size);

	return genes;
}

node* initialize_genes_viciado(node* initial_array, int& size)
{
	node* genes = new node[size];

	memcpy(genes, initial_array, size);

	shuffle(initial_array, size, sizeof(node));

	return genes;
}

int sub_graph_sum(node* sub_graph, int& size, int& k, int length, int partition_start, NodeArray<NodeArray<int>>& dist)
{
	int sum = 0;
	int final_node = partition_start + length;
	for(int n = partition_start; n < final_node; n++)
	{
		int i;
		for(i = final_node; i < size; i++)
		{
			sum += dist[sub_graph[n]][sub_graph[i]];
		}
	}

	return sum;
}

int partition_sum(node* nodes, int& size, int& k, int& size1, int& len1, int& len2, NodeArray<NodeArray<int>>& dist)
{
	int sum = 0;
	int i;
	for(i = 0; i < size1; ++i)
	{
		//printf("P[%d]: ", i);
		//printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[i * len1], k, len1));
		sum += sub_graph_sum(nodes, size, k, len1, i * len1, dist);
	}
	int starter = size1 * len1;
	int counter = 0;
	for(; i < k - 1; ++i)
	{
		//printf("P[%d]: ", i);
		//printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[starter + counter * len2], k, len2));
		sum += sub_graph_sum(nodes, size, k, len2, starter + counter * len2, dist);
		counter++;
	}

	return sum;
}

void fitness_cromossomo_somatorio(Cromossomo& cromossomo, int& size, int& k, int& size1, int& len1, int& len2, NodeArray<NodeArray<int>>& dist)
{
	cromossomo.fitness = partition_sum(cromossomo.genes, size, k, size1, len1, len2, dist) / (double)10;
}

void fitness_cromossomo(Cromossomo& cromossomo, Graph& g, int& k, int& size1, int& len1, int& len2)
{
	int bfs_times = 0;
	cromossomo.fitness = partition_density_array(g, cromossomo.genes, k, size1, len1, len2, bfs_times);
	printf("\nBfs_Times: %d\n", bfs_times);
}

//penalize partitions with more disconected nodes
void fitness_cromossomo_consider_bfs_times(Cromossomo& cromossomo, Graph& g, int& k, int& size1, int& len1, int& len2)
{
	int bfs_times = 0;
	cromossomo.fitness = partition_density_array(g, cromossomo.genes, k, size1, len1, len2, bfs_times);

	cromossomo.fitness -= cromossomo.fitness * (bfs_times / (double)100);
}

int roleta(Cromossomo* cromossomos, double& sum_rate)
{
	if(sum_rate < EPSILON)
	{
		//printf("entrou na roleta zerada");
		return rand() % QUANTIDADE_CROMOSSOMOS;
	}
	float rand_number = RANDOM * sum_rate;
	int i;
	for(i = 0; rand_number > 0 && i < QUANTIDADE_CROMOSSOMOS; ++i)
	{
		rand_number -= cromossomos[i].fitness;
	}

	if( i == QUANTIDADE_CROMOSSOMOS )
	{
		return i - 1;
	}

	return i;
}

void trocar_filhos(Cromossomo* cromossomos, Cromossomo* filhos, double& sum_rate, int& size, int filhos_troca)
{
	int worst_pai;
	int best_filho;

	for(int j = 0; j < filhos_troca; ++j)
	{
		worst_pai = j;
		best_filho = j;
		for(int i = j + 1; i < QUANTIDADE_CROMOSSOMOS; ++i)
		{
			if(cromossomos[i].fitness < cromossomos[worst_pai].fitness)
			{
				worst_pai = i;
			}
		}
		for(int i = j + 1; i < QUANT_FILHOS; ++i)
		{
			if(filhos[i].fitness > filhos[best_filho].fitness)
			{
				best_filho = i;
			}
		}

		//printf("Troca %d: Pior Pai %d por Melhor Filho %d", j, worst_pai, best_filho);

		//manipulate sum_rate in sons trade
		sum_rate -= cromossomos[worst_pai].fitness;
		sum_rate += cromossomos[best_filho].fitness;

		Cromossomo temp;

		memcpy(cromossomos[worst_pai].genes, filhos[best_filho].genes, sizeof(node) * size);
		cromossomos[worst_pai].fitness = filhos[best_filho].fitness;

		//first one now is best son
		temp = cromossomos[j];
		cromossomos[j] = cromossomos[worst_pai];
		cromossomos[worst_pai] = temp;

		//trade it in sons too
		temp = filhos[j];
		filhos[j] = filhos[best_filho];
		filhos[best_filho] = temp;
	}

	return;
}

void create_filho_particao_order_based(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	//rand from minimal to maximal crossing rate( random interval is [0,1] )


	int* selecionados = new int[QUANTIDADE_PARA_PERMUTAR];
	int* selecionados_cluster = new int[QUANTIDADE_PARA_PERMUTAR];
	int* mom_pos = new int[QUANTIDADE_PARA_PERMUTAR];
	int* dad_pos = new int[QUANTIDADE_PARA_PERMUTAR];

	int counter = 0;

	while(counter < QUANTIDADE_PARA_PERMUTAR)
	{
		int new_rand = rand() % size;

		int new_rand_cluster;
		if(new_rand < size1 * len1)
		{
			new_rand_cluster = new_rand / len1;
		}
		else
		{
			new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
		}

		int i;
		for(i = 0; i < counter && selecionados_cluster[i] != new_rand_cluster; i++);

		//printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

		if(i == counter)
		{
			//printf("selecionados[%d]: Pos %d\n", counter, new_rand);
			selecionados[i] = new_rand;
			selecionados_cluster[i] = new_rand_cluster;
			counter++;
		}
	}

	for(int i = 0; i < QUANTIDADE_PARA_PERMUTAR; i++)
	{
		int found = 0;
		for(int j = 0; j < size; j++)
		{
			if(pai.genes[j]->index() == mae.genes[selecionados[i]]->index())
			//if(pai.genes[j] == mae.genes[selecionados[i]])
			{
				//printf("dad_pos[%d]: %d\n", i, j);
				dad_pos[i] = j;
				found += 1;
			}
			if(mae.genes[j]->index() == pai.genes[selecionados[i]]->index())
			//if(mae.genes[j] == pai.genes[selecionados[i]])
			{
				//printf("mom_pos[%d]: %d\n", i, j);
				mom_pos[i] = j;
				found += 1;
			}

			if(found == 2)
			{
				break;
			}
		}
	}

	memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	//memcpy(filho1.genes, pai.genes, sizeof(int) * size);
	//memcpy(filho2.genes, mae.genes, sizeof(int) * size);

	for(int i = 0; i < QUANTIDADE_PARA_PERMUTAR; i++)
	{
		int menor_mae = 0;
		int menor_pai = 0;
		int menor_selecionado = 0;
		for(int j = 1; j < QUANTIDADE_PARA_PERMUTAR; j++)
		{
			if(mom_pos[j] < mom_pos[menor_mae])
			{
				menor_mae = j;
			}
			if(dad_pos[j] < dad_pos[menor_pai])
			{
				menor_pai = j;
			}
			if(selecionados[j] < selecionados[menor_selecionado])
			{
				menor_selecionado = j;
			}
		}

		//printf("Menor mae: %d, Menor pai: %d\n", menor_mae, menor_pai);

		filho1.genes[selecionados[menor_selecionado]] = mae.genes[mom_pos[menor_mae]];
		filho2.genes[selecionados[menor_selecionado]] = pai.genes[dad_pos[menor_pai]];

		dad_pos[menor_pai] = size;
		mom_pos[menor_mae] = size;
		selecionados[menor_selecionado] = size;
	}

  /*  printf("Filho1: ");
	for(int i = 0; i < size; i++)
	{
		printf("%d ", filho1.genes[i]);
	}
	printf("\nFilho2: ");
	for(int i = 0; i < size; i++)
	{
		printf("%d ", filho2.genes[i]);
	}*/

	return;
}

//after test, change it to rand from 1 to k, then get x different from the fixed cluster;
void create_filho_spread_same_clusters_rand(Cromossomo& filho, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int max_quantity = 10;
	int quantity = 2 + rand() % (max_quantity - 2);
	int* selecionados = new int[quantity];
	int* previously_selected = new int[quantity];
	int selecionados_cluster = -1;

	int counter = 0;

	while(counter < quantity)
	{
		int new_rand;
		int new_rand_cluster;
		do
		{
			int new_rand = rand() % size;
			if(new_rand < size1 * len1)
			{
			   new_rand_cluster = new_rand / len1;
			}
			else
			{
			   new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
			}
		}while(new_rand_cluster == selecionados_cluster);

		//printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

		//printf("selecionados[%d]: Pos %d\n", counter, new_rand);

		selecionados[counter] = new_rand;
		selecionados_cluster = new_rand_cluster;
		counter++;
   }

	memcpy(filho.genes, pai.genes, sizeof(node) * size);

	int starter = size1 * len1;

	int j;

	for(int i = 0; i < quantity; i++)
	{
		int clusterLength = len1;
		int clusterInitial = 0;

		int selected = 0;
		int counter = 0;
		int new_rand;
		int new_rand_cluster;

		do
		{
			new_rand = rand() % size;

			if(new_rand < size1 * len1)
			{
				new_rand_cluster = new_rand / len1;
				clusterLength = len1;
				clusterInitial = (new_rand/len1) * len1;
			}
			else
			{
				new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
				clusterLength = len2;
				clusterInitial = starter + (((new_rand - starter)/len2) * len2);
			}

			for(j = 0; j < counter && previously_selected[j] != new_rand && new_rand_cluster != selecionados_cluster; j++);
		}while(j != counter);

		previously_selected[j] = new_rand;
		counter++;

		node temp = filho.genes[selecionados[i]];
		filho.genes[selecionados[i]] = filho.genes[new_rand];
		filho.genes[new_rand] = temp;
	}
}

//after test, change it to rand from 1 to k, then get x different from the fixed cluster;
void create_filho_spread_same_clusters_block_rand(Cromossomo& filho, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int cluster_from = rand() % k;

	int quantity = 1 + rand() % (len2 - 1);
	int selected_from = rand() % (len2 - quantity);
	int selected_to = rand() % (len2 - quantity);

	memcpy(filho.genes, pai.genes, sizeof(node) * size);
	int starter = size1 * len1;

	int cluster_to;
	do
	{
		cluster_to = rand() % k;
	}while(cluster_from == cluster_to);

	int initial_from;
	int initial_to;

	if(cluster_from < size1)
	{
		initial_from = cluster_from * len1;
	}
	else
	{
		initial_from = starter + ((cluster_from - size1) * len2);
	}

	if(cluster_to < size1)
	{
		initial_to = cluster_to * len1;
	}
	else
	{
		initial_to = starter + ((cluster_to - size1) * len2);
	}

	int counter = 0;

	while(counter < quantity)
	{
		node temp = filho.genes[selected_from + initial_from];
		filho.genes[selected_from + initial_from] = filho.genes[selected_to + initial_to];
		filho.genes[selected_to + initial_to] = temp;

		selected_to++;
		selected_from++;
		counter++;
   }

}

void create_filho_shuffle_between_clusters_rand(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int quantity = 1 + rand() % (k - 1);
	int* selecionados = new int[quantity];
	int* selecionados_cluster = new int[quantity];
	int* new_selecionados = new int[quantity];

	int counter = 0;

	while(counter < quantity)
	{
	   int new_rand = rand() % size;

	   int new_rand_cluster;
	   if(new_rand < size1 * len1)
	   {
		   new_rand_cluster = new_rand / len1;
	   }
	   else
	   {
		   new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
	   }

	   int i;
	   for(i = 0; i < counter && selecionados_cluster[i] != new_rand_cluster; i++);

	   //printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

	   if(i == counter)
	   {
		   //printf("selecionados[%d]: Pos %d\n", counter, new_rand);
		   selecionados[i] = new_rand;
		   new_selecionados[i] = new_rand;
		   selecionados_cluster[i] = new_rand_cluster;
		   counter++;
	   }
	}

	memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	shuffle(new_selecionados, quantity, sizeof(int));

	for(int i = 0; i < quantity; i++)
	{
	   node temp = filho1.genes[new_selecionados[i]];
	   filho1.genes[new_selecionados[i]] = filho1.genes[selecionados[i]];
	   filho1.genes[selecionados[i]] = temp;

	   temp = filho2.genes[new_selecionados[i]];
	   filho2.genes[new_selecionados[i]] = filho2.genes[selecionados[i]];
	   filho2.genes[selecionados[i]] = temp;
	}
}

void create_filho_shuffle_between_clusters_fixed(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int* selecionados = new int[QUANTIDADE_PARA_PERMUTAR];
	int* selecionados_cluster = new int[QUANTIDADE_PARA_PERMUTAR];
	int* new_selecionados = new int[QUANTIDADE_PARA_PERMUTAR];

	int counter = 0;

	while(counter < QUANTIDADE_PARA_PERMUTAR)
	{
	   int new_rand = rand() % size;

	   int new_rand_cluster;
	   if(new_rand < size1 * len1)
	   {
		   new_rand_cluster = new_rand / len1;
	   }
	   else
	   {
		   new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
	   }

	   int i;
	   for(i = 0; i < counter && selecionados_cluster[i] != new_rand_cluster; i++);

	   //printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

	   if(i == counter)
	   {
		   //printf("selecionados[%d]: Pos %d\n", counter, new_rand);
		   selecionados[i] = new_rand;
		   new_selecionados[i] = new_rand;
		   selecionados_cluster[i] = new_rand_cluster;
		   counter++;
	   }
	}

	memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	shuffle(new_selecionados, QUANTIDADE_PARA_PERMUTAR, sizeof(int));

	for(int i = 0; i < QUANTIDADE_PARA_PERMUTAR; i++)
	{
	   node temp = filho1.genes[new_selecionados[i]];
	   filho1.genes[new_selecionados[i]] = filho1.genes[selecionados[i]];
	   filho1.genes[selecionados[i]] = temp;

	   temp = filho2.genes[new_selecionados[i]];
	   filho2.genes[new_selecionados[i]] = filho2.genes[selecionados[i]];
	   filho2.genes[selecionados[i]] = temp;
	}
}

void create_filho_mom_cluster_rand(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int quantity = 1 + rand() % (k - 1);
	int* selecionados = new int[quantity];
	int* selecionados_cluster = new int[quantity];

	int counter = 0;

	//printf("\n");
	while(counter < quantity)
	{
	   int new_rand = rand() % size;

	   int new_rand_cluster;
	   if(new_rand < size1 * len1)
	   {
		   new_rand_cluster = new_rand / len1;
	   }
	   else
	   {
		   new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
	   }

	   int i;
	   for(i = 0; i < counter && selecionados_cluster[i] != new_rand_cluster; i++);

	   //printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

	   if(i == counter)
	   {
		   //printf("selecionados[%d]: Pos %d\n", counter, new_rand);
		   selecionados[i] = new_rand;
		   selecionados_cluster[i] = new_rand_cluster;
		   counter++;
	   }
	}

	memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	int starter = size1 * len1;

	int j;

	for(int i = 0; i < quantity; i++)
	{
		//find dad and mom cluster in graph
		int dad_cluster = 0;
		int mom_cluster = 0;
		int clusterLength = len1;
		int clusterInitial1 = 0;
		int clusterInitial2 = 0;

		int selected = 0;
		int counter = 0;
		for(; selected < size1; ++selected)
		{
			for(j = 0; j < len1; ++j)
			{
				if(pai.genes[selecionados[selected]]->index() == mae.genes[counter]->index())
				{
					mom_cluster = selected;
					clusterLength = len1;
					clusterInitial1 = len1 * selected;
				}
				if(mae.genes[selecionados[selected]]->index() == pai.genes[counter]->index())
				{
					dad_cluster = selected;
					clusterLength = len1;
					clusterInitial2 = len1 * selected;
				}
				counter++;
			}
		}

		int cluster_counter = 0;
		for(; selected < k; ++selected)
		{
			for(j = 0; j < len2; ++j)
			{
				if(pai.genes[selecionados[selected]]->index() == mae.genes[counter]->index())
				{
					mom_cluster = selected;
					clusterLength = len2;
					clusterInitial1 = starter + len2 * cluster_counter;
				}
				if(mae.genes[selecionados[i]]->index() == pai.genes[counter]->index())
				{
					dad_cluster = selected;
					clusterLength = len2;
					clusterInitial2 = starter + len2 * cluster_counter;
				}
				counter++;
			}
			cluster_counter++;
		}

		int rand_from_cluster = rand() % clusterLength;

		node temp = filho1.genes[selecionados[i]];
		filho1.genes[selecionados[i]] = filho1.genes[clusterInitial1 + rand_from_cluster];
		filho1.genes[clusterInitial1 + rand_from_cluster] = temp;

		temp = filho2.genes[selecionados[i]];
		filho2.genes[selecionados[i]] = filho2.genes[clusterInitial2 + rand_from_cluster];
		filho2.genes[clusterInitial2 + rand_from_cluster] = temp;
	}
}

void create_filho_mom_cluster_fixed(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int& size, int& k, int& size1, int& len1, int& len2)
{
	int* selecionados = new int[QUANTIDADE_PARA_PERMUTAR];
	int* selecionados_cluster = new int[QUANTIDADE_PARA_PERMUTAR];

	int counter = 0;

	//printf("\n");
	while(counter < QUANTIDADE_PARA_PERMUTAR)
	{
	   int new_rand = rand() % size;

	   int new_rand_cluster;
	   if(new_rand < size1 * len1)
	   {
		   new_rand_cluster = new_rand / len1;
	   }
	   else
	   {
		   new_rand_cluster = ((new_rand - (size1 * len1))/len2) + size1;
	   }

	   int i;
	   for(i = 0; i < counter && selecionados_cluster[i] != new_rand_cluster; i++);

	   //printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

	   if(i == counter)
	   {
		   //printf("selecionados[%d]: Pos %d\n", counter, new_rand);
		   selecionados[i] = new_rand;
		   selecionados_cluster[i] = new_rand_cluster;
		   counter++;
	   }
	}

	memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	int starter = size1 * len1;

	int j;

	for(int i = 0; i < QUANTIDADE_PARA_PERMUTAR; i++)
	{
		//find dad and mom cluster in graph
		int dad_cluster = 0;
		int mom_cluster = 0;
		int clusterLength = len1;
		int clusterInitial1 = 0;
		int clusterInitial2 = 0;

		int selected = 0;
		int counter = 0;
		for(; selected < size1; ++selected)
		{
			for(j = 0; j < len1; ++j)
			{
				if(pai.genes[selecionados[selected]]->index() == mae.genes[counter]->index())
				{
					mom_cluster = selected;
					clusterLength = len1;
					clusterInitial1 = len1 * selected;
				}
				if(mae.genes[selecionados[selected]]->index() == pai.genes[counter]->index())
				{
					dad_cluster = selected;
					clusterLength = len1;
					clusterInitial2 = len1 * selected;
				}
				counter++;
			}
		}

		int cluster_counter = 0;
		for(; selected < k; ++selected)
		{
			for(j = 0; j < len2; ++j)
			{
				if(pai.genes[selecionados[selected]]->index() == mae.genes[counter]->index())
				{
					mom_cluster = selected;
					clusterLength = len2;
					clusterInitial1 = starter + len2 * cluster_counter;
				}
				if(mae.genes[selecionados[i]]->index() == pai.genes[counter]->index())
				{
					dad_cluster = selected;
					clusterLength = len2;
					clusterInitial2 = starter + len2 * cluster_counter;
				}
				counter++;
			}
			cluster_counter++;
		}

		int rand_from_cluster = rand() % clusterLength;

		node temp = filho1.genes[selecionados[i]];
		filho1.genes[selecionados[i]] = filho1.genes[clusterInitial1 + rand_from_cluster];
		filho1.genes[clusterInitial1 + rand_from_cluster] = temp;

		temp = filho2.genes[selecionados[i]];
		filho2.genes[selecionados[i]] = filho2.genes[clusterInitial2 + rand_from_cluster];
		filho2.genes[clusterInitial2 + rand_from_cluster] = temp;
	}
}

Cromossomo* genetic(Graph& graph, node* initial_array, int& size, int& k, int& size1, int& len1, int& len2, NodeArray<NodeArray<int>>& dist)
{
	double sum_rate = 0;
	int best_initial = 0;

	Cromossomo* best_cromossomo = new Cromossomo;
	best_cromossomo->genes = new node[size];
	best_cromossomo->fitness = 0;

	//initialize cromossomo
	Cromossomo* cromossomos = new Cromossomo[QUANTIDADE_CROMOSSOMOS];
	Cromossomo* filhos = new Cromossomo[QUANT_FILHOS];
	for(int i = 0; i < QUANTIDADE_CROMOSSOMOS; ++i)
	{
		cromossomos[i].genes = new node[size];

		//shuffle(initial_array, size, sizeof(node));

		memcpy(cromossomos[i].genes, initial_array, size * sizeof(node));

		shuffle(initial_array, size, sizeof(node));

		//alterar o fitness
		//fitness_cromossomo(cromossomos[i], graph, k, size1, len1, len2);
		fitness_cromossomo_somatorio(cromossomos[i], size, k, size1, len1, len2, dist);

		//printf("---Rate:%.3f\t", cromossomos[i].rate);
		sum_rate += cromossomos[i].fitness;
		if(cromossomos[i].fitness > best_cromossomo->fitness)
		{
			best_cromossomo->fitness = cromossomos[i].fitness;
			best_initial = i;
		}
	}

	//printf("melhor fitness inicial: %.4lf\n", best_cromossomo->rate);
	memcpy(best_cromossomo->genes, cromossomos[best_initial].genes, sizeof(node) * size);

	//initialize filhos, QUANTIDADE PAR
	for(int i = 0; i < QUANT_FILHOS; ++i)
	{
		filhos[i].genes = new node[size];
	}

	int iter = 0;
	int pai_index;
	int mae_index;

	//parte um, onde o crossing é mais violento
	while(iter < MAX_ITER/2)
	{
		//printf("\nIteracao: %d\n", iter);

		pai_index = roleta(cromossomos, sum_rate);
		//printf("\nPai da roleta: %d\t", pai_index);

		//jogar o pai para o inicio da lista, para tira-lo da roleta
		sum_rate -= cromossomos[pai_index].fitness;
		Cromossomo temp = cromossomos[pai_index];
		cromossomos[pai_index] = cromossomos[0];
		cromossomos[0] = temp;


		mae_index = roleta(&cromossomos[1], sum_rate);
		//printf("Mae da roleta: %d\n", mae_index);
		//retornar ao sum_rate original
		sum_rate += cromossomos[0].fitness;

		//printf("\n\n---Filhos---\n");

		//cria os filhos
		for(int i = 0; i < QUANT_FILHOS; i++)
		{
			//printf("\nFilho: %d\t", i);
			//create_filho_spread_same_clusters_rand(filhos[i], filhos[i], cromossomos[0], cromossomos[mae_index], size, k, size1, len1, len2);
			create_filho_spread_same_clusters_block_rand(filhos[i], filhos[i], cromossomos[0], cromossomos[mae_index], size, k, size1, len1, len2);

			//mutacao
			if(RANDOM < 0.1)
			{
				int new_rand1 = rand() % size;
				int new_rand1_cluster;

				if(new_rand1 < size1 * len1)
				{
					 new_rand1_cluster = new_rand1 / len1;
				}
				else
				{
					new_rand1_cluster = ((new_rand1 - (size1 * len1))/len2) + size1;
				}

				int new_rand2;
				int new_rand2_cluster;

				do
				{
					new_rand2 = rand() % size;
					if(new_rand2 < size1 * len1)
					{
						new_rand2_cluster = new_rand2 / len1;
					}
					else
					{
						new_rand2_cluster = ((new_rand2 - (size1 * len1))/len2) + size1;
					}
				}while(new_rand1_cluster == new_rand2_cluster);

				node temp = filhos[i].genes[new_rand1];
				filhos[i].genes[new_rand1] = filhos[i].genes[new_rand2];
				filhos[i].genes[new_rand2] = temp;

				/*temp = filhos[i+1].genes[new_rand1];
				filhos[i+1].genes[new_rand1] = filhos[i+1].genes[new_rand2];
				filhos[i+1].genes[new_rand2] = temp;*/
			}

			fitness_cromossomo_somatorio(filhos[i], size, k, size1, len1, len2, dist);
			//fitness_cromossomo(filhos[i], graph, k, size1, len1, len2);
			//fitness_cromossomo(filhos[i+1], graph, k, size1, len1, len2);
			//printf("Fitness: %.3lf\n", filhos[i].rate);

			if(filhos[i].fitness > best_cromossomo->fitness)
			{
				memcpy(best_cromossomo->genes, filhos[i].genes, sizeof(node) * size);
				best_cromossomo->fitness = filhos[i].fitness;
				//printf("achou filho melhor\n");
			}
			/*if(filhos[i+1].fitness > best_cromossomo->fitness)
			{
				memcpy(best_cromossomo->genes, filhos[i+1].genes, sizeof(node) * size);
				best_cromossomo->fitness = filhos[i+1].fitness;
				//printf("achou filho melhor\n");
			}*/
		}

		//printf("Troca de Filhos\n");
		trocar_filhos(cromossomos, filhos, sum_rate, size, FILHOS_TROCA);

		//printf("\n\n");

		++iter;
	}

	while(iter < MAX_ITER)
	{
		//printf("\nIteracao: %d\n", iter);

		pai_index = roleta(cromossomos, sum_rate);
		//printf("\nPai da roleta: %d\t", pai_index);

		//jogar o pai para o inicio da lista, para tira-lo da roleta
		sum_rate -= cromossomos[pai_index].fitness;
		Cromossomo temp = cromossomos[pai_index];
		cromossomos[pai_index] = cromossomos[0];
		cromossomos[0] = temp;


		mae_index = roleta(&cromossomos[1], sum_rate);
		//printf("Mae da roleta: %d\n", mae_index);
		//retornar ao sum_rate original
		sum_rate += cromossomos[0].fitness;

		//printf("\n\n---Filhos---\n");

		//cria os filhos
		for(int i = 0; i < QUANT_FILHOS; i+=2)
		{
			//printf("\nFilho: %d\t", i);
			create_filho_shuffle_between_clusters_rand(filhos[i], filhos[i+1], cromossomos[0], cromossomos[mae_index], size, k, size1, len1, len2);


			//mutacao
			if(RANDOM < 0.1)
			{
				int new_rand1 = rand() % size;
				int new_rand1_cluster;

				if(new_rand1 < size1 * len1)
				{
					 new_rand1_cluster = new_rand1 / len1;
				}
				else
				{
					new_rand1_cluster = ((new_rand1 - (size1 * len1))/len2) + size1;
				}

				int new_rand2;
				int new_rand2_cluster;

				do
				{
					new_rand2 = rand() % size;
					if(new_rand2 < size1 * len1)
					{
						new_rand2_cluster = new_rand2 / len1;
					}
					else
					{
						new_rand2_cluster = ((new_rand2 - (size1 * len1))/len2) + size1;
					}
				}while(new_rand1_cluster == new_rand2_cluster);

				node temp = filhos[i].genes[new_rand1];
				filhos[i].genes[new_rand1] = filhos[i].genes[new_rand2];
				filhos[i].genes[new_rand2] = temp;

				temp = filhos[i+1].genes[new_rand1];
				filhos[i+1].genes[new_rand1] = filhos[i+1].genes[new_rand2];
				filhos[i+1].genes[new_rand2] = temp;
			}


			//fitness_cromossomo(filhos[i], graph, k, size1, len1, len2);
			fitness_cromossomo_somatorio(filhos[i], size, k, size1, len1, len2, dist);
			//fitness_cromossomo(filhos[i+1], graph, k, size1, len1, len2);
			fitness_cromossomo_somatorio(filhos[i], size, k, size1, len1, len2, dist);
			//printf("Fitness: %.3lf\n", filhos[i].rate);

			if(filhos[i].fitness > best_cromossomo->fitness)
			{
				memcpy(best_cromossomo->genes, filhos[i].genes, sizeof(node) * size);
				best_cromossomo->fitness = filhos[i].fitness;
				//printf("achou filho melhor\n");
			}
			if(filhos[i+1].fitness > best_cromossomo->fitness)
			{
				memcpy(best_cromossomo->genes, filhos[i+1].genes, sizeof(node) * size);
				best_cromossomo->fitness = filhos[i+1].fitness;
				//printf("achou filho melhor\n");
			}
		}

		//printf("Troca de Filhos\n");
		trocar_filhos(cromossomos, filhos, sum_rate, size, FILHOS_TROCA);

		//printf("\n\n");

		++iter;
	}

	//printf("Iter: %d\n", iter);

	return best_cromossomo;
}


//				FIM DA PARTE DO GENETICO				//



void aCluster(ClusterGraph& CG, ClusterGraphAttributes& CGA, SList<node>& vc, cluster& parent, const Color::Name col) {
	cluster c = CG.createCluster(vc, parent);
	CGA.strokeColor( c ) = Color( col );
	CGA.fillColor( c ) = Color( col );
	CGA.fillBgColor( c ) = Color( col );

	CGA.height( c ) = 15000.0;
	CGA.width( c ) = 15000.0;
	CGA.strokeWidth( c ) = 350.0;

	for(node v : c->nodes) {
		CGA.fillColor( v ) = Color( col );
		CGA.fillBgColor( v ) = Color( col );

		CGA.height( v ) = 1500.0;
		CGA.width( v ) = 1500.0;
		CGA.strokeWidth( v ) = 1.0;
		CGA.shape(v) = ogdf::Shape::Ellipse;

//			CGA.label( v ) = std::to_string(v->index()).c_str();
	}
	vc.clear();
}

void theClusters(ClusterGraph& CG, SList<node> vc[], int nclus) {
    int i = 0, c = 0;
    int nnodes = CG.constGraph().numberOfNodes();
    for(node v : CG.constGraph().nodes) {
    	if (i++ == (nnodes * (c+1))/nclus)
    		c++;
    	vc[c].pushBack(v);
    }
}

void bfs_pred_matrix(Graph& G, NodeArray<NodeArray<node>>& pred, NodeArray<NodeArray<int>>& dist)
{
	NodeArray<bool>visitados = NodeArray<bool>(G);

	for(node iter_node = G.firstNode(); iter_node; iter_node = iter_node->succ())
	{
		visitados.fill(false);

	    Queue<node> node_q;  // init q for store currently discovered nodes.
	    Queue<int> depth_q;  // init q for store currently discovered nodes.
	    /*
	    q.append(initial_node);
	    visitados[initial_node] = true;
	    */

	    node_q.append(iter_node);
	    depth_q.append(1);
	    visitados[iter_node] = true;

	    //pred[initial_node] = NULL;
	    pred[iter_node][iter_node] = NULL;
	    dist[iter_node][iter_node] = 0;

	    while(node_q.size()!=0)
	    {
	        node currentNode = node_q.pop();    //pop elements from the queue
	        int currentDepth = depth_q.pop();    //pop elements from the queue
	        adjEntry adj = currentNode->firstAdj();

	        for(; adj; adj = adj->succ())
			{
				node adjNode = adj->twinNode();
	            if(visitados[adjNode] == false)
	            {
	                node_q.append(adjNode);
	                depth_q.append(currentDepth + 1);
	                visitados[adjNode] = true;
	                //pred[adjNode] = currentNode;
	                pred[iter_node][adjNode] = currentNode;
	                dist[iter_node][adjNode] = currentDepth;
	            }
	        }

	        //printf("%d ", currentNode->index());
	    }
	}
}

void bfs_pred_array(Graph& G, node initial_node, NodeArray<node>& pred)
{
	NodeArray<bool>visitados = NodeArray<bool>(G);
	visitados.fill(false);

    Queue<node> q;  // init q for store currently discovered nodes.
    q.append(initial_node);
    visitados[initial_node] = true;
    pred[initial_node] = NULL;

    while(q.size()!=0)
    {
        node currentNode = q.pop();    //pop elements from the q
        adjEntry adj = currentNode->firstAdj();

        for(; adj; adj = adj->succ())
		{
			node adjNode = adj->twinNode();
            if(visitados[adjNode] == false)
            {
                q.append(adjNode);
				visitados[adjNode] = true;
                pred[adjNode] = currentNode;
            }
        }

        //printf("%d ", currentNode->index());
    }
}

int partitionEdges(Graph& G, ClusterGraph& cg, SList<node>& sub_graph)
{
	int edge_count = 0;
	NodeArray<bool>visitados = NodeArray<bool>(G);
    visitados.fill(false);
    NodeArray<bool>visitados_total = NodeArray<bool>(G);
	visitados_total.fill(false);

    SListIterator<node> iter_node = sub_graph.begin();
	visitados[*iter_node] = true;
    int visited_nodes = 1;
    int list_size = sub_graph.size();
    while(visited_nodes < list_size)
    {
        visitados.fill(false);
        Queue<node> q;
        q.append(*iter_node);
        while(q.size()!=0)
        {
            node currentNode = q.pop();    //pop elements from the q

			if(visitados[currentNode] == true)
			{
				continue;
			}

			visitados[currentNode] = true;
			visitados_total[currentNode] = true;
			visited_nodes++;
            adjEntry adj = currentNode->firstAdj();

            for(; adj; adj = adj->succ())
    		{
    			node adjNode = adj->twinNode();
				//printf("Current: %d, Twin: %d\n", currentNode->index(), adjNode->index());
                if(cg.clusterOf(*iter_node) == cg.clusterOf(adjNode))
                {
                    if(visitados[adjNode] == false)
                    {
						++edge_count;
                        q.append(adjNode);
                    }
                }
            }
        }

        iter_node = sub_graph.begin();
		while(visited_nodes < list_size && iter_node.valid() && visitados_total[*iter_node] == true)
        {
        	iter_node = iter_node.succ();
        }
	}

	return edge_count;
}

double sub_graph_density(Graph& g, ClusterGraph& cg, SList<node>& subGraph)
{
	int pEdgeCount = partitionEdges(g, cg, subGraph);
	int nodeQuant = subGraph.size();

    printf("Edges: %d - Nodes: %d - ", pEdgeCount, nodeQuant);

	return ((double)(2 * pEdgeCount)) / ( (nodeQuant) * (nodeQuant - 1) );
}

void partition_density(Graph& g, ClusterGraph& cg, SList<node>* partition, int& nClusters)
{
    printf("\n");
	for(int i = 0; i < nClusters; ++i)
	{
		printf("P[%d]: ", i);
        printf("Group density: %lf\n", i, sub_graph_density(g, cg, partition[i]));
	}

    return;
}

int path_length(ClusterGraph cg, node source, node target, NodeArray<NodeArray<node>> pred, NodeArray<NodeArray<int>> dist, int n)
{
	node actual_node = pred[source][target];

	while(actual_node->index() != source->index())
	{
		if(cg.clusterOf(actual_node) != cg.clusterOf(source))
		{
			return n;
		}

		//fazer o que precisar com o path
		actual_node = pred[source][actual_node];
	}

	return dist[source][target];
}

int main(int argc, char *argv[])
{
	int k = 3;

	if(argc > 2)
	{
		k = atoi(argv[2]);
	}

	if(k <= 1)
	{
		printf("Valor Inválido para K\n");
		return 0;
	}
	//////////// Graphs and clusters
	clock_t begin = clock();

	srand(time(NULL));

	Graph G;
	//FILE* file = fopen("random/teste.clq", "r");

	readJostleGraph(G, fopen(argv[1], "r"));
	//readDimacsGraph(G, file);

	printf("K: %d\n", k);

	int n = G.numberOfNodes();
	int m = G.numberOfEdges();

	printf("N: %d, M: %d", n, m);

	/*List<node> list = List<node>();
	G.allNodes(list);

	int size1 = n % k;
	int size2 = k - size1;

	int len1 = ceil(n / (double)k);
	int len2 = floor(n / (double)k);

	node* initial_node_array = new node[n];

	printf("Len1: %d, Len2: %d, Size1: %d, Size2: %d\n", len1, len2, size1, size2);

	//NodeArray<node> pred = NodeArray<node>(G);
	NodeArray<NodeArray<node>> pred = NodeArray<NodeArray<node>>(G);
	NodeArray<NodeArray<int>> dist = NodeArray<NodeArray<int>>(G);
	int i;
	int j;
	i = 0;
	for(node iter_node = G.firstNode(); iter_node; iter_node = iter_node->succ())
	{
		initial_node_array[i] = iter_node;
		pred[iter_node] = NodeArray<node>(G);
		dist[iter_node] = NodeArray<int>(G);
		i++;
	}

	bfs_pred_matrix(G, pred, dist);*/

	/*int starter = size1 * len1;
	ListIterator<node> highest_degree_node;
	i = 0;
	for(; i < size1; ++i)
	{
		highest_degree_node = list.begin();

		for(ListIterator<node> iter = highest_degree_node.succ(); iter.valid(); iter = iter.succ())
		{
			//printf("iter: %d-", (*iter)->index());
			if((*iter)->degree() > (*highest_degree_node)->degree())
			{
				highest_degree_node = iter;
			}
		}
		initial_node_array[i * len1] = *highest_degree_node;
		list.del(highest_degree_node);
	}
	int second_part_counter = 0;
	for(; i < k; ++i)
	{
		highest_degree_node = list.begin();

		for(ListIterator<node> iter = highest_degree_node.succ(); iter.valid(); iter = iter.succ())
		{
			//printf("iter: %d-", (*iter)->index());
			if((*iter)->degree() > (*highest_degree_node)->degree())
			{
				highest_degree_node = iter;
			}
		}
		initial_node_array[starter + second_part_counter * len2] = *highest_degree_node;
		list.del(highest_degree_node);
		second_part_counter++;
	}*/

	//para viciar o primeiro cromossomo
	/*int counter_initial = 1;
	i = 0;
	for(; i < size1; ++i)
	{
		//verif
		node partition_first = initial_node_array[i * len1];

		for(j = 1; j < len1; ++j)
		{
			ListIterator<node> new_node = list.begin();
			int closer_distance = dist[partition_first][*new_node];

			for(ListIterator<node> iter = new_node.succ(); iter.valid(); iter = iter.succ())
			{
				if(dist[partition_first][*iter] < closer_distance)
				{
					new_node = iter;
					closer_distance = dist[partition_first][*iter];
				}
			}

			initial_node_array[counter_initial] = *new_node;
			list.del(new_node);
			++counter_initial;
		}
		//jump highest degree nodes pre selected
		++counter_initial;
	}

	second_part_counter = 0;
	for(; i < k; ++i)
	{
		//verif
		node partition_first = initial_node_array[starter + second_part_counter * len2];

		for(j = 1; j < len2; ++j)
		{
			ListIterator<node> new_node = list.begin();
			int closer_distance = dist[partition_first][*new_node];

			for(ListIterator<node> iter = new_node.succ(); iter.valid(); iter = iter.succ())
			{
				if(dist[partition_first][*iter] < closer_distance)
				{
					new_node = iter;
					closer_distance = dist[partition_first][*iter];
				}
			}

			initial_node_array[counter_initial] = *new_node;
			list.del(new_node);
			++counter_initial;
		}

		++counter_initial;
		++second_part_counter;
	}

	Cromossomo* best_cromo = new Cromossomo;
	best_cromo->genes = new node[n];
	memcpy(best_cromo->genes, initial_node_array, sizeof(node) * n);*/

	/*Cromossomo* best_cromo = genetic(G, initial_node_array, n, k, size1, len1, len2, dist);
	clock_t end = clock();
	fitness_cromossomo_somatorio(*best_cromo, n, k, size1, len1, len2, dist);
	printf("best fitness: %.6lf\n", best_cromo->fitness);

	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time Spent: %.3lf segundos\n", time_spent);

	//printf("\n");

	list = List<node>();

	G.allNodes(list);

	//node initial_node = G.firstNode()->succ()->succ()->succ()->succ();

	SList<node> partition[k];

	int counter = 0;

	i = 0;
	for(; i < size1; ++i)
	{
		for(j = 0; j < len1; ++j)
		{
			partition[i].pushBack(best_cromo->genes[counter]);
			counter++;
		}
	}

	for(; i < k; ++i)
	{
		for(j = 0; j < len2; ++j)
		{
			partition[i].pushBack(best_cromo->genes[counter]);
			counter++;
		}
	}

	//print partitions
	for(i = 0; i < k; i++)
	{
		printf("Part %d -> ( ", i);
		for(SListIterator<node> iter = partition[i].begin(); iter.valid(); iter = iter.succ())
		{
			printf("%d ", (*iter)->index());
		}

		printf(")\n");
	}*/

	ClusterGraph cg = ClusterGraph(G);

	//////////// Attributes

	ClusterGraphAttributes CGA(cg);
	CGA.addAttributes( GraphAttributes::nodeLabel |
	        GraphAttributes::nodeStyle |
	        GraphAttributes::edgeType |
	        GraphAttributes::edgeArrow |
	        GraphAttributes::edgeStyle);

	//create cluster
	/*for(int i = 0; i < k; i++)
	{
		cluster c = cg.createCluster(partition[i]);

		CGA.strokeColor( c ) = Color( clustCol[i] );
		CGA.fillColor( c ) = Color( clustCol[i] );
		CGA.fillBgColor( c ) = Color( clustCol[i] );

		CGA.height( c ) = 15000.0;
		CGA.width( c ) = 15000.0;
		CGA.strokeWidth( c ) = 350.0;

		for(node v : c->nodes) {
			CGA.fillColor( v ) = Color( clustCol[i] );
			CGA.fillBgColor( v ) = Color( clustCol[i] );

			CGA.height( v ) = 1500.0;
			CGA.width( v ) = 1500.0;
			CGA.strokeWidth( v ) = 1.0;
			CGA.shape(v) = ogdf::Shape::Ellipse;
			CGA.label( v ) = std::to_string(v->index()).c_str();
		}
	}*/

    printf("Total Graph Density: %lf\n", ((double)(2 * G.numberOfEdges() )) / ( G.numberOfNodes() * (G.numberOfNodes() - 1)));

	//partition_density(G, cg, partition, k);

	//fclose(file);


	//int nclusters = 2; // must be <= nmaxClusters
	//ClusterGraph CG(G);
    //SList<node> vc[nclusters];
    //theClusters(CG, vc, nclusters);

	//////////// Attributes

	/*ClusterGraphAttributes CGA(CG);
	CGA.addAttributes( GraphAttributes::nodeLabel |
	        GraphAttributes::nodeStyle |
	        GraphAttributes::edgeType |
	        GraphAttributes::edgeArrow |
	        GraphAttributes::edgeStyle);

	for (int c = 0; c < k; c++)
	{
		aCluster()
	}*/

    for(edge e : G.edges)
    {
        CGA.bends(e);
        CGA.strokeWidth(e) = 100.0;
        CGA.arrowType(e) = ogdf::EdgeArrow::None;
        CGA.strokeColor(e) = Color( Color::Name::Gray );

    }

	GraphAttributes GA( G );
	//GraphAttributes GA( CGA );

	//////////// Layouts

	/*SugiyamaLayout SL;
	SL.setRanking( new LongestPathRanking );
    SL.setCrossMin( new MedianHeuristic );
    FastSimpleHierarchyLayout *fshl = new FastSimpleHierarchyLayout;
    SL.setLayout( fshl );
    OptimalHierarchyClusterLayout *ohl = new OptimalHierarchyClusterLayout;
    SL.setClusterLayout( ohl );*/

    ClusterPlanarizationLayout cpl;

//    GEMLayout GEM;

    FMMMLayout fmmm;
    fmmm.useHighLevelOptions(true);
    fmmm.unitEdgeLength(15.0);
    fmmm.newInitialPlacement(true);
    fmmm.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

	//////////// Output file

	string gname(argv[1]);
	int b = gname.find_last_of('/');
	int l = gname.find_last_of(".");
	std::filebuf fb;

	//fmmm.call( GA );

	//////////// Graph drawing

//    SL.call( GA );

    fmmm.call( GA );
//    GEM.call(GA);

    GA.scale(8,false);
	fb.open (string("output-")+gname.substr(b+1,l-b-1)+string(".svg"),std::ios::out);
	std::ostream os(&fb);
	GraphIO::drawSVG(GA, os);
	fb.close();
	cout << string("output-")+gname.substr(b+1,l-b-1)+string(".svg") << endl;

	//////////// Clustered graph drawing

    //fmmm.call( CGA );

    /*cpl.call( G, CGA, cg );

	//SL.call( CGA );

    CGA.scale(4,false);
	fb.open (string("output-clustered-")+gname.substr(b+1,l-b-1)+string(".svg"),std::ios::out);
	std::ostream osc(&fb);
	GraphIO::drawSVG(CGA, osc);
	fb.close();
	cout << string("output-clustered-")+gname.substr(b+1,l-b-1)+string(".svg") << endl;*/

	return 0;
}
