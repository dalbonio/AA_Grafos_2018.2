#include "genetic_algorithm_partition.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
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


#define EPSILON 0.0000000000000000001
#define QUANTIDADE_CROMOSSOMOS 400
#define QUANT_FILHOS 10
#define FILHOS_TROCA 2
#define MAX_ITER 20000
#define MIN_CROSSING_OVER_RATE 20/(double)100
#define MAX_CROSSING_OVER_RATE 80/(double)100

#define RANDOM ((double)rand() / RAND_MAX)

namespace Genetics
{
	//ok
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

	double sub_graph_density_array(Graph& g, node* subGraph, int k, int length)
	{
		int pEdgeCount = partitionEdgesArray(g, subGraph, length);
		int nodeQuant = length;

	    printf("Edges: %d - Nodes: %d - ", pEdgeCount, nodeQuant);

		return ((double)(2 * pEdgeCount)) / ( (nodeQuant) * (nodeQuant - 1) );
	}

	//genetic one
	double partition_density_array(Graph& g, node* partition, int& k, int size1, int len1, int len2)
	{
	    printf("\n");
		double sum = 0;
		int i;
		for(i = 0; i < size1; ++i)
		{
			//printf("P[%d]: ", i);
	        //printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[i * len1], k, len1));
			sum += sub_graph_density_array(g, &partition[i * len1], k, len1)
		}
		int starter = size1 * len1;
		int counter = 0;
		for(; i < k; ++i)
		{
			//printf("P[%d]: ", i);
	        //printf("Group density: %lf\n", i, sub_graph_density_array(g, &partition[starter + counter * len2], k, len2));
			sum += sub_graph_density_array(g, &partition[starter + counter * len2], k, len2);
			counter++;
		}

	    return sum;
	}

	double* initialize_genes(node* initial_array, int size)
	{
		node* genes = new node[size];

		shuffle(initial_array, size, sizeof(node));

		memcpy(genes, initial_array, size);

		return genes;
	}

	void fitness_cromossomo(Cromossomo& cromossomo, Graph& g, int& k, int size1, int len1, int len2)
	{
		cromossomo.fitness = partition_density_array(g, cromossomo.genes, k, size1, len1, len2);
	}

	int roleta(Cromossomo* cromossomos, double sum_rate)
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

	void trocar_filhos(Cromossomo* cromossomos, Cromossomo* filhos, double& sum_rate, int size)
	{
		int worst_pai;
		int best_filho;

		for(int j = 0; j < FILHOS_TROCA; ++j)
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

	void create_filho_shuffle_between_clusters(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int size, int k, int size1, int len1, int len2)
	{
		int* selecionados = new int[QUANTIDADE_PARA_PERMUTAR];
   		int* selecionados_cluster = new int[QUANTIDADE_PARA_PERMUTAR];
		int* new_selecionados = new int[QUANTIDADE_PARA_PERMUTAR];

		int counter = 0;

		printf("\n");
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

		   printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

		   if(i == counter)
		   {
			   printf("selecionados[%d]: Pos %d\n", counter, new_rand);
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


	void create_filho_particao_order_based(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int size, int k, int size1, int len1, int len2)
	{
	    //rand from minimal to maximal crossing rate( random interval is [0,1] )


	    int* selecionados = new int[QUANTIDADE_PARA_PERMUTAR];
		int* selecionados_cluster = new int[QUANTIDADE_PARA_PERMUTAR];
	    int* mom_pos = new int[QUANTIDADE_PARA_PERMUTAR];
	    int* dad_pos = new int[QUANTIDADE_PARA_PERMUTAR];

	    int counter = 0;

	    printf("\n");
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

			printf("new_rand: %d, new_rand_cluster: %d\n", new_rand, new_rand_cluster);

	        if(i == counter)
	        {
	            printf("selecionados[%d]: Pos %d\n", counter, new_rand);
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
	            //if(pai.genes[j]->index() == mae.genes[selecionados[i]]->index())
	            if(pai.genes[j] == mae.genes[selecionados[i]])
	            {
	                printf("dad_pos[%d]: %d\n", i, j);
	                dad_pos[i] = j;
	                found += 1;
	            }
	            //if(mae.genes[j]->index() == pai.genes[selecionados[i]]->index())
	            if(mae.genes[j] == pai.genes[selecionados[i]])
	            {
	                printf("mom_pos[%d]: %d\n", i, j);
	                mom_pos[i] = j;
	                found += 1;
	            }

	            if(found == 2)
	            {
	                break;
	            }
	        }
	    }

	    //memcpy(filho1.genes, pai.genes, sizeof(node) * size);
	    //memcpy(filho2.genes, mae.genes, sizeof(node) * size);

	    memcpy(filho1.genes, pai.genes, sizeof(int) * size);
	    memcpy(filho2.genes, mae.genes, sizeof(int) * size);

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

	        printf("Menor mae: %d, Menor pai: %d\n", menor_mae, menor_pai);

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

	//USAR ESTRATEGIA 1 PARA SEPARAR OS NODES NA HORA DE FAZER O FITNESS COM AS DENSIDADES

	Cromossomo* genetic(Graph& graph, int* initialArray, int size, int k, int size1, int len1, int len2)
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
			cromossomo[i].genes = new node[size];

			shuffle(initial_array, size, sizeof(node));
			shuffle(initial_array, size, sizeof(node));

			memcpy(cromossomo[i].genes, initial_array, size);

			//alterar o fitness
			fitness_cromossomo(cromossomo[i], graph, k, size1, len1, len2);

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
			for(int i = 0; i < QUANT_FILHOS; i += 2)
			{
	            //printf("\nFilho: %d\t", i);
				create_filho_particao_order_based(filhos[i], filhos[i+1], cromossomos[0], cromossomos[mae_index], size, k, size1, len1, len2);

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

                    int temp = filhos[i].genes[new_rand1];
					filhos[i].genes[new_rand1] = filhos[i].genes[new_rand2];
                    filhos[i].genes[new_rand2] = temp;

                    temp = filhos[i+1].genes[new_rand1];
                    filhos[i+1].genes[new_rand1] = filhos[i+1].genes[new_rand2];
                    filhos[i+1].genes[new_rand2]; = temp;
				}

				fitness_cromossomo(filhos[i], graph, k, size1, len1, len2);
                fitness_cromossomo(filhos[i+1], graph, k, size1, len1, len2);
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
			trocar_filhos(cromossomos, filhos, sum_rate, size);

			//printf("\n\n");

			++iter;
		}

		//printf("Iter: %d\n", iter);

		return best_cromossomo;
	}
}
