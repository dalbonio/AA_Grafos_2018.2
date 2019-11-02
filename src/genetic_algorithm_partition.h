#ifndef GENETIC_ALGORITHM_PARTITION_H
#define GENETIC_ALGORITHM_PARTITION_H

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

namespace Genetics
{
    struct Cromossomo
    {
        double fitness;

        node* genes;
    };

    void shuffle(void *array, size_t n, size_t size);

    double sub_graph_density_array(Graph& g, node* subGraph, int k, int length);

    double partition_density_array(Graph& g, node* partition, int& k, int size1, int len1, int len2);

    double* initialize_genes(node* initial_array, int size);

    void fitness_cromossomo(Cromossomo& cromossomo, Graph& g, int& k, int size1, int len1, int len2);

    int roleta(Cromossomo* cromossomos, double sum_rate);

    void trocar_filhos(Cromossomo* cromossomos, Cromossomo* filhos, double& sum_rate, int size);

    void create_filho_particao_order_based(Cromossomo& filho1, Cromossomo& filho2, Cromossomo& pai, Cromossomo& mae, int size, int k, int size1, int len1, int len2);

    Cromossomo* genetic(Graph& graph, int* initialArray, int size, int k, int size1, int len1, int len2);

}


#endif
