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
	//////////// Graphs and clusters

	Graph G;
	//FILE* file = fopen("random/teste.clq", "r");

	readDimacsGraph(G, fopen(argv[1], "r"));
	//readDimacsGraph(G, file);

	int k = 3;

	printf("K: %d\n", k);

	if(argc > 2)
	{
		k = atoi(argv[2]);
	}

	int n = G.numberOfNodes();

	int size1 = n % k;
	int size2 = k - size1;

	int len1 = ceil(n / (double)k);
	int len2 = floor(n / (double)k);

	node* initial_node_array = new node[n];

	printf("Len1: %d, Len2: %d, Size1: %d, Size2: %d\n", len1, len2, size1, size2);

	//NodeArray<node> pred = NodeArray<node>(G);
	NodeArray<NodeArray<node>> pred = NodeArray<NodeArray<node>>(G);
	NodeArray<NodeArray<int>> dist = NodeArray<NodeArray<int>>(G);

	int i = 0;
	for(node iter_node = G.firstNode(); iter_node; iter_node = iter_node->succ())
	{
		initial_node_array[i] = iter_node;
		pred[iter_node] = NodeArray<node>(G);
		dist[iter_node] = NodeArray<int>(G);
		i++;
	}

	bfs_pred_matrix(G, pred, dist);

	//printf("\n");

	List<node> list = List<node>();

	G.allNodes(list);

	//node initial_node = G.firstNode()->succ()->succ()->succ()->succ();

	SList<node> partition[k];

	int j;
	/*for(i = 0; i < k; i++)
	{
		ListIterator<node> highest_degree_node = list.begin();

		for(ListIterator<node> iter = highest_degree_node.succ(); iter.valid(); iter = iter.succ())
		{
			//printf("iter: %d-", (*iter)->index());
			if((*iter)->degree() > (*highest_degree_node)->degree())
			{
				highest_degree_node = iter;
			}
		}

		printf("\n");

		partition[i] = SList<node>();
		partition[i].pushFront(*highest_degree_node);

		//printf("%d\n", *highest_degree->index);

		list.del(highest_degree_node);
	}

	i = 0;
	for(; i < size1; ++i)
	{
		//verif
		node partition_first = partition[i].front();
		//printf("P_First: %d\n", partition_first->index());

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

			printf("new_node: %d\n", (*new_node)->index());
			partition[i].pushBack(*new_node);
			//printf("%d ", (*new_node)->index());
			list.del(new_node);
		}
	}

	for(; i < k; ++i)
	{
		//verif
		node partition_first = partition[i].front();
		//printf("P_First: %d\n", partition_first->index());

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

			printf("new_node: %d\n", (*new_node)->index());
			partition[i].pushBack(*new_node);
			//printf("%d ", (*new_node)->index());
			list.del(new_node);
		}
	}*/

	//dividir de um jeito diferente

	i = 0;
    ListIterator<node> new_node = list.begin();
    do
    {
        node partition_first = partition[i].front();
		int less_distant = dist[partition_first][*new_node];

        for(ListIterator<node> iter = new_node.succ(); iter.valid(); iter = iter.succ())
        {
            if(dist[partition_first][*iter] < less_distant)
            {
                new_node = iter;
            }
        }

        partition[i].pushBack(*new_node);

        list.del(new_node);

        //reset partition when reach end of iters
        if(++i == k)
        {
            i = 0;
        }

        new_node = list.begin();
    }while( new_node.valid() );
	//fim da divisao

	//print partitions
	for(i = 0; i < k; i++)
	{
		printf("Part %d -> ( ", i);
		for(SListIterator<node> iter = partition[i].begin(); iter.valid(); iter = iter.succ())
		{
			printf("%d ", (*iter)->index());
		}

		printf(")\n");
	}

	ClusterGraph cg = ClusterGraph(G);

	//////////// Attributes

	ClusterGraphAttributes CGA(cg);
	CGA.addAttributes( GraphAttributes::nodeLabel |
	        GraphAttributes::nodeStyle |
	        GraphAttributes::edgeType |
	        GraphAttributes::edgeArrow |
	        GraphAttributes::edgeStyle);

	//create cluster
	for(int i = 0; i < k; i++)
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
	}

    printf("Total Graph Density: %lf\n", ((double)(2 * G.numberOfEdges() )) / ( G.numberOfNodes() * (G.numberOfNodes() - 1)));

    partition_density(G, cg, partition, k);

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

	GraphAttributes GA( CGA );

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

    GA.scale(16,false);
	fb.open (string("output-")+gname.substr(b+1,l-b-1)+string(".svg"),std::ios::out);
	std::ostream os(&fb);
	GraphIO::drawSVG(GA, os);
	fb.close();
	cout << string("output-")+gname.substr(b+1,l-b-1)+string(".svg") << endl;

	//////////// Clustered graph drawing

    //fmmm.call( CGA );

    //cpl.call( G, CGA, cg );

	//SL.call( CGA );

    /*CGA.scale(4,false);
	fb.open (string("output-clustered-")+gname.substr(b+1,l-b-1)+string(".svg"),std::ios::out);
	std::ostream osc(&fb);
	GraphIO::drawSVG(CGA, osc);
	fb.close();
	cout << string("output-clustered-")+gname.substr(b+1,l-b-1)+string(".svg") << endl;*/

	return 0;
}
