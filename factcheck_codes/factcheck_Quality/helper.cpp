# include <iostream>
# include <fstream>
# include <vector>
# include <cmath>
# include <queue>
# include <cstdlib>
# include <ctime>
# include <set>
# include <algorithm>
//# include <random>
using namespace std;

//int HEPT_seed[] = {287, 639, 131, 267, 100, 559, 608, 124, 3683, 359, 80, 1775, 8, 15, 535, 1292, 247, 76, 274};

struct gnode{
	std::vector<int> outnode;
	std::vector<int> innode;
	std::vector<float> outweight;
	std::vector<float> inweight;
	int indeg;
	int outdeg;
};

struct vosnode{
	std::set<int> vlist;
};

struct vovnode{
	std::vector<int> vlist;
};

typedef struct vosnode VoV;
typedef struct vovnode VoV1;

typedef struct gnode Gnode;
typedef std::vector<Gnode> Graph;

Graph G; int edges=0;

void erase_v3(std::vector<int> &vec, int value)
{
    auto pr = std::equal_range(std::begin(vec), std::end(vec), value);
    vec.erase(pr.first, pr.second);
}

void new_edge(Graph &G, int src, int dest, float val)
{
	if(G.size()<(1+std::max(src, dest)))
		G.resize(1+std::max(src, dest));

	G[src].outnode.push_back(dest);
	G[dest].innode.push_back(src);
	G[src].outweight.push_back(val);
	G[dest].inweight.push_back(val);
	(G[src].outdeg)++;
	(G[dest].indeg)++;
}

void print_outneighbors(Graph &G, int src)
{
	int i;
	Gnode v = G[src];
	for(i=0; i< v.outnode.size(); ++i)
	{
		std::cout<<v.outnode[i]<<" "<<v.outweight[i]<<std::endl;
	}
}

void print_inneighbors(Graph &G, int dest)
{
	int i;
	Gnode v = G[dest];
	for(i=0; i< v.innode.size(); ++i)
	{
		std::cout<<v.innode[i]<<" "<<v.inweight[i]<<std::endl;
	}
}

void load_graph(Graph &G, char *fname, float p, int start_idx, float pval = 0.01)
{
	int src, dest;
	float val;
  	ifstream Gfile (fname);
  	if (Gfile.is_open())
  	{
    	while (!Gfile.eof())
    	{
      		if(p < 0)
      		{
      			Gfile>>src>>dest>>val;
      			if(src==dest)
      				continue;

      			new_edge(G, dest-start_idx, src-start_idx, -val*p);
      		}
      		else if(p > 0 && p<=1)
      		{
      			Gfile>>src>>dest>>val;
      			if(src==dest)
      				continue;

      			new_edge(G, src-start_idx, dest-start_idx, val*p);
      		}
      		else
      		{
      			Gfile>>src>>dest;
      			if(src==dest)
      				continue;
      			new_edge(G, src-start_idx, dest-start_idx, pval);
      		}

      		edges++;
    	}
 	   Gfile.close();
	}
	else
		std::cout<<"Failed to open "<<fname<<std::endl;
}

Graph flip_edges(Graph &G)
{
	Graph G1;
	G1.resize(G.size());
	for(int i=0; i<G.size(); i++)
	{
		for(int j=0; j<G[i].outnode.size(); ++j)
			new_edge(G1, G[i].outnode[j], i, G[i].outweight[j]);
	}
	return G1;
}

Graph reverse_sign(Graph &G)
{
	Graph G1;
	G1.resize(G.size());
	for(int i=0; i<G.size(); i++)
	{
		for(int j=0; j<G[i].outnode.size(); ++j)
			new_edge(G1, i, G[i].outnode[j], -G[i].outweight[j]);
	}
	return G1;
}

void remove_node(Graph &G, int nidx)
{
	int thisi, thiso;
	for(int ins=0; ins<G[nidx].innode.size(); ins++)
	{
		G[nidx].inweight[ins] = 0;
		thisi = G[nidx].innode[ins];

		for(int i=0; i<G[thisi].outnode.size(); i++)
			if(G[thisi].outnode[i]==nidx)
			{
				G[thisi].outweight[i] = 0;
				break;
			}
	}
	for(int outs=0; outs<G[nidx].outnode.size(); outs++)
	{
		G[nidx].outweight[outs] = 0;
		thiso = G[nidx].outnode[outs];

		for(int i=0; i<G[thiso].innode.size(); i++)
			if(G[thiso].innode[i]==nidx)
			{
				G[thiso].inweight[i] = 0;
				break;
			}
	}

}
