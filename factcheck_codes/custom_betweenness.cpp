# include "ICM_unsigned.cpp"
# define MAXFLOAT 100000
# include <functional> 

using namespace std;

void Dijkstra(Graph &G1, int src, vector<int> &pred, vector<double> &dist)
{
	priority_queue<std::pair<double, int>, std::vector<pair<double, int>>, std::greater<pair<double, int>>> q;
	
	int n = G1.size();
	pred.resize(n);
	fill(pred.begin(), pred.end(), -1);
	dist.resize(n);
	vector<int> visit(n);
	for(int i=0; i<n; i++)
	{
		if(i==src)
			dist[i] = 0;
		else
			dist[i] = MAXFLOAT;
		q.push(pair<double, int>(dist[i], i));
	}

	while(!q.empty())
	{
		auto tt = q.top(); int u = tt.second;
		q.pop();
		if(visit[u]==1)
			continue;

		visit[u] = 1;
		for(int j=0; j<G1[u].outnode.size(); j++)
		{
			int v = G1[u].outnode[j];
			if(visit[v]==1)
				continue;
			double new_dist = dist[u]+G1[u].outweight[j];
			if(dist[v] > new_dist)
			{
				dist[v] = new_dist;
				pred[v] = u;
		
				q.push(pair<double, int>(dist[v], v));
			}
		}
	}	
}

vector<double> fixed_BC(Graph G, vector<int> init, vector<double> &dist)
{
	Graph G1;
	int i, j;
	
	for(i=0; i<G.size(); i++)
        {
                for(j=0; j<G[i].outnode.size(); ++j)
                        new_edge(G1, i, G[i].outnode[j], -log(G[i].outweight[j]));
        }
	int newsrc = G.size();
	for(i=0; i<init.size(); i++)
		new_edge(G1, newsrc, init[i], 0);

	vector<int> pred;
//	vector<double> dist;

	Dijkstra(G1, newsrc, pred, dist);
	
	vector<int> leaves(G1.size());
	fill(leaves.begin(), leaves.end(), 1);
	for(i=0; i<pred.size(); i++)
	{
		if(pred[i]==-1)
		{
			leaves[i] = 0;	//No predecessor means not a leaf
			continue;
		}

		leaves[pred[i]] = 0; // Someones predecessor means not a leaf
	}

	vector<double> BC(G1.size()); fill(BC.begin(), BC.end(), 1);
	set<int> prevleaves, newleaves; 
	for(i=0; i<leaves.size(); i++)
	{
		if(leaves[i]==1)
		{
			newleaves.insert(i);
		}
	}

	do{
		prevleaves = newleaves;
		newleaves.clear();
		for(auto t: prevleaves)
		{
			if(pred[t]==-1)
				continue;
			BC[pred[t]] += BC[t];
			newleaves.insert(pred[t]);
		}
	}while(newleaves.size()>0);
	
	for(i=0; i<init.size(); i++)
		BC[init[i]] = 0;
	return BC;
}

void write_BC(char *fname, vector<double> BC)
{
	priority_queue<pair<double, int>> q;

	for(int i=0; i<BC.size()-1; i++)
	        q.push(std::pair<double, int>(BC[i], i));
        
        for (int i=0; i<50; ++i)
        {
                int xx = q.top().second;
                cout<<xx<<endl;
                q.pop();
        }
	cout<<endl;
}

int main(int argc, char *argv[])
{
/*	Graph G;
	new_edge(G, 0, 1, 0.3);
	new_edge(G, 1, 3, 0.2);
	new_edge(G, 0, 2, 1);
	new_edge(G, 2, 1, 0.5);
	new_edge(G, 4, 2, 0.5);
	new_edge(G, 1, 5, 0.9);
	vector<double> dist, BC;
	vector<int> pred;
	vector<int> init; init.push_back(0);
	BC = fixed_BC(G, init, dist);
	for(int i=0; i<G.size(); i++)
		cout<<i<<" "<<dist[i]<<" "<<BC[i]<<endl;*/

        int numseeds=50;
        load_graph(G, argv[1], 1, 0, 0.01);
        std::vector<double> newB, prevB;
        std::cout<<"Loaded graph with "<<G.size()<<" nodes and "<<edges<<" edges"<<std::endl;
        std::vector<int> selseeds;
        vector<double> dist, BC;
        vector<int> pred;

        selseeds = degseeds(G, numseeds);
        std::vector<int> infseeds(numseeds);
        
	BC = fixed_BC(G, selseeds, dist);

	write_BC("", BC);
	return 0;
}

