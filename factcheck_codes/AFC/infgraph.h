typedef pair<double,int> dipair;


#include "iheap.h"

class InfGraph: public Graph
{
private:
//    vector<bool> visit;
//    vector<int> visit_mark;
public:
    vector<vector<int>> hyperG;
    vector<bool> visit;
    vector<int> visit_mark;
    vector<vector<int>> hyperGT;
	vector<set<int>> StoU;
	int StoUsize;
//	void Buildstou(vector<double> stou, int NUMSIM);
	InfGraph():Graph()
	{
		sfmt_init_gen_rand(&sfmtSeed , 95082);

	}
    InfGraph(string folder, string graph_file): Graph(folder, graph_file)
    {
        sfmt_init_gen_rand(&sfmtSeed , 95082);
        init_hyper_graph();
        visit = vector<bool> (n);
        visit_mark = vector<int> (n);
	StoUsize = 100000;
    }


    void init_hyper_graph(){
        hyperG.clear();
        for (int i = 0; i < n; i++)
	{
            hyperG.push_back(vector<int>());
//		StoU.push_back(set<int>());
	}
        hyperGT.clear();

    }

	void init_StoU(InfGraph &gred)
	{
#ifdef DISCRETE
		StoUsize = 2*log(n)/(0.1);
		Timer t(0, "Step 0");
	    int R = StoUsize;
		int sexp = 0;
	    StoU.resize(R);
		vector<int> num_visited(n);

	        std::random_device rd;
        	std::mt19937 gen(rd());
	        std::geometric_distribution<> gd(0.01);

	    for(int i=0; i<R; i++)
	    {
                vector<int> visit1(n);
                vector<int> visit_mark1(n);
                int n_visit_mark1 = 0;
                q.clear();
		selseeds.clear();

                int thisnode = 0;
                while(1)
                {
                        thisnode = thisnode + gd(gen);
                        if(thisnode >= n) break;
                        selseeds.push_back(thisnode);
                }

                for(int j=0; j< (int)selseeds.size(); j++)
                {
                        q.push_back(selseeds[j]);
                        visit_mark1[n_visit_mark1++] = selseeds[j];
                        visit1[selseeds[j]] = 1;
                }


		while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
                        for(int j=0; j< (int)sG[expand].size(); j++)
                        {
                                int v = sG[expand][j];
                                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                                if (randDouble > probs[expand][j])
                                    continue;
                                if (visit1[v]==1)
                                    continue;
                                if (visit1[v]==0)
                                {
                                    ASSERT(n_visit_mark1 < n);
                                    visit_mark1[n_visit_mark1++] = v;
                                    visit1[v] = 1;
                                  
//                                    StoU[i].insert(v);
                                  }
                                q.push_back(v);
                        }
                }
		sexp += n_visit_mark1;
                for (int i = 0; i < n_visit_mark1; i++)
		{
                        visit1[visit_mark1[i]] = 0;
			(num_visited[visit_mark1[i]])++;
//			redgraph[visit_mark1[i]] = 1;
//			gmap.insert(visit_mark1[i]);
		}
	    }

		for(int i=0; i<n; i++)
		{
		    if(num_visited[i]>7)
		    {
			redgraph[i] = 1;
			gmap.insert(i);
		    }
		}

		int cnt = 1;
		for(int i=0; i<n; i++)
		{
			if(redgraph[i]>0)
			{
				gmapvec.push_back(i);
				redgraph[i]=cnt;
				cnt++;
			}
		}
		cout<<"Expected spread from source:: "<<(double)sexp/R<<". Total nodes reachable:: "<<gmap.size()<<endl;
		int idx = 0;
		int xx = gmap.size(), x, y;

		gred.m = m; gred.init_update(xx); gred.StoU.resize(R); gred.StoUsize = R;
		gred.visit = vector<bool> (xx); gred.visit_mark = vector<int> (xx);
		gred.selseeds = selseeds; gred.init.resize(xx);
		gred.hyperG.resize(gred.n);// gred.hyperGT.resize(R);
		for(int i=0; i<(int)selseeds.size(); i++)
		{
			gred.selseeds[i] = redgraph[selseeds[i]]-1;
			gred.init[gred.selseeds[i]] = 1;
		}
		
		for(int i=0; i< (int)gmapvec.size(); i++)
		{
			x = gmapvec[i];
			for(int j=0; j<(int) sG[x].size(); j++)
			{
				y = redgraph[sG[x][j]]-1;
				if(redgraph[sG[x][j]]>0)
				{
					gred.add_edge(i, y, probs[x][j]); //cout<<x<<" "<<y<<" "<<probs[x][j]<<endl;
				}
			}
		}
		for(int i=0; i<R; i++)
		{
			for(set<int>::iterator it=StoU[i].begin(); it!=StoU[i].end(); it++)
				gred.StoU[i].insert(redgraph[*it]-1);
		}
#endif
	}
	void build_hyper_graph_r(int64 R, const Argument & arg)
	{
		if( R > INT_MAX ){
			cout<<"Error:R too large"<<endl;
			exit(1);
		}
		INFO("build_hyper_graph_r", R);



        int prevSize = hyperGT.size();
        while ((int)hyperGT.size() <= R)
            hyperGT.push_back( vector<int>() );



        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            random_number.push_back(  sfmt_genrand_uint32(&sfmtSeed) % n);
        }

        //trying BFS start from same node
        
	double avg = 0;
        for (int i = prevSize; i < R; i++)
        {
#ifdef CONTINUOUS
            BuildHypergraphNode(random_number[i], i, arg );
#endif
#ifdef DISCRETE
            newBuildHypergraphNode(random_number[i], i );
#endif
        }

        int totAddedElement = 0;
        for (int i = prevSize; i < R; i++)
        {
            for (int t : hyperGT[i])
            {
                hyperG[t].push_back(i);
                //hyperG.addElement(t, i);
                totAddedElement++;
            }
        }
    }

    void build_hyper_graph_ts(int64 R, const Argument & arg)
    {
        if( R > INT_MAX ){
            cout<<"Error:R too large"<<endl;
            exit(1);
        }
        INFO("build_hyper_graph_r", R);

	hyperGT.clear(); hyperG.clear();

        int prevSize = 0;
        while ((int)hyperGT.size() <= R)
            hyperGT.push_back( vector<int>() );

	while ((int)hyperG.size()<=n)
	{
		hyperG.push_back(vector<int>());
		StoU.push_back(set<int>());
	}

        vector<int> random_number;
        for (int i = 0; i < R; i++)
        {
            random_number.push_back(  sfmt_genrand_uint32(&sfmtSeed) % n);
        }

        //trying BFS start from same node


        for (int i = prevSize; i < R; i++)
        {
#ifdef DISCRETE
            newBuildHypergraphNode(random_number[i], i );
#endif
        }
    }

#ifdef DISCRETE
#include "discrete_rrset.h"
#endif
#ifdef CONTINUOUS
#include "continuous_rrset.h"
#endif

    //return the number of edges visited
    deque<int> q, q1;
    sfmt_t sfmtSeed;
    set<int> seedSet;
	vector<int> seedloc;

	int pick_best_seed()
	{
		vector<int> degree;
		int id = 0;
#ifdef DISCRETE
		for (int i = 0; i < n; i++)
	        {
         	    degree.push_back( hyperG[i].size());
        	}
		auto t = max_element(degree.begin(), degree.end());
		id = t-degree.begin();
		seedSet.insert(id); cout<<id<<" spreads "<<(double)degree[id]/hyperGT.size()<<" seedloc "<<seedloc[id];
		seedloc[id] = 1;
             // sG[id].clear(); gT[id].clear();
#endif
		return (id);
	}

    void build_seedset(int k, vector<double> stou)
    {
        Counter cnt(1);
        vector< int > degree;
        vector< bool> visit_local(hyperGT.size());

        seedSet.clear();
       
/*	for (int i=0; i<(int)selseeds.size(); i++)
	{
		if(hyperG[selseeds[i]].size() > 0)
		{cout<<"ERROR!!!!";
		exit(1);
		hyperG[selseeds[i]].clear();
		}
	}*/
	for (int i = 0; i < n; i++)
        {
            degree.push_back( hyperG[i].size());
		if(degree[i]<0)
			cout<<"ERROR!";
        }
        ASSERT(k > 0);
        ASSERT(k < (int)degree.size());

        for (int i = 0; i < k; i++)
        {
            auto t = max_element(degree.begin(), degree.end());
            int id = t - degree.begin(); cout<<id<<" "; // cout<<gmapvec[id]<<" ";
            seedSet.insert(id);
            degree[id] = 0;
            for (int t : hyperG[id])
            {
                if (!visit_local[t])
                {
                    visit_local[t] = true;
                    for (int node : hyperGT[t])
                    {
                        degree[node]--;
                    }
                }
            }
        }
        TRACE(seedSet);
    }
    double InfluenceHyperGraph()
    {

        set<int> s;
        TRACE(seedSet);
        for (auto t : seedSet)
        {
            for (auto tt : hyperG[t])
            {
                s.insert(tt);
            }
        }
        double inf = (double)(n-selseeds.size()) * s.size() / hyperGT.size();
        return inf;
    }

    void Buildstou(vector<double> &stou, int NUMSIM)
    {
#ifdef DISCRETE
        stou.resize(n);
	fill(stou.begin(), stou.end(), 1.0);
/*        vector<double> counts(n);
	double exp_spread = 0;
        for(int i=0; i<NUMSIM; i++)
        {
                int n_visit_mark = 0;
                q.clear();
                for(int j=0; j< (int)selseeds.size(); j++)
                {
                        q.push_back(selseeds[j]);
                        visit_mark[n_visit_mark++] = selseeds[j];
                        visit[selseeds[j]] = true;
                }

                while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
                        for(int j=0; j< (int)sG[expand].size(); j++)
                        {
                                int v = sG[expand][j];
                                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                                if (randDouble > probs[expand][j])
                                    continue;
                                if (visit[v])
                                    continue;
                                if (!visit[v])
                                {
					if(n_visit_mark >=n )
					{
						cout<<n_visit_mark<<" "<<n<<endl;
						exit(0);
					}
                                    ASSERT(n_visit_mark < n);
                                    visit_mark[n_visit_mark++] = v;
                                        counts[v]=counts[v]+1;
                                }
                                q.push_back(v);
                        }
                }
                for (int j = 0; j < n_visit_mark; j++)
                        visit[visit_mark[j]] = false;
		n_visit_mark=0;

        }
        for(int j=0; j<n; j++)
	{
		stou[j] = counts[j]/NUMSIM;
		exp_spread = exp_spread+stou[j];
	}

	cout<<endl<<exp_spread<<endl;
*/
#endif
        return;
}


};





