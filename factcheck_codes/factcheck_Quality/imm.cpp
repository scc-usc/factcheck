# include <algorithm>
# include <queue>
# include <deque>
# include "ICM_unsigned.cpp"
# define NUMSIM 1000
# define MAXT 50
# include <random>
# include <string.h>
int pick_rand(vector<double> pdist)
{
	double r;

	r= (double)rand()/RAND_MAX;
	int low = 0, high = pdist.size(), mid = (low+high)/2;

	if(r<pdist[0]) return 0;

	while(low < high)
	{
		if(r<=pdist[mid])
			high = mid;
		else if(r>pdist[mid])
			low = mid;
		if(mid == (low+high)/2)
			return mid+1;

		mid = (low+high)/2;
	}
}

vector<int> readseeds(char* fname, int Gsize, int bin, int k=1000)
{
	vector<int> seedvec(Gsize), seedlocs;
	int xx;
	ifstream fi(fname);
	if (fi.is_open())
  	{
    	while (!fi.eof() && k>0)
    	{
    		fi>>xx;
    		seedvec[xx] = 1;
		seedlocs.push_back(xx);
    		k--;
    	}

    }
    else
    {
    	cout<<"Failed to read seedfile"<<endl;
    }
    if(bin==1)
	return seedvec;
    else
	return seedlocs;
}


double simulationICM(char *fname, Graph &G, std::vector<int> init, std::vector<int> immune1, int maxt, int num_sim, std::vector<double> &newB, int MCSP)
{

        int n = G.size();
	int conv = 0;
	deque<int> q, q1;
	vector<int> immune(n);
	vector<set<int>> liveG(n);
	
	if(immune1.size()<G.size())
        {
                for(int i1 = 0; i1<immune1.size(); i1++)
                        immune[immune1[i1]] = 1;
        }
        else
                immune = immune1;

 //       vector<int> setloc;

//        for(int i=0; i<G.size(); i++)
//                if(init[i]==1) setloc.push_back(i);

	vector<int> visit(n), visit1(n);
        vector<int> visit_mark(n), visit_mark1(n);

	for(int simnum = 0; simnum<NUMSIM; simnum++)
	{
	int success = 0;
        vector<int> nodeclear;

                int n_visit_mark = 0;
                q.clear(); q1.clear();
                
                n_visit_mark = 0;

		for(int i=0; i<init.size(); i++)
		{
			q.push_back(init[i]);
                	visit_mark[n_visit_mark++] = init[i];
                	visit[init[i]] = true;
		}
                while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
                        if(immune[expand] == 1) {q1.push_back(expand); success=1;}
                        for(int j=0; j< (int)G[expand].outnode.size(); j++)
                        {
                                int v = G[expand].outnode[j];
                                double randDouble = (double)rand()/RAND_MAX;;
                                if (randDouble > G[expand].outweight[j])
                                    continue;


                                if (visit[v]==true)
                                    continue;

                                liveG[expand].insert(v);
				nodeclear.push_back(expand);

                                if (visit[v]==false)
                                {
                                    visit_mark[n_visit_mark++] = v;
                                    visit[v] = true;
                                }
                                q.push_back(v);
                        }
                }

	              for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;


                if(success==1)
                {
                n_visit_mark = 0;

                for(auto t: q1)
                {
                        visit_mark[n_visit_mark++] = t;
                        visit[t] = true;
			conv++;
                }

                int expand;
                while(!q1.empty())
                {
                        expand = q1.front();
                        q1.pop_front();
                        for(int j=0; j< (int) G[expand].outnode.size(); j++)
                        {
                                int v = G[expand].outnode[j];
  				 if(liveG[expand].find(v)== liveG[expand].end())
                                        continue;

                                if(visit[v]==true)
                                        continue;

                                visit_mark[n_visit_mark++] = v;
                                visit[v] = true;
                                conv++;

                                q1.push_back(v);

                        }
                }
                for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;
                }

                for(int i=0; i< (int)nodeclear.size(); i++)
                        liveG[nodeclear[i]].clear();
                nodeclear.clear();


		if(simnum%(NUMSIM/10)==0)
			cout<<"simulation number "<<simnum<<" value is"<<(double)conv/(simnum+1)<<endl;
		}
		return (double)conv/NUMSIM;

}

double simulationICM1(char *fname, Graph &G, std::vector<int> init, std::vector<int> immune1, int maxt, int num_sim, std::vector<double> &newB, int MCSP)
{
        int n = G.size();
        int conv = 0;
        deque<int> q, q1;
        vector<int> immune(n);
        vector<set<int>> liveG(n);
	vector<vector<int>> hyperG(n);

	double SPROB = 0.01;
	Graph G1 = flip_edges(G);

        if(immune1.size()<G.size())
        {
                for(int i1 = 0; i1<immune1.size(); i1++)
                        immune[immune1[i1]] = 1;
        }
        else
                immune = immune1;
 
         vector<int> binit(n);
	vector<int> immunelocs;

	for(int i=0; i<n; i++)
		if(immune[i]==1) immunelocs.push_back(i);

//        for(int i=0; i<init.size(); i++)
//                binit[init[i]] = 1;

        vector<int> visit(n), visit1(n);
        vector<int> visit_mark(n), visit_mark1(n);

        for(int simnum = 0; simnum<NUMSIM; simnum++)
        {
        int success = 0;
        vector<int> nodeclear;

                int n_visit_mark = 0;
                q.clear(); q1.clear();

                n_visit_mark = 0;
//		int uStart = rand()%(n-init.size());
//		if(binit[uStart]>0) {simnum--; continue;}
		int uStart = rand()%n;
                q.push_back(uStart);
                visit_mark[n_visit_mark++] = uStart;
                visit[uStart] = true;
                while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
			double rr = (double)rand()/RAND_MAX;
                        if(rr < SPROB) {q1.push_back(expand); success=1;}
			for(int j=0; j< (int)G1[expand].outnode.size(); j++)
                        {
                                int v = G1[expand].outnode[j];
                                double randDouble = (double)rand()/RAND_MAX;;
                                if (randDouble > G1[expand].outweight[j])
                                    continue;


                                if (visit[v]==true)
                                    continue;

                                liveG[expand].insert(v);
                                nodeclear.push_back(expand);

                                if (visit[v]==false)
                                {
                                    visit_mark[n_visit_mark++] = v;
                                    visit[v] = true;
                                }
                                q.push_back(v);
                        }
                }

                      for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;


                if(success==1)
                {
                n_visit_mark = 0;

                int expand;
		int init_reached = 0;
                while(!q1.empty())
                {
                        expand = q1.front();
			if(visit[expand]==false)
                        {
                                visit[expand]=true;
                                visit_mark[n_visit_mark++] = expand;
                                hyperG[expand].push_back(simnum);
                        }

                        q1.pop_front();
                        for(int j=0; j< (int) G[expand].outnode.size(); j++)
                        {
                                int v = G[expand].outnode[j];
                                 if(liveG[v].find(expand)== liveG[v].end())
                                        continue;

                                if(visit[v]==true)
                                        continue;

                                visit_mark[n_visit_mark++] = v;
                                visit[v] = true;
				hyperG[v].push_back(simnum);
                                q1.push_back(v);

                        }
                }
                for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;
                }

                for(int i=0; i< (int)nodeclear.size(); i++)
                        liveG[nodeclear[i]].clear();
                nodeclear.clear();


//                if(simnum%(NUMSIM/10)==0)
//                        cout<<"simulation number "<<simnum<<" value is"<<(double)n*conv/(simnum+1)<<endl;
                }
		set<int> infset;
		for(int i=0; i<immunelocs.size(); i++)
		{
			for(auto tt:hyperG[immunelocs[i]])
				infset.insert(tt);
		}
                return (double)n*infset.size()/NUMSIM;

}


double simulationICM(char *fname, Graph &G, std::vector<int> init, std::vector<int> immune, int maxt, int num_sim, std::vector<double> &newB)
{
	int n, i, j, t;
	std::vector<int> prevI(G.size()), newI(G.size()), infec(G.size()), tinfec;
	newB.resize(G.size());
	std::vector<double> spread(maxt);
	int pcount, ncount, temp;
	double r;

	srand(time(NULL));
	for(n=0; n<num_sim; ++n)
	{
		pcount = 0; ncount = 0;
//		fill(infec.begin(), infec.end(), 0);
		tinfec.clear(); tinfec.reserve(1000);
		for(i=0; i<init.size(); i++)
		{
			prevI[i] = init[i];
			infec[init[i]] = 1;
			pcount++;
			newB[init[i]] += 1;
		}
		spread[0] += pcount;
		for(t=1; t<maxt; ++t)
		{
			for(i=0; i<pcount; ++i)
			{
				for(j=0; j<G[prevI[i]].outnode.size(); ++j)
				{
					r = (double)rand()/RAND_MAX;
					if(immune[G[prevI[i]].outnode[j]]==1) continue;

					if(r < G[prevI[i]].outweight[j] && infec[G[prevI[i]].outnode[j]]==0)
					{
						newI[ncount] = G[prevI[i]].outnode[j];
						infec[newI[ncount]] = 1;
						tinfec.push_back(newI[ncount]);
						newB[newI[ncount]] += 1.0;
						ncount++;
					}
				}
			}
			spread[t] += ncount;
//			prevI = newI;
			pcount = ncount;
			for(j=0; j<pcount; j++)
				prevI[j] = newI[j];
			ncount = 0;

			if(pcount<1)
			{
			//	cout<<"EARLY EXIT";
				break;
			}
		}

		for(i=0; i<tinfec.size(); i++)
			infec[tinfec[i]] = 0;

		if(n % (num_sim/10) == 0)
			cout<<"Simulation "<<n<<endl;
	}

	std::transform(newB.begin(), newB.end(), newB.begin(), std::bind1st(std::multiplies<double>(),1.0/num_sim));
	for(i=0; i<maxt; ++i)
	{
		spread[i] = spread[i]/num_sim;
		if(i>0)
			spread[i] += spread[i-1];

	//	cout<<spread[i]<<endl;
	}
	double num_infec = 0;
	for(int i = 0; i<G.size(); i++)
		num_infec += newB[i];

//	return spread[maxt-1];
	return num_infec;
}

int pick_best_imm(Graph &G, vector<double> newB1, std::vector<int> init, vector<int> immnodes, int low, int high)
{
	int i, j, t;
	double r;
	std::vector<double> spread(MAXT), cdist(newB1.size()), newB(G.size());
	std::vector<int> prevI(G.size()), newI(G.size()), infec(G.size());
	int pcount, ncount, temp;

/*	cdist[0] = newB1[0];
	for(i=1; i<newB1.size(); i++)
	{
		cdist[i] = cdist[i-1]+newB1[i];
		if (cdist[i] > G.size())
			cout<<"ERROR";
	}
//	std::transform(cdist.begin(), cdist.end(), cdist.begin(), std::bind1st(std::multiplies<double>(),1.0/cdist[cdist.size()-1]));
	for(i=0; i<newB1.size(); i++)
	{
		cdist[i] = cdist[i]/cdist[cdist.size()-1];
	}*/


	for(int n=0; n<NUMSIM; n++)
	{
		//int src = pick_rand(cdist);
		int src = rand() % G.size();

		pcount = 0; ncount = 0;
		fill(infec.begin(), infec.end(), 0);

		prevI[0] = src;
		infec[src] = 1;
		pcount++;
		//newB[src] += newB1[src];
		spread[0] += pcount;
		for(t=1; t<MAXT; ++t)
		{
			for(i=0; i<pcount; ++i)
			{
				for(j=0; j<G[prevI[i]].outnode.size(); ++j)
				{
					r = (double)rand()/RAND_MAX;

					if(r < G[prevI[i]].outweight[j] && infec[G[prevI[i]].outnode[j]]==0)
					{
						newI[ncount] = G[prevI[i]].outnode[j];
						infec[newI[ncount]] = 1;
						if(newI[ncount]>=low)
							newB[newI[ncount]] += newB1[newI[ncount]-low];
						ncount++;
					}
				}
			}
			spread[t] += ncount;
			prevI = newI;
			pcount = ncount;
			ncount = 0;
			if(pcount<1)
			{
//				cout<<"EARLY EXIT";
				break;
			}
		}

		if(n % (NUMSIM/10) == 0)
			cout<<"Simulation "<<n<<endl;
	}

	for(int i=0; i<init.size(); i++)
	{
		newB[init[i]-low] = 0;
	}

	double maxel = newB[low]; int idx = low;
	for(int i=low; i<=high; i++)
	{
		if(immnodes[i-low]==1) newB[i]=0;
		if(newB[i]>maxel){maxel = newB[i]; idx = i;}
	}
	return idx;

}



void all_imm(Graph &G, std::vector<int> init, std::vector<int> &immune, int res)//using top k
{
	double xx;
	int i, j;
        vector<double> newB, newB1, newB2(G.size());
	Graph G1 = flip_edges(G);
	xx = simulationICM("", G, init, 20, 5000, newB);
	xx = simulationICM("", G1, G1.size(),20, NUMSIM, newB1);

	for(i=0; i<G.size(); i++)
		newB2[i] = newB[i]*newB1[i];

	for(int k=0; k< res; k++)
	{
		double maxel = newB2[0]; int idx = 0;
		for(int i=0; i<newB2.size(); i++)
		{
			if(newB2[i]>maxel){maxel = newB2[i]; idx = i;}
		}
		immune[idx] = 1;
		newB2[idx] = 0;
	}
}
void AR_imm(Graph &G, std::vector<int> init, std::vector<int> &immune, int res)
{
	int i, I, thisnode;
	double xx;
	vector<double> newB, newB1, cdist, probs(G.size());
	vector<int> immnodes;
	fill(immune.begin(), immune.end(), 0);
	Graph G1 = G;
	for(i=0; i<res; i++)
	{
		
/*		I = G.size();

		G1.resize(1+2*G.size());
		for(int j=0; j<immnodes.size(); j++)
			new_edge(G1, immnodes[j], I, 1);

		for(int j=0; j<G.size(); j++)
		{
			new_edge(G1, j, I+j+1, 1);
			new_edge(G1, I, I+j+1, 1);
		}

		xx = simulationICM("", G1, init, 20, 1000, newB1);
		for(int j=0; j< G.size(); j++)
		{
			probs[j] = newB1[j+I+1];
			if(immune[j]==1)
				probs[j] = 0;
		}
		for(int j=0; j<init.size(); j++)
			probs[init[j]] = 0;

		Graph G2 = flip_edges(G);
		G2.resize(1+2*G.size());

		for(int j=0; j<immnodes.size(); j++)
			new_edge(G2, immnodes[j], I, 1);

		for(int j=0; j<G.size(); j++)
		{
			if(immune[j]==0) new_edge(G2, j, I+j+1, 1);
			new_edge(G2, I, I+j+1, 1);
		}
		thisnode = pick_best_imm(G2, probs, init, immune, I+1, G2.size()-1);
		thisnode = thisnode - I - 1; */

		xx = simulationICM("", G1, init, 20, 5000, newB1);
		Graph G2 = flip_edges(G1);
		thisnode = pick_best_imm(G2, newB1, init, immune, 0, G2.size()-1);

		immune[thisnode] = 1;
		immnodes.push_back(thisnode);
		remove_node(G1, thisnode);

		cout<<"Iteration "<<i<<": Picked node # "<<thisnode<<endl;
	}
	cout<<"--------------------------"<<endl;
	for(i=0; i<res; i++)
		cout<<immnodes[i]<<endl;
	cout<<"--------------------------"<<endl;
}

void varyseeds(char *fname, Graph G, vector<int> init, vector<int> immune1)
{
        int n = G.size();
        int conv = 0;
        deque<int> q, q1;
        vector<int> immune(n);
        vector<set<int>> liveG(n);
	double SPROB = 0.01;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::geometric_distribution<> gd(0.01);

        ofstream fid;
        fid.open(fname);

        if(!fid.is_open())
        {
                cout<<"Couldn't open file "<<fname<<endl;
                return;
        }


        if(immune1.size()<G.size())
        {
                for(int i1 = 0; i1<immune1.size(); i1++)
                        immune[immune1[i1]] = 1;
        }
        else
                immune = immune1;

//        for(int i=0; i<init.size(); i++)
 //               if(immune[init[i]]==1) immune[init[i]] = 0;

        vector<int> visit(n), visit1(n);
        vector<int> visit_mark(n), visit_mark1(n);
        Graph G1 = flip_edges(G);

        vector<vector<pair<int, int>>> hyperGT(n);
        for(int simnum = 0; simnum<NUMSIM; simnum++)
        {
                int success = 0; int isrepeat = 0;
                vector<int> nodeclear, activeimm;

                int n_visit_mark = 0;
                q.clear(); q1.clear();
		init.clear();

		int thisnode = 0;
		while(1)
		{
			thisnode = thisnode + gd(gen);

			if(thisnode >= n) break;

			init.push_back(thisnode);
		}
		if(simnum==1)
			cout<<"Inserted "<<init.size()<<" nodes"<<endl; 

                n_visit_mark = 0;

                for(int i=0; i<init.size(); i++)
                {
                        isrepeat = 0;
                        int newrand = init[i];
                        q.push_back(newrand);
                        visit_mark[n_visit_mark++] = newrand;
                        visit[newrand] = true;
                }

                while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
                        if(immune[expand] == 1) {activeimm.push_back(expand); success=1;}
                        for(int j=0; j< (int)G[expand].outnode.size(); j++)
                        {
                                int v = G[expand].outnode[j];
                                double randDouble = (double)rand()/RAND_MAX;;
                                if (randDouble > G[expand].outweight[j])
                                    continue;


                                if (visit[v]==true)
                                    continue;

                                liveG[expand].insert(v);
                                nodeclear.push_back(v);

                                if (visit[v]==false)
                                {
                                    visit_mark[n_visit_mark++] = v;
                                    visit[v] = true;
                                }
                                q.push_back(v);
                        }
                }

                      for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;

                for(int iiter=0; iiter < activeimm.size(); iiter++)
//              if(success == 1)
                {
                        n_visit_mark = 0;
                        int t = activeimm[iiter];
//                      int t=rand()%n;
                        visit_mark[n_visit_mark++] = t;
                        visit[t] = true;

                        q1.push_back(t);

                        hyperGT[t].push_back(make_pair(simnum, t));
                        while(!q1.empty())
                        {
                         int expand = q1.front(); //{if(immune[expand]==1) hyperGT[expand].push_back(simnum);}
                         q1.pop_front();
                         for(int j=0; j< (int)G[expand].outnode.size(); j++)
                         {
                                int v = G[expand].outnode[j];

                                if(liveG[expand].find(v)== liveG[expand].end())
                                    continue;


                                if (visit[v]==true)
                                    continue;


                                visit_mark[n_visit_mark++] = v;
                                visit[v] = true;
                                //if(immune[v]==1)
                                        hyperGT[t].push_back(make_pair(simnum, v));
                                q1.push_back(v);
                         }
                        }
                        for (int i = 0; i < n_visit_mark; i++)
                                visit[visit_mark[i]] = false;

                }
                 for(int  i = 0; i < nodeclear.size(); i++)
                        liveG[nodeclear[i]].clear();

                nodeclear.clear();

                if(simnum%(NUMSIM/10)==0)
                        cout<<"Simulation no. "<<simnum<<endl;

        }
        set<pair<int, int>> inflist; float spread;
        for(int i=0; i<immune1.size(); i++)
        {
                for(auto tt:hyperGT[immune1[i]])
                        inflist.insert(tt);

                spread = (double)inflist.size()/NUMSIM;
                cout<<spread<<endl;
                fid<<spread<<endl;
        }
        fid.close();
}

double ICM_innoc(char *fname, Graph G, vector<int> init, vector<int> immune1)
{
        int n = G.size();
        int conv = 0;
        deque<int> q, q1;
        vector<int> immune(n);
   
        ofstream fid;
        fid.open(fname);

        if(!fid.is_open())
        {
                cout<<"Couldn't open file "<<fname<<endl;
                return 0;
        }

        for(int i1 = 0; i1<immune1.size(); i1++)
	        immune[immune1[i1]] = 1;
        
        vector<int> visit(n), visit1(n);
        vector<int> visit_mark(n), visit_mark1(n);
        
        for(int simnum = 0; simnum<NUMSIM; simnum++)
        {
                int success = 0; int isrepeat = 0;
                vector<int> nodeclear, activeimm;

                int n_visit_mark = 0;
                q.clear(); q1.clear();

                n_visit_mark = 0;

                for(int i=0; i<init.size(); i++)
                {
                        int newrand = init[i];
                        q.push_back(newrand);
                        visit_mark[n_visit_mark++] = newrand;
                        visit[newrand] = true;
                }

                while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
                        if(immune[expand] == 1) {continue;}
                        for(int j=0; j< (int)G[expand].outnode.size(); j++)
                        {
                                int v = G[expand].outnode[j];
                                double randDouble = (double)rand()/RAND_MAX;;
                                if (randDouble > G[expand].outweight[j])
                                    continue;


                                if (visit[v]==true)
                                    continue;

                                if (visit[v]==false)
                                {
                                    visit_mark[n_visit_mark++] = v;
                                    visit[v] = true;
                                }
				if(immune[v]==0) conv++;

                                q.push_back(v);
                        }
                }

                      for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;

                if(simnum%(NUMSIM/10)==0)
                        cout<<"simulation number "<<simnum<<" value is"<<(double)conv/(simnum+1)<<endl;
                
		
	}
	double retval =  (double)conv/NUMSIM;
	fid<<retval;
	fid.close();

	return retval;

}

int main(int argc, char *argv[])
{
	srand(time(NULL));

	int numseeds=50, K = 50;
	load_graph(G, argv[1], 1, 0, 0.01);
	std::vector<double> newB, prevB;
	std::vector<int> seeds, dist(G.size());
	std::cout<<"Loaded graph with "<<G.size()<<" nodes and "<<edges<<" edges"<<std::endl;
	std::vector<int> selseeds, dimmune(G.size()), ARimmune;

	selseeds = degseeds(G, numseeds+K);
	std::vector<int> infseeds(numseeds);
	for(int i=0; i<K; i++)
	{
		dimmune[selseeds[i]] = 1;
	}

	if(strcmp(argv[2], "_")==0)
		ARimmune.assign(selseeds.begin(), selseeds.begin()+K);
	else
		ARimmune = readseeds(argv[2], G.size(), 0, K);

//	cout<<"With Degree Immunization"<<endl;
//	cout<<simulationICM("", G, infseeds, dimmune, 20, 10000, newB, 0)<<endl;
//	cout<<simulationSI("", G, infseeds, 50, 5000, dimmune)<<endl;
	cout<<"With AR Immunization"<<endl;
//	cout<<ICM_innoc(argv[3], G, infseeds, ARimmune)<<endl;
//	cout<<simulationICM1("", G, infseeds, ARimmune, 50, 10000, newB, 0.01);

	varyseeds(argv[3], G, infseeds, ARimmune);
	return 0;
}
