# include <math.h>
# include "helper.cpp"

std::vector<int> degseeds(Graph &G, int num_seeds, int vec = 0)
{
	std::priority_queue<std::pair<double, int> > q;
	std::vector<int> selseeds(num_seeds), dist(G.size());

	for (int i = 0; i < G.size(); ++i) 
	{
		double wdeg = 0;
		for(int j=0; j<G[i].outnode.size(); j++)
			wdeg += G[i].outweight[j];

    	q.push(std::pair<double, int>(wdeg, i));
  	}
  	for (int i=0; i<num_seeds; ++i)
  	{
  		selseeds[i] = q.top().second;
 		cout<<selseeds[i]<<'\t';
  		q.pop();
  	}
  	cout<<endl;
  	if(vec==0)
  		return selseeds;
  	else
  	{
  		for(int i=0; i<num_seeds; ++i)
  		{
  			(dist[selseeds[i]])++;
  		}
  		return dist;
  	}
}

double utilityf(int x)
{
	//return pow(0.5, x);
	return x>0?0.0:1.0;
}

double my_formula_ICM(char *fname, Graph &G, std::vector<int> init, int maxt, std::vector<double> &newB)
{
	ofstream outfile(fname);
	std::vector<double> prevB2(G.size()), prevB(G.size());
	newB.resize((G.size()));
	int i, num_nodes = G.size(), node_idx;
	double rt=1, num_infec, ex;
	fill(newB.begin(), newB.end(), 0);

	for(i=0; i<init.size(); ++i)
	{
		prevB[init[i]] = 1;
	}
	

	for(int t=1; t<maxt; ++t)
	{
		num_infec = 0.0;
		for(i=0; i<num_nodes; i++)
		{
			rt = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt *= (1 - G[i].inweight[j]*(prevB[node_idx]-prevB2[node_idx]));
			} 
			rt = 1.0 - rt;
			if(rt > 1.0)
			{
				cout<<"ERROR"<<endl;
				print_inneighbors(G, i);
				return 0;
			}
			newB[i] = prevB[i] + rt*(1-prevB[i]);
			num_infec += newB[i];
		}
		prevB2 = prevB;
		prevB = newB;
		if(t%1==0)
			cout<<num_infec<<endl;
	}

	return num_infec;
}

double my_formula_ICM(char *fname, Graph &G, std::vector<double> prevB, int maxt, std::vector<double> &newB)
{
	ofstream outfile(fname);
	std::vector<double> prevB2(G.size());
	newB.resize((G.size()));
	int i, num_nodes = G.size(), node_idx;
	double rt=1, num_infec, ex;

	if (prevB.size() < G.size())
	{
		cout<<"ERROR!"<<endl;
		exit(0);
	}

	for(int t=1; t<maxt; ++t)
	{
		num_infec = 0.0;
		for(i=0; i<num_nodes; i++)
		{
			rt = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt *= (1 - G[i].inweight[j]*(prevB[node_idx]-prevB2[node_idx]));
			} 
			rt = 1.0 - rt;
			if(rt > 1.0)
			{
				cout<<"ERROR"<<endl;
				print_inneighbors(G, i);
				return 0;
			}
			newB[i] = prevB[i] + rt*(1-prevB[i]);
			num_infec += newB[i];
		}
		prevB2 = prevB;
		prevB = newB;
	}

	return num_infec;
}

double my_formula_CCM(char *fname, Graph &G, std::vector<int> init, int maxt)
{
// Add Complex Contagion Model code here 
	;
}

double simulationICM(char *fname, Graph &G, int up_bound, int maxt, int num_sim, std::vector<double> &newB)
{
	int n, i, j, t;
	std::vector<int> prevI(G.size()), newI(G.size()), infec(G.size()), tinfec;
	newB.resize(G.size());
	std::vector<double> spread(maxt);
	int pcount, ncount, temp;
	double r;
	int init;
	for(n=0; n<num_sim; ++n)
	{
		pcount = 0; ncount = 0;
		fill(infec.begin(), infec.end(), 0);
//		tinfec.clear(); tinfec.reserve(1000);
		init = rand() % up_bound;
		
		prevI[0] = init;
		infec[init] = 1;
		pcount++;
		newB[init] += 1;
		
		spread[0] += pcount;
		for(t=1; t<maxt; ++t)
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
			pcount = ncount;
			ncount = 0;
		}

//		if(n % (num_sim/10) == 0)
//			cout<<"Simulation "<<n<<endl;
	}

	std::transform(newB.begin(), newB.end(), newB.begin(), std::bind1st(std::multiplies<double>(),1.0/num_sim));
	for(i=0; i<maxt; ++i)
	{
		spread[i] = spread[i]/num_sim;
		if(i>0)
			spread[i] += spread[i-1];

//		cout<<spread[i]<<endl;
	}
	double num_infec = 0;
	for(int i = 0; i<G.size(); i++)
		num_infec += newB[i];

	cout<<num_infec<<endl;
	return num_infec;
}

double simulationICM(char *fname, Graph &G, std::vector<int> init, int maxt, int num_sim, std::vector<double> &newB)
{
	int n, i, j, t;
	std::vector<int> prevI(G.size()), newI(G.size()), infec(G.size()), tinfec;
	fill(newB.begin(), newB.end(), 0);
	newB.resize(G.size());
	std::vector<double> spread(maxt);
	int pcount, ncount, temp;
	double r;

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

//		cout<<spread[i]<<endl;
	}
	double num_infec = 0;
	for(int i = 0; i<G.size(); i++)
		num_infec += newB[i];

//	return spread[maxt-1];
	return num_infec;
}

double simulationICM(char *fname, Graph &G, std::vector<double> dinit, int maxt, int num_sim, std::vector<double> &newB)
{
	int n, i, j, t;
	std::vector<int> prevI(G.size()), newI(G.size()), infec(G.size()), tinfec;
	newB.resize(G.size());
	std::vector<double> spread(maxt);
	int pcount, ncount, temp;
	double r;

	for(n=0; n<num_sim; ++n)
	{
		pcount = 0; ncount = 0;
//		fill(infec.begin(), infec.end(), 0);
		tinfec.clear(); tinfec.reserve(1000);
		std::vector<int> init;
		for(i=0; i<dinit.size(); i++)
		{
			r = (double)rand()/RAND_MAX;
			if(r <= dinit[i])
				init.push_back(i);
		}
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
					if(r <= G[prevI[i]].outweight[j] && infec[G[prevI[i]].outnode[j]]==0)
					{
						newI[ncount] = G[prevI[i]].outnode[j];
						tinfec.push_back(newI[ncount]);
						infec[newI[ncount]] = 1;
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

//		cout<<spread[i]<<endl;
	}
	double num_infec = 0;
	for(int i = 0; i<G.size(); i++)
		num_infec += newB[i];

//	return spread[maxt-1];
	return num_infec;
}

double simulationSI(char *fname, Graph &G, std::vector<int> init, int maxt, int num_sim, std::vector<int> immune)
{
	int n, i, j, t, icount;
	double r;
	std::vector<int> prevI(G.size()), infec(G.size());
	std::vector<double> spread(maxt);
	if(immune.size() < G.size())
		immune.resize(G.size());

	for(n=0; n<num_sim; n++)
	{
		icount = 0;
		fill(infec.begin(), infec.end(), 0);
		fill(prevI.begin(), prevI.end(), 0);

		for(i=0; i<init.size(); i++)
		{
			infec[init[i]] = 1;
			icount++;
			prevI[i] = init[i];
		}

		spread[0] += icount;

		for(t=1; t<maxt; t++)
		{
			for(i=0; i<icount; i++)
			{
				for(j=0; j<G[prevI[i]].outnode.size(); ++j)
				{
					if(immune[G[prevI[i]].outnode[j]]==1) continue;
					r = (double)rand()/RAND_MAX;
					if((r < G[prevI[i]].outweight[j]) && infec[G[prevI[i]].outnode[j]]==0)
					{
						prevI[icount] = G[prevI[i]].outnode[j];
						infec[prevI[icount]] = 1;
						icount++;
					}
				}
			}
			spread[t] += icount;
		}
		if(n % (num_sim/20) == 0)
			cout<<"Simulation "<<n<<endl;
	}
	
	for(i=0; i<maxt; ++i)
	{
		spread[i] = spread[i]/num_sim;
		

		if(i%1==0)
			cout<<spread[i]<<endl;
	}

	return spread[maxt-1];
}

double my_formula_SI(char *fname, Graph &G, std::vector<int> init, std::vector<int> dist, int maxt)
{
	std::vector<double> prevB(G.size()), newB(G.size());
	int i, num_nodes = G.size(), node_idx;
	double rt1, rt2, num_infec, ex;


	if (init.size()==num_nodes)
	{
		for(i=0; i<init.size(); ++i)
			prevB[i] = init[i];
	}
	else
	{
		for(i=0; i<init.size(); ++i)
			prevB[init[i]] = 1;
	}

	for(int t=1; t<maxt; ++t)
	{
		num_infec = (double)init.size();
		for(i=0; i<num_nodes; i++)
		{
			rt1 = 1.0;
			rt2 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]);
			} 
			
			if(prevB[i]==1) {rt1 = 0; rt2 = 0;}
			newB[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			num_infec += newB[i];
		}
		prevB = newB;

		if(t%10==0)
			cout<<"Finished round "<<t<<" with "<<num_infec<<" infections "<<endl;
	}
	return num_infec;
}


double my_formula_SI(char *fname, Graph &G, std::vector<double> init, std::vector<int> dist, int maxt) //Float
{
	std::vector<double> prevB(G.size()), newB(G.size());
	int i, num_nodes = G.size(), node_idx;
	double rt1, rt2, num_infec, ex;


	if (init.size()==num_nodes)
	{
		for(i=0; i<init.size(); ++i)
			prevB[i] = init[i];
	}
	else
		return -1.0;

	for(int t=1; t<maxt; ++t)
	{
		num_infec = (double)init.size();
		for(i=0; i<num_nodes; i++)
		{
			rt1 = 1.0;
			rt2 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]);
			}

			if(prevB[i]==1) {rt1 = 0; rt2 = 0;}

			newB[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			num_infec += newB[i];
		}
		prevB = newB;

		if(t%10==0)
			cout<<"Finished round "<<t<<" with "<<num_infec<<" infections "<<endl;
	}
	return num_infec;
}


std::vector<int> iOSSUM_SI(char *fname, Graph &G, std::vector<int> init, int num_seeds, int r=1)
{
	std::vector<double> prevB(G.size()), newB(G.size()), diffB(G.size()), prevB1(G.size()), newB1(G.size()), prevB2(G.size()), newB2(G.size());
	std::vector<int> dist(G.size()), seeds(num_seeds);
	int i, num_nodes = G.size(), node_idx, max_idx;
	double rt1, rt2, rt11, rt22, num_infec, ex, maxdiff=0, xx1, xx2, zz1, zz2;
	int maxt = num_seeds*r + 1;

	for(i=0; i<init.size(); ++i)
	{
		prevB[init[i]] = 1;
	}

	for(int t=1; t<maxt; ++t)
	{
		maxdiff = 0;
		max_idx = 0;
		num_infec = (double)init.size();
		
		for(i=0; i<num_nodes; i++)
		{
			diffB[i] = 0;

			rt1 = 1.0; rt11 = 1.0;
			rt2 = 1.0; rt22 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]); rt11 *= (1 - utilityf(dist[i]+1)*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]); rt22 *= ((1 - prevB[node_idx])*G[i].inweight[j]*utilityf(dist[i]+1));
			}

			if(prevB[i]==1) {rt1 = 0; rt2 = 0; rt11 = 0; rt22 = 0;}

			newB1[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			newB2[i] = 1 - (1-prevB[i])*rt11 - prevB[i]*rt22;

			diffB[i] = newB1[i] - newB2[i];

			rt1 = 1; rt2 = 1;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]);
			} 
			if(prevB[i]==1) {rt1 = 0; rt2 = 0;}
			
			newB[i] = 1 - (1-newB[i])*rt1 - prevB[i]*rt2;
			

			for(int j=0; j<G[i].outnode.size(); j++)
			{
				rt1 = 1.0; rt11 = 1.0;
				rt2 = 1.0; rt22 = 1.0;
				int i1 = G[i].outnode[j];
				for(int k=0; k<G[i1].innode.size(); k++)
				{
					node_idx = G[i1].innode[k];

					if(node_idx == i)
					{
						zz1 = newB1[i]; zz2 = newB2[i];
					}
					else
					{
						zz1 = prevB[node_idx]; zz2 = prevB[node_idx];
					}

					rt1 *= (1 - utilityf(dist[i1])*G[i1].inweight[k]*zz1); rt11 *= (1 - utilityf(dist[i1])*G[i1].inweight[k]*zz2);
					rt2 *= ((1 - zz1)*G[i1].inweight[k])*utilityf(dist[i1]); rt22 *= ((1 - zz2)*G[i1].inweight[k]*utilityf(dist[i1]));
				}

				if(prevB[i1]==1) {rt1 = 0; rt2 = 0; rt11 = 0; rt22 = 0;}

				xx1 = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
				xx2 = 1 - (1-prevB[i])*rt11 - prevB[i]*rt22;
				
				if(i==474 && xx1!=xx2)
				{
					cout<<endl;
				}


			//	diffB[i] += (xx1-xx2);

			}

			
			if(diffB[i]>maxdiff)
			{
				maxdiff = diffB[i];
				max_idx = i;
			}
		}

		(dist[max_idx])++;

		for(i=0; i<num_nodes; i++)
		{
			rt1 = 1.0;
			rt2 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]);
			} 
			
			if(prevB[i]==1) {rt1 = 0; rt2 = 0;}

			newB[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			num_infec += newB[i];
		}
		prevB = newB;
		cout<<max_idx<<'\t';
//		if(t%2==0)
//			cout<<"Finished round "<<t<<" with "<<num_infec<<" infections "<<endl;
	}
	cout<<endl<<"-------"<<endl;
	return dist;
}


std::vector<int> iapprox_SI(char *fname, Graph &G, std::vector<int> init, int num_seeds, int r=1)
{
	std::vector<double> prevB(G.size()), newB(G.size()), diffB(G.size()), prevB1(G.size()), newB1(G.size()), prevB2(G.size()), newB2(G.size());
	std::vector<int> dist(G.size()), seeds(num_seeds);
	int i, num_nodes = G.size(), node_idx, max_idx;
	double rt1, rt2, rt11, rt22, num_infec, ex, maxdiff=0, xx1, xx2, zz1, zz2;
	int maxt = num_seeds*r + 1;

	for(i=0; i<init.size(); ++i)
	{
		prevB[init[i]] = 1;
	}

	for(int t=1; t<maxt; ++t)
	{
		maxdiff = 0;
		max_idx = 0;
		num_infec = (double)init.size();
		double thisdiff = 0;
		for(i=0; i<num_nodes; i++)
		{
			diffB[i] = 0;

			rt1 = 1.0; rt11 = 1.0;
			rt2 = 1.0; rt22 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]); rt11 *= (1 - utilityf(dist[i]+1)*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]); rt22 *= ((1 - prevB[node_idx])*G[i].inweight[j]*utilityf(dist[i]+1));
			}

			if(prevB[i]==1) {rt1 = 0; rt2 = 0; rt11 = 0; rt22 = 0;}

			newB1[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			newB2[i] = 1 - (1-prevB[i])*rt11 - prevB[i]*rt22;

			diffB[i] = newB1[i] - newB2[i];
			thisdiff = diffB[i];
			for(int j=0; j<G[i].outnode.size(); j++)
			{
				node_idx = G[i].outnode[j];
				//thisdiff += (G[i].outweight[j]*(1-prevB[node_idx])+ prevB[node_idx])*(dist[i]==0);//*(diffB[i]);
				thisdiff += G[i].outweight[j]*(1-prevB[node_idx])*(diffB[i]);
				//thisdiff += utilityf(dist[i]);
			}

			if(thisdiff>maxdiff)
			{
				maxdiff = thisdiff;
				max_idx = i;
			}

		}

		(dist[max_idx])++;

		/*
		for(i=0; i<num_nodes; i++)
		{
			rt1 = 1.0;
			rt2 = 1.0;
			for(int j=0; j<G[i].innode.size(); j++)
			{
				node_idx = G[i].innode[j];
				rt1 *= (1 - utilityf(dist[i])*G[i].inweight[j]*prevB[node_idx]);
				rt2 *= ((1 - prevB[node_idx])*G[i].inweight[j])*utilityf(dist[i]);
			} 
			
			if(prevB[i]==1) {rt1 = 0; rt2 = 0;}

			newB[i] = 1 - (1-prevB[i])*rt1 - prevB[i]*rt2;
			num_infec += newB[i];
		}
		prevB = newB;
		*/
		cout<<"("<<max_idx<<", "<<num_infec<<")"<<'\t';
//		if(t%2==0)
//			cout<<"Finished round "<<t<<" with "<<num_infec<<" infections "<<endl;
	}
	cout<<endl<<"-------"<<endl;
	return dist;
}

/*std::vector<int> degseeds(Graph &G, int num_seeds)
{
	std::priority_queue<std::pair<int, int> > q;
	std::vector<int> selseeds(num_seeds);

	for (int i = 0; i < G.size(); ++i) 
	{
    	q.push(std::pair<double, int>(G[i].outnode.size(), i));
  	}
  	for (int i=0; i<num_seeds; ++i)
  	{
  		selseeds[i] = q.top().second;
 //		cout<<selseeds[i]<<" "<<q.top().first<<endl;
  		q.pop();
  	}

  	return selseeds;
}*/

/*std::vector<int> greedy_seeds(Graph &G, int num_seeds)
{
	int thisseed, num_nodes = G.size();
	std::vector<int> sel_seeds;
	double spread = 0, maxspread = 0;
	for(int k=0; k<num_seeds; ++k)
	{
		sel_seeds.push_back(0);
		for(int i=0; i<num_nodes; i++)
		{
			sel_seeds[k] = i;
			spread = my_formula_ICM("", G, sel_seeds, 20);
//			spread = my_formula_SI("", G, sel_seeds, 2);
			if(maxspread < spread)
			{
				maxspread = spread;
				thisseed = i;
			}
		}
		sel_seeds[k] = thisseed;
		cout<<"Added seed "<<thisseed<<", total spread "<<maxspread<<endl;
	}
	return sel_seeds;
}*/

/*int main(int argc, char *argv[])
{
	load_graph(G, argv[1], 0, 1, 0.01);
	std::vector<int> seeds, dist(G.size());
	std::cout<<"Loaded graph with "<<G.size()<<" nodes and "<<edges<<" edges"<<std::endl;
//  cout<<my_formula_ICM("theory_ICM_HEPT", G, degseeds(G, 50), 20);

//	cout<<simulationSI("", G, degseeds(G, 20), 200, 20);
//	seeds = greedy_seeds(G, 10);
	seeds = degseeds(G, 20, 0);
	cout<<my_formula_SI("", G, seeds, dist, 50)<<endl<<endl;
	dist = iapprox_SI("", G, seeds, 200, 1);
//	dist = degseeds(G, 200, 1);
//	fill(dist.begin(), dist.end(), 1);
	cout<<my_formula_SI("", G, seeds, dist, 50)<<endl<<endl;
	return 0;
}*/
