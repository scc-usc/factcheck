/*
 * BFS stf (influModel == IC)
        {
            int i = expand;
            for (int j = 0; j < (int)gT[i].size(); j++)
            {
                //int u=expand;
arting from one node
 */
# define SPROB 1.1

int BuildHypergraphNode(int uStart, int hyperiiid)
{

//	vector<int> visit_mark(n);
//	vector<bool> visit(n);
    int n_visit_edge = 1;
        ASSERT((int)hyperGT.size() > hyperiiid);
        hyperGT[hyperiiid].push_back(uStart);

    int n_visit_mark = 0;
	
    int rnum =sfmt_genrand_uint32(&sfmtSeed) % StoUsize; 
	q.clear();
    q.push_back(uStart);
    ASSERT(n_visit_mark < n);
    visit_mark[n_visit_mark++] = uStart;
    visit[uStart] = true;
    while (!q.empty())
    {

        int expand = q.front();
        q.pop_front();
        if (influModel == IC)
        {
            int i = expand;
            for (int j = 0; j < (int)gT[i].size(); j++)
            {
                //int u=expand;
                int v = gT[i][j];
                n_visit_edge++;
                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                if (randDouble > probT[i][j])
                    continue;
                if (visit[v])
                    continue;
                if (!visit[v])
                {
                    ASSERT(n_visit_mark < n);
                    visit_mark[n_visit_mark++] = v;
                    visit[v] = true;
                }
                q.push_back(v);
		if(StoU[rnum].find(v) != StoU[rnum].end())
		{ 
                    ASSERT((int)hyperGT.size() > hyperiiid);
                    hyperGT[hyperiiid].push_back(v);
		}
            }

        }
        else if (influModel == LT)
        {
            if (gT[expand].size() == 0)
                continue;
            ASSERT(gT[expand].size() > 0);
            n_visit_edge += gT[expand].size();
            double randDouble = sfmt_genrand_real1(&sfmtSeed);
            for (int i = 0; i < (int)gT[expand].size(); i++)
            {
                ASSERT( i < (int)probT[expand].size());
                randDouble -= probT[expand][i];
                if (randDouble > 0)
                    continue;
                //int u=expand;
                int v = gT[expand][i];

                if (visit[v])
                    break;
                if (!visit[v])
                {
                    visit_mark[n_visit_mark++] = v;
                    visit[v] = true;
                }
                q.push_back(v);
//                    ASSERT((int)hyperGT.size() > hyperiiid);
//                   hyperGT[hyperiiid].push_back(v);
                break;
            }
        }
        else
            ASSERT(false);

    }
    

                for (int j = 0; j < n_visit_mark; j++)
                        visit[visit_mark[j]] = false;

		
 //   		for (int i = 0; i < n_visit_mark1; i++)
//		        visit1[visit_mark1[i]] = 0;

	    return n_visit_edge;
}

int newBuildHypergraphNode(int uStart, int hyperiiid)
{
//    int n_visit_edge = 1;
	 int success = 0;
        ASSERT((int)hyperGT.size() > hyperiiid);
//        hyperGT[hyperiiid].push_back(uStart);
//        vector<int> visit_mark(n);
//        vector<bool> visit(n);
	vector<int> nodeclear;

//                vector<int> visit1(n);
//                vector<int> visit_mark1(n);
                int n_visit_mark = 0;
                q.clear(); q1.clear();
        	q.push_back(uStart);       
		n_visit_mark = 0;
                visit_mark[n_visit_mark++] = uStart;
                visit[uStart] = true;

		while(!q.empty())
                {
                        int expand = q.front();
                        q.pop_front();
			double rr =  sfmt_genrand_real1(&sfmtSeed);
			if(rr < SPROB) {q1.push_back(expand); success=1;}
                        for(int j=0; j< (int)gT[expand].size(); j++)
                        {
                                int v = gT[expand][j];
                                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                                if (randDouble > probT[expand][j])
                                    continue;

				liveG[v].insert(expand);
				
				nodeclear.push_back(v);

                                if (visit[v]==true)
                                    continue;
				

                                if (visit[v]==false)
                                {
                                    ASSERT(n_visit_mark < n);
                                    visit_mark[n_visit_mark++] = v;
                                    visit[v] = true;
                                }
                                q.push_back(v);
                        }
                }
		


              for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;

//		visit_mark.clear();

		if(success==1)
		{
		n_visit_mark = 0;

/*		for(auto t: q1)
		{
			visit_mark[n_visit_mark++] = t;
			visit[t] = true;
		}
*/		
		int expand;
		while(!q1.empty())
		{
			expand = q1.front();
			if(visit[expand]==false)
			{
				visit[expand]=true;
				visit_mark[n_visit_mark++] = expand;
				hyperGT[hyperiiid].push_back(expand);
			}

			q1.pop_front();
		
			for(int j=0; j< (int) sG[expand].size(); j++)
			{
				int v = sG[expand][j];
				if(liveG[expand].find(v)== liveG[expand].end())
					continue;

				if(visit[v]==true)
					continue;
				

				ASSERT(n_visit_mark <n );
				visit_mark[n_visit_mark++] = v;
                                visit[v] = true;
//				if(init[v]==1) continue;

				ASSERT((int)hyperGT.size() > hyperiiid);
		                hyperGT[hyperiiid].push_back(v);
		
				q1.push_back(v);

			}
		}
		for (int i = 0; i < n_visit_mark; i++)
                      visit[visit_mark[i]] = false;
		}

		for(int i=0; i< (int)nodeclear.size(); i++)
			liveG[nodeclear[i]].clear();
		nodeclear.clear();

		return hyperGT[hyperiiid].size();
}
