#define HEAD_INFO
//#define HEAD_TRACE
#include "sfmt/SFMT.h"
#include "head.h"
using namespace std;
typedef double (*pf)(int, int);
vector<int> redgraph;
set<int> gmap;
vector<int> gmapvec;
class Graph
{
public:
    int n, m, k;
    vector<int> inDeg, outDeg;
    vector<vector<int>> gT;
    vector<vector<int>> sG;
    vector<set<int>> liveG;
	vector<int> selseeds;
	vector<int> init;
#ifdef CONTINUOUS
    vector<vector<dpair>> probT;
#endif
#ifdef DISCRETE
    vector<vector<double>> probT;
    vector<vector<double>> probs;
#endif


    enum InfluModel {IC, LT, CONT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
        TRACE(influModel == IC);
        TRACE(influModel == LT);
        TRACE(influModel == CONT);
    }

    string folder;
    string graph_file;
    void readNM()
    {
        ifstream cin((folder + "attribute.txt").c_str());
        ASSERT(!cin == false);
        string s;
        while (cin >> s)
        {
            if (s.substr(0, 2) == "n=")
            {
                n = atoi(s.substr(2).c_str());
                continue;
            }
            if (s.substr(0, 2) == "m=")
            {
                m = atoi(s.substr(2).c_str());
                continue;
            }
            ASSERT(false);
        }
        TRACE(n, m );
        cin.close();
    }
#ifdef DISCRETE
    void add_edge(int a, int b, double p)
    {
        probT[b].push_back(p);
	probs[a].push_back(p);
        gT[b].push_back(a);
	sG[a].push_back(b);

//        inDeg[b]++;
//	outDeg[a]++;
    }
#endif
#ifdef CONTINUOUS
    void add_edge(int a, int b, double p1, double p2)
    {
        probT[b].push_back(MP(p1, p2));
        gT[b].push_back(a);
        inDeg[b]++;
    }
#endif
    vector<bool> hasnode;
    void readGraph()
    {
        FILE *fin = fopen((graph_file).c_str(), "r");
        ASSERT(fin != false);
        int readCnt = 0;
        for (int i = 0; i < m; i++)
        {
            readCnt ++;
            int a, b;
#ifdef DISCRETE
            double p;
            int c = fscanf(fin, "%d%d%lf", &a, &b, &p);
            ASSERTT(c == 3, a, b, p, c);
#endif
#ifdef CONTINUOUS
            double p1, p2;
            int c = fscanf(fin, "%d%d%lf%lf", &a, &b, &p1, &p2);
            ASSERT(c == 4);
#endif

            //TRACE_LINE(a, b);
            ASSERT( a < n );
            ASSERT( b < n );
            hasnode[a] = true;
            hasnode[b] = true;
#ifdef DISCRETE
            add_edge(a, b, p);
#endif
#ifdef CONTINUOUS
            add_edge(a, b, p1, p2);
#endif
        }
        TRACE_LINE_END();
        int s = 0;
        for (int i = 0; i < n; i++)
            if (hasnode[i])
                s++;
        INFO(s);
        ASSERT(readCnt == m);

#ifdef DISCRETE
	selseeds.resize(50);
	std::priority_queue<std::pair<double, int> > q11;

	for (int i = 0; i < n; ++i) 
	{
		double wdeg = 0;
		for(int j=0; j<(int)probs[i].size(); j++)
			wdeg += (double)(probs[i][j]);

    		q11.push(std::pair<double, int>(wdeg, i));
  	}
	init.resize(n);
	redgraph.resize(n);
  	for (int i=0; i<(int)selseeds.size(); ++i)
  	{
  		selseeds[i] = q11.top().second;
 		cout<<selseeds[i]<<'\t';
		init[selseeds[i]]=1;
  		q11.pop();
  	}
  	cout<<endl;
#endif
    }
#ifdef DISCRETE
    void readGraphBin()
    {
        string graph_file_bin = graph_file.substr(0, graph_file.size() - 3) + "bin";
        INFO(graph_file_bin);
        FILE *fin = fopen(graph_file_bin.c_str(), "rb");
        //fread(fin);
        struct stat filestatus;
        stat( graph_file_bin.c_str(), &filestatus );
        int64 sz = filestatus.st_size;
        char *buf = new char[sz];
        int64 sz2 = fread(buf, 1, sz, fin);
        INFO("fread finish", sz, sz2);
        ASSERT(sz == sz2);
        for (int64 i = 0; i < sz / 12; i++)
        {
            int a = ((int *)buf)[i * 3 + 0];
            int b = ((int *)buf)[i * 3 + 1];
            float p = ((float *)buf)[i * 3 + 2];
            //INFO(a,b,p);
            add_edge(a, b, p);
        }
        delete []buf;
        fclose(fin);
    }
#endif

	Graph()
	{
		;
	}
	void init_update(int xx)
	{
		n  = xx;
		            gT.resize(n);
                sG.resize(n);
        	liveG.resize(n);   
#ifdef DISCRETE
            probT.resize(n);
                probs.resize(n);
#endif

            //hyperGT.push_back(vector<int>());
            inDeg.resize(n);
            outDeg.resize(n);
        
	}
    Graph(string folder, string graph_file): folder(folder), graph_file(graph_file)
    {
        readNM();

        //init vector
        FOR(i, n)
        {
            gT.push_back(vector<int>());
		sG.push_back(vector<int>());
		liveG.push_back(set<int>());
            hasnode.push_back(false);
#ifdef DISCRETE
            probT.push_back(vector<double>());
		probs.push_back(vector<double>());
#endif
#ifdef CONTINUOUS
            probT.push_back(vector<dpair>());
#endif
            //hyperGT.push_back(vector<int>());
            inDeg.push_back(0);
		outDeg.push_back(0);
        }
        readGraph();
        //system("sleep 10000")i;
    }

};
double sqr(double t)
{
    return t * t;
}



