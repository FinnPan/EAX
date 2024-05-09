#include "util.h"
#include "eax.h"
#include "clk.h"

using namespace thu;

static const char *g_progname;
static Evaluator* g_eval;
static int g_maxIter;
static int g_bestCost;
static double g_avgCost;

enum HeuristicKey {
    HK_2opt,
    HK_eax,
    HK_clk,
    HK_NUM
};
static const char* HeuristicKey2Str (int k);
static int HeuristicStr2Key (const char* str);

static void Usage ();
static bool Proc2OPT (int argc, char **argv);
static bool ProcEAX (int argc, char **argv);
static bool ProcCLK (int argc, char **argv);
static bool Run2OPT (int times);
static bool RunEAX (bool verbose);
static bool RunCLK (int times, bool verbose);

int main(int argc, char* argv[])
{
    g_progname = argv[0];
    if (argc < 3) {
        Usage();
        return 1;
    }

    int hk = HeuristicStr2Key(argv[1]);
    if (hk < 0) {
        Usage();
        return 1;
    }

    Evaluator eval;
    if (!eval.Init(argv[2])) {
        return 1;
    }
    g_eval = &eval;

    RUsage ru;
    g_bestCost = INT_MAX;
    bool res = false;
    switch (hk) {
        case HK_2opt:
            res = Proc2OPT(argc, argv);
            break;
        case HK_eax:
            res = ProcEAX(argc, argv);
            break;
        case HK_clk:
            res = ProcCLK(argc, argv);
            break;
        default:
            Usage();
    };

    if (res) {
        ru.Report(HeuristicKey2Str(hk));
        double gap = g_eval->ComputeGap(g_bestCost);
        printf("[%s results] iter = %d; best = %d; avg = %0.3f; gap = %.3f%%\n",
                HeuristicKey2Str(hk), g_maxIter, g_bestCost, g_avgCost, 100*gap);
    }
    return 0;
}

const char* HeuristicKey2Str (int k)
{
    switch (k) {
        case HK_2opt:
            return "2opt";
        case HK_eax:
            return "eax";
        case HK_clk:
            return "clk";
        default:
            return "unknown";
    };
}

int HeuristicStr2Key (const char* str)
{
    int hk = -1;
    for (int i = 0; i < HK_NUM; ++i) {
        if (!strcmp(str, HeuristicKey2Str(i))) {
            hk = i;
            break;
        }
    }
    return hk;
}

void Usage ()
{
    fprintf(stderr, "%s: the command line utility of the TSP heuristics\n", g_progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "usage:\n");
    fprintf(stderr, "  %s 2opt  instance_file [-times t]\n", g_progname);
    fprintf(stderr, "  %s eax   instance_file [-verbose]\n", g_progname);
    fprintf(stderr, "  %s clk   instance_file [-times t] [-verbose]\n", g_progname);
}

bool Proc2OPT (int argc, char **argv)
{
    int times = 20;
    for (int i = 3; i < argc; i++){
        if(!strcmp(argv[i], "-times")){
            if (++i >= argc) {
                Usage();
                return false;
            }
            times = atoi(argv[i]);
            if (times <= 0) {
                fprintf(stderr, "ERROR: invalid -times\n");
                return false;
            }
        } else {
            Usage();
            return false;
        }
    }

    return Run2OPT(times);
}

bool ProcEAX (int argc, char **argv)
{
    bool verbose = false;
    for (int i = 3; i < argc; i++){
        if (!strcmp(argv[i], "-verbose")){
            verbose = true;
        } else {
            Usage();
            return false;
        }
    }

    return RunEAX(verbose);
}

bool ProcCLK (int argc, char **argv)
{
    int times = 1;
    bool verbose = false;
        for (int i = 3; i < argc; i++){
        if (!strcmp(argv[i], "-times")){
            if(++i >= argc) {
                Usage();
                return false;
            }
            times = atoi(argv[i]);
            if (times <= 0) {
                fprintf(stderr, "ERROR: invalid -times\n");
                return false;
            }
        } else if (!strcmp(argv[i], "-verbose")){
            verbose = true;
        } else {
            Usage();
            return false;
        }
    }

    return RunCLK(times, verbose);
}

bool Run2OPT (int times)
{
    TwoOpt opt2(g_eval);
    double totalCost = 0;
    g_maxIter = times;
    for (int i = 0; i < g_maxIter; i++) {
        opt2.DoIt();
        int cost = g_eval->DoIt(opt2.GetFlipper());
        totalCost += cost;
        if (cost < g_bestCost) {
            g_bestCost = cost;
        }
    }
    g_avgCost = totalCost / g_maxIter;
    return true;
}

bool RunEAX (bool verbose)
{
    EAXGA eax(g_eval);
    eax.SetVerbose(verbose);
    eax.DoIt();
    g_maxIter = eax.GetGenNum();
    g_bestCost = eax.GetBestIndi()._cost;
    g_avgCost = eax.GetAvgCost();
    return true;
}

bool RunCLK (int times, bool verbose)
{
    const int num = g_eval->GetNumCity();

    /* generate edge by k-nearest */
    const int max_edge_num = 10;
    int edgeNum = 0;
    for (int ci = 0; ci < num; ++ci) {
        int eCnt = 0;
        for (int ni = 0; ni < g_eval->GetMaxNumNear() && eCnt < max_edge_num; ++ni) {
            int cj = g_eval->GetNear(ci, ni);
            if (cj > ci) {
                eCnt++;
                edgeNum++;
            }
        }
    }
    assert(edgeNum > 0);
    int *edgeList = new int[2*edgeNum];
    int listIdx = 0;
    for (int ci = 0; ci < num; ++ci) {
        int eCnt = 0;
        for (int ni = 0; ni < g_eval->GetMaxNumNear() && eCnt < max_edge_num; ++ni) {
            int cj = g_eval->GetNear(ci, ni);
            if (cj > ci) {
                edgeList[listIdx] = ci;
                edgeList[listIdx+1] = cj;
                eCnt++;
                listIdx += 2;
            }
        }
    }
    assert(listIdx == 2*edgeNum);

    srand((unsigned int)time(0));
    TspLib* dat = const_cast<TspLib*>(g_eval->GetTspLib());
    int in_repeater = num;
    double totalCost = 0;
    g_maxIter = times;
    for (int i = 0;  i< g_maxIter; ++i) {
        if (verbose) {
            printf ("Starting Run %d\n", i);
        }
        double cost;
        if (CClinkern_tour (num, dat, edgeNum, edgeList, 100000000,
                in_repeater, nullptr, nullptr, &cost, !verbose)) {
            fprintf (stderr, "CClinkern_tour failed\n");
        }
        totalCost += cost;
        if (cost < g_bestCost) {
            g_bestCost = (int)cost;
        }
    }
    g_avgCost = totalCost / g_maxIter;

    delete[] edgeList;
    return true;
}
