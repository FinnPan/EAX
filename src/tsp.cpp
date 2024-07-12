#include "util.h"
#include "eax.h"
#include "clk.h"

static const char *g_progname;
static thu::Evaluator *g_eval;
static int g_bestCost;
static int g_maxIter;
static double g_avgCost;

enum HeuristicKey {
    HK_2opt,
    HK_eax,
    HK_clk,
    HK_NUM
};
static const char *HeuristicKey2Str (int k);
static int HeuristicStr2Key (const char *str);

static void Usage ();
static bool Proc2OPT (int argc, char* argv[]);
static bool ProcEAX (int argc, char* argv[]);
static bool ProcCLK (int argc, char* argv[]);
static bool Run2OPT (int times);
static bool RunEAX (bool verbose);
static bool RunCLK (int times, bool verbose);

int main (int argc, char* argv[])
{
    g_progname = argv[0];
    if (argc < 3) {
        Usage();
        return 1;
    }

    const char* hkStr = argv[1];
    const int hk = HeuristicStr2Key(hkStr);
    if (hk < 0) {
        Usage();
        return 1;
    }

    const char* instFileName = argv[2];
    thu::Evaluator eval;
    if (!eval.Init(instFileName)) {
        return 1;
    }

    g_eval = &eval;
    g_bestCost = INT_MAX;

    bool res = false;
    thu::RUsage ru;
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
        ru.Report(hkStr);
        double gap = g_eval->ComputeGap(g_bestCost);
        printf("[%s results] iter = %d; avg = %0.3f; best = %d; gap = %.3f%%\n",
                hkStr, g_maxIter, g_avgCost, g_bestCost, 100 * gap);
        return 0;
    }
    return 1;
}

const char *HeuristicKey2Str (int k)
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

int HeuristicStr2Key (const char *str)
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
    fprintf(stderr, "  %s 2opt instance_file [-times t]\n", g_progname);
    fprintf(stderr, "  %s eax  instance_file [-verbose]\n", g_progname);
    fprintf(stderr, "  %s clk  instance_file [-times t] [-verbose]\n", g_progname);
}

bool Proc2OPT (int argc, char* argv[])
{
    int times = 20;
    for (int i = 3; i < argc; i++) {
        if(!strcmp(argv[i], "-times")) {
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

bool ProcEAX (int argc, char* argv[])
{
    bool verbose = false;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "-verbose")) {
            verbose = true;
        } else {
            Usage();
            return false;
        }
    }
    return RunEAX(verbose);
}

bool ProcCLK (int argc, char* argv[])
{
    int times = 1;
    bool verbose = false;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "-times")) {
            if(++i >= argc) {
                Usage();
                return false;
            }
            times = atoi(argv[i]);
            if (times <= 0) {
                fprintf(stderr, "ERROR: invalid -times\n");
                return false;
            }
        } else if (!strcmp(argv[i], "-verbose")) {
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
    double totalCost = 0;
    thu::TwoOpt opt2(g_eval);
    for (int i = 0; i < times; i++) {
        opt2.DoIt();
        int cost = g_eval->ComputeCost(opt2.GetFlipper());
        totalCost += cost;
        if (cost < g_bestCost) {
            g_bestCost = cost;
        }
    }
    g_maxIter = times;
    g_avgCost = totalCost / times;
    return true;
}

bool RunEAX (bool verbose)
{
    thu::GA_EAX eax(g_eval);
    eax.SetVerbose(verbose);
    eax.DoIt();
    g_bestCost = eax.GetBestCost();
    g_maxIter = eax.GetGenNum();
    g_avgCost = eax.GetAvgCost();
    return true;
}

bool RunCLK (int times, bool verbose)
{
    double totalCost = 0;
    thu::ChainedLK clk(g_eval);
    clk.SetVerbose(verbose);
    for (int i = 0; i < times; ++i) {
        double cost;
        if (!clk.DoIt(cost)) {
            fprintf (stderr, "ChainedLK failed\n");
        }
        totalCost += cost;
        if (cost < g_bestCost) {
            g_bestCost = static_cast<int>(cost);
        }
    }
    g_maxIter = times;
    g_avgCost = totalCost / times;
    return true;
}
