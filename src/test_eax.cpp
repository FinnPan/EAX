#include "eax.h"

using namespace th;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "ERROR: please input instance file\n");
        return 0;
    }

    Evaluator eval;
    if (!eval.Init(argv[1])) {
        return 1;
    }

    EAXGA eax(&eval);
    eax.SetSilent(true);

    RUsage ru;
    eax.DoIt();
    printf("iter = %d, best = %d, avg = %0.3f\n",
           eax.GetGenNum(), eax.GetBestIndi()._cost, eax.GetAvgCost());
    ru.Report("GA-EAX");

    return 0;
}
