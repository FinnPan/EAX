#include "util.h"

using namespace th;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "ERROR: please input instance file\n");
        return 1;
    }

    Evaluator eval;
    if (!eval.Init(argv[1])) {
        return 1;
    }

    constexpr int max_iter = 20;
    TwoOpt opt2(&eval);
    int totalCost = 0;

    RUsage ru;
    for (int i = 0; i < max_iter; i++) {
        opt2.DoIt();
        totalCost += eval.DoIt(opt2.GetFlipper());
    }
    printf("iter = %d, avg = %0.3lf\n", max_iter, (double)totalCost / max_iter);
    ru.Report("2opt");

    return 0;
}
