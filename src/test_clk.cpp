/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*               CODE FOR TESTING ITERATED LIN-KERNIGHAN                    */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 17, 1995                                                    */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/

#include "clk.h"

using namespace th;

int main (int ac, char **av)
{
    const int run_silently = 1;
    const int number_runs = 0;
    const int max_edge_num = 10;

    int rval = 0;
    double val, best;
    int edgeNum, *edgeList = (int *) NULL;
    int ncount, in_repeater, eCnt, listIdx, iter;
    Evaluator eval;
    TspLib* dat = const_cast<TspLib*>(eval.GetTspLib());
    unsigned int seed = (unsigned int) time (0);
    RUsage ru;

    srand(seed);

    if (ac != 2) {
        printf ("please input tsp file\n");
        return 1;
    }
    if (!eval.Init(av[1])) {
        fprintf (stderr, "could not read the TSPLIB file\n");
        rval = 1;
        goto CLEANUP;
    }
    ncount = eval.GetNumCity();
    in_repeater = ncount;

    /* generate edge by k-nearest */
    edgeNum = 0;
    for (int ci = 0; ci < ncount; ++ci) {
        eCnt = 0;
        for (int ni = 0; ni < eval.GetMaxNumNear() && eCnt < max_edge_num; ++ni) {
            int cj = eval.GetNear(ci, ni);
            if (cj > ci) {
                eCnt++;
                edgeNum++;
            }
        }
    }
    assert(edgeNum > 0);
    edgeList = new int[2*edgeNum];
    listIdx = 0;
    for (int ci = 0; ci < ncount; ++ci) {
        eCnt = 0;
        for (int ni = 0; ni < eval.GetMaxNumNear() && eCnt < max_edge_num; ++ni) {
            int cj = eval.GetNear(ci, ni);
            if (cj > ci) {
                edgeList[listIdx] = ci;
                edgeList[listIdx+1] = cj;
                eCnt++;
                listIdx += 2;
            }
        }
    }
    assert(listIdx == 2*edgeNum);

    ru.Reset();
    iter = 0;
    best = DBL_MAX;
    do {
        printf ("Starting Run %d\n", iter);
        if (CClinkern_tour (ncount, dat, edgeNum, edgeList, 100000000,
                in_repeater, nullptr, nullptr, &val, run_silently)) {
            fprintf (stderr, "CClinkern_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (val < best) {
            best = val;
        }
    } while (++iter < number_runs);
    printf ("Overall Best Cycle: %.0f\n", val);
    ru.Report("CLK");

CLEANUP:
    delete[] edgeList;

    return rval;
}
