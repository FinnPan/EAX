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
/*  correctness or usefulness of curr code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*  int CClinkern_tour (int ncount, TspLib *dat, int ecount,                */
/*      int *elist, int stallcount, int repeatcount, int *incycle,          */
/*      int *outcycle, double *val, int silent)                             */
/*    RUNS Chained Lin-Kernighan.                                           */
/*    -ncount (the number of nodes in the graph)                            */
/*    -dat (coordinate dat)                                                 */
/*    -ecount (the number of good edges - should not be 0)                  */
/*    -elist (the good edges in end1 end2 format)                           */
/*    -stallcount (the max number of 4-swaps without progress               */
/*    -repeatcount (the number of 4-swap kicks)                             */
/*    -incycle (a starting cycle, in node node node format - can be NULL)   */
/*    -outcycle (returns the cycle - can be NULL)                           */
/*    -run_slightly (if nonzero, then very little info will be printed)     */
/*                                                                          */
/*    NOTES: If incycle is NULL, then a random starting cycle is used. If   */
/*     outcycle is not NULL, then it should point to an array of length     */
/*     at least ncount.                                                     */
/*                                                                          */
/****************************************************************************/

#include <cstdlib>
#include <list>
#include "clk.h"

namespace th { /* tsp heuristics */

#define MY_MALLOC(num, type)                                              \
    (type*)malloc((size_t)(num) * sizeof(type))

#define MY_FREE(object, type) {                                            \
    free(object);                                                          \
    object = nullptr;                                                      \
}

#define MAXDEPTH       25   /* Shouldn't really be less than 2.             */
#define KICK_MAXDEPTH  50
#define IMPROVE_SWITCH -1   /* When to start using IMPROVE_KICKS (-1 never) */
#define MARK_LEVEL 10       /* Number of tour neighbors after 4-swap kick   */
#define BACKTRACK   4
#define MAX_BACK   12       /* Upper bound on the XXX_count entries         */
static const int backtrack_count[BACKTRACK] = {4, 3, 3, 2};
static const int weird_backtrack_count[3] = {4, 3, 3};

#define BIGINT 2000000000
#define Edgelen(n1, n2, D)  D->EdgeLen (n1, n2)

#define FLIP(aprev, a, b, bnext, f, x) {                                   \
    (x)->Flip((a), (b));                                 \
    (f)->stack[(f)->counter].first = (a);                                  \
    (f)->stack[(f)->counter++].last = (b);                                 \
}

#define UNFLIP(aprev, a, b, bnext, f, x) {                                 \
    (x)->Flip((b), (a));                                \
    (f)->counter--;                                                        \
}

#define MARK(xn, xQ, xF, xD, xG)  turn ((xn), (xQ))

#define markedge_add(n1, n2, E)    E->add_edges[n1 ^ n2] = 1
#define markedge_del(n1, n2, E)    E->del_edges[n1 ^ n2] = 1
#define unmarkedge_add(n1, n2, E)  E->add_edges[n1 ^ n2] = 0
#define unmarkedge_del(n1, n2, E)  E->del_edges[n1 ^ n2] = 0
#define is_it_added(n1, n2, E)     E->add_edges[n1 ^ n2]
#define is_it_deleted(n1, n2, E)   E->del_edges[n1 ^ n2]

typedef struct edge {
    int other;
    int weight;
} edge;

typedef struct edgelook {
    int other;
    int diff;
    int over;
    int seq;
    int side;
    int mm;
} edgelook;

typedef struct flippair {
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack {
    flippair *stack{nullptr};
    int counter;
    int max;
} flipstack;

typedef struct graph {
    edge **goodlist;
    edge *edgespace;
    int  *degree;
    int  *weirdmark;
    int   weirdmagic;
    int   ncount;
} graph;

typedef struct adddel {
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct aqueue {
    char *active;
    std::list<int> qu;
} aqueue;


static void
   lin_kernighan (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       double *val, int *win_cycle, flipstack *w, flipstack *fstack),
   look_ahead_noback (graph *G, TspLib *D, adddel *E, Flipper *F,
       int first, int last, int gain, edgelook *winner),
   turn (int n, aqueue *Q),
   kickturn (int n, aqueue *Q, TspLib *D, graph *G, Flipper *F),
   bigturn (graph *G, int n, int tonext, aqueue *Q, Flipper *F,
        TspLib *D),
   first_kicker (graph *G, TspLib *D, Flipper *F, int *t1, int *t2),
   find_walk_four (graph *G, TspLib *D, Flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
   insertedge (graph *G, int n1, int n2, int w),
   initgraph (graph *G),
   freegraph (graph *G),
   init_adddel (adddel *E),
   free_adddel (adddel *E),
   init_aqueue (aqueue *Q),
   free_aqueue (aqueue *Q),
   add_to_active_queue (int n, aqueue *Q),
   free_flipstack (flipstack *f),
   rand_cycle (int ncount, int *cyc);

static int
   buildgraph (graph *G, int ncount, int ecount, int *elist, TspLib *D),
   repeated_lin_kernighan (graph *G, TspLib *D, int *cyc,
       int stallcount, int repeatcount, double *val, int silent),
   weird_second_step (graph *G, TspLib *D, adddel *E, aqueue *Q,
       Flipper *F, int gain, int t1, int t2, flipstack *fstack),
   step (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       int level, int gain, int *Gstar, int first, int last, flipstack *fstack),
   step_noback (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       int level, int gain, int *Gstar, int first, int last,
       flipstack *fstack),
   kick_step_noback (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       int level, int gain, int *Gstar, int first, int last, flipstack *win,
       flipstack *fstack),
   random_four_swap (graph *G, TspLib *D, aqueue *Q, Flipper *F,
       int *delta, flipstack *win, flipstack *fstack),
   build_adddel (adddel *E, int ncount),
   build_aqueue (aqueue *Q, int ncount),
   pop_from_active_queue (aqueue *Q),
   init_flipstack (flipstack *f, int total, int single);

static double
   improve_tour (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       int start, flipstack *fstack),
   kick_improve (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
       flipstack *win, flipstack *fstack),
   cycle_length (int ncount, int *cyc, TspLib *D);

static std::vector<edgelook>
   look_ahead (graph *G, TspLib *D, adddel *E, Flipper *F, int first,
       int last, int gain, int level),
   weird_look_ahead (graph *G, TspLib *D, Flipper *F, int gain, int t1,
       int t2),
   weird_look_ahead2 (graph *G, TspLib *D, Flipper *F, int gain, int t2,
       int t3, int t4),
   weird_look_ahead3 (graph *G, TspLib *D, Flipper *F, int gain, int t2,
       int t3, int t6);


int CClinkern_tour (int ncount, TspLib *dat, int ecount,
        int *elist, int stallcount, int repeatcount, int *incycle,
        int *outcycle, double *val, int silent)
{
    int rval = 0;
    int i;
    int *tcyc = (int *) NULL;
    graph G;

    if (silent != 1) {
        printf ("linkern ...\n"); fflush (stdout);
    }

    initgraph (&G);

    if (ncount < 10 && repeatcount > 0) {
        printf ("Less than 10 nodes, setting repeatcount to 0\n");
        fflush (stdout);
        repeatcount = 0;
    }

    tcyc = MY_MALLOC (ncount, int);
    if (tcyc == (int *) NULL) {
        fprintf (stderr, "out of memory in linkern\n");
        rval = 1; goto CLEANUP;
    }

    rval = buildgraph (&G, ncount, ecount, elist, dat);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    if (incycle) {
        for (i = 0; i < ncount; i++) tcyc[i] = incycle[i];
    } else {
        rand_cycle (ncount, tcyc);
    }
    *val = cycle_length (ncount, tcyc, dat);
    if (silent != 1) {
        printf ("Starting Cycle: %.0f\n", *val); fflush (stdout);
    }

    rval = repeated_lin_kernighan (&G, dat, tcyc, stallcount, repeatcount,
                 val, silent);
    if (rval) {
        fprintf (stderr, "repeated_lin_kernighan failed\n"); goto CLEANUP;
    }

    if (outcycle) {
        for (i = 0; i < ncount; i++) outcycle[i] = tcyc[i];
    }

CLEANUP:

    MY_FREE (tcyc, int);
    freegraph (&G);

    return rval;
}

static int repeated_lin_kernighan (graph *G, TspLib *D, int *cyc,
        int stallcount, int count, double *val, int silent)
{
    int rval    = 0;
    int round   = 0;
    int quitcount, hit, delta;
    int *win_cycle = (int *) NULL;
    flipstack winstack, fstack;
    double t, best = *val, oldbest = *val;
    RUsage ru;
    int ncount = G->ncount;
    adddel E;
    Flipper* F;
    aqueue Q;

    (void)oldbest;

    init_aqueue (&Q);
    init_adddel (&E);
    rval = build_aqueue (&Q, ncount);
    if (rval) {
        fprintf (stderr, "build_aqueue failed\n"); goto CLEANUP;
    }
    rval = build_adddel (&E, ncount);
    if (rval) {
        fprintf (stderr, "build_adddel failed\n"); goto CLEANUP;
    }

    hit = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    rval = init_flipstack (&fstack, hit, 0);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }
    rval = init_flipstack (&winstack, 500 + ncount / 50, hit);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }

    win_cycle = MY_MALLOC (ncount, int);
    if (win_cycle == (int *) NULL) {
        fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
        rval = 1; goto CLEANUP;
    }
    win_cycle[0] = -1;

    quitcount = stallcount;
    if (quitcount > count) quitcount = count;

    F = new Flipper(ncount, cyc);
    fstack.counter = 0;
    winstack.counter = 0;
    win_cycle[0] = -1;

    {
        int *tcyc = (int *) NULL;
        int i;

        tcyc = MY_MALLOC (ncount, int);
        if (tcyc == (int *) NULL) {
            fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
            rval = 1; goto CLEANUP;
        }
        /* init active_queue with random order */
        rand_cycle (ncount, tcyc);
        for (i = 0; i < ncount; i++) {
            add_to_active_queue (tcyc[i], &Q);
        }
        MY_FREE (tcyc, int);
    }

    lin_kernighan (G, D, &E, &Q, F, &best, win_cycle, &winstack, &fstack);

    winstack.counter = 0;
    win_cycle[0] = -1;

    if (silent != 1) {
        if (quitcount > 0) {
            printf ("%4d Steps   Best: %.0f\n", round, best);
            ru.Report("quitcount");
        } else {
            printf ("LK Cycle: %.0f\n", best);
        }
        fflush (stdout);
    }

    while (round < quitcount) {
        hit = 0;
        fstack.counter = 0;

        if (IMPROVE_SWITCH == -1 || round < IMPROVE_SWITCH) {
            rval = random_four_swap (G, D, &Q, F, &delta,
                                     &winstack, &fstack);
            if (rval) {
                fprintf (stderr, "random_four_swap failed\n"); goto CLEANUP;
            }
        } else {
            delta = (int)kick_improve (G, D, &E, &Q, F, &winstack, &fstack);
        }

        fstack.counter = 0;
        t = best + delta;
        lin_kernighan (G, D, &E, &Q, F, &t, win_cycle, &winstack, &fstack);

        if (t <= best) {
            winstack.counter = 0;
            win_cycle[0] = -1;
            if (t < best) {
                best = t;
                quitcount = round + stallcount;
                if (quitcount > count)
                    quitcount = count;
                hit++;
            }
        } else {
            if (win_cycle[0] == -1) {
                while (winstack.counter) {
                    winstack.counter--;
                    F->Flip(winstack.stack[winstack.counter].last,
                           winstack.stack[winstack.counter].first);
                }
            } else {
                delete F;
                F = new Flipper(ncount, win_cycle);
                while (winstack.counter) {
                    winstack.counter--;
                    F->Flip(winstack.stack[winstack.counter].last,
                            winstack.stack[winstack.counter].first);
                }
                win_cycle[0] = -1;
            }
        }

        round++;
        if (silent != 1 && (hit || (round % 1000 == 999))) {
            printf ("%4d Steps   Best: %.0f\n", round, best);
            ru.Report("hit");
        }
    }
    if (silent != 1 && round > 0) {
        printf ("%4d Total Steps.\n", round); fflush (stdout);
    }

    F->GetCycle(cyc);
    delete F;
    F = nullptr;

    t = cycle_length (ncount, cyc, D);
    if (t != best) {
        printf ("WARNING: LK incremental counter was off by %.0f\n", t-best);
        fflush (stdout);
        best = t;
    }
    *val = best;

CLEANUP:

    free_aqueue (&Q);
    free_adddel (&E);
    free_flipstack (&fstack);
    free_flipstack (&winstack);
    MY_FREE (win_cycle, int);
    return rval;
}

static void lin_kernighan (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, double *val, int *win_cycle, flipstack *win,
        flipstack *fstack)
{
    int start, i;
    double delta, totalwin = 0.0;

    while (1) {
        start = pop_from_active_queue (Q);
        if (start == -1) break;

        delta = improve_tour (G, D, E, Q, F, start, fstack);
        if (delta > 0.0) {
            totalwin += delta;
            if (win->counter < win->max) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->stack[win->counter].firstprev =
                             fstack->stack[i].firstprev;
                    win->stack[win->counter].lastnext =
                             fstack->stack[i].lastnext;
                    win->counter++;
                }
            } else if (win_cycle[0] == -1) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->counter++;
                }
                F->GetCycle(win_cycle);
            }
            fstack->counter = 0;
        }
    }

    if (win_cycle[0] == -1) {
        for (i = 0; i < fstack->counter; i++) {
            win->stack[win->counter].first = fstack->stack[i].first;
            win->stack[win->counter].last  = fstack->stack[i].last;
            win->counter++;
        }
    }
    (*val) -= totalwin;
}

static double improve_tour (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, int t1, flipstack *fstack)
{
    int t2 = F->Next(t1);
    int gain, Gstar = 0;

    gain = Edgelen (t1, t2, D);
    markedge_del (t1, t2, E);

    if (step (G, D, E, Q, F, 0, gain, &Gstar, t1, t2, fstack)
        == 0) {
        Gstar = weird_second_step (G, D, E, Q, F, gain, t1, t2, fstack);
    }
    unmarkedge_del (t1, t2, E);

    if (Gstar) {
        MARK (t1, Q, F, D, G);
        MARK(t2, Q, F, D, G);
    }
    return (double) Gstar;
}

static int step (graph *G, TspLib *D, adddel *E, aqueue *Q, Flipper *F,
        int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack)
{
    int val, curr, newlast, hit = 0, oldG = gain;

    if (level >= BACKTRACK) {
        return step_noback (G, D, E, Q, F, level, gain, Gstar, first, last,
                            fstack);
    }

    auto list = look_ahead (G, D, E, F, first, last, gain, level);
    for (const auto& e : list) {
        {
            curr = e.other;
            newlast = e.over;

            gain = oldG - e.diff;
            val = gain - Edgelen (newlast, first, D);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }

            FLIP (first, last, newlast, curr, fstack, F);

            if (level < MAXDEPTH) {
                markedge_add (last, curr, E);
                markedge_del (curr, newlast, E);
                hit += step (G, D, E, Q, F, level + 1, gain, Gstar, first,
                             newlast, fstack);
                unmarkedge_add (last, curr, E);
                unmarkedge_del (curr, newlast, E);
            }

            if (!hit) {
                UNFLIP (first, last, newlast, curr, fstack, F);
            } else {
                MARK (curr, Q, F, D, G);
                MARK (newlast, Q, F, D, G);
                return 1;
            }
        }
    }
    return 0;
}

static int step_noback (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack)
{
    edgelook e;

    look_ahead_noback (G, D, E, F, first, last, gain - *Gstar - level, &e);

    if (e.diff < BIGINT) {
        {
            if (e.mm) {
                int hit = 0;
                int curr = e.other;
                int newfirst = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newfirst, last, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }
                FLIP (curr, newfirst, first, last, fstack, F);

                if (level < MAXDEPTH) {
                    markedge_add (first, curr, E);
                    markedge_del (curr, newfirst, E);
                    hit += step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                        newfirst, last, fstack);
                    unmarkedge_add (first, curr, E);
                    unmarkedge_del (curr, newfirst, E);
                }

                if (!hit) {
                    UNFLIP (curr, newfirst, first, last, fstack, F);
                    return 0;
                } else {
                    MARK (curr, Q, F, D, G);
                    MARK (newfirst, Q, F, D, G);
                    return 1;
                }
            } else
            {
                int hit = 0;
                int curr = e.other;
                int newlast = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newlast, first, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }

                FLIP (first, last, newlast, curr, fstack, F);

                if (level < MAXDEPTH) {
                    markedge_add (last, curr, E);
                    markedge_del (curr, newlast, E);
                    hit += step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                        first, newlast, fstack);
                    unmarkedge_add (last, curr, E);
                    unmarkedge_del (curr, newlast, E);
                }

                if (!hit) {
                    UNFLIP (first, last, newlast, curr, fstack, F);
                    return 0;
                } else {
                    MARK (curr, Q, F, D, G);
                    MARK (newlast, Q, F, D, G);
                    return 1;
                }
            }
        }
    } else {
        return 0;
    }
}

static double kick_improve (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, flipstack *win, flipstack *fstack)
{
    int t1, t2;
    int gain, Gstar = 0;
    int hit = 0;

    do {
        first_kicker (G, D, F, &t1, &t2);
        gain = Edgelen (t1, t2, D);
        markedge_del (t1, t2, E);
        hit = kick_step_noback (G, D, E, Q, F, 0, gain, &Gstar, t1, t2, win, fstack);
        unmarkedge_del (t1, t2, E);
    } while (!hit);

    kickturn (t1, Q, D, G, F);
    kickturn (t2, Q, D, G, F);

    return (double) -Gstar;
}

#define G_MULT 1.5

static int kick_step_noback (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *win, flipstack *fstack)
{
    edgelook winner;
    int val;
    int curr, prev, newlast;
    int lastnext = F->Next(last);
    int i;
    int cutoff = (int) (G_MULT * (double) gain);
    edge **goodlist = G->goodlist;

    winner.diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < cutoff; i++) {
        curr = goodlist[last][i].other;
        if (!is_it_deleted (last, curr, E) && curr != first &&
                                              curr != lastnext) {
            prev = F->Prev(curr);
            if (!is_it_added (curr, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (curr, prev, D);
                if (val < winner.diff) {
                    winner.diff = val;
                    winner.other = curr;
                    winner.over = prev;
                }
            }
        }
    }

    if (winner.diff < BIGINT) {
        curr = winner.other;
        newlast = winner.over;
        gain -= winner.diff;
        *Gstar = gain - Edgelen (newlast, first, D);

        FLIP (first, last, newlast, curr, fstack, F);
        kickturn (curr, Q, D, G, F);
        kickturn (newlast, Q, D, G, F);
        if (win->counter < win->max) {
            win->stack[win->counter].first = last;
            win->stack[win->counter].last = newlast;
            win->counter++;
        }

        if (level < KICK_MAXDEPTH) {
            markedge_add (last, curr, E);
            markedge_del (curr, newlast, E);
            kick_step_noback (G, D, E, Q, F, level+1, gain, Gstar, first,
                              newlast, win, fstack);
            unmarkedge_add (last, curr, E);
            unmarkedge_del (curr, newlast, E);
        }
        return 1;
    } else {
        return 0;
    }
}

static int weird_second_step (graph *G, TspLib *D, adddel *E, aqueue *Q,
        Flipper *F, int len_t1_t2, int t1, int t2, flipstack *fstack)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit;
    int t3prev, t4next;

    (void)t3prev;

    auto list = weird_look_ahead (G, D, F, len_t1_t2, t1, t2);
    for (const auto& h : list) {
        t3 = h.other;
        t4 = h.over;

        oldG = len_t1_t2 - h.diff;

        t3prev = F->Prev(t3);
        t4next = F->Next(t4);

        markedge_add (t2, t3, E);
        markedge_del (t3, t4, E);
        G->weirdmagic++;
        G->weirdmark[t1] = G->weirdmagic;
        G->weirdmark[t2] = G->weirdmagic;
        G->weirdmark[t3] = G->weirdmagic;
        G->weirdmark[t4next] = G->weirdmagic;

        auto list2 = weird_look_ahead2 (G, D, F, oldG, t2, t3, t4);
        for (const auto& e : list2) {
            t5 = e.other;
            t6 = e.over;

            markedge_add (t4, t5, E);
            if (e.seq) {
                if (!e.side) {
                    gain = oldG - e.diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP (t1, t2, t6, t5, fstack, F);
                    FLIP (t2, t5, t3, t4, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6, fstack);
                    unmarkedge_del (t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit) {
                        UNFLIP (t2, t5, t3, t4, fstack, F);
                        UNFLIP (t1, t2, t6, t5, fstack, F);
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G);
                        MARK (t4, Q, F, D, G);
                        MARK (t5, Q, F, D, G);
                        MARK (t6, Q, F, D, G);
                        return Gstar;
                    }
                } else {
                    gain = oldG - e.diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP (t1, t2, t3, t4, fstack, F);
                    FLIP (t6, t5, t2, t4, fstack, F);
                    FLIP (t1, t3, t6, t2, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6, fstack);
                    unmarkedge_del (t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit) {
                        UNFLIP (t1, t3, t6, t2, fstack, F);
                        UNFLIP (t6, t5, t2, t4, fstack, F);
                        UNFLIP (t1, t2, t3, t4, fstack, F);
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G);
                        MARK (t4, Q, F, D, G);
                        MARK (t5, Q, F, D, G);
                        MARK (t6, Q, F, D, G);
                        return Gstar;
                    }
                }
            } else {
                tG = oldG - e.diff;
                markedge_del (t5, t6, E);
                auto list3 = weird_look_ahead3 (G, D, F, tG, t2, t3, t6);
                for (const auto& f : list3) {
                    t7 = f.other;
                    t8 = f.over;
                    gain = tG - f.diff;
                    if (!f.side) {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP (t1, t2, t8, t7, fstack, F);
                        FLIP (t2, t7, t3, t4, fstack, F);
                        FLIP (t7, t4, t6, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8, fstack);
                        unmarkedge_del (t6, t7, E);
                        unmarkedge_del (t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit) {
                            UNFLIP (t7, t4, t6, t5, fstack, F);
                            UNFLIP (t2, t7, t3, t4, fstack, F);
                            UNFLIP (t1, t2, t8, t7, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G);
                            MARK (t4, Q, F, D, G);
                            MARK (t5, Q, F, D, G);
                            MARK (t6, Q, F, D, G);
                            MARK (t7, Q, F, D, G);
                            MARK (t8, Q, F, D, G);
                            return Gstar;
                        }
                    } else {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP (t1, t2, t6, t5, fstack, F);
                        FLIP (t1, t6, t8, t7, fstack, F);
                        FLIP (t3, t4, t2, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8, fstack);
                        unmarkedge_add (t6, t7, E);
                        unmarkedge_del (t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit) {
                            UNFLIP (t3, t4, t2, t5, fstack, F);
                            UNFLIP (t1, t6, t8, t7, fstack, F);
                            UNFLIP (t1, t2, t6, t5, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G);
                            MARK (t4, Q, F, D, G);
                            MARK (t5, Q, F, D, G);
                            MARK (t6, Q, F, D, G);
                            MARK (t7, Q, F, D, G);
                            MARK (t8, Q, F, D, G);
                            return Gstar;
                        }
                    }
                }
                unmarkedge_del (t5, t6, E);
            }
            unmarkedge_add (t4, t5, E);
        }
        unmarkedge_add (t2, t3, E);
        unmarkedge_del (t3, t4, E);
    }
    return 0;
}

static std::vector<edgelook> look_ahead (graph *G, TspLib *D, adddel *E,
    Flipper *F, int first, int last, int gain, int level)
{
    std::vector<edgelook> list;
    int i, val;
    int curr, prev;
    int lastnext = F->Next(last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, ahead = backtrack_count[level];
    edge **goodlist = G->goodlist;

    for (i = 0; i < ahead; i++) {
        value[i] = BIGINT;
    }
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[last][i].weight <= gain; i++) {
        curr = goodlist[last][i].other;
        if (!is_it_deleted (last, curr, E) && curr != first &&
                                              curr != lastnext) {
            prev = F->Prev(curr);
            if (!is_it_added (curr, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (curr, prev, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                    }
                    value[k] = val;
                    other[k] = curr;
                    save[k] = prev;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            list.push_back(edgelook{});
            list.back().diff = value[i];
            list.back().other = other[i];
            list.back().over = save[i];
        }
    }
    std::reverse(list.begin(), list.end());

    return list;
}

static void look_ahead_noback (graph *G, TspLib *D, adddel *E, Flipper *F,
        int first, int last, int gain, edgelook *winner)
{
    int val;
    int curr, prev;
    int lastnext = F->Next(last);
    int i;
    int next;
    edge **goodlist = G->goodlist;

    winner->diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < gain; i++) {
        curr = goodlist[last][i].other;
        if (!is_it_deleted (last, curr, E) && curr != first &&
                                              curr != lastnext) {
            prev = F->Prev(curr);
            if (!is_it_added (curr, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (curr, prev, D);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = curr;
                    winner->over = prev;
                    winner->mm = 0;
                }
            }
        }
    }
    {
        int firstprev = F->Prev(first);

        for (i = 0; goodlist[first][i].weight < gain; i++) {
            curr = goodlist[first][i].other;
            if (!is_it_deleted (first, curr, E) && curr != last &&
                                                   curr != firstprev) {
                next = F->Next(curr);
                if (!is_it_added (curr, next, E)) {
                    val = goodlist[first][i].weight - Edgelen (curr, next, D);
                    if (val < winner->diff) {
                        winner->diff = val;
                        winner->other = curr;
                        winner->over = next;
                        winner->mm = 1;
                    }
                }
            }
        }
    }
}

static std::vector<edgelook> weird_look_ahead (graph *G, TspLib *D, Flipper *F,
        int gain, int t1, int t2)
{
    std::vector<edgelook> list;
    int i, curr, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val, ahead;
    edge **goodlist = G->goodlist;

    ahead = weird_backtrack_count[0];
    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t2][i].weight <= gain; i++) {
        curr = goodlist[t2][i].other;
        if (curr != t1) {
            next = F->Next(curr);
            val = goodlist[t2][i].weight - Edgelen (curr, next, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                }
                value[k] = val;
                other[k] = curr;
            save[k] = next;
            }
        }
    }
    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            list.push_back(edgelook{});
            list.back().diff = value[i];
            list.back().other = other[i];
            list.back().over = save[i];
        }
    }
    std::reverse(list.begin(), list.end());
    return list;
}

static std::vector<edgelook> weird_look_ahead2 (graph *G, TspLib *D, Flipper *F,
       int gain, int t2, int t3, int t4)
{
    std::vector<edgelook> list;
    int i, t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[1];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t4][i].weight <= gain; i++) {
        t5 = goodlist[t4][i].other;
        if (weirdmark[t5] != weirdmagic) {
            if (F->Sequence(t2, t5, t3)) {
                t6 = F->Prev(t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 0;
                }
                t6 = F->Next(t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 1;
                }
            } else {
                t6 = F->Prev(t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 0;
                    side[k] = 0;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            list.push_back(edgelook{});
            list.back().diff = value[i];
            list.back().other = other[i];
            list.back().over = save[i];
            list.back().seq = seq[i];
            list.back().side = side[i];
        }
    }
    std::reverse(list.begin(), list.end());
    return list;
}

static std::vector<edgelook> weird_look_ahead3 (graph *G, TspLib *D, Flipper *F,
        int gain, int t2, int t3, int t6)
{
    std::vector<edgelook> list;
    int i, t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[2];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

    for (i = 0; goodlist[t6][i].weight <= gain; i++) {
        t7 = goodlist[t6][i].other;   /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdmark[t7] != weirdmagic && F->Sequence(t2, t7, t3)) {
            t8 = F->Prev(t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 0;
            }
            t8 = F->Next(t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 1;
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            list.push_back(edgelook{});
            list.back().diff = value[i];
            list.back().other = other[i];
            list.back().over = save[i];
            list.back().side = side[i];
        }
    }
    std::reverse(list.begin(), list.end());
    return list;
}

static double cycle_length (int ncount, int *cyc, TspLib *D)
{
    int i;
    double val = 0.0;

    for (i = 1; i < ncount; i++) {
        val += (double) Edgelen (cyc[i - 1], cyc[i], D);
    }
    val += (double) Edgelen (cyc[0], cyc[ncount - 1], D);

    return val;
}

static int random_four_swap (graph *G, TspLib *D, aqueue *Q, Flipper *F,
       int *delta, flipstack *win, flipstack *fstack)
{
    int t1, t2, t3, t4, t5, t6, t7, t8;

    find_walk_four (G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);

    if (!F->Sequence(t1, t3, t5)) {
        std::swap (t3, t5);
        std::swap (t4, t6);
    }
    if (!F->Sequence(t1, t5, t7)) {
        std::swap (t5, t7);
        std::swap (t6, t8);
        if (!F->Sequence(t1, t3, t5)) {
            std::swap (t3, t5);
            std::swap (t4, t6);
        }
    }
    FLIP (t1, t2, t5, t6, fstack, F);
    FLIP (t4, t3, t7, t8, fstack, F);
    FLIP (t1, t5, t6, t2, fstack, F);

    if (win->counter < win->max) {
        win->stack[win->counter].first = t2;
        win->stack[win->counter].last = t5;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t3;
        win->stack[win->counter].last = t7;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t5;
        win->stack[win->counter].last = t6;
        win->counter++;
    }

    bigturn (G ,t1, 0, Q, F, D);
    bigturn (G, t2, 1, Q, F, D);
    bigturn (G, t3, 0, Q, F, D);
    bigturn (G, t4, 1, Q, F, D);
    bigturn (G, t5, 0, Q, F, D);
    bigturn (G, t6, 1, Q, F, D);
    bigturn (G, t7, 0, Q, F, D);
    bigturn (G, t8, 1, Q, F, D);

    *delta =
           Edgelen (t1, t6, D) + Edgelen (t2, t5, D) +
           Edgelen (t3, t8, D) + Edgelen (t4, t7, D) -
           Edgelen (t1, t2, D) - Edgelen (t3, t4, D) -
           Edgelen (t5, t6, D) - Edgelen (t7, t8, D);
    return 0;
}

#define HUNT_PORTION_LONG 0.001

static void first_kicker (graph *G, TspLib *D, Flipper *F, int *t1,
        int *t2)
{
    int longcount = (int) ((double) G->ncount * HUNT_PORTION_LONG) + 10;
    int i, best, try1, len, next, prev, nextl, prevl;
    int ncount = G->ncount;
    edge **goodlist = G->goodlist;

    try1 = rand() % ncount;
    next = F->Next(try1);
    prev = F->Prev(try1);
    nextl = Edgelen (try1, next, D);
    prevl = Edgelen (try1, prev, D);
    if (nextl >= prevl) {
        *t1 = try1;
        *t2 = next;
        best = nextl - goodlist[*t1][0].weight;
    } else {
        *t1 = prev;
        *t2 = try1;
        best = prevl - goodlist[*t1][0].weight;
    }

    for (i = 0; i < longcount; i++) {
        try1 = rand() % ncount;
        next = F->Next(try1);
        prev = F->Prev(try1);
        nextl = Edgelen (try1, next, D);
        prevl = Edgelen (try1, prev, D);
        if (nextl >= prevl) {
            len = nextl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = try1;
                *t2 = next;
            }
        } else {
            len = prevl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
}


#define WALK_STEPS 50

static void find_walk_four (graph *G, TspLib *D, Flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, n, i, j;

    first_kicker (G, D, F, &s1, &s2);

    do {
        old = -1;
        n = s2;

        for (i = 0;  i < WALK_STEPS; i++) {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s3 = n;
        s4 = F->Next(s3);

        n = s4;
        for (i = 0; i < WALK_STEPS; i++) {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s5 = n;
        s6 = F->Next(s5);

        n = s6;
        for (i = 0; i < WALK_STEPS; i++) {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s7 = n;
        s8 = F->Next(s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7 ||
             s1 == s8 ||
             s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6 || s2 == s7 ||
             s2 == s8 ||
             s3 == s5 || s3 == s6 || s3 == s7 || s3 == s8 ||
             s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8 ||
             s5 == s7 || s5 == s8 ||
             s6 == s7 || s6 == s8);

    *t1 = s1;  *t2 = s2;  *t3 = s3;  *t4 = s4;
    *t5 = s5;  *t6 = s6;  *t7 = s7;  *t8 = s8;
}

static void turn (int n, aqueue *Q)
{
    add_to_active_queue (n, Q);
}

static void kickturn (int n, aqueue *Q, TspLib *D,
        graph *G, Flipper *F)
{
    (void)D;
    (void)G;
    add_to_active_queue (n, Q);
    {
        int k;
        k = F->Next(n);
        add_to_active_queue (k, Q);
        k = F->Next(k);
        add_to_active_queue (k, Q);
        k = F->Prev(n);
        add_to_active_queue (k, Q);
        k = F->Prev(k);
        add_to_active_queue (k, Q);
    }
}

static void bigturn (graph *G, int n, int tonext, aqueue *Q, Flipper *F,
        TspLib *D)
{
    int i, k;

    (void)D;

    add_to_active_queue (n, Q);
    if (tonext) {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = F->Next(k);
            add_to_active_queue (k, Q);
        }
    } else {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = F->Prev(k);
            add_to_active_queue (k, Q);
        }
    }

    for (i = 0; i < G->degree[n]; i++) {
        add_to_active_queue (G->goodlist[n][i].other, Q);
    }
}

static void initgraph (graph *G)
{
    G->goodlist   = (edge **) NULL;
    G->edgespace  = (edge *) NULL;
    G->degree     = (int *) NULL;
    G->weirdmark  = (int *) NULL;
    G->weirdmagic = 0;
    G->ncount     = 0;
}

static void freegraph (graph *G)
{
    if (G) {
        MY_FREE (G->goodlist, edge *);
        MY_FREE (G->edgespace, edge);
        MY_FREE (G->degree, int);
        MY_FREE (G->weirdmark, int);
        G->weirdmagic = 0;
        G->ncount     = 0;
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist,
        TspLib *D)
{
    int rval = 0;
    int n1, n2, w, i;
    edge *p;

    G->goodlist  = MY_MALLOC (ncount, edge *);
    G->degree    = MY_MALLOC (ncount, int);
    G->weirdmark = MY_MALLOC (ncount, int);
    G->edgespace = MY_MALLOC ((2 * ecount) + ncount, edge);
    if (G->goodlist == (edge **) NULL || G->degree == (int *) NULL ||
        G->edgespace == (edge *) NULL)  {
        fprintf (stderr, "out of memory in buildgraph\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        G->degree[i] = 1;
        G->weirdmark[i] = 0;
    }
    for (i = ecount - 1; i >= 0; i--) {
        G->degree[elist[2 * i]]++;
        G->degree[elist[(2 * i) + 1]]++;
    }

    for (i = 0, p = G->edgespace; i < ncount; i++) {
        G->goodlist[i] = p;
        p += (G->degree[i]);
        G->goodlist[i][G->degree[i] - 1].weight = BIGINT;
        G->degree[i] = 0;
    }

    for (i = ecount - 1; i >= 0; i--) {
        n1 = elist[2 * i];
        n2 = elist[(2 * i) + 1];
        w = Edgelen (n1, n2, D);
        insertedge (G, n1, n2, w);
        insertedge (G, n2, n1, w);
    }
    G->ncount     = ncount;
    G->weirdmagic = 0;

CLEANUP:

    if (rval) freegraph (G);
    return rval;
}

static void insertedge (graph *G, int n1, int n2, int w)
{
    int i;
    edge *e = G->goodlist[n1];

    for (i = G->degree[n1] - 1; i >= 0 && e[i].weight >= w; i--) {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other  = e[i].other;
    }
    e[i + 1].weight = w;
    e[i + 1].other  = n2;
    G->degree[n1]++;
}

static int init_flipstack (flipstack *f, int total, int single)
{
    f->counter = 0;
    f->max     = 0;
    f->stack   = (flippair *) NULL;

    f->stack = MY_MALLOC (total + single, flippair);
    if (f->stack == (flippair *) NULL) {
        fprintf (stderr, "out of memory in init_flipstack\n"); return 1;
    }
    f->max = total;

    return 0;
}

static void free_flipstack (flipstack *f)
{
    f->counter = 0;
    f->max     = 0;
    MY_FREE (f->stack, flippair);
}

static void init_adddel (adddel *E)
{
    E->add_edges = (char *) NULL;
    E->del_edges = (char *) NULL;
}

static void free_adddel (adddel *E)
{
    if (E) {
        MY_FREE (E->add_edges, char);
        MY_FREE (E->del_edges, char);
    }
}

static int build_adddel (adddel *E, int ncount)
{
    int rval = 0;
    int i, M;

    i = 0;
    while ((1 << i) < ncount)
        i++;
    M = (1 << i);

    E->add_edges = MY_MALLOC (M, char);
    E->del_edges = MY_MALLOC (M, char);
    if (E->add_edges == (char *) NULL || E->del_edges == (char *) NULL) {
        fprintf (stderr, "out of memory in build_adddel\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < M; i++) {
        E->add_edges[i] = 0;
        E->del_edges[i] = 0;
    }

CLEANUP:

    if (rval) {
        free_adddel (E);
    }
    return rval;
}

static void init_aqueue (aqueue *Q)
{
    Q->active = (char *) NULL;
    Q->qu.clear();
}

static void free_aqueue (aqueue *Q)
{
    if (Q) {
        MY_FREE (Q->active, char);
        Q->qu.clear();
    }
}

static int build_aqueue (aqueue *Q, int ncount)
{
    int rval = 0;
    int i;

    init_aqueue (Q);

    Q->active = MY_MALLOC (ncount, char);
    if (Q->active == (char *) NULL) {
        fprintf (stderr, "out of memory in build_aqueue\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) Q->active[i] = 0;

CLEANUP:

    if (rval) {
        free_aqueue (Q);
    }
    return rval;
}

static void add_to_active_queue (int n, aqueue *Q)
{
    if (Q->active[n] == 0) {
        Q->active[n] = 1;
        Q->qu.push_back(n);
    }
}

static int pop_from_active_queue (aqueue *Q)
{
    int n = -1;
    if (!Q->qu.empty()) {
        n = Q->qu.front();
        Q->qu.pop_front();
        Q->active[n] = 0;
    }
    return n;
}

static void rand_cycle (int ncount, int *cyc)
{
    int i, k;

    for (i = 0; i < ncount; i++) cyc[i] = i;
    for (i = ncount; i > 1; i--) {
        k = rand() % i;
        std::swap(cyc[i - 1], cyc[k]);
    }
}

} /* namespace th */
