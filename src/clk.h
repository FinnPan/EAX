#ifndef __THU_CLK_H
#define __THU_CLK_H

#include <list>
#include "util.h"

namespace thu { /* tsp heuristics */

class ChainedLK {
    struct FlipStack {
        struct Pair {
            int firstPrev;
            int first;
            int last;
            int lastNext;
        };

        Pair *stack;
        int cnt;
        int num, num0;
        Flipper *flipper;

        FlipStack (int total, int single, Flipper *f);
        ~FlipStack ();
        void Reset ();
        void Flip (int aPrev, int a, int b, int bNext);
        void Unflip (int aPrev, int a, int b, int bNext);
    };

    struct Graph {
        struct Edge {
            int other;
            int weight;
        };

        struct EdgeLook {
            int other;
            int diff;
            int over;
            int seq;
            int side;
            int mm;
        };

        int n;
        Edge **goodList;
        Edge *edgeSpace;
        int *degree;
        int *weirdMark;
        int weirdMagic;

        Graph (int num, int edgeNum);
        ~Graph ();
        void Reset (const Evaluator *eval, int edgeNum, int *edgeList);
        void Insertedge (int c1, int c2, int w);
    };

    struct AddDel {
        int m;
        char *addEdges;
        char *delEdges;

        explicit AddDel (int num);
        ~AddDel ();
        void Reset ();
        void MarkEdgeAdd (int c1, int c2) { addEdges[c1 ^ c2] = 1; }
        void MarkEdgeDel (int c1, int c2) { delEdges[c1 ^ c2] = 1; }
        void UnmarkEdgeAdd (int c1, int c2) { addEdges[c1 ^ c2] = 0; }
        void UnmarkEdgeDel (int c1, int c2) { delEdges[c1 ^ c2] = 0; }
        bool IsAdded (int c1, int c2) const { return addEdges[c1 ^ c2]; }
        bool IsDeleted (int c1, int c2) const { return delEdges[c1 ^ c2]; }
    };

    struct Queue {
        int n;
        char *active;
        std::list<int> q;

        explicit Queue (int num);
        ~Queue ();
        void Reset ();
        void AddToActiveQueue (int c);
        int PopFromActiveQueue ();
    };

public:
    explicit ChainedLK (const Evaluator *eval);
    ~ChainedLK ();
    void SetVerbose (bool v) { _verbose = v; }
    bool DoIt (double &cost);

private:
    void InitEdgeList ();

    void LinKernighan (double *val);
    void LookAheadNoBack (int first, int last, int gain, Graph::EdgeLook* winner);
    void KickTurn (int c);
    void BigTurn (int c, int toNext);
    void FirstKicker (int *t1, int *t2);
    void FindWalkTour (int *t1, int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8);

    int WeirdSecondStep (int lenT1T2, int t1, int t2);
    int Step (int level, int gain, int *Gstar, int first, int last);
    int StepNoBack (int level, int gain, int *Gstar, int first, int last);
    int KickStepNoBack (int level, int gain, int *Gstar, int first, int last);
    int RandomFourSwap (int *delta);

    double ImproveTour (int start);
    double KickImprove ();

    void LookAhead (int first, int last, int gain, int level, std::vector<Graph::EdgeLook>&);
    void WeirdLookAhead (int gain, int t1, int t2, std::vector<Graph::EdgeLook>&);
    void WeirdLookAhead2 (int gain, int t2, int t3, int t4, std::vector<Graph::EdgeLook>&);
    void WeirdLookAhead3 (int gain, int t2, int t3, int t6, std::vector<Graph::EdgeLook>&);

private:
    const Evaluator *_eval;
    int *_edgeList;
    int _edgeNum;
    int *_winCycle;
    Flipper *_flipper;
    FlipStack *_winStack, *_stack;
    Graph *_g;
    AddDel *_ad;
    Queue *_q;
    bool _verbose;
};

} /* namespace thu */

#endif /* __THU_CLK_H */
