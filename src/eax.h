#ifndef __THU_EAX_H
#define __THU_EAX_H

#include "util.h"

namespace thu { /* tsp heuristics */

/* Genetic Algorithm powered by Edge Assembly Crossover.
 * Some codes are adapted from GA-EAX by Yuichi Nagata. */
class GA_EAX {
public:
    GA_EAX (const Evaluator* eval, int nPop = 100, int nKid = 30);
    ~GA_EAX ();
    void DoIt ();
    void SetVerbose (bool s) { _verbose = s; }
    int GetPopNum () const { return _numPop; }
    int GetKidNum () const { return _numKid; }
    int GetBestCost () const { return GetBestIndi().GetCost(); }
    int GetGenNum () const { return _numGen; }
    double GetAvgCost () const { return _avgCost; }

private:
    /* An individual based on double-linked list. */
    class Indi {
    public:
        Indi ();
        ~Indi ();
        Indi (const Indi&) = delete;
        Indi& operator= (const Indi&);
        Indi (Indi&&) = delete;
        Indi& operator= (Indi&&) = delete;
        void Init (int n);
        int GetCost () const { return _cost; }
        void SetCost (int c) { _cost = c; }
        int GetPrev (int i) const { return _link[i][0]; }
        int GetNext (int i) const { return _link[i][1]; }
        void SetPrev (int i, int c) { _link[i][0] = c; }
        void SetNext (int i, int c) { _link[i][1] = c; }
        void ComputeCost (const Evaluator* e);
        void FromArray (const Evaluator* e, const int* route);
        void FromFlipper (const Evaluator* e, const Flipper* f);
    private:
        int   _n;
        int **_link;
        int   _cost;
    };

    /* Edge Assembly Crossover.
     * For simplicity, Diversity preservation by Entropy and E-sets by Block2
     * are removed, even if they help achieve the optimal solution. */
    class Cross {
    public:
        Cross (const Evaluator* e);
        ~Cross ();
        void DoIt (Indi& kid, Indi& pa2, int nKid);
    private:
        void ToArray (const Indi& kid, int* arr, int* arrInv) const;
        void BuildABcycle (const Indi& pa1, const Indi& pa2, int nKid);
        void BuildABcycle_0 (int stAppear, int& posiCurr);
        void ChangeSol (Indi& kid, int idx, bool reverse, bool updateSeg = true);
        void MakeUnit ();
        int MakeCompleteSol (Indi& kid);
        void BackToPa1 (Indi& kid, int appliedCycle);
        void GoToBest (Indi& kid, int bestAppliedCycle);
    private:
        const Evaluator *_eval;
        const int        _numCity;
        const int        _maxNumABcycle;

        int *_pa1Route, *_pa1RouteInv;

        int _numABcycle;
        int** _ABcycleList;
        int* _permuABCycle;
        int* _gainABcycle;

        int** _overlapEdges;
        int *_cycBuf1, *_cycBuf1Inv;
        int *_cycBuf2, *_cycBuf2Inv;
        int _cycBuf1Num, _cycBuf2Num;
        int* _cycRoute;
        int* _ABCycle;

        int _numModiEdge;
        int** _modiEdge;
        int _numBestModiEdge;
        int** _bestModiEdge;

        int** _segment;
        int* _segUnit;
        int _numUnit;
        int _numSeg;
        int* _segPosiList;
        int _numSPL;
        int* _linkAPosi;
        int** _linkBPosi;
        int* _posiSeg;
        int* _numElementInUnit;
        int* _centerUnit;
        int* _listCenterUnit;

        int* _routeBuf;
    };

private:
    void SelectBest ();
    bool ShouldTerminate ();
    Indi& GetBestIndi () const { return _pop[_numPop]; }
private:
    const Evaluator *_eval;
    const int        _numPop;
    const int        _numKid;
    TwoOpt          *_2opt;
    Cross           *_cross;
    Indi            *_pop;
    int             *_matingSeq;
    bool             _verbose;
    int              _numGen;
    double           _avgCost;
    int              _stagnateGen;
};

} /* namespace thu */

#endif /* __THU_EAX_H */
