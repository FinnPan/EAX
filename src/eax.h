#ifndef  __EAX_H
#define  __EAX_H

#include "util.h"

namespace th { /* tsp heuristics */

class Indi {
public:
    Indi ();
    ~Indi ();

    Indi (const Indi&) = delete;
    Indi& operator= (const Indi&);
    Indi (Indi&&) = delete;
    Indi& operator= (Indi&&) = delete;

    void Define (int n);
    void ToArray (int* arr, int* arrInv = nullptr) const;
    void DoneBy (const Evaluator* e);
    void FromFlipper (const Evaluator* e, const Flipper* f);
    void FromArray (const Evaluator* e, const int* route);
    void MadeRand (const Evaluator* e);

    int _numCity;
    int** _link;
    int _cost;
};


/* TODO:
 * 1. Diversity preservation: Entropy
 * 2. E-sets Type: Block2
 */
class Cross {
public:
    Cross (const Evaluator* e, int nPop);
    ~Cross ();

    void DoIt (Indi& kid, Indi& pa2, int nKids);

private:
    void BuildABcycle (const Indi& pa1, const Indi& pa2, int nKids);
    void BuildABcycle_0 (int stAppear, int& posiCurr);
    void ChangeSol (Indi& kid, int idx, bool reverse, bool updateSeg = true);
    void MakeUnit ();
    int MakeCompleteSol (Indi& kid);
    void BackToPa1 (Indi& kid, int appliedCycle);
    void GoToBest (Indi& kid, int bestAppliedCycle);

private:
    const Evaluator* const _eval;
    const int _numCity;
    const int _numPop;
    const int _maxNumABcycle;

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


class EAXGA {
public:
    EAXGA (const Evaluator* eval, int nPop = 100, int nKid = 30);
    ~EAXGA ();

    void DoIt ();

    const Indi& GetBestIndi () const { return _best; }
    int GetPopNum () const { return _numPop; }
    int GetKidsNum () const { return _numKids; }
    int GetGenNum () const { return _numGen; }
    double GetAvgCost () const { return _avgCost; }
    void SetSilent (bool s) { _silent = s; }

private:
    bool Init ();
    void SelectBest ();
    bool ShouldTerminate ();
    void SelectForMating ();

    const Evaluator* _eval;
    int* const _matingSeq;
    TwoOpt* _opt2;
    Cross* _cross;
    Indi* _pop;
    Indi _best;
    const int _numPop;
    const int _numKids;
    bool _silent;
    int _numGen;
    double _avgCost;
    int _stagnGen;
};

} /* namespace th */

#endif /*__EAX_H */
