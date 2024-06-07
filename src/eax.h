#ifndef __THU_EAX_H
#define __THU_EAX_H

#include "util.h"

namespace thu { /* tsp heuristics */

/* Genetic Algorithm powered by Edge Assembly Crossover.
 * Some codes are adapted from GA-EAX. */
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
        /* A cycle, such that edges of E_A and edges of E_B are alternately linked. */
        class ABcycle {
        public:
            explicit ABcycle (int n);
            ~ABcycle ();
            int GetCyc (int i) const { return _cyc[i]; }
            void SetCyc (int i, int v) { _cyc[i] = v; }
            int GetGain () const { return _gain; }
            void SetGain (int g) { _gain = g; }
            void ChangeIndi (ABcycle* buf, bool reverse, Indi& pa1) const;
        private:
            int* _cyc;
            int _gain;
        };
        /* Segment representation of an intermediate solution. */
        struct Segment {
            int segId;
            int beginPos, endPos;
            int prevPos, nextPos;
            int tourId;
        };
    public:
        Cross (const Evaluator* e);
        ~Cross ();
        void DoIt (Indi& pa1, Indi& pa2, int nKid);
    private:
        void InitPa1CityPos (const Indi& pa1) const;
        void BuildABcycle (const Indi& pa1, const Indi& pa2, int nKid);
        void BuildABcycle_0 (int stAppear, int& posiCurr);
        void UpdateSeg ();
        void MakeUnit ();
        int MakeCompleteSol (Indi& pa1);
        void BackToPa1 (Indi& pa1, const ABcycle* abc);
        void GoToBest (Indi& pa1, const ABcycle* abc);
    private:
        const Evaluator *_eval;
        const int _numCity;
        const int _maxNumABcycle;

        int *_pa1City, *_pa1Pos;

        ABcycle** _ABcycleList;
        int _numABcycle;

        ABcycle* _ABCycle;
        int** _overlapEdges;
        int *_cycBuf1, *_cycBuf1Inv;
        int *_cycBuf2, *_cycBuf2Inv;
        int _cycBuf1Num, _cycBuf2Num;
        int* _cycRoute;
        int* _checkCycBuf1;

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
