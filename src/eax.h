#ifndef __THU_EAX_H
#define __THU_EAX_H

#include <tuple>
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
        void SetGain (int g) { _cost -= g; }
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

    /* A cycle, such that A edges and B edges are alternately linked. */
    class ABcycle {
    public:
        /* Iterate all B+A+B or A+B+A edge-tuples in ABcycle. */
        class Iterator {
        public:
            Iterator () = default;
            ~Iterator () = default;
            void Begin (const ABcycle* abc, bool reverse);
            /* return edge-tuples: b1--r1, r1--r2, r2--b2 */
            bool End (int& b1, int& r1, int& r2, int& b2) const;
            void operator++ (int) { ++_i; }
        private:
            const ABcycle* _abc;
            int _len, _i, _start, _step, _offset[4];
        };

        explicit ABcycle (int n);
        ~ABcycle () { delete _cyc; }
        int GetMaxLen () const { return _maxLen; }
        int GetLen () const { return _len; }
        void SetLen (int l) { _len = l; }
        int GetCyc (int i) const { return _cyc[i]; }
        void SetCyc (int i, int v) { _cyc[i] = v; }
        int GetGain () const { return _gain; }
        void SetGain (int g) { _gain = g; }
    private:
        const int _maxLen;
        int _len;
        /* gain when replacing A edges with B edges */
        int _gain;
        /* (2*i)--(2*i+1) are A edges, (2*i+1)--(2*i+2) are B edges. */
        int* _cyc;
    };

    using EdgeTuple = std::tuple<int, int, int, int>;
    using EdgeTuples = std::vector<EdgeTuple>;

    /* Build all AB-cycles from pa1 and pa2. */
    class ABcyleMgr {
    public:
        explicit ABcyleMgr (const Evaluator *eval);
        ~ABcyleMgr ();
        void Build (const Indi& pa1, const Indi& pa2, int nKid);
        const ABcycle* GetABCycle (int idx) const { return _ABcycleList[idx]; }
        int GetNumABCycle () const { return _numABcycle; }
        void ChangeIndi (int idx, bool reverse, Indi& pa1) const;
        void BackToPa1 (int idx, const EdgeTuples& me, Indi& pa1) const;
        void GoToBest (int idx, const EdgeTuples& me, Indi& pa1) const;
        static void BreakEdge (int b1, int r1, int r2, int b2, Indi& pa1);
    private:
        bool Build_0 (const int stAppear, const int nKid, int& posiCurr);
    private:
        const Evaluator *_eval;
        const int _numCity;
        const int _maxNumABcycle;
        ABcycle** _ABcycleList;
        int _numABcycle;
        int** _overlapEdges;
        int *_cycBuf1, *_cycBuf1Inv;
        int *_cycBuf2, *_cycBuf2Inv;
        int _cycBuf1Num, _cycBuf2Num;
        int* _cycRoute;
        int* _checkCycBuf1;
    };

    /* Segment representation of an intermediate solution. */
    struct Segment {
        int segId;
        int beginPos, endPos;
        int prevPos, nextPos;
        int tourId;
    };

    /* Edge Assembly Crossover.
     * For simplicity, Diversity preservation by Entropy and E-sets by Block2
     * are removed, even if they help achieve the optimal solution. */
    class Cross {
    public:
        explicit Cross (const Evaluator* e);
        ~Cross ();
        void DoIt (Indi& pa1, Indi& pa2, int nKid);
    private:
        void InitPa1CityPos (const Indi& pa1) const;
        void UpdateSeg (int idx);
        void MakeUnit ();
        int MakeCompleteSol (Indi& pa1);
    private:
        const Evaluator *_eval;
        const int _numCity;
        int *_pa1City, *_pa1Pos;
        ABcyleMgr* _abcMgr;
        EdgeTuples _modiEdges, _bestModiEdges;

        int** _segment;
        int _numSeg;
        int* _segUnit;
        int _numUnit;
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
