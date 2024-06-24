#ifndef __THU_EAX_H
#define __THU_EAX_H

#include "util.h"

namespace thu { /* tsp heuristics */

/* Genetic Algorithm powered by Edge Assembly Crossover.
 * Some codes are adapted from GA-EAX. */
class GA_EAX {
public:
    GA_EAX (const Evaluator* eval, int numPop = 100, int numKid = 30);
    ~GA_EAX ();
    void DoIt ();
    void SetVerbose (bool s) { _verbose = s; }
    int GetPopNum () const { return _numPop; }
    int GetKidNum () const { return _numKid; }
    int GetBestCost () const { return GetBestTour().GetCost(); }
    int GetGenNum () const { return _numGen; }
    double GetAvgCost () const { return _avgCost; }

private:
    /* A tour based on double-linked list. */
    class Tour {
    public:
        Tour ();
        ~Tour ();
        Tour (const Tour&) = delete;
        Tour& operator= (const Tour&);
        Tour (Tour&&) = delete;
        Tour& operator= (Tour&&) = delete;
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

    /* 3 edges: _b1--_r1, _r1--_r2, _r2--_b2. Maybe B+A+B or A+B+A. */
    class EdgeTriple {
    public:
        explicit EdgeTriple (int b1 = -1, int r1 = -1, int r2 = -1, int b2 = -1) {
            Set(b1, r1, r2, b2);
        }
        ~EdgeTriple () = default;
        void Get (int& b1, int& r1, int& r2, int& b2) const {
            b1 = _b1; r1 = _r1; r2 = _r2; b2 = _b2;
        }
        void Set (int b1, int r1, int r2, int b2) {
            _b1 = b1; _r1 = r1; _r2 = r2; _b2 = b2;
        }
        /* _r1--_r2 will be broken, _b1--_r1 and _r2--_b2 will be connected. */
        void Reconnect (Tour& pa) const;
    private:
        int _b1, _r1, _r2, _b2;
    };

    /* A cycle, such that A edges and B edges are alternately linked. */
    class ABcycle {
    public:
        /* Iterate all edge-triples in ABcycle. */
        class Iterator {
        public:
            Iterator () = default;
            ~Iterator () = default;
            void Begin (const ABcycle* abc, bool reverse);
            bool End (EdgeTriple& et) const;
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
        int GetCity (int i) const { return _cyc[i]; }
        void SetCity (int i, int v) { _cyc[i] = v; }
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

    /* Build all AB-cycles from pa and pb. */
    class ABcyleMgr {
    public:
        explicit ABcyleMgr (const Evaluator *eval);
        ~ABcyleMgr ();
        void Build (const Tour& pa, const Tour& pb, int numKid);
        const ABcycle* GetCycle (int idx) const { return _ABcycleList[idx]; }
        int GetNumCycle () const { return _numABcycle; }
        int Apply (int idx, bool reverse, Tour& pa) const;
    private:
        bool Build_0 (const int stAppear, const int numKid, int& posiCurr);
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

    /* Edge Assembly Crossover.
     * For simplicity, Diversity preservation by Entropy and E-sets by Block2
     * are removed, even if they help achieve the optimal solution. */
    class Cross {
    public:
        explicit Cross (const Evaluator* e);
        ~Cross ();
        void DoIt (Tour& pa, Tour& pb, int numKid);
    private:
        void InitCityPosi (const Tour& pa) const;
        void UndoApply (int idx, const std::vector<EdgeTriple>&, Tour& pa) const;
        void DoApply (int idx, const std::vector<EdgeTriple>&, Tour& pa) const;
        void MakeSegment (int idx);
        void MakeUnit ();
        int MakeCompleteTour (Tour& pa);
    private:
        const Evaluator *_eval;
        const int _numCity;
        const int _maxNumNear;
        int *_city, *_posi;
        ABcyleMgr* _abcMgr;
        std::vector<EdgeTriple> _modiEdges, _bestModiEdges;

        /* link-posi and segment at a poistion */
        struct PosiInfo {
            int linkPosiB1, linkPosiB2;
            int linkPosiA;
            int segId;
            void SetLinkPosiB1 (int b) { linkPosiB2 = linkPosiB1; linkPosiB1 = b; }
            void SetSegIdAndLinkPosiA (int i, int a) { segId = i; linkPosiA = a; }
        };

        /* Segment representation for tour parts of a sub-tour (unit) */
        struct Segment {
            int unitId;
            int begPosi, endPosi;
        };

        int _numSeg;
        Segment* _segment;
        PosiInfo* _posiInfo;
        int _numUnit;
        int* _numEleInUnit;
        int* _cityInCU;
        int* _routeOfCU;
    };

private:
    void SelectBest ();
    bool ShouldTerminate ();
    Tour& GetBestTour () const { return _pop[_numPop]; }
private:
    const Evaluator *_eval;
    const int        _numPop;
    const int        _numKid;
    TwoOpt          *_2opt;
    Cross           *_cross;
    Tour            *_pop;
    int             *_matingSeq;
    bool             _verbose;
    int              _numGen;
    double           _avgCost;
    int              _stagnateGen;
};

} /* namespace thu */

#endif /* __THU_EAX_H */
