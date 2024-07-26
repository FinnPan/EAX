#ifndef __THU_UTIL_H
#define __THU_UTIL_H

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <climits>
#include <random>
#include <chrono>
#include <algorithm>
#include <vector>
#include <queue>

namespace thu { /* tsp heuristics */

/* Report elapsed time, cpu time, memory usage. */
class RUsage {
    using Clock = std::chrono::system_clock;
    using TimePoint = std::chrono::time_point<Clock>;
public:
    RUsage ();
    ~RUsage ();
    void Reset ();
    void Report (const char *tag) const;
private:
    static TimePoint GetTimeOfDay () { return Clock::now(); }
    static double GetTimeOfCPU ();
    static void GetProcessMem (double &physMem, double &virtMem);
private:
    TimePoint _dayTime0;
    double    _cpuTime0;
};

/* Parse .tsp file, adapted from Concorde. */
class TspLib {
public:
    TspLib ();
    ~TspLib ();
    bool Init (const char *fileName);
    int GetDimension () const { return _dimension; }
    int HasOptimal () const { return (_optimal > 0); }
    int GetOptimal () const { return _optimal; }
    int EdgeLen (int i, int j) const {
        return (_edgeLen)(_x[i]-_x[j], _y[i]-_y[j]);
    }
private:
    static int EdgeLen_euclid (double xd, double yd);
    static int EdgeLen_att (double xd, double yd);
    static int EdgeLen_euclid_ceiling (double xd, double yd);
private:
    int     _dimension;
    int     _optimal;
    double *_x;
    double *_y;
    int   (*_edgeLen)(double i, double j);
};

/* Two-level tree, adapted from Concorde.
   This is described in the paper "Data structures for traveling
   salesman" by Fredman, Johnson, McGeoch, and Ostheimer in 1995.
   Uses the "group_size" approach described in the paper. */
class Flipper {
    struct ChildNode;
    struct ParentNode {
        ParentNode *adj[2];
        ChildNode  *ends[2];
        int         size;
        int         id;
        int         rev;
    };
    struct ChildNode {
        ParentNode *parent;
        ChildNode  *adj[2];
        int         id;
        int         name;
    };
public:
    explicit Flipper (int count);
    explicit Flipper (int count, const int *cyc);
    ~Flipper ();
    void SetCycle (const int *cyc);
    void GetCycle (int *cyc) const;
    /* returns the successor of x in the current cycle */
    int Next (int x) const;
    /* returns the predecessor of x in the current cycle */
    int Prev (int x) const;
    /* returns true only if xyz occur as an increasing subsequence */
    bool Sequence (int x, int y, int z) const;
    /* flips the portion of the cycle from x to y (inclusive) */
    void Flip (int x, int y);
private:
    void Init_0 ();
    void Init_1 ();
    void Flip_0 (ChildNode *xc,  ChildNode *yc);
    void Reverse () { _reversed ^= 1; }
    bool IsBackward (ParentNode *p) const { return (_reversed^(p->rev)); }
    bool IsForward (ParentNode *p) const { return !IsBackward(p); }
    bool SameSegment (ChildNode *a, ChildNode *b) const;
    void SameSegmentFlip (ChildNode *a, ChildNode *b) const;
    void ConsecutiveSegmentFlip (ParentNode *a, ParentNode *b) const;
    /* split between a and aPrev */
    void SegmentSplit (ParentNode *p, ChildNode *aPrev, ChildNode *a,
            int left_or_right) const;
private:
    int         _numCity;
    int         _groupSize;
    int         _numSegments;
    int         _splitCutoff;
    int         _reversed;
    ParentNode *_parents;
    ChildNode  *_children;
};

/* Utilities: cost, random, neighbor-lists. */
class Evaluator {
public:
    using RandEngine = std::mt19937;
    using RandType = RandEngine::result_type;
    Evaluator ();
    ~Evaluator ();
    bool Init (const char *filename);
    int ComputeCost (const int *arr) const;
    int ComputeCost (const Flipper *f) const;
    const int *MakeRand () const;
    RandEngine &GetRandEngine () const { return *_randEng; }
    RandType GetRand () const { return GetRandEngine()(); }
    int GetMaxNumNear () const { return _maxNumNear; }
    /* return jth-nearest neighbor of i
     * j is valid from 0 to GetMaxNumNear()-1. */
    int GetNear (int i, int j) const { return _nearTbl[i][j]; }
    const TspLib *GetTspLib () const { return &_tspLib; }
    int GetNumCity () const { return _tspLib.GetDimension(); }
    int GetCost (int i, int j) const { return _tspLib.EdgeLen(i, j); }
    double ComputeGap (int cost) const;
private:
    void BuildNeighborLists ();
private:
    const int   _maxNumNear;
    RandEngine *_randEng;
    int        *_randBuf;
    int*       *_nearTbl;
    TspLib      _tspLib;
};

/* 2-exchange speeded up by neighbor-lists.
   This is described in Section 3.3 of "The traveling salesman problem:
   a case study in local optimization" by Johnson, David, Lyle in 1997. */
class TwoOpt {
public:
    explicit TwoOpt (const Evaluator* e) :
        _eval(e), _flipper(e->GetNumCity()) {}
    ~TwoOpt () = default;
    const Flipper* GetFlipper () const { return &_flipper; }
    void DoIt ();
private:
    void TwoExchange ();
private:
    const Evaluator *_eval;
    Flipper          _flipper;
};

} /* namespace thu */

#endif /* __THU_UTIL_H */
