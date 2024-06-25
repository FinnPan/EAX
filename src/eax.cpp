#include "eax.h"

namespace thu { /* tsp heuristics */

GA_EAX::GA_EAX (const Evaluator* eval, int numPop, int numKid)
    : _eval(eval), _numPop(numPop), _numKid(numKid),
      _2opt(nullptr), _cross(nullptr), _pop(nullptr), _matingSeq(nullptr),
      _verbose(false), _numGen(0), _avgCost(0), _stagnateGen(0)
{
    const int n = eval->GetNumCity();

    _2opt = new TwoOpt(eval);
    _cross = new Cross(eval);
    _pop = new Tour[numPop+1];
    for (int i = 0; i < numPop+1; ++i) {
        _pop[i].Init(n);
    }
    _matingSeq = new int[numPop];

    for (int i = 0; i < numPop; ++i) {
        _matingSeq[i] = i;
    }
}

GA_EAX::~GA_EAX ()
{
    delete _2opt;
    delete _cross;
    delete[] _pop;
    delete _matingSeq;
}

void GA_EAX::DoIt ()
{
    _numGen = 0;
    _stagnateGen = 0;

    for (int i = 0; i < _numPop; ++i) {
        _2opt->DoIt();
        _pop[i].FromFlipper(_eval, _2opt->GetFlipper());
    }

    do {
        SelectBest();
        if (_verbose) {
            printf("=%d: %d %.3f\n", _numGen, GetBestTour().GetCost(), _avgCost);
        }

        if (ShouldTerminate()) {
            break;
        }

        std::shuffle(_matingSeq, _matingSeq+_numPop, _eval->GetRandEngine());
        for (int i = 0; i < _numPop; ++i) {
            int idx0 = _matingSeq[i];
            int idx1 = _matingSeq[(i+1)%_numPop];
            _cross->DoIt(_pop[idx0], _pop[idx1], _numKid);
        }

        ++_numGen;
    } while (1);
}

void GA_EAX::SelectBest ()
{
    _avgCost = 0.0;

    const int stockBestCost = GetBestTour().GetCost();
    int bestIndex = 0;

    for (int i = 0; i < _numPop; ++i) {
        _avgCost += _pop[i].GetCost();
        if (_pop[i].GetCost() < _pop[bestIndex].GetCost()) {
            bestIndex = i;
        }
    }

    _avgCost /= _numPop;
    GetBestTour() = _pop[bestIndex];

    if (GetBestTour().GetCost() < stockBestCost) {
        _stagnateGen = 0;
    } else {
        _stagnateGen++;
    }
}

bool GA_EAX::ShouldTerminate ()
{
    if (_avgCost - GetBestTour().GetCost() < 0.001) {
        return true;
    }

    if (_stagnateGen > (1500 / _numKid)) {
        return true;
    }

    return false;
}

GA_EAX::Tour::Iterator::Iterator (const Tour& pa, int st)
{
    _pa = &pa;
    _start = st;
    _curr = -1;
    _next = st;
    _cnt = -1;
}

void GA_EAX::Tour::Iterator::operator++ (int)
{
    int prev = _curr;
    _curr = _next;
    if (_pa->GetPrev(_curr) != prev) {
        _next = _pa->GetPrev(_curr);
    } else {
        _next = _pa->GetNext(_curr);
    }
    ++_cnt;
}

GA_EAX::Tour::Tour ()
    : _n(0), _cost(std::numeric_limits<int>::max()), _link(nullptr)
{}

GA_EAX::Tour::~Tour ()
{
    for (int i = 0; i < _n; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void GA_EAX::Tour::Init (int n)
{
    _n = n;
    _link = new int*[n];
    for (int i = 0; i < n; ++i) {
        _link[i] = new int[2];
    }
}

GA_EAX::Tour& GA_EAX::Tour::operator= (const Tour& rhs)
{
    if (this != &rhs) {
        _n = rhs._n;
        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < 2; ++j) {
                _link[i][j] = rhs._link[i][j];
            }
        }
        _cost = rhs._cost;
    }
    return *this;
}

void GA_EAX::Tour::ComputeCost (const Evaluator* e)
{
    _cost = 0;
    for (int i = 0; i < _n; ++i) {
        _cost += e->GetCost(i, _link[i][1]);
    }
}

void GA_EAX::Tour::FromArray (const Evaluator* e, const int* route)
{
    const int n = _n;
    for (int i = 1; i < n - 1; ++i) {
        _link[route[i]][0] = route[i - 1];
        _link[route[i]][1] = route[i + 1];
    }

    _link[route[0]][0] = route[n - 1];
    _link[route[0]][1] = route[1];
    _link[route[n - 1]][0] = route[n - 2];
    _link[route[n - 1]][1] = route[0];

    ComputeCost(e);
}

void GA_EAX::Tour::FromFlipper (const Evaluator* e, const Flipper* f)
{
    for (int i = 0; i < _n; ++i) {
        _link[i][0] = f->Prev(i);
        _link[i][1] = f->Next(i);
    }
    ComputeCost(e);
}

void GA_EAX::EdgeTriple::Reconnect (Tour& pa) const
{
    if (pa.GetPrev(_r1) == _r2) {
        pa.SetPrev(_r1, _b1);
    } else {
        pa.SetNext(_r1, _b1);
    }
    if (pa.GetPrev(_r2) == _r1) {
        pa.SetPrev(_r2, _b2);
    } else {
        pa.SetNext(_r2, _b2);
    }
}

void GA_EAX::ABcycle::Iterator::Begin (const ABcycle* abc, bool reverse)
{
    _abc = abc;
    _len = _abc->GetLen();
    _i = 0;

    if (reverse) {
        /* start from pos_1, backward */
        _start = _len + 1;
        _step = -2;
        for (int k = 0; k < 4; ++k) {
            _offset[k] = -k;
        }
    } else {
        /* start from pos_(len - 1), forward */
        _start = _len - 1;
        _step = 2;
        for (int k = 0; k < 4; ++k) {
            _offset[k] = k;
        }
    }
}

bool GA_EAX::ABcycle::Iterator::End (EdgeTriple& et) const
{
    if (_i < (_len / 2)) {
        int i = _start + _step * _i;
        int b1 = _abc->GetCity((i + _offset[0]) % _len);
        int r1 = _abc->GetCity((i + _offset[1]) % _len);
        int r2 = _abc->GetCity((i + _offset[2]) % _len);
        int b2 = _abc->GetCity((i + _offset[3]) % _len);
        et.Set(b1, r1, r2, b2);
        return true;
    }
    return false;
}

GA_EAX::ABcycle::ABcycle (int n) :
    _maxLen{2 * n}, _len{0}, _gain{0}
{
    _cyc = new int[_maxLen];
}

GA_EAX::ABcyleMgr::ABcyleMgr (const Evaluator *e) :
    _eval(e), _numCity(e->GetNumCity()), _maxNumABcycle(2000)
{
    const int n = _numCity;

    _ABcycleList = new ABcycle*[_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _ABcycleList[j] = new ABcycle(n);
    }

    _overlapEdges = new int*[n];
    for (int j = 0; j < n; ++j) {
        _overlapEdges[j] = new int[5];
    }
    _cycBuf1 = new int[n];
    _cycBuf2 = new int[n];
    _cycBuf1Inv = new int[n];
    _cycBuf2Inv = new int[n];
    _cycRoute = new int[2 * n + 1];
    _checkCycBuf1 = new int[n];
}

GA_EAX::ABcyleMgr::~ABcyleMgr ()
{
    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete _ABcycleList[j];
    }
    delete[] _ABcycleList;

    for (int j = 0; j < _numCity; ++j) {
        delete[] _overlapEdges[j];
    }
    delete[] _overlapEdges;
    delete[] _cycBuf1;
    delete[] _cycBuf2;
    delete[] _cycBuf1Inv;
    delete[] _cycBuf2Inv;
    delete[] _cycRoute;
    delete[] _checkCycBuf1;
}

void GA_EAX::ABcyleMgr::Build (const Tour& pa, const Tour& pb, int numKid)
{
    _numABcycle = 0;

    _cycBuf1Num = 0;
    _cycBuf2Num = 0;
    for (int j = 0; j < _numCity; ++j) {
        _overlapEdges[j][1] = pa.GetPrev(j);
        _overlapEdges[j][3] = pa.GetNext(j);
        _overlapEdges[j][0] = 2;
        _overlapEdges[j][2] = pb.GetPrev(j);
        _overlapEdges[j][4] = pb.GetNext(j);

        _cycBuf1[_cycBuf1Num++] = j;
        _cycBuf1Inv[_cycBuf1[j]] = j;
        _checkCycBuf1[j] = -1;
    }

    /**************************************************/

    int flagSt = 1;
    int prType = 2;
    int flagCircle = 0;
    int posiCurr = 0;
    int r = 0, pr = 0, st = 0, ci = 0;
    while (_cycBuf1Num != 0) {
        if (flagSt == 1) {
            posiCurr = 0;
            r = _eval->GetRand() % _cycBuf1Num;
            st = _cycBuf1[r];
            _checkCycBuf1[st] = posiCurr;
            _cycRoute[posiCurr] = st;
            ci = st;
            prType = 2;
        } else if (flagSt == 0) {
            ci = _cycRoute[posiCurr];
        }

        flagCircle = 0;
        while (flagCircle == 0) {
            posiCurr++;
            pr = ci;

            switch (prType) {
                case 1:
                    ci = _overlapEdges[pr][posiCurr % 2 + 1];
                    break;
                case 2:
                    r = _eval->GetRand() % 2;
                    ci = _overlapEdges[pr][posiCurr % 2 + 1 + 2 * r];
                    if (r == 0) {
                        std::swap(_overlapEdges[pr][posiCurr % 2 + 1],
                                  _overlapEdges[pr][posiCurr % 2 + 3]);
                    }
                    break;
                case 3:
                    ci = _overlapEdges[pr][posiCurr % 2 + 3];
            }

            _cycRoute[posiCurr] = ci;

            if (_overlapEdges[ci][0] == 2) {
                if (ci == st) {
                    if (_checkCycBuf1[st] == 0) {
                        if ((posiCurr - _checkCycBuf1[st]) % 2 == 0) {
                            if (_overlapEdges[st][posiCurr % 2 + 1] == pr) {
                                std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                                          _overlapEdges[ci][posiCurr % 2 + 3]);
                            }
                            if (!Build_0(1, numKid, posiCurr)) {
                                goto LLL;
                            }

                            flagSt = 0;
                            flagCircle = 1;
                            prType = 1;
                        } else {
                            std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                                      _overlapEdges[ci][posiCurr % 2 + 3]);
                            prType = 2;
                        }
                        _checkCycBuf1[st] = posiCurr;
                    } else {
                        if (!Build_0(2, numKid, posiCurr)) {
                            goto LLL;
                        }

                        flagSt = 1;
                        flagCircle = 1;
                    }
                } else if (_checkCycBuf1[ci] == -1) {
                    _checkCycBuf1[ci] = posiCurr;
                    if (_overlapEdges[ci][posiCurr % 2 + 1] == pr) {
                        std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                                  _overlapEdges[ci][posiCurr % 2 + 3]);
                    }
                    prType = 2;
                } else if (_checkCycBuf1[ci] > 0) {
                    std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                              _overlapEdges[ci][posiCurr % 2 + 3]);
                    if ((posiCurr - _checkCycBuf1[ci]) % 2 == 0) {
                        if (!Build_0(1, numKid, posiCurr)) {
                            goto LLL;
                        }

                        flagSt = 0;
                        flagCircle = 1;
                        prType = 1;
                    } else {
                        std::swap(_overlapEdges[ci][(posiCurr + 1) % 2 + 1],
                                  _overlapEdges[ci][(posiCurr + 1) % 2 + 3]);
                        prType = 3;
                    }
                }
            } else if (_overlapEdges[ci][0] == 1) {
                if (ci == st) {
                    if (!Build_0(1, numKid, posiCurr)) {
                        goto LLL;
                    }

                    flagSt = 1;
                    flagCircle = 1;
                } else {
                    prType = 1;
                }
            }
        }
    }

    while (_cycBuf2Num != 0) {
        posiCurr = 0;
        r = _eval->GetRand() % _cycBuf2Num;
        st = _cycBuf2[r];
        _cycRoute[posiCurr] = st;
        ci = st;

        flagCircle = 0;
        while (flagCircle == 0) {
            pr = ci;
            posiCurr++;
            ci = _overlapEdges[pr][posiCurr % 2 + 1];
            _cycRoute[posiCurr] = ci;
            if (ci == st) {
                if (!Build_0(1, numKid, posiCurr)) {
                    goto LLL;
                }

                flagCircle = 1;
            }
        }
    }

LLL:
    std::shuffle(_ABcycleList, _ABcycleList + _numABcycle, _eval->GetRandEngine());
}

bool GA_EAX::ABcyleMgr::Build_0 (const int stAppear, const int numKid, int& posiCurr)
{
    ABcycle* abc = _ABcycleList[_numABcycle];
    const int st = _cycRoute[posiCurr];
    int len = 0;
    int st_count = 0;

    abc->SetCity(len++, st);
    while (1) {
        posiCurr--;
        int ci = _cycRoute[posiCurr];
        if (_overlapEdges[ci][0] == 2) {
            _cycBuf1[_cycBuf1Inv[ci]] = _cycBuf1[_cycBuf1Num - 1];
            _cycBuf1Inv[_cycBuf1[_cycBuf1Num - 1]] = _cycBuf1Inv[ci];
            _cycBuf1Num--;
            _cycBuf2[_cycBuf2Num] = ci;
            _cycBuf2Inv[ci] = _cycBuf2Num;
            _cycBuf2Num++;
        } else if (_overlapEdges[ci][0] == 1) {
            _cycBuf2[_cycBuf2Inv[ci]] = _cycBuf2[_cycBuf2Num - 1];
            _cycBuf2Inv[_cycBuf2[_cycBuf2Num - 1]] = _cycBuf2Inv[ci];
            _cycBuf2Num--;
        }

        _overlapEdges[ci][0]--;
        if (ci == st) {
            st_count++;
        }
        if (st_count == stAppear) {
            break;
        }

        assert(len < abc->GetMaxLen());
        abc->SetCity(len++, ci);
    }

    if (len <= 2) {
        return true;
    }

    abc->SetLen(len);
    if (posiCurr % 2 != 0) {
        int stock = abc->GetCity(0);
        for (int j = 0; j < len - 1; j++) {
            abc->SetCity(j, abc->GetCity(j + 1));
        }
        abc->SetCity(len - 1, stock);
    }

    int diff = 0;
    for (int j = 0; j < len / 2; ++j) {
        int a = abc->GetCity(2 * j);
        int b = abc->GetCity((2 * j + 1) % len);
        int c = abc->GetCity((2 * j + 2) % len);
        diff += _eval->GetCost(a, b) - _eval->GetCost(b, c);
    }
    abc->SetGain(diff);

    _numABcycle++;

    if (_numABcycle >= _maxNumABcycle) {
        fprintf(stderr, "WARNING: maximum number of AB-cycles (%d) must be increased\n",
                _maxNumABcycle);
        return false;
    }

    if (_numABcycle >= numKid) {
        return false;
    }

    return true;
}

int GA_EAX::ABcyleMgr::Apply (int idx, bool reverse, Tour& pa) const
{
    const ABcycle* abc = _ABcycleList[idx];
    ABcycle::Iterator iter;
    EdgeTriple et;
    for (iter.Begin(abc, reverse); iter.End(et); iter++) {
        et.Reconnect(pa);
    }
    return abc->GetGain();
}

GA_EAX::Cross::Cross (const Evaluator* e)
    : _eval(e), _numCity(e->GetNumCity()), _maxNumNear(10)
{
    assert(_maxNumNear <= _eval->GetMaxNumNear());

    const int n = _numCity;

    _city = new int[n];
    _posi = new int[n];
    _abcMgr = new ABcyleMgr(_eval);

    _posiInfo = new PosiInfo[n];
    _numEleInUnit = new int[n];
    _cityInCU = new int[n];
    _routeOfCU = new int[n];

    memset(_cityInCU, 0, sizeof(_cityInCU[0]) * n);
}

GA_EAX::Cross::~Cross ()
{
    delete[] _city;
    delete[] _posi;
    delete _abcMgr;

    delete[] _posiInfo;
    delete[] _numEleInUnit;
    delete[] _cityInCU;
    delete[] _routeOfCU;
}

void GA_EAX::Cross::DoIt (Tour& pa, Tour& pb, int numKid)
{
    /* init _city and _posi */
    Tour::Iterator iter(pa, 0);
    do {
        iter++;
        if (iter.GetCnt() >= _numCity) {
            break;
        }
        int posi = iter.GetCnt();
        int city = iter.GetCurr();
        _city[posi] = city;
        _posi[city] = posi;
    } while (1);

    _abcMgr->Build(pa, pb, numKid);

    /* main loop to generate kids */
    const int numAbc = std::min(numKid, _abcMgr->GetNumCycle());
    int bestIdx = -1, bestGain = 0;
    int aa, bb, a1, b1;
    for (int idx = 0; idx < numAbc; ++idx) {
        int gain = _abcMgr->Apply(idx, false /*reverse*/, pa);

        MakeSegment(idx);
        MakeUnit();
        gain += MakeCompleteTour(pa);
        pa.SetGain(gain);

        if (bestGain < gain) {
            bestIdx = idx;
            _bestModiEdges = _modiEdges;
            bestGain = gain;
        }

        /* back to pa */
        for (auto it = _modiEdges.rbegin(); it != _modiEdges.rend(); ++it) {
            it->Get(aa, a1, bb, b1);
            EdgeTriple(a1, aa, bb, b1).Reconnect(pa);
            EdgeTriple(bb, b1, a1, aa).Reconnect(pa);
        }
        _abcMgr->Apply(idx, true /*reverse*/, pa);
        pa.SetGain(-gain);
    }

    if (bestIdx != -1) {
        /* go to best */
        _abcMgr->Apply(bestIdx, false /*reverse*/, pa);
        for (auto it = _bestModiEdges.begin(); it != _bestModiEdges.end(); ++it) {
            it->Get(aa, bb, a1, b1);
            EdgeTriple(a1, aa, bb, b1).Reconnect(pa);
            EdgeTriple(aa, a1, b1, bb).Reconnect(pa);
        }
        pa.SetGain(bestGain);
    }
}

void GA_EAX::Cross::MakeSegment (int idx)
{
    _segments.clear();

    const ABcycle* abc = _abcMgr->GetCycle(idx);
    const int n = _numCity;
    ABcycle::Iterator iter;
    EdgeTriple et;
    int b1, r1, r2, b2;
    int p1, p2, p0{0};
    bool hasZeroPosi = false;

    for (iter.Begin(abc, false /* reverse */); iter.End(et); iter++) {
        et.Get(b1, r1, r2, b2);
        p1 = _posi[r1];
        p2 = _posi[r2];

        if (p1 == 0 && p2 == n - 1) {
            p0 = p1;
        } else if (p1 == n - 1 && p2 == 0) {
            p0 = p2;
        } else if (p1 < p2) {
            p0 = p2;
        } else if (p2 < p1) {
            p0 = p1;
        } else {
            assert(0);
        }

        _segments.push_back(Segment(p0));
        _posiInfo[p1].SetLinkPosiB1(_posi[b1]);
        _posiInfo[p2].SetLinkPosiB1(_posi[b2]);

        if (p0 == 0) {
            hasZeroPosi = true;
        }
    }

    if (!hasZeroPosi) {
        _segments.push_back(Segment(0));
        _posiInfo[n-1].SetLinkPosiB1(0);
        _posiInfo[0].SetLinkPosiB1(n - 1);
    }

    const int numSeg = _segments.size();
    if (numSeg > n) {
        fprintf(stderr, "ERROR: #segments reach maximum (%d)\n", numSeg);
        exit(1);
    }

    std::sort(_segments.begin(), _segments.end(),
        [](const Segment& l, const Segment& r) {
            return l.begPosi < r.begPosi;
        }
    );
    for (int s = 0; s < numSeg - 1; ++s) {
        _segments[s].endPosi = _segments[s + 1].begPosi - 1;
    }
    _segments[numSeg - 1].endPosi = n - 1;
}

void GA_EAX::Cross::MakeUnit ()
{
    _numUnit = 0;

    int numSeg = 0;
    for (auto& seg : _segments) {
        seg.unitId = -1;
        int p1 = seg.begPosi;
        int p2 = seg.endPosi;
        _posiInfo[p1].SetSegIdAndLinkPosiA(numSeg, p2);
        _posiInfo[p2].SetSegIdAndLinkPosiA(numSeg, p1);
        numSeg++;
    }

    int p_st, p1, p2, p_pre, p_next;
    while (1) {
        bool foundUninit = false;
        for (const auto& seg : _segments) {
            if (seg.unitId == -1) {
                p_st = seg.begPosi;
                p1 = p_st;
                p_pre = -1;
                foundUninit = true;
                break;
            }
        }
        if (!foundUninit) {
            break;
        }

        while (1) {
            int segId = _posiInfo[p1].segId;
            _segments[segId].unitId = _numUnit;

            p2 = _posiInfo[p1].linkPosiA;
            p_next = _posiInfo[p2].linkPosiB1;
            if (p1 == p2) {
                if (p_next == p_pre) {
                    p_next = _posiInfo[p2].linkPosiB2;
                }
            }

            if (p_next == p_st) {
                ++_numUnit;
                break;
            }

            p_pre = p2;
            p1 = p_next;
        }
    }

    int unitId = -1, segId = -1;
    for (auto& seg : _segments) {
        if (seg.unitId != unitId) {
            ++segId;
            _segments[segId] = seg;
            unitId = seg.unitId;
        } else {
            _segments[segId].endPosi = seg.endPosi;
        }
    }
    _segments.resize(segId + 1);

    for (int s = 0; s < _numUnit; ++s) {
        _numEleInUnit[s] = 0;
    }
    for (const auto& seg : _segments) {
        unitId = seg.unitId;
        _numEleInUnit[unitId] += seg.endPosi - seg.begPosi + 1;
    }
}

int GA_EAX::Cross::MakeCompleteTour (Tour& pa)
{
    _modiEdges.clear();

    int gain = 0;

    while (_numUnit != 1) {
        int cu = -1;
        int minNumEle = std::numeric_limits<int>::max();
        for (int u = 0; u < _numUnit; ++u) {
            if (_numEleInUnit[u] < minNumEle) {
                cu = u;
                minNumEle = _numEleInUnit[u];
            }
        }

        int st = -1;
        for (const auto& seg : _segments) {
            if (seg.unitId == cu) {
                int posi = seg.begPosi;
                st = _city[posi];
            }
        }
        assert(st != -1);

        int numEleInCU = 0;
        Tour::Iterator iter(pa, st);
        do {
            iter++;
            int curr = iter.GetCurr();
            _cityInCU[curr] = 1;
            _routeOfCU[numEleInCU++] = curr;
            if (iter.GetNext() == iter.GetStart()) {
                break;
            }
        } while(1);
        assert(numEleInCU == _numEleInUnit[cu]);

        int numNear = _maxNumNear;
        int aa = -1, bb = -1, a1 = -1, b1 = -1;
        int maxDiff = std::numeric_limits<int>::min();

    RESTART:;
        for (int s = 1; s <= numEleInCU; ++s) {
            int a = _routeOfCU[s % numEleInCU];
            for (int n = 0; n < numNear; ++n) {
                int c = _eval->GetNear(a, n);
                if (_cityInCU[c] == 0) {
                    for (int j1 = 0; j1 < 2; ++j1) {
                        int b = _routeOfCU[(s - 1 + 2 * j1) % numEleInCU];
                        for (int j2 = 0; j2 < 2; ++j2) {
                            int d = (j2 == 0)? pa.GetPrev(c) : pa.GetNext(c);
                            int diff = _eval->GetCost(a, b) + _eval->GetCost(c, d)
                                     - _eval->GetCost(a, c) - _eval->GetCost(b, d);
                            if (diff > maxDiff) {
                                aa = a;
                                bb = b;
                                a1 = c;
                                b1 = d;
                                maxDiff = diff;
                            }
                            diff = _eval->GetCost(a, b) + _eval->GetCost(d, c)
                                 - _eval->GetCost(a, d) - _eval->GetCost(b, c);
                            if (diff > maxDiff) {
                                aa = a;
                                bb = b;
                                a1 = d;
                                b1 = c;
                                maxDiff = diff;
                            }
                        }
                    }
                }
            }
        }

        if (a1 == -1 && numNear == _maxNumNear) {
            numNear = _eval->GetMaxNumNear();
            goto RESTART;
        } else if (a1 == -1 && numNear == _eval->GetMaxNumNear()) {
            int r = _eval->GetRand() % (numEleInCU - 1);
            int a = _routeOfCU[r];
            int b = _routeOfCU[r + 1];
            for (int j = 0; j < _numCity; ++j) {
                if (_cityInCU[j] == 0) {
                    aa = a;
                    bb = b;
                    a1 = j;
                    b1 = pa.GetPrev(j);
                    break;
                }
            }
            maxDiff = _eval->GetCost(aa, bb) + _eval->GetCost(a1, b1)
                    - _eval->GetCost(a, a1) - _eval->GetCost(b, b1);
        }

        EdgeTriple(a1, aa, bb, b1).Reconnect(pa);
        EdgeTriple(aa, a1, b1, bb).Reconnect(pa);
        _modiEdges.push_back(EdgeTriple(aa, bb, a1, b1));
        gain += maxDiff;

        int selectUnit = -1;
        for (const auto& seg : _segments) {
            if (seg.begPosi <= _posi[a1] && _posi[a1] <= seg.endPosi) {
                selectUnit = seg.unitId;
                break;
            }
        }
        assert(selectUnit != -1);

        for (auto& seg : _segments) {
            if (seg.unitId == selectUnit) {
                seg.unitId = cu;
            }
        }
        _numEleInUnit[cu] += _numEleInUnit[selectUnit];

        for (auto& seg : _segments) {
            if (seg.unitId == _numUnit - 1) {
                seg.unitId = selectUnit;
            }
        }
        _numEleInUnit[selectUnit] = _numEleInUnit[_numUnit - 1];
        --_numUnit;

        for (int s = 0; s < numEleInCU; ++s) {
            int c = _routeOfCU[s];
            _cityInCU[c] = 0;
        }
    }

    return gain;
}

} /* namespace thu */
