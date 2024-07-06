#include "eax.h"

namespace thu { /* tsp heuristics */

GA_EAX::GA_EAX (const Evaluator *eval, int numPop, int numKid)
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

GA_EAX::Tour::Iterator::Iterator (const Tour& pa, int start)
{
    _pa = &pa;
    _start = start;
    _curr = -1;
    _next = start;
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

void GA_EAX::Tour::ComputeCost (const Evaluator *e)
{
    _cost = 0;
    for (int i = 0; i < _n; ++i) {
        _cost += e->GetCost(i, _link[i][1]);
    }
}

void GA_EAX::Tour::FromArray (const Evaluator *e, const int *arr)
{
    const int n = _n;
    for (int i = 1; i < n - 1; ++i) {
        _link[arr[i]][0] = arr[i - 1];
        _link[arr[i]][1] = arr[i + 1];
    }

    _link[arr[0]][0] = arr[n - 1];
    _link[arr[0]][1] = arr[1];
    _link[arr[n - 1]][0] = arr[n - 2];
    _link[arr[n - 1]][1] = arr[0];

    ComputeCost(e);
}

void GA_EAX::Tour::FromFlipper (const Evaluator *e, const Flipper *f)
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

void GA_EAX::ABcycle::Iterator::Begin (const ABcycle *abc, bool reverse)
{
    _abc = abc;
    _len = _abc->GetLength();
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
    _capacity{2 * n}, _length{0}, _offset(0), _gain{0}
{
    _cyc = new int[_capacity];
}

int GA_EAX::ABcycle::Apply (bool reverse, Tour& pa) const
{
    Iterator iter;
    EdgeTriple et;
    for (iter.Begin(this, reverse); iter.End(et); iter++) {
        et.Reconnect(pa);
    }
    return GetGain();
}

GA_EAX::ABcycleMgr::ABcycleMgr (const Evaluator *e) :
    _eval(e), _numCity(e->GetNumCity()), _maxNumABcycle(2000)
{
    const int n = _numCity;

    _ABcycleList = new ABcycle*[_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _ABcycleList[j] = new ABcycle(n);
    }

    _cycRank2City = new int[n];
    _cycRank2Posi = new int[n];
    _cycRank1City = new int[n];
    _cycRank1Posi = new int[n];
    _cycCity = new int[2 * n];
    _cycPosi = new int[n];
    _quadEdge = new QuadEdge[n];
}

GA_EAX::ABcycleMgr::~ABcycleMgr ()
{
    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete _ABcycleList[j];
    }
    delete[] _ABcycleList;

    delete[] _cycRank2City;
    delete[] _cycRank2Posi;
    delete[] _cycRank1City;
    delete[] _cycRank1Posi;
    delete[] _cycCity;
    delete[] _cycPosi;
    delete[] _quadEdge;
}

void GA_EAX::ABcycleMgr::QuadEdge::Init (int a0, int a1, int b0, int b1)
{
    _edge[1] = a0;
    _edge[3] = a1;
    _edge[0] = b0;
    _edge[2] = b1;
    _remain = 2;
}

void GA_EAX::ABcycleMgr::Build (const Tour& pa, const Tour& pb, int numKid)
{
    _numABcycle = 0;
    _numCycRank2 = 0;
    _numCycRank1 = 0;
    for (int j = 0; j < _numCity; ++j) {
        _quadEdge[j].Init(pa.GetPrev(j), pa.GetNext(j), pb.GetPrev(j), pb.GetNext(j));
        _cycRank2City[_numCycRank2] = j;
        _cycRank2Posi[j] = _numCycRank2;
        _numCycRank2++;
        _cycPosi[j] = -1;
    }

    bool restart = true;
    int currIdx = 0;
    EdgeType select = ET_Rand;
    int start = -1, curr = -1, prev = -1;
    while (_numCycRank2 != 0) {
        if (restart) {
            currIdx = 0;
            int r = _eval->GetRand() % _numCycRank2;
            start = _cycRank2City[r];
            _cycCity[currIdx] = start;
            _cycPosi[start] = currIdx;
            curr = start;
            select = ET_Rand;
        } else {
            curr = _cycCity[currIdx];
        }

        do {
            currIdx++;
            prev = curr;
            if (select == ET_Rand) {
                EdgeType et = static_cast<EdgeType>(_eval->GetRand() % 2);
                curr = _quadEdge[prev].GetEdge(currIdx, et);
                if (et == ET_First) {
                    _quadEdge[prev].SwapEdge(currIdx);
                }
            } else {
                curr = _quadEdge[prev].GetEdge(currIdx, select);
            }
            _cycCity[currIdx] = curr;

            if (_quadEdge[curr].GetRemain() == 2) {
                if (curr == start) {
                    if (_cycPosi[start] == 0) {
                        if ((currIdx - _cycPosi[start]) % 2 == 0) {
                            if (_quadEdge[start].GetEdge(currIdx, ET_First) == prev) {
                                _quadEdge[curr].SwapEdge(currIdx);
                            }
                            if (!Build_0(1, numKid, currIdx)) {
                                goto END;
                            }
                            restart = false;
                            select = ET_First;
                            break;
                        } else {
                            _quadEdge[curr].SwapEdge(currIdx);
                            select = ET_Rand;
                        }
                        _cycPosi[start] = currIdx;
                    } else {
                        if (!Build_0(2, numKid, currIdx)) {
                            goto END;
                        }
                        restart = true;
                        break;
                    }
                } else if (_cycPosi[curr] == -1) {
                    _cycPosi[curr] = currIdx;
                    if (_quadEdge[curr].GetEdge(currIdx, ET_First) == prev) {
                        _quadEdge[curr].SwapEdge(currIdx);
                    }
                    select = ET_Rand;
                } else if (_cycPosi[curr] > 0) {
                    _quadEdge[curr].SwapEdge(currIdx);
                    if ((currIdx - _cycPosi[curr]) % 2 == 0) {
                        if (!Build_0(1, numKid, currIdx)) {
                            goto END;
                        }
                        restart = false;
                        select = ET_First;
                        break;
                    } else {
                        _quadEdge[curr].SwapEdge(currIdx + 1);
                        select = ET_Second;
                    }
                }
            } else if (_quadEdge[curr].GetRemain() == 1) {
                if (curr == start) {
                    if (!Build_0(1, numKid, currIdx)) {
                        goto END;
                    }
                    restart = true;
                    break;
                } else {
                    select = ET_First;
                }
            }
        } while (1);
    }

    while (_numCycRank1 != 0) {
        currIdx = 0;
        int r = _eval->GetRand() % _numCycRank1;
        start = _cycRank1City[r];
        _cycCity[currIdx] = start;
        curr = start;
        do {
            prev = curr;
            currIdx++;
            curr = _quadEdge[prev].GetEdge(currIdx, ET_First);
            _cycCity[currIdx] = curr;
            if (curr == start) {
                if (!Build_0(1, numKid, currIdx)) {
                    goto END;
                }
                break;
            }
        } while (1);
    }

END:
    std::shuffle(_ABcycleList, _ABcycleList + _numABcycle, _eval->GetRandEngine());
}

bool GA_EAX::ABcycleMgr::Build_0 (const int numStart, const int numKid, int& currIdx)
{
    const int start = _cycCity[currIdx];
    int len = 0, cntStart = 0, curr = start;
    ABcycle *abc = _ABcycleList[_numABcycle];

    do {
        assert(len < abc->GetCapacity());
        abc->SetCity(len++, curr);
        curr = _cycCity[--currIdx];

        if (_quadEdge[curr].GetRemain() == 2) {
            _cycRank2City[_cycRank2Posi[curr]] = _cycRank2City[_numCycRank2 - 1];
            _cycRank2Posi[_cycRank2City[_numCycRank2 - 1]] = _cycRank2Posi[curr];
            _numCycRank2--;
            _cycRank1City[_numCycRank1] = curr;
            _cycRank1Posi[curr] = _numCycRank1;
            _numCycRank1++;
        } else if (_quadEdge[curr].GetRemain() == 1) {
            _cycRank1City[_cycRank1Posi[curr]] = _cycRank1City[_numCycRank1 - 1];
            _cycRank1Posi[_cycRank1City[_numCycRank1 - 1]] = _cycRank1Posi[curr];
            _numCycRank1--;
        }
        _quadEdge[curr].DecrRemain();

        if (curr == start) {
            cntStart++;
        }
    } while (cntStart != numStart);

    if (len <= 2) {
        return true;
    }

    abc->SetLength(len);
    abc->SetOffset(0);
    if (currIdx % 2 != 0) {
        abc->SetOffset(1);
    }

    int gain = 0;
    for (int j = 0; j < len / 2; ++j) {
        int a = abc->GetCity(2 * j);
        int b = abc->GetCity(2 * j + 1);
        int c = abc->GetCity((2 * j + 2) % len);
        gain += _eval->GetCost(a, b) - _eval->GetCost(b, c);
    }
    abc->SetGain(gain);

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

GA_EAX::Cross::Cross (const Evaluator *e)
    : _eval(e), _numCity(e->GetNumCity()), _maxNumNear(10)
{
    assert(_maxNumNear <= _eval->GetMaxNumNear());

    const int n = _numCity;

    _city = new int[n];
    _posi = new int[n];
    _abcMgr = new ABcycleMgr(_eval);

    _posiInfo = new PosiInfo[n];
    _numEleInUnit = new int[n];
    _cuFlag = new int[n];
    _cuCity = new int[n];

    memset(_cuFlag, 0, sizeof(_cuFlag[0]) * n);
}

GA_EAX::Cross::~Cross ()
{
    delete[] _city;
    delete[] _posi;
    delete _abcMgr;

    delete[] _posiInfo;
    delete[] _numEleInUnit;
    delete[] _cuFlag;
    delete[] _cuCity;
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
    const int numAbc = std::min(numKid, _abcMgr->GetCycleNum());
    int bestIdx = -1, bestGain = 0;
    int aa, bb, a1, b1;
    for (int idx = 0; idx < numAbc; ++idx) {
        int gain = _abcMgr->GetCycle(idx)->Apply(false /*reverse*/, pa);

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
        _abcMgr->GetCycle(idx)->Apply(true /*reverse*/, pa);
        pa.SetGain(-gain);
    }

    if (bestIdx != -1) {
        /* go to best */
        _abcMgr->GetCycle(bestIdx)->Apply(false /*reverse*/, pa);
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

    const ABcycle *abc = _abcMgr->GetCycle(idx);
    const int n = _numCity;
    ABcycle::Iterator iter;
    EdgeTriple et;
    int b1, r1, r2, b2;
    int p1, p2, p0{0};
    bool hasZeroPosi = false;

    for (iter.Begin(abc, false /*reverse*/); iter.End(et); iter++) {
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

    for (int s = 0; s < numSeg; ++s) {
        _segments[s].unitId = -1;
        int p1 = _segments[s].begPosi;
        int p2 = _segments[s].endPosi;
        _posiInfo[p1].SetSegIdAndLinkPosiA(s, p2);
        _posiInfo[p2].SetSegIdAndLinkPosiA(s, p1);
    }
}

void GA_EAX::Cross::MakeUnit ()
{
    _numUnit = 0;

    int p_st, p1, p2, p_pre, p_next;
    while (1) {
        bool found = false;
        for (const auto& seg : _segments) {
            if (seg.unitId == -1) {
                p_st = seg.begPosi;
                p1 = p_st;
                p_pre = -1;
                found = true;
                break;
            }
        }
        if (!found) {
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

        int start = -1;
        for (const auto& seg : _segments) {
            if (seg.unitId == cu) {
                int posi = seg.begPosi;
                start = _city[posi];
            }
        }
        assert(start != -1);

        int numEleInCU = 0;
        Tour::Iterator iter(pa, start);
        do {
            iter++;
            int curr = iter.GetCurr();
            _cuFlag[curr] = 1;
            _cuCity[numEleInCU++] = curr;
            if (iter.GetNext() == iter.GetStart()) {
                break;
            }
        } while(1);
        assert(numEleInCU == _numEleInUnit[cu]);

        int numNear = _maxNumNear;
        int aa = -1, bb = -1, a1 = -1, b1 = -1;
        int maxDiff = std::numeric_limits<int>::min();

    RESTART:
        for (int s = 1; s <= numEleInCU; ++s) {
            int a = _cuCity[s % numEleInCU];
            for (int n = 0; n < numNear; ++n) {
                int c = _eval->GetNear(a, n);
                if (_cuFlag[c] == 0) {
                    for (int j1 = 0; j1 < 2; ++j1) {
                        int b = _cuCity[(s - 1 + 2 * j1) % numEleInCU];
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
            int a = _cuCity[r];
            int b = _cuCity[r + 1];
            for (int j = 0; j < _numCity; ++j) {
                if (_cuFlag[j] == 0) {
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
            int c = _cuCity[s];
            _cuFlag[c] = 0;
        }
    }

    return gain;
}

} /* namespace thu */
