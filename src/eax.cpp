#include "eax.h"

namespace thu { /* tsp heuristics */

GA_EAX::GA_EAX (const Evaluator* eval, int nPop, int nKid)
    : _eval(eval), _numPop(nPop), _numKid(nKid),
      _2opt(nullptr), _cross(nullptr), _pop(nullptr), _matingSeq(nullptr),
      _verbose(false), _numGen(0), _avgCost(0), _stagnateGen(0)
{
    const int n = eval->GetNumCity();

    _2opt = new TwoOpt(eval);
    _cross = new Cross(eval);
    _pop = new Indi[nPop+1];
    for (int i = 0; i < nPop+1; ++i) {
        _pop[i].Init(n);
    }
    _matingSeq = new int[nPop];

    for (int i = 0; i < nPop; ++i) {
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
            printf("=%d: %d %.3f\n", _numGen, GetBestIndi().GetCost(), _avgCost);
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

    const int stockBestCost = GetBestIndi().GetCost();
    int bestIndex = 0;

    for (int i = 0; i < _numPop; ++i) {
        _avgCost += _pop[i].GetCost();
        if (_pop[i].GetCost() < _pop[bestIndex].GetCost()) {
            bestIndex = i;
        }
    }

    _avgCost /= _numPop;
    GetBestIndi() = _pop[bestIndex];

    if (GetBestIndi().GetCost() < stockBestCost) {
        _stagnateGen = 0;
    } else {
        _stagnateGen++;
    }
}

bool GA_EAX::ShouldTerminate ()
{
    if (_avgCost - GetBestIndi().GetCost() < 0.001) {
        return true;
    }

    if (_stagnateGen > (1500 / _numKid)) {
        return true;
    }

    return false;
}

GA_EAX::Indi::Indi ()
    : _n(0), _link(nullptr), _cost(std::numeric_limits<int>::max())
{}

GA_EAX::Indi::~Indi ()
{
    for (int i = 0; i < _n; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void GA_EAX::Indi::Init (int n)
{
    _n = n;
    _link = new int*[n];
    for (int i = 0; i < n; ++i) {
        _link[i] = new int[2];
    }
}

GA_EAX::Indi& GA_EAX::Indi::operator= (const Indi& rhs)
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

void GA_EAX::Indi::ComputeCost (const Evaluator* e)
{
    _cost = 0;
    for (int i = 0; i < _n; ++i) {
        _cost += e->GetCost(i, _link[i][1]);
    }
}

void GA_EAX::Indi::FromArray (const Evaluator* e, const int* route)
{
    const int n = _n;
    for (int i = 1; i < n - 1; ++i) {
        _link[route[i]][0] = route[i - 1];
        _link[route[i]][1] = route[i + 1];
    }

    if (n <= 1) {
        fprintf(stderr, "ERROR: invalid city number (%d)\n", n);
        exit(1);
    }

    _link[route[0]][0] = route[n - 1];
    _link[route[0]][1] = route[1];
    _link[route[n - 1]][0] = route[n - 2];
    _link[route[n - 1]][1] = route[0];

    ComputeCost(e);
}

void GA_EAX::Indi::FromFlipper (const Evaluator* e, const Flipper* f)
{
    for (int i = 0; i < _n; ++i) {
        _link[i][0] = f->Prev(i);
        _link[i][1] = f->Next(i);
    }
    ComputeCost(e);
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

bool GA_EAX::ABcycle::Iterator::End (int& b1, int& r1, int& r2, int& b2) const
{
    if (_i < (_len / 2)) {
        int i = _start + _step * _i;
        b1 = _abc->GetCyc((i + _offset[0]) % _len);
        r1 = _abc->GetCyc((i + _offset[1]) % _len);
        r2 = _abc->GetCyc((i + _offset[2]) % _len);
        b2 = _abc->GetCyc((i + _offset[3]) % _len);
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
    const int n = _numCity;

    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete _ABcycleList[j];
    }
    delete[] _ABcycleList;

    for (int j = 0; j < n; ++j) {
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

void GA_EAX::ABcyleMgr::Build (const Indi& pa1, const Indi& pa2, int nKid)
{
    const int n = _numCity;

    _numABcycle = 0;

    _cycBuf1Num = 0;
    _cycBuf2Num = 0;

    for (int j = 0; j < n; ++j) {
        _overlapEdges[j][1] = pa1.GetPrev(j);
        _overlapEdges[j][3] = pa1.GetNext(j);
        _overlapEdges[j][0] = 2;
        _overlapEdges[j][2] = pa2.GetPrev(j);
        _overlapEdges[j][4] = pa2.GetNext(j);

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
                            if (!Build_0(1, nKid, posiCurr)) {
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
                        if (!Build_0(2, nKid, posiCurr)) {
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
                        if (!Build_0(1, nKid, posiCurr)) {
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
                    if (!Build_0(1, nKid, posiCurr)) {
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
                if (!Build_0(1, nKid, posiCurr)) {
                    goto LLL;
                }

                flagCircle = 1;
            }
        }
    }

LLL:
    std::shuffle(_ABcycleList, _ABcycleList + _numABcycle, _eval->GetRandEngine());
}

bool GA_EAX::ABcyleMgr::Build_0 (const int stAppear, const int nKid, int& posiCurr)
{
    ABcycle* abc = _ABcycleList[_numABcycle];
    const int st = _cycRoute[posiCurr];
    int len = 0;
    int st_count = 0;

    abc->SetCyc(len++, st);
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

        if (len >= abc->GetMaxLen()) {
            fprintf(stderr, "ERROR: AB-cycle length reach max (%d)\n", abc->GetMaxLen());
            exit(1);
        }
        abc->SetCyc(len++, ci);
    }

    if (len <= 2) {
        return true;
    }

    abc->SetLen(len);
    if (posiCurr % 2 != 0) {
        int stock = abc->GetCyc(0);
        for (int j = 0; j < len - 1; j++) {
            abc->SetCyc(j, abc->GetCyc(j + 1));
        }
        abc->SetCyc(len - 1, stock);
    }

    int diff = 0;
    for (int j = 0; j < len / 2; ++j) {
        int a = abc->GetCyc(2 * j);
        int b = abc->GetCyc((2 * j + 1) % len);
        int c = abc->GetCyc((2 * j + 2) % len);
        diff += _eval->GetCost(a, b) - _eval->GetCost(b, c);
    }
    abc->SetGain(diff);

    _numABcycle++;

    if (_numABcycle >= _maxNumABcycle) {
        fprintf(stderr, "WARNING: maximum number of AB-cycles (%d) must be increased\n",
                _maxNumABcycle);
        return false;
    }

    if (_numABcycle >= nKid) {
        return false;
    }

    return true;
}

void GA_EAX::ABcyleMgr::ChangeIndi (int idx, bool reverse, Indi& pa1) const
{
    const ABcycle* abc = _ABcycleList[idx];
    ABcycle::Iterator iter;
    int b1, r1, r2, b2;
    for (iter.Begin(abc, reverse); iter.End(b1, r1, r2, b2); iter++) {
        BreakEdge(b1, r1, r2, b2, pa1);
    }
}

void GA_EAX::ABcyleMgr::BackToPa1 (int idx, const EdgeTuples& me, Indi& pa1) const
{
    for (auto it = me.rbegin(); it != me.rend(); ++it) {
        int aa = std::get<0>(*it);
        int a1 = std::get<1>(*it);
        int bb = std::get<2>(*it);
        int b1 = std::get<3>(*it);

    #if 1
        BreakEdge(a1, aa, bb, b1, pa1);
        BreakEdge(bb, b1, a1, aa, pa1);
    #else
        if (pa1.GetPrev(aa) == bb) {
            pa1.SetPrev(aa, a1);
        } else {
            pa1.SetNext(aa, a1);
        }
        if (pa1.GetPrev(b1) == a1) {
            pa1.SetPrev(b1, bb);
        } else {
            pa1.SetNext(b1, bb);
        }
        if (pa1.GetPrev(bb) == aa) {
            pa1.SetPrev(bb, b1);
        } else {
            pa1.SetNext(bb, b1);
        }
        if (pa1.GetPrev(a1) == b1) {
            pa1.SetPrev(a1, aa);
        } else {
            pa1.SetNext(a1, aa);
        }
    #endif
    }

    ChangeIndi(idx, true /*reverse*/, pa1);
}

void GA_EAX::ABcyleMgr::GoToBest (int idx, const EdgeTuples& me, Indi& pa1) const
{
    ChangeIndi(idx, false /*reverse*/, pa1);

    for (auto it = me.begin(); it != me.end(); ++it) {
        int aa = std::get<0>(*it);
        int bb = std::get<1>(*it);
        int a1 = std::get<2>(*it);
        int b1 = std::get<3>(*it);

    #if 1
        BreakEdge(a1, aa, bb, b1, pa1);
        BreakEdge(aa, a1, b1, bb, pa1);
    #else
        if (pa1.GetPrev(aa) == bb) {
            pa1.SetPrev(aa, a1);
        } else {
            pa1.SetNext(aa, a1);
        }
        if (pa1.GetPrev(bb) == aa) {
            pa1.SetPrev(bb, b1);
        } else {
            pa1.SetNext(bb, b1);
        }
        if (pa1.GetPrev(a1) == b1) {
            pa1.SetPrev(a1, aa);
        } else {
            pa1.SetNext(a1, aa);
        }
        if (pa1.GetPrev(b1) == a1) {
            pa1.SetPrev(b1, bb);
        } else {
            pa1.SetNext(b1, bb);
        }
    #endif
    }
}

/* r1--r2 will be broken, b1--r1 and r2--b2 will be connected. */
void GA_EAX::ABcyleMgr::BreakEdge (int b1, int r1, int r2, int b2, Indi& pa1)
{
    if (pa1.GetPrev(r1) == r2) {
        pa1.SetPrev(r1, b1);
    } else {
        pa1.SetNext(r1, b1);
    }
    if (pa1.GetPrev(r2) == r1) {
        pa1.SetPrev(r2, b2);
    } else {
        pa1.SetNext(r2, b2);
    }
}

GA_EAX::Cross::Cross (const Evaluator* e)
    : _eval(e), _numCity(e->GetNumCity())
{
    const int n = _numCity;

    _pa1City = new int[n];
    _pa1Pos = new int[n];
    _abcMgr = new ABcyleMgr(_eval);
    _modiEdges.resize(n);
    _bestModiEdges.resize(n);

    _segment = new int*[n];
    for (int j = 0; j < n; ++j) {
        _segment[j] = new int[2];
    }
    _segUnit = new int[n];
    _segPosiList = new int[n];
    _linkAPosi = new int[n];
    _linkBPosi = new int*[n];
    for (int j = 0; j < n; ++j) {
        _linkBPosi[j] = new int[2];
    }
    _posiSeg = new int[n];
    _numElementInUnit = new int[n];
    _centerUnit = new int[n];
    for (int j = 0; j < n; ++j) {
        _centerUnit[j] = 0;
    }
    _listCenterUnit = new int[n + 2];
}

GA_EAX::Cross::~Cross ()
{
    const int n = _numCity;

    delete[] _pa1City;
    delete[] _pa1Pos;
    delete _abcMgr;

    for (int j = 0; j < n; ++j) {
        delete[] _segment[j];
    }
    delete[] _segment;
    delete[] _segUnit;
    delete[] _segPosiList;
    delete[] _linkAPosi;
    for (int j = 0; j < n; ++j) {
        delete[] _linkBPosi[j];
    }
    delete[] _linkBPosi;
    delete[] _posiSeg;
    delete[] _numElementInUnit;
    delete[] _centerUnit;
    delete[] _listCenterUnit;
}

void GA_EAX::Cross::DoIt (Indi& pa1, Indi& pa2, int nKid)
{
    InitPa1CityPos(pa1);
    _abcMgr->Build(pa1, pa2, nKid);

    /* main loop to generate nKid kids */
    int bestGain = 0, bestAbcIdx = -1, gain;
    nKid = std::min(nKid, _abcMgr->GetNumABCycle());
    for (int j = 0; j < nKid; ++j) {
        gain = 0;
        _numSPL = 0;

        _abcMgr->ChangeIndi(j, false /*reverse*/, pa1);
        UpdateSeg(j);
        gain += _abcMgr->GetABCycle(j)->GetGain();

        MakeUnit();
        gain += MakeCompleteSol(pa1);
        pa1.SetGain(gain);

        if (bestGain < gain) {
            bestAbcIdx = j;
            bestGain = gain;
            _bestModiEdges.clear();
            for (const auto& me : _modiEdges) {
                _bestModiEdges.push_back(me);
            }
        }

        _abcMgr->BackToPa1(j, _modiEdges, pa1);
        pa1.SetGain(-gain);
    }

    if (bestAbcIdx != -1) {
        _abcMgr->GoToBest(bestAbcIdx, _bestModiEdges, pa1);
        pa1.SetGain(bestGain);
    }
}

void GA_EAX::Cross::InitPa1CityPos (const Indi& pa1) const
{
    const int n = _numCity;
    int prev = -1, curr = -1, next = 0;
    for (int i = 0; i < n; ++i) {
        prev = curr;
        curr = next;
        if (pa1.GetPrev(curr) != prev) {
            next = pa1.GetPrev(curr);
        } else {
            next = pa1.GetNext(curr);
        }
        _pa1City[i] = curr;
        _pa1Pos[curr] = i;
    }
}

void GA_EAX::Cross::UpdateSeg (int idx)
{
    const ABcycle* abc = _abcMgr->GetABCycle(idx);
    const int n = _numCity;
    ABcycle::Iterator iter;
    int b1, r1, r2, b2;

    for (iter.Begin(abc, false /* reverse */); iter.End(b1, r1, r2, b2); iter++) {
        if (_numSPL >= n) {
            fprintf(stderr, "ERROR: numSPL reach max (%d) in UpdateSeg\n", n);
            exit(1);
        }

        if (_pa1Pos[r1] == 0 && _pa1Pos[r2] == n - 1) {
            _segPosiList[_numSPL++] = _pa1Pos[r1];
        } else if (_pa1Pos[r1] == n - 1 && _pa1Pos[r2] == 0) {
            _segPosiList[_numSPL++] = _pa1Pos[r2];
        } else if (_pa1Pos[r1] < _pa1Pos[r2]) {
            _segPosiList[_numSPL++] = _pa1Pos[r2];
        } else if (_pa1Pos[r2] < _pa1Pos[r1]) {
            _segPosiList[_numSPL++] = _pa1Pos[r1];
        } else {
            fprintf(stderr, "ERROR: invalid else branch in UpdateSeg\n");
            exit(1);
        }

        _linkBPosi[_pa1Pos[r1]][1] = _linkBPosi[_pa1Pos[r1]][0];
        _linkBPosi[_pa1Pos[r2]][1] = _linkBPosi[_pa1Pos[r2]][0];
        _linkBPosi[_pa1Pos[r1]][0] = _pa1Pos[b1];
        _linkBPosi[_pa1Pos[r2]][0] = _pa1Pos[b2];
    }
}

void GA_EAX::Cross::MakeUnit ()
{
    const int n = _numCity;
    int flag = 1;
    for (int s = 0; s < _numSPL; ++s) {
        if (_segPosiList[s] == 0) {
            flag = 0;
            break;
        }
    }
    if (flag == 1) {
        if (_numSPL >= n) {
            fprintf(stderr, "ERROR: numSPL reach max (%d) in MakeUnit", n);
            exit(1);
        }
        _segPosiList[_numSPL++] = 0;

        _linkBPosi[n - 1][1] = _linkBPosi[n - 1][0];
        _linkBPosi[0][1] = _linkBPosi[0][0];
        _linkBPosi[n - 1][0] = 0;
        _linkBPosi[0][0] = n - 1;
    }

    std::sort(_segPosiList, _segPosiList + _numSPL);

    _numSeg = _numSPL;
    for (int s = 0; s < _numSeg - 1; ++s) {
        _segment[s][0] = _segPosiList[s];
        _segment[s][1] = _segPosiList[s + 1] - 1;
    }

    _segment[_numSeg - 1][0] = _segPosiList[_numSeg - 1];
    _segment[_numSeg - 1][1] = n - 1;

    for (int s = 0; s < _numSeg; ++s) {
        _linkAPosi[_segment[s][0]] = _segment[s][1];
        _linkAPosi[_segment[s][1]] = _segment[s][0];
        _posiSeg[_segment[s][0]] = s;
        _posiSeg[_segment[s][1]] = s;
    }

    for (int s = 0; s < _numSeg; ++s) {
        _segUnit[s] = -1;
    }
    _numUnit = 0;

    int p_st, p1, p2, p_next, p_pre;
    int segNum;

    while (1) {
        flag = 0;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == -1) {
                p_st = _segment[s][0];
                p_pre = -1;
                p1 = p_st;
                flag = 1;
                break;
            }
        }
        if (flag == 0) {
            break;
        }

        while (1) {
            segNum = _posiSeg[p1];
            _segUnit[segNum] = _numUnit;

            p2 = _linkAPosi[p1];
            p_next = _linkBPosi[p2][0];
            if (p1 == p2) {
                if (p_next == p_pre) {
                    p_next = _linkBPosi[p2][1];
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

    for (int s = 0; s < _numUnit; ++s) {
        _numElementInUnit[s] = 0;
    }

    int unitNum = -1;
    int tmpNumSeg = -1;
    for (int s = 0; s < _numSeg; ++s) {
        if (_segUnit[s] != unitNum) {
            ++tmpNumSeg;
            _segment[tmpNumSeg][0] = _segment[s][0];
            _segment[tmpNumSeg][1] = _segment[s][1];
            unitNum = _segUnit[s];
            _segUnit[tmpNumSeg] = unitNum;
            _numElementInUnit[unitNum] += _segment[s][1] - _segment[s][0] + 1;
        } else {
            _segment[tmpNumSeg][1] = _segment[s][1];
            _numElementInUnit[unitNum] += _segment[s][1] - _segment[s][0] + 1;
        }
    }
    _numSeg = tmpNumSeg + 1;
}

int GA_EAX::Cross::MakeCompleteSol (Indi& pa1)
{
    _modiEdges.clear();

    int gainModi = 0;
    constexpr int NearMaxDef = 10;
    int center_un = 0;
    int numEleInCU = 0;

    int st, prev, curr, next, a, b, c, d, aa = 0, bb = 0, a1 = 0, b1 = 0;

    while (_numUnit != 1) {
        int min_unit_city = _numCity + 12345;
        for (int u = 0; u < _numUnit; ++u) {
            if (_numElementInUnit[u] < min_unit_city) {
                center_un = u;
                min_unit_city = _numElementInUnit[u];
            }
        }

        st = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == center_un) {
                int posi = _segment[s][0];
                st = _pa1City[posi];
            }
        }
        if (st == -1) {
            fprintf(stderr, "ERROR: invalid st\n");
            exit(1);
        }

        curr = -1;
        next = st;
        numEleInCU = 0;
        while (1) {
            prev = curr;
            curr = next;
            _centerUnit[curr] = 1;
            _listCenterUnit[numEleInCU] = curr;
            ++numEleInCU;

            if (pa1.GetPrev(curr) != prev) {
                next = pa1.GetPrev(curr);
            } else {
                next = pa1.GetNext(curr);
            }

            if (next == st) {
                break;
            }
        }
        _listCenterUnit[numEleInCU] = _listCenterUnit[0];
        _listCenterUnit[numEleInCU + 1] = _listCenterUnit[1];

        if (numEleInCU != _numElementInUnit[center_un]) {
            fprintf(stderr, "ERROR: invalid numEleInCU (%d)\n", numEleInCU);
            exit(1);
        }

        int max_diff = std::numeric_limits<int>::min();
        int diff = 0;
        a1 = -1;
        b1 = -1;
        /* N_near (see Step 5.3 in Section 2.2 of the Online Supplement)
         * nearMax must be smaller than or equal to eva->_maxNumNear (Evaluator) */
        if (NearMaxDef > _eval->GetMaxNumNear()) {
            fprintf(stderr, "ERROR: invalid NearMaxDef (%d)\n", NearMaxDef);
            exit(1);
        }
        int nearMax = NearMaxDef;

    RESTART:;
        for (int s = 1; s <= numEleInCU; ++s) {
            a = _listCenterUnit[s];
            for (int near_num = 0; near_num < nearMax; ++near_num) {
                c = _eval->GetNear(a, near_num);
                if (_centerUnit[c] == 0) {
                    for (int j1 = 0; j1 < 2; ++j1) {
                        b = _listCenterUnit[s - 1 + 2 * j1];
                        for (int j2 = 0; j2 < 2; ++j2) {
                            d = (j2 == 0)? pa1.GetPrev(c) : pa1.GetNext(c);
                            diff = _eval->GetCost(a, b) + _eval->GetCost(c, d)
                                   - _eval->GetCost(a, c) - _eval->GetCost(b, d);
                            if (diff > max_diff) {
                                aa = a;
                                bb = b;
                                a1 = c;
                                b1 = d;
                                max_diff = diff;
                            }
                            diff = _eval->GetCost(a, b) + _eval->GetCost(d, c)
                                   - _eval->GetCost(a, d) - _eval->GetCost(b, c);
                            if (diff > max_diff) {
                                aa = a;
                                bb = b;
                                a1 = d;
                                b1 = c;
                                max_diff = diff;
                            }
                        }
                    }
                }
            }
        }

        /* This value must also be changed if nearMax is chenged above */
        if (a1 == -1 && nearMax == NearMaxDef) {
            nearMax = _eval->GetMaxNumNear();
            goto RESTART;
        } else if (a1 == -1 && nearMax == _eval->GetMaxNumNear()) {
            int r = _eval->GetRand() % (numEleInCU - 1);
            a = _listCenterUnit[r];
            b = _listCenterUnit[r + 1];
            for (int j = 0; j < _numCity; ++j) {
                if (_centerUnit[j] == 0) {
                    aa = a;
                    bb = b;
                    a1 = j;
                    b1 = pa1.GetPrev(j);
                    break;
                }
            }
            max_diff = _eval->GetCost(aa, bb) + _eval->GetCost(a1, b1)
                       - _eval->GetCost(a, a1) - _eval->GetCost(b, b1);
        }

    #if 1
        ABcyleMgr::BreakEdge(a1, aa, bb, b1, pa1);
        ABcyleMgr::BreakEdge(aa, a1, b1, bb, pa1);
    #else
        if (pa1.GetPrev(aa) == bb) {
            pa1.SetPrev(aa, a1);
        } else {
            pa1.SetNext(aa, a1);
        }
        if (pa1.GetPrev(bb) == aa) {
            pa1.SetPrev(bb, b1);
        } else {
            pa1.SetNext(bb, b1);
        }
        if (pa1.GetPrev(a1) == b1) {
            pa1.SetPrev(a1, aa);
        } else {
            pa1.SetNext(a1, aa);
        }
        if (pa1.GetPrev(b1) == a1) {
            pa1.SetPrev(b1, bb);
        } else {
            pa1.SetNext(b1, bb);
        }
    #endif

        _modiEdges.push_back(std::make_tuple(aa, bb, a1, b1));
        gainModi += max_diff;

        int posi_a1 = _pa1Pos[a1];
        int select_un = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segment[s][0] <= posi_a1 && posi_a1 <= _segment[s][1]) {
                select_un = _segUnit[s];
                break;
            }
        }
        if (select_un == -1) {
            fprintf(stderr, "ERROR: invalid select_un\n");
            exit(1);
        }

        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == select_un) {
                _segUnit[s] = center_un;
            }
        }
        _numElementInUnit[center_un] += _numElementInUnit[select_un];

        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == _numUnit - 1) {
                _segUnit[s] = select_un;
            }
        }
        _numElementInUnit[select_un] = _numElementInUnit[_numUnit - 1];
        --_numUnit;

        for (int s = 0; s < numEleInCU; ++s) {
            c = _listCenterUnit[s];
            _centerUnit[c] = 0;
        }
    }

    return gainModi;
}

} /* namespace thu */
