#include "eax.h"

namespace thu { /* tsp heuristics */

Indi::Indi ()
    : _numCity(0), _link(nullptr), _cost(std::numeric_limits<int>::max())
{
}

Indi::~Indi ()
{
    for (int i = 0; i < _numCity; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void Indi::Define (int n)
{
    _numCity = n;
    _link = new int*[n];
    for (int i = 0; i < n; ++i) {
        _link[i] = new int[2];
    }
}

Indi& Indi::operator= (const Indi& rhs)
{
    if (this != &rhs) {
        _numCity = rhs._numCity;
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < 2; ++j) {
                _link[i][j] = rhs._link[i][j];
            }
        }
        _cost = rhs._cost;
    }
    return *this;
}

void Indi::ToArray (int* arr, int* arrInv) const
{
    int curr = -1, next = 0, prev = -1;
    for (int i = 0; i < _numCity; ++i) {
        prev = curr;
        curr = next;
        if (_link[curr][0] != prev) {
            next = _link[curr][0];
        } else {
            next = _link[curr][1];
        }
        arr[i] = curr;
        if (arrInv) {
            arrInv[curr] = i;
        }
    }
}

void Indi::DoneBy (const Evaluator* e)
{
    _cost = 0;
    for (int i = 0; i < e->GetNumCity(); ++i) {
        _cost += e->GetCost(i, _link[i][1]);
    }
}

void Indi::FromArray (const Evaluator* e, const int* route)
{
    const int n = _numCity;
    for (int i = 1; i < n - 1; ++i) {
        _link[route[i]][0] = route[i - 1];
        _link[route[i]][1] = route[i + 1];
    }
    _link[route[0]][0] = route[n - 1];
    _link[route[0]][1] = route[1];
    _link[route[n - 1]][0] = route[n - 2];
    _link[route[n - 1]][1] = route[0];

    DoneBy(e);
}

void Indi::MadeRand (const Evaluator* e)
{
    const int* route = e->MakeRand();
    FromArray(e, route);
}

void Indi::FromFlipper (const Evaluator* e, const Flipper* f)
{
    for (int i = 0; i < _numCity; ++i) {
        _link[i][0] = f->Prev(i);
        _link[i][1] = f->Next(i);
    }
    DoneBy(e);
}


EAXGA::EAXGA (const Evaluator* eval, int nPop, int nKid)
    : _eval(eval), _matingSeq(new int[nPop + 1]),
      _opt2(nullptr), _cross(nullptr), _pop(nullptr),
      _numPop(nPop), _numKids(nKid), _verbose(false),
      _numGen(0), _avgCost(0), _stagnGen(0)
{
    Init();
}

EAXGA::~EAXGA ()
{
    delete[] _matingSeq;
    delete _opt2;
    delete _cross;
    delete[] _pop;
}

void EAXGA::Init ()
{
    const int n = _eval->GetNumCity();
    if (n <= 0) {
        fprintf(stderr, "ERROR: failed to init EAXGA, #cities = %d\n", n);
        return;
    }

    _best.Define(n);
    _opt2 = new TwoOpt(_eval);
    _cross = new Cross(_eval, _numPop);
    _pop = new Indi[_numPop];
    for (int i = 0; i < _numPop; ++i) {
        _pop[i].Define(n);
    }
}

void EAXGA::DoIt ()
{
    _numGen = 0;
    _stagnGen = 0;
    _best.MadeRand(_eval);

    for (int i = 0; i < _numPop; ++i) {
        _opt2->DoIt();
        _pop[i].FromFlipper(_eval, _opt2->GetFlipper());
    }

    while (1) {
        SelectBest();
        if (_verbose) {
            printf("%d: %lld %lf\n", _numGen, (long long)_best._cost, _avgCost);
        }

        if (ShouldTerminate()) {
            break;
        }

        SelectForMating();
        for (int s = 0; s < _numPop; ++s) {
            _cross->DoIt(_pop[_matingSeq[s]], _pop[_matingSeq[s + 1]],
                         _numKids);
        }

        ++_numGen;
    }
}

void EAXGA::SelectBest ()
{
    int stockBest = _best._cost;

    _avgCost = 0.0;
    int bestIndex = 0;
    int bestCost = _pop[0]._cost;

    for (int i = 0; i < _numPop; ++i) {
        _avgCost += _pop[i]._cost;
        if (_pop[i]._cost < bestCost) {
            bestIndex = i;
            bestCost = _pop[i]._cost;
        }
    }

    _best = _pop[bestIndex];
    _avgCost /= _numPop;

    if (_best._cost < stockBest) {
        _stagnGen = 0;
    } else {
        _stagnGen++;
    }
}

bool EAXGA::ShouldTerminate ()
{
    if (_avgCost - _best._cost < 0.001) {
        return true;
    }

    if (_stagnGen > (1500 / _numKids)) {
        return true;
    }

    return false;
}

void EAXGA::SelectForMating ()
{
    for (int i = 0; i < _numPop; i++) {
        _matingSeq[i] = i;
    }
    std::shuffle(_matingSeq, _matingSeq + _numPop, _eval->GetRandEnginge());
    _matingSeq[_numPop] = _matingSeq[0];
}

Cross::Cross (const Evaluator* e, int nPop)
    : _eval(e), _numCity(e->GetNumCity()), _numPop(nPop), _maxNumABcycle(2000)
{
    const int n = _numCity;
    _ABcycleList = new int*[_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _ABcycleList[j] = new int[2 * n + 4];
    }
    _permuABCycle = new int[_maxNumABcycle];
    _gainABcycle = new int[_maxNumABcycle];

    _pa1Route = new int[n];
    _pa1RouteInv = new int[n];

    _overlapEdges = new int*[n];
    for (int j = 0; j < n; ++j) {
        _overlapEdges[j] = new int[5];
    }
    _cycBuf1 = new int[n];
    _cycBuf2 = new int[n];
    _cycBuf1Inv = new int[n];
    _cycBuf2Inv = new int[n];
    _cycRoute = new int[2 * n + 1];
    _ABCycle = new int[2 * n + 4];

    _modiEdge = new int*[n];
    for (int j = 0; j < n; ++j) {
        _modiEdge[j] = new int[4];
    }
    _bestModiEdge = new int*[n];
    for (int j = 0; j < n; ++j) {
        _bestModiEdge[j] = new int[4];
    }

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

    _routeBuf = new int[n];
}

Cross::~Cross ()
{
    const int n = _numCity;
    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete[] _ABcycleList[j];
    }
    delete[] _ABcycleList;
    delete[] _permuABCycle;
    delete[] _gainABcycle;

    for (int j = 0; j < n; ++j) {
        delete[] _overlapEdges[j];
    }
    delete[] _overlapEdges;
    delete[] _cycBuf1;
    delete[] _cycBuf2;
    delete[] _cycBuf1Inv;
    delete[] _cycBuf2Inv;
    delete[] _cycRoute;
    delete[] _ABCycle;

    delete[] _pa1Route;
    delete[] _pa1RouteInv;

    for (int j = 0; j < n; ++j) {
        delete[] _modiEdge[j];
    }
    delete[] _modiEdge;
    for (int j = 0; j < n; ++j) {
        delete[] _bestModiEdge[j];
    }
    delete[] _bestModiEdge;

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

    delete[] _routeBuf;
}

void Cross::DoIt (Indi& kid, Indi& pa2, int nKids)
{
    int bestGain = 0, gain;
    int bestAppliedCycle, appliedCycle;

    /* init _pa1Route and _pa1RouteInv for pa1 */
    kid.ToArray(_pa1Route, _pa1RouteInv);

    BuildABcycle(kid, pa2, nKids);

    if (nKids > _numABcycle) {
        nKids = _numABcycle;
    }
    for (int i = 0; i < _numABcycle; i++) {
        _permuABCycle[i] = i;
    }
    std::shuffle(_permuABCycle, _permuABCycle + _numABcycle, _eval->GetRandEnginge());

    /* main loop to generate nKids kids */
    for (int j = 0; j < nKids; ++j) {
        gain = 0;
        _numModiEdge = 0;
        _numSPL = 0;

        appliedCycle = _permuABCycle[j];
        ChangeSol(kid, appliedCycle, false /*reverse*/);
        gain += _gainABcycle[appliedCycle];

        MakeUnit();
        gain += MakeCompleteSol(kid);
        kid._cost = kid._cost - gain;

        if (bestGain < gain) {
            bestGain = gain;
            bestAppliedCycle = appliedCycle;
            _numBestModiEdge = _numModiEdge;
            for (int s = 0; s < _numBestModiEdge; ++s) {
                memcpy(_bestModiEdge[s], _modiEdge[s], sizeof(int) * 4);
            }
        }

        BackToPa1(kid, appliedCycle);
        kid._cost = kid._cost + gain;
    }

    if (bestGain != 0) {
        GoToBest(kid, bestAppliedCycle);
        kid._cost = kid._cost - bestGain;
    }
}

void Cross::BuildABcycle (const Indi& pa1, const Indi& pa2, int nKids)
{
    const int n = _numCity;
    int* checkCycBuf1 = _routeBuf;
    _cycBuf2Num = 0;
    _cycBuf1Num = 0;

    for (int j = 0; j < n; ++j) {
        _overlapEdges[j][1] = pa1._link[j][0];
        _overlapEdges[j][3] = pa1._link[j][1];
        _overlapEdges[j][0] = 2;
        _overlapEdges[j][2] = pa2._link[j][0];
        _overlapEdges[j][4] = pa2._link[j][1];

        _cycBuf1[_cycBuf1Num++] = j;
        _cycBuf1Inv[_cycBuf1[j]] = j;
        checkCycBuf1[j] = -1;
    }

    /**************************************************/

    _numABcycle = 0;
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
            checkCycBuf1[st] = posiCurr;
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
                    if (checkCycBuf1[st] == 0) {
                        if ((posiCurr - checkCycBuf1[st]) % 2 == 0) {
                            if (_overlapEdges[st][posiCurr % 2 + 1] == pr) {
                                std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                                          _overlapEdges[ci][posiCurr % 2 + 3]);
                            }
                            BuildABcycle_0(1, posiCurr);
                            if (_numABcycle == nKids) {
                                goto LLL;
                            }
                            if (_numABcycle == _maxNumABcycle) {
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
                        checkCycBuf1[st] = posiCurr;
                    } else {
                        BuildABcycle_0(2, posiCurr);
                        if (_numABcycle == nKids) {
                            goto LLL;
                        }
                        if (_numABcycle == _maxNumABcycle) {
                            goto LLL;
                        }

                        flagSt = 1;
                        flagCircle = 1;
                    }
                } else if (checkCycBuf1[ci] == -1) {
                    checkCycBuf1[ci] = posiCurr;
                    if (_overlapEdges[ci][posiCurr % 2 + 1] == pr) {
                        std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                                  _overlapEdges[ci][posiCurr % 2 + 3]);
                    }
                    prType = 2;
                } else if (checkCycBuf1[ci] > 0) {
                    std::swap(_overlapEdges[ci][posiCurr % 2 + 1],
                              _overlapEdges[ci][posiCurr % 2 + 3]);
                    if ((posiCurr - checkCycBuf1[ci]) % 2 == 0) {
                        BuildABcycle_0(1, posiCurr);
                        if (_numABcycle == nKids) {
                            goto LLL;
                        }
                        if (_numABcycle == _maxNumABcycle) {
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
                    BuildABcycle_0(1, posiCurr);
                    if (_numABcycle == nKids) {
                        goto LLL;
                    }
                    if (_numABcycle == _maxNumABcycle) {
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
                BuildABcycle_0(1, posiCurr);
                if (_numABcycle == nKids) {
                    goto LLL;
                }
                if (_numABcycle == _maxNumABcycle) {
                    goto LLL;
                }
                flagCircle = 1;
            }
        }
    }

LLL:;

    if (_numABcycle == _maxNumABcycle) {
        fprintf(stderr, "WARNING: _maxNumABcycle(%d) must be increased\n",
               _maxNumABcycle);
    }
}

void Cross::BuildABcycle_0 (int stAppear, int& posiCurr)
{
    const int st = _cycRoute[posiCurr];
    int st_count = 0;
    int cem = 0;
    _ABCycle[cem] = st;

    while (1) {
        cem++;
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
        _ABCycle[cem] = ci;
    }

    if (cem == 2) {
        return;
    }

    _ABcycleList[_numABcycle][0] = cem;

    if (posiCurr % 2 != 0) {
        int stock = _ABCycle[0];
        for (int j = 0; j < cem - 1; j++) {
            _ABCycle[j] = _ABCycle[j + 1];
        }
        _ABCycle[cem - 1] = stock;
    }

    for (int j = 0; j < cem; j++) {
        _ABcycleList[_numABcycle][j + 2] = _ABCycle[j];
    }
    _ABcycleList[_numABcycle][1] = _ABCycle[cem - 1];
    _ABcycleList[_numABcycle][cem + 2] = _ABCycle[0];
    _ABcycleList[_numABcycle][cem + 3] = _ABCycle[1];

    _ABCycle[cem] = _ABCycle[0];
    _ABCycle[cem + 1] = _ABCycle[1];
    int diff = 0;
    for (int j = 0; j < cem / 2; ++j) {
        diff += _eval->GetCost(_ABCycle[2 * j], _ABCycle[1 + 2 * j])
                - _eval->GetCost(_ABCycle[1 + 2 * j], _ABCycle[2 + 2 * j]);
    }
    _gainABcycle[_numABcycle] = diff;
    ++_numABcycle;
}

void Cross::ChangeSol (Indi& kid, int idx, bool reverse, bool updateSeg)
{
    const int n = _numCity;
    int cem, r1, r2, b1, b2;

    cem = _ABcycleList[idx][0];
    _ABCycle[0] = _ABcycleList[idx][0];

    if (reverse) {
        for (int j = 0; j < cem + 3; j++) {
            _ABCycle[cem + 3 - j] = _ABcycleList[idx][j + 1];
        }
    } else {
        for (int j = 1; j <= cem + 3; j++) {
            _ABCycle[j] = _ABcycleList[idx][j];
        }
    }

    for (int j = 0; j < cem / 2; j++) {
        r1 = _ABCycle[2 + 2 * j];
        r2 = _ABCycle[3 + 2 * j];
        b1 = _ABCycle[1 + 2 * j];
        b2 = _ABCycle[4 + 2 * j];

        if (kid._link[r1][0] == r2) {
            kid._link[r1][0] = b1;
        } else {
            kid._link[r1][1] = b1;
        }
        if (kid._link[r2][0] == r1) {
            kid._link[r2][0] = b2;
        } else {
            kid._link[r2][1] = b2;
        }

        if (updateSeg) {
            if (_numSPL >= n) {
                fprintf(stderr, "ERROR: numSPL reach max(%d) in ChangeSol\n", n);
                exit(1);
            }
            if (_pa1RouteInv[r1] == 0 && _pa1RouteInv[r2] == n - 1) {
                _segPosiList[_numSPL++] = _pa1RouteInv[r1];
            } else if (_pa1RouteInv[r1] == n - 1 && _pa1RouteInv[r2] == 0) {
                _segPosiList[_numSPL++] = _pa1RouteInv[r2];
            } else if (_pa1RouteInv[r1] < _pa1RouteInv[r2]) {
                _segPosiList[_numSPL++] = _pa1RouteInv[r2];
            } else if (_pa1RouteInv[r2] < _pa1RouteInv[r1]) {
                _segPosiList[_numSPL++] = _pa1RouteInv[r1];
            } else {
                assert(0);
            }

            _linkBPosi[_pa1RouteInv[r1]][1] = _linkBPosi[_pa1RouteInv[r1]][0];
            _linkBPosi[_pa1RouteInv[r2]][1] = _linkBPosi[_pa1RouteInv[r2]][0];
            _linkBPosi[_pa1RouteInv[r1]][0] = _pa1RouteInv[b1];
            _linkBPosi[_pa1RouteInv[r2]][0] = _pa1RouteInv[b2];
        }
    }
}

void Cross::MakeUnit ()
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
            fprintf(stderr, "ERROR: numSPL reach max(%d) in MakeUnit", n);
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

int Cross::MakeCompleteSol (Indi& kid)
{
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
                st = _pa1Route[posi];
            }
        }
        assert(st != -1);

        curr = -1;
        next = st;
        numEleInCU = 0;
        while (1) {
            prev = curr;
            curr = next;
            _centerUnit[curr] = 1;
            _listCenterUnit[numEleInCU] = curr;
            ++numEleInCU;

            if (kid._link[curr][0] != prev) {
                next = kid._link[curr][0];
            } else {
                next = kid._link[curr][1];
            }

            if (next == st) {
                break;
            }
        }
        _listCenterUnit[numEleInCU] = _listCenterUnit[0];
        _listCenterUnit[numEleInCU + 1] = _listCenterUnit[1];

        assert(numEleInCU == _numElementInUnit[center_un]);

        int max_diff = std::numeric_limits<int>::min();
        int diff = 0;
        a1 = -1;
        b1 = -1;
        /*
         *N_near (see Step 5.3 in Section 2.2 of the Online Supplement)
         * nearMax must be smaller than or equal to eva->_maxNumNear (Evaluator)
         */
        assert(NearMaxDef <= _eval->GetMaxNumNear());
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
                            d = kid._link[c][j2];
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
                    b1 = kid._link[j][0];
                    break;
                }
            }
            max_diff = _eval->GetCost(aa, bb) + _eval->GetCost(a1, b1)
                       - _eval->GetCost(a, a1) - _eval->GetCost(b, b1);
        }

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }

        _modiEdge[_numModiEdge][0] = aa;
        _modiEdge[_numModiEdge][1] = bb;
        _modiEdge[_numModiEdge][2] = a1;
        _modiEdge[_numModiEdge][3] = b1;
        ++_numModiEdge;

        gainModi += max_diff;

        int posi_a1 = _pa1RouteInv[a1];
        int select_un = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segment[s][0] <= posi_a1 && posi_a1 <= _segment[s][1]) {
                select_un = _segUnit[s];
                break;
            }
        }
        assert(select_un != -1);

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

void Cross::BackToPa1 (Indi& kid, int appliedCycle)
{
    int aa, bb, a1, b1;
    for (int s = _numModiEdge - 1; s >= 0; --s) {
        aa = _modiEdge[s][0];
        a1 = _modiEdge[s][1];
        bb = _modiEdge[s][2];
        b1 = _modiEdge[s][3];

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
    }

    ChangeSol(kid, appliedCycle, true /*reverse*/, false /*updateSeg*/);
}

void Cross::GoToBest (Indi& kid, int bestAppliedCycle)
{
    int aa, bb, a1, b1;

    ChangeSol(kid, bestAppliedCycle, false /*reverse*/, false /*updateSeg*/);

    for (int s = 0; s < _numBestModiEdge; ++s) {
        aa = _bestModiEdge[s][0];
        bb = _bestModiEdge[s][1];
        a1 = _bestModiEdge[s][2];
        b1 = _bestModiEdge[s][3];

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }
    }
}

} /* namespace thu */
