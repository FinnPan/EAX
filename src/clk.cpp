#include <cstdlib>
#include "clk.h"

namespace thu { /* tsp heuristics */

#define MAXDEPTH       25   /* Shouldn't really be less than 2.             */
#define KICK_MAXDEPTH  50
#define IMPROVE_SWITCH -1   /* When to start using IMPROVE_KICKS (-1 never) */
#define MARK_LEVEL 10       /* Number of tour neighbors after 4-swap kick   */
#define BACKTRACK   4
#define MAX_BACK   12       /* Upper bound on the XXX_count entries         */
#define G_MULT 1.5
#define HUNT_PORTION_LONG 0.001
#define WALK_STEPS 50

ChainedLK::FlipStack::FlipStack (int total, int single, Flipper *f)
{
    stack = new Pair [total + single];
    num0 = total;
    flipper = f;
}

ChainedLK::FlipStack::~FlipStack ()
{
    delete[] stack;
    stack = nullptr;
}

void ChainedLK::FlipStack::Reset ()
{
    cnt = 0;
    num = num0;
}

void ChainedLK::FlipStack::Flip (int aPrev, int a, int b, int bNext)
{
    flipper->Flip(a, b);
    stack[cnt].first = a;
    stack[cnt].last = b;
    cnt++;
}

void ChainedLK::FlipStack::Unflip (int aPrev, int a, int b, int bNext)
{
    flipper->Flip(b, a);
    cnt--;
}

ChainedLK::Graph::Graph (int num, int edgeNum)
{
    n = num;
    goodList = new Edge* [n];
    edgeSpace = new Edge [(2 * edgeNum) + n];
    degree = new int [n];
    weirdMark = new int [n];
    weirdMagic = 0;
}

ChainedLK::Graph::~Graph ()
{
    delete[] goodList;
    goodList = nullptr;

    delete[] degree;
    degree = nullptr;

    delete[] weirdMark;
    weirdMark = nullptr;

    delete[] edgeSpace;
    edgeSpace = nullptr;
}

void ChainedLK::Graph::Reset (const Evaluator *eval, int edgeNum, int *edgeList)
{
    int c1, c2, w;
    Edge *p;

    weirdMagic = 0;

    for (int i = 0; i < n; i++) {
        degree[i] = 1;
        weirdMark[i] = 0;
    }
    for (int i = edgeNum - 1; i >= 0; i--) {
        degree[edgeList[2 * i]]++;
        degree[edgeList[(2 * i) + 1]]++;
    }

    p = edgeSpace;
    for (int i = 0; i < n; i++) {
        goodList[i] = p;
        p += (degree[i]);
        goodList[i][degree[i] - 1].weight = std::numeric_limits<int>::max();
        degree[i] = 0;
    }

    for (int i = edgeNum - 1; i >= 0; i--) {
        c1 = edgeList[2 * i];
        c2 = edgeList[(2 * i) + 1];
        w = eval->GetCost(c1, c2);
        Insertedge(c1, c2, w);
        Insertedge(c2, c1, w);
    }
}

void ChainedLK::Graph::Insertedge (int c1, int c2, int w)
{
    int i;
    auto *e = goodList[c1];
    for (i = degree[c1] - 1; i >= 0 && e[i].weight >= w; i--) {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other = e[i].other;
    }
    e[i + 1].weight = w;
    e[i + 1].other = c2;
    degree[c1]++;
}

ChainedLK::AddDel::AddDel (int num)
{
    int i = 0;
    while ((1 << i) < num) {
        i++;
    }
    m = (1 << i);

    addEdges = new char [m];
    delEdges = new char [m];
}

ChainedLK::AddDel::~AddDel ()
{
    delete[] addEdges;
    addEdges = nullptr;

    delete[] delEdges;
    delEdges = nullptr;
}

void ChainedLK::AddDel::Reset ()
{
    for (int i = 0; i < m; i++) {
        addEdges[i] = 0;
        delEdges[i] = 0;
    }
}

ChainedLK::Queue::Queue (int num)
{
    n = num;
    active = new char [n];
}

ChainedLK::Queue::~Queue ()
{
    delete[] active;
    active = nullptr;
}

void ChainedLK::Queue::Reset ()
{
    for (int i = 0; i < n; i++) {
        active[i] = 0;
    }
    q.clear();
}

void ChainedLK::Queue::AddToActiveQueue (int c)
{
    if (active[c] == 0) {
        active[c] = 1;
        q.push_back(c);
    }
}

int ChainedLK::Queue::PopFromActiveQueue ()
{
    int c = -1;
    if (!q.empty()) {
        c = q.front();
        q.pop_front();
        active[c] = 0;
    }
    return c;
}

void ChainedLK::LinKernighan (double *val)
{
    int start, cnt;
    double delta, totalWin = 0.0;

    while (1) {
        start = _q->PopFromActiveQueue();
        if (start == -1) break;

        delta = ImproveTour(start);
        if (delta > 0.0) {
            totalWin += delta;
            if (_winStack->cnt < _winStack->num) {
                for (int i = 0; i < _stack->cnt; i++) {
                    cnt = _winStack->cnt;
                    _winStack->stack[cnt].first = _stack->stack[i].first;
                    _winStack->stack[cnt].last = _stack->stack[i].last;
                    _winStack->stack[cnt].firstPrev = _stack->stack[i].firstPrev;
                    _winStack->stack[cnt].lastNext = _stack->stack[i].lastNext;
                    _winStack->cnt++;
                }
            } else if (_winCycle[0] == -1) {
                for (int i = 0; i < _stack->cnt; i++) {
                    cnt = _winStack->cnt;
                    _winStack->stack[cnt].first = _stack->stack[i].first;
                    _winStack->stack[cnt].last = _stack->stack[i].last;
                    _winStack->cnt++;
                }
                _flipper->GetCycle(_winCycle);
            }
            _stack->cnt = 0;
        }
    }

    if (_winCycle[0] == -1) {
        for (int i = 0; i < _stack->cnt; i++) {
            _winStack->stack[_winStack->cnt].first = _stack->stack[i].first;
            _winStack->stack[_winStack->cnt].last = _stack->stack[i].last;
            _winStack->cnt++;
        }
    }
    (*val) -= totalWin;
}

double ChainedLK::ImproveTour (int t1)
{
    int t2 = _flipper->Next(t1);
    int gain, Gstar = 0;

    gain = _eval->GetCost(t1, t2);
    _ad->MarkEdgeDel(t1, t2);

    if (Step(0, gain, &Gstar, t1, t2) == 0) {
        Gstar = WeirdSecondStep(gain, t1, t2);
    }
    _ad->UnmarkEdgeDel(t1, t2);

    if (Gstar) {
        _q->AddToActiveQueue(t1);
        _q->AddToActiveQueue(t2);
    }

    return (double) Gstar;
}

int ChainedLK::Step (int level, int gain, int *Gstar, int first, int last)
{
    int val, curr, newLast, hit = 0, oldG = gain;

    if (level >= BACKTRACK) {
        return StepNoBack(level, gain, Gstar, first, last);
    }

    std::vector<Graph::EdgeLook> lookList;
    LookAhead(first, last, gain, level, lookList);

    for (const auto &e : lookList) {
        curr = e.other;
        newLast = e.over;

        gain = oldG - e.diff;
        val = gain - _eval->GetCost(newLast, first);
        if (val > *Gstar) {
            *Gstar = val;
            hit++;
        }

        _stack->Flip(first, last, newLast, curr);

        if (level < MAXDEPTH) {
            _ad->MarkEdgeAdd(last, curr);
            _ad->MarkEdgeDel(curr, newLast);
            hit += Step(level + 1, gain, Gstar, first, newLast);
            _ad->UnmarkEdgeAdd(last, curr);
            _ad->UnmarkEdgeDel(curr, newLast);
        }

        if (!hit) {
            _stack->Unflip(first, last, newLast, curr);
        } else {
            _q->AddToActiveQueue(curr);
            _q->AddToActiveQueue(newLast);
            return 1;
        }
    }

    return 0;
}

int ChainedLK::StepNoBack (int level, int gain, int *Gstar, int first, int last)
{
    Graph::EdgeLook e;
    LookAheadNoBack(first, last, gain - *Gstar - level, &e);

    if (e.diff < std::numeric_limits<int>::max()) {
        if (e.mm) {
            int hit = 0;
            int curr = e.other;
            int newFirst = e.over;
            int val;

            gain -= e.diff;
            val = gain - _eval->GetCost(newFirst, last);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }

            _stack->Flip(curr, newFirst, first, last);

            if (level < MAXDEPTH) {
                _ad->MarkEdgeAdd(first, curr);
                _ad->MarkEdgeDel(curr, newFirst);
                hit += StepNoBack(level + 1, gain, Gstar, newFirst, last);
                _ad->UnmarkEdgeAdd(first, curr);
                _ad->UnmarkEdgeDel(curr, newFirst);
            }

            if (!hit) {
                _stack->Unflip(curr, newFirst, first, last);
                return 0;
            } else {
                _q->AddToActiveQueue(curr);
                _q->AddToActiveQueue(newFirst);
                return 1;
            }
        } else {
            int hit = 0;
            int curr = e.other;
            int newLast = e.over;
            int val;

            gain -= e.diff;
            val = gain - _eval->GetCost(newLast, first);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }

            _stack->Flip(first, last, newLast, curr);

            if (level < MAXDEPTH) {
                _ad->MarkEdgeAdd(last, curr);
                _ad->MarkEdgeDel(curr, newLast);
                hit += StepNoBack(level + 1, gain, Gstar, first, newLast);
                _ad->UnmarkEdgeAdd(last, curr);
                _ad->UnmarkEdgeDel(curr, newLast);
            }

            if (!hit) {
                _stack->Unflip(first, last, newLast, curr);
                return 0;
            } else {
                _q->AddToActiveQueue(curr);
                _q->AddToActiveQueue(newLast);
                return 1;
            }
        }
    } else {
        return 0;
    }
}

double ChainedLK::KickImprove ()
{
    int t1, t2;
    int gain, Gstar = 0;
    int hit = 0;

    do {
        FirstKicker(&t1, &t2);
        gain = _eval->GetCost(t1, t2);
        _ad->MarkEdgeDel(t1, t2);
        hit = KickStepNoBack(0, gain, &Gstar, t1, t2);
        _ad->UnmarkEdgeDel(t1, t2);
    } while (!hit);

    KickTurn(t1);
    KickTurn(t2);

    return (double)(-Gstar);
}

int ChainedLK::KickStepNoBack (int level, int gain, int *Gstar, int first, int last)
{
    Graph::EdgeLook winner;
    int val;
    int curr, prev, newLast;
    int lastNext = _flipper->Next(last);
    int cutoff = (int) (G_MULT * (double) gain);
    auto **goodList = _g->goodList;

    winner.diff = std::numeric_limits<int>::max();
    for (int i = 0; goodList[last][i].weight < cutoff; i++) {
        curr = goodList[last][i].other;
        if (!_ad->IsDeleted(last, curr) && curr != first && curr != lastNext) {
            prev = _flipper->Prev(curr);
            if (!_ad->IsAdded(curr, prev)) {
                val = goodList[last][i].weight - _eval->GetCost(curr, prev);
                if (val < winner.diff) {
                    winner.diff = val;
                    winner.other = curr;
                    winner.over = prev;
                }
            }
        }
    }

    if (winner.diff < std::numeric_limits<int>::max()) {
        curr = winner.other;
        newLast = winner.over;
        gain -= winner.diff;
        *Gstar = gain - _eval->GetCost(newLast, first);

        _stack->Flip(first, last, newLast, curr);
        KickTurn(curr);
        KickTurn(newLast);
        if (_winStack->cnt < _winStack->num) {
            _winStack->stack[_winStack->cnt].first = last;
            _winStack->stack[_winStack->cnt].last = newLast;
            _winStack->cnt++;
        }

        if (level < KICK_MAXDEPTH) {
            _ad->MarkEdgeAdd(last, curr);
            _ad->MarkEdgeDel(curr, newLast);
            KickStepNoBack(level+1, gain, Gstar, first, newLast);
            _ad->UnmarkEdgeAdd(last, curr);
            _ad->UnmarkEdgeDel(curr, newLast);
        }
        return 1;
    } else {
        return 0;
    }
}

int ChainedLK::WeirdSecondStep (int lenT1T2, int t1, int t2)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit;
    int t4Next;

    std::vector<Graph::EdgeLook> lookList1, lookList2, lookList3;
    WeirdLookAhead(lenT1T2, t1, t2, lookList1);

    for (const auto &h : lookList1) {
        t3 = h.other;
        t4 = h.over;

        oldG = lenT1T2 - h.diff;

        t4Next = _flipper->Next(t4);

        _ad->MarkEdgeAdd(t2, t3);
        _ad->MarkEdgeDel(t3, t4);
        _g->weirdMagic++;
        _g->weirdMark[t1] = _g->weirdMagic;
        _g->weirdMark[t2] = _g->weirdMagic;
        _g->weirdMark[t3] = _g->weirdMagic;
        _g->weirdMark[t4Next] = _g->weirdMagic;

        WeirdLookAhead2(oldG, t2, t3, t4, lookList2);

        for (const auto &e : lookList2) {
            t5 = e.other;
            t6 = e.over;

            _ad->MarkEdgeAdd(t4, t5);
            if (e.seq) {
                if (!e.side) {
                    gain = oldG - e.diff;
                    val = gain - _eval->GetCost(t6, t1);
                    if (val > Gstar) {
                        Gstar = val;
                    }
                    _stack->Flip(t1, t2, t6, t5);
                    _stack->Flip(t2, t5, t3, t4);

                    _ad->MarkEdgeDel(t5, t6);
                    hit = Step(2, gain, &Gstar, t1, t6);
                    _ad->UnmarkEdgeDel(t5, t6);

                    if (!hit && Gstar) {
                        hit = 1;
                    }

                    if (!hit) {
                        _stack->Unflip(t2, t5, t3, t4);
                        _stack->Unflip(t1, t2, t6, t5);
                    } else {
                        _ad->UnmarkEdgeAdd(t2, t3);
                        _ad->UnmarkEdgeDel(t3, t4);
                        _ad->UnmarkEdgeAdd(t4, t5);
                        _q->AddToActiveQueue(t3);
                        _q->AddToActiveQueue(t4);
                        _q->AddToActiveQueue(t5);
                        _q->AddToActiveQueue(t6);
                        return Gstar;
                    }
                } else {
                    gain = oldG - e.diff;
                    val = gain - _eval->GetCost(t6, t1);
                    if (val > Gstar) {
                        Gstar = val;
                    }
                    _stack->Flip(t1, t2, t3, t4);
                    _stack->Flip(t6, t5, t2, t4);
                    _stack->Flip(t1, t3, t6, t2);

                    _ad->MarkEdgeDel(t5, t6);
                    hit = Step(2, gain, &Gstar, t1, t6);
                    _ad->UnmarkEdgeDel(t5, t6);

                    if (!hit && Gstar) {
                        hit = 1;
                    }

                    if (!hit) {
                        _stack->Unflip(t1, t3, t6, t2);
                        _stack->Unflip(t6, t5, t2, t4);
                        _stack->Unflip(t1, t2, t3, t4);
                    } else {
                        _ad->UnmarkEdgeAdd(t2, t3);
                        _ad->UnmarkEdgeDel(t3, t4);
                        _ad->UnmarkEdgeAdd(t4, t5);
                        _q->AddToActiveQueue(t3);
                        _q->AddToActiveQueue(t4);
                        _q->AddToActiveQueue(t5);
                        _q->AddToActiveQueue(t6);
                        return Gstar;
                    }
                }
            } else {
                tG = oldG - e.diff;
                _ad->MarkEdgeDel(t5, t6);

                WeirdLookAhead3(tG, t2, t3, t6, lookList3);

                for (const auto &f : lookList3) {
                    t7 = f.other;
                    t8 = f.over;
                    gain = tG - f.diff;
                    if (!f.side) {
                        val = gain - _eval->GetCost(t8, t1);
                        if (val > Gstar) {
                            Gstar = val;
                        }
                        _stack->Flip(t1, t2, t8, t7);
                        _stack->Flip(t2, t7, t3, t4);
                        _stack->Flip(t7, t4, t6, t5);

                        _ad->MarkEdgeAdd(t6, t7);
                        _ad->MarkEdgeDel(t7, t8);
                        hit = Step(3, gain, &Gstar, t1, t8);
                        _ad->UnmarkEdgeDel(t6, t7);
                        _ad->UnmarkEdgeDel(t7, t8);

                        if (!hit && Gstar) {
                            hit = 1;
                        }

                        if (!hit) {
                            _stack->Unflip(t7, t4, t6, t5);
                            _stack->Unflip(t2, t7, t3, t4);
                            _stack->Unflip(t1, t2, t8, t7);
                        } else {
                            _ad->UnmarkEdgeAdd(t2, t3);
                            _ad->UnmarkEdgeDel(t3, t4);
                            _ad->UnmarkEdgeAdd(t4, t5);
                            _ad->UnmarkEdgeDel(t5, t6);
                            _q->AddToActiveQueue(t3);
                            _q->AddToActiveQueue(t4);
                            _q->AddToActiveQueue(t5);
                            _q->AddToActiveQueue(t6);
                            _q->AddToActiveQueue(t7);
                            _q->AddToActiveQueue(t8);
                            return Gstar;
                        }
                    } else {
                        val = gain - _eval->GetCost(t8, t1);
                        if (val > Gstar) {
                            Gstar = val;
                        }
                        _stack->Flip(t1, t2, t6, t5);
                        _stack->Flip(t1, t6, t8, t7);
                        _stack->Flip(t3, t4, t2, t5);

                        _ad->MarkEdgeAdd(t6, t7);
                        _ad->MarkEdgeDel(t7, t8);
                        hit = Step(3, gain, &Gstar, t1, t8);
                        _ad->UnmarkEdgeAdd(t6, t7);
                        _ad->UnmarkEdgeDel(t7, t8);

                        if (!hit && Gstar) {
                            hit = 1;
                        }

                        if (!hit) {
                            _stack->Unflip(t3, t4, t2, t5);
                            _stack->Unflip(t1, t6, t8, t7);
                            _stack->Unflip(t1, t2, t6, t5);
                        } else {
                            _ad->UnmarkEdgeAdd(t2, t3);
                            _ad->UnmarkEdgeDel(t3, t4);
                            _ad->UnmarkEdgeAdd(t4, t5);
                            _ad->UnmarkEdgeDel(t5, t6);
                            _q->AddToActiveQueue(t3);
                            _q->AddToActiveQueue(t4);
                            _q->AddToActiveQueue(t5);
                            _q->AddToActiveQueue(t6);
                            _q->AddToActiveQueue(t7);
                            _q->AddToActiveQueue(t8);
                            return Gstar;
                        }
                    }
                }
                _ad->UnmarkEdgeDel(t5, t6);
            }
            _ad->UnmarkEdgeAdd(t4, t5);
        }
        _ad->UnmarkEdgeAdd(t2, t3);
        _ad->UnmarkEdgeDel(t3, t4);
    }

    return 0;
}

void ChainedLK::LookAhead (int first, int last, int gain, int level,
        std::vector<Graph::EdgeLook>& lookList)
{
    static const int backtrack_count[BACKTRACK] = {4, 3, 3, 2};

    lookList.clear();

    int val;
    int curr, prev;
    int lastNext = _flipper->Next(last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, ahead = backtrack_count[level];
    auto **goodList = _g->goodList;

    for (int i = 0; i < ahead; i++) {
        value[i] = std::numeric_limits<int>::max();
    }
    value[ahead] = std::numeric_limits<int>::min();

    for (int i = 0; goodList[last][i].weight <= gain; i++) {
        curr = goodList[last][i].other;
        if (!_ad->IsDeleted(last, curr) && curr != first && curr != lastNext) {
            prev = _flipper->Prev(curr);
            if (!_ad->IsAdded(curr, prev)) {
                val = goodList[last][i].weight - _eval->GetCost(curr, prev);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                    }
                    value[k] = val;
                    other[k] = curr;
                    save[k] = prev;
                }
            }
        }
    }

    for (int i = 0; i < ahead; i++) {
        if (value[i] < std::numeric_limits<int>::max()) {
            lookList.push_back(Graph::EdgeLook{});
            lookList.back().diff = value[i];
            lookList.back().other = other[i];
            lookList.back().over = save[i];
        }
    }
    std::reverse(lookList.begin(), lookList.end());
}

void ChainedLK::LookAheadNoBack (int first, int last, int gain, Graph::EdgeLook* winner)
{
    int val;
    int curr, prev;
    int firstPrev, lastNext = _flipper->Next(last);
    int next;
    auto **goodList = _g->goodList;

    winner->diff = std::numeric_limits<int>::max();
    for (int i = 0; goodList[last][i].weight < gain; i++) {
        curr = goodList[last][i].other;
        if (!_ad->IsDeleted(last, curr) && curr != first && curr != lastNext) {
            prev = _flipper->Prev(curr);
            if (!_ad->IsAdded(curr, prev)) {
                val = goodList[last][i].weight - _eval->GetCost(curr, prev);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = curr;
                    winner->over = prev;
                    winner->mm = 0;
                }
            }
        }
    }

    firstPrev = _flipper->Prev(first);
    for (int i = 0; goodList[first][i].weight < gain; i++) {
        curr = goodList[first][i].other;
        if (!_ad->IsDeleted(first, curr) && curr != last && curr != firstPrev) {
            next = _flipper->Next(curr);
            if (!_ad->IsAdded(curr, next)) {
                val = goodList[first][i].weight - _eval->GetCost(curr, next);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = curr;
                    winner->over = next;
                    winner->mm = 1;
                }
            }
        }
    }
}

void ChainedLK::WeirdLookAhead (int gain, int t1, int t2,
        std::vector<Graph::EdgeLook>& lookList)
{
    const int weird_backtrack_count_1 = 4;

    lookList.clear();

    int curr, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val, ahead;
    auto **goodList = _g->goodList;

    ahead = weird_backtrack_count_1;
    for (int i = 0; i < ahead; i++) {
        value[i] = std::numeric_limits<int>::max();
    }
    value[ahead] = std::numeric_limits<int>::min();

    for (int i = 0; goodList[t2][i].weight <= gain; i++) {
        curr = goodList[t2][i].other;
        if (curr != t1) {
            next = _flipper->Next(curr);
            val = goodList[t2][i].weight - _eval->GetCost(curr, next);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                }
                value[k] = val;
                other[k] = curr;
                save[k] = next;
            }
        }
    }

    for (int i = 0; i < ahead; i++) {
        if (value[i] < std::numeric_limits<int>::max()) {
            lookList.push_back(Graph::EdgeLook{});
            lookList.back().diff = value[i];
            lookList.back().other = other[i];
            lookList.back().over = save[i];
        }
    }
    std::reverse(lookList.begin(), lookList.end());
}

void ChainedLK::WeirdLookAhead2 (int gain, int t2, int t3, int t4,
        std::vector<Graph::EdgeLook>& lookList)
{
    const int weird_backtrack_count_2 = 3;

    lookList.clear();

    int t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count_2;
    auto **goodList = _g->goodList;
    int  *weirdMark = _g->weirdMark;
    int  weirdMagic = _g->weirdMagic;

    for (int i = 0; i < ahead; i++) {
        value[i] = std::numeric_limits<int>::max();
    }
    value[ahead] = std::numeric_limits<int>::min();

    for (int i = 0; goodList[t4][i].weight <= gain; i++) {
        t5 = goodList[t4][i].other;
        if (weirdMark[t5] != weirdMagic) {
            if (_flipper->Sequence(t2, t5, t3)) {
                t6 = _flipper->Prev(t5);
                val = goodList[t4][i].weight - _eval->GetCost(t5, t6);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 0;
                }
                t6 = _flipper->Next(t5);
                val = goodList[t4][i].weight - _eval->GetCost(t5, t6);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 1;
                }
            } else {
                t6 = _flipper->Prev(t5);
                val = goodList[t4][i].weight - _eval->GetCost(t5, t6);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 0;
                    side[k] = 0;
                }
            }
        }
    }

    for (int i = 0; i < ahead; i++) {
        if (value[i] < std::numeric_limits<int>::max()) {
            lookList.push_back(Graph::EdgeLook{});
            lookList.back().diff = value[i];
            lookList.back().other = other[i];
            lookList.back().over = save[i];
            lookList.back().seq = seq[i];
            lookList.back().side = side[i];
        }
    }
    std::reverse(lookList.begin(), lookList.end());
}

void ChainedLK::WeirdLookAhead3 (int gain, int t2, int t3, int t6,
        std::vector<Graph::EdgeLook>& lookList)
{
    const int weird_backtrack_count_3 = 3;

    lookList.clear();

    int t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count_3;
    auto **goodList = _g->goodList;
    int  *weirdMark = _g->weirdMark;
    int  weirdMagic = _g->weirdMagic;

    for (int i = 0; i < ahead; i++) {
        value[i] = std::numeric_limits<int>::max();
    }
    value[ahead] = std::numeric_limits<int>::min();

    for (int i = 0; goodList[t6][i].weight <= gain; i++) {
        t7 = goodList[t6][i].other;   /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdMark[t7] != weirdMagic && _flipper->Sequence(t2, t7, t3)) {
            t8 = _flipper->Prev(t7);
            val = goodList[t6][i].weight - _eval->GetCost(t7, t8);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 0;
            }
            t8 = _flipper->Next(t7);
            val = goodList[t6][i].weight - _eval->GetCost(t7, t8);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 1;
            }
        }
    }

    for (int i = 0; i < ahead; i++) {
        if (value[i] < std::numeric_limits<int>::max()) {
            lookList.push_back(Graph::EdgeLook{});
            lookList.back().diff = value[i];
            lookList.back().other = other[i];
            lookList.back().over = save[i];
            lookList.back().side = side[i];
        }
    }
    std::reverse(lookList.begin(), lookList.end());
}

int ChainedLK::RandomFourSwap (int *delta)
{
    int t1, t2, t3, t4, t5, t6, t7, t8;

    FindWalkTour(&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);

    if (!_flipper->Sequence(t1, t3, t5)) {
        std::swap(t3, t5);
        std::swap(t4, t6);
    }
    if (!_flipper->Sequence(t1, t5, t7)) {
        std::swap(t5, t7);
        std::swap(t6, t8);
        if (!_flipper->Sequence(t1, t3, t5)) {
            std::swap(t3, t5);
            std::swap(t4, t6);
        }
    }
    _stack->Flip(t1, t2, t5, t6);
    _stack->Flip(t4, t3, t7, t8);
    _stack->Flip(t1, t5, t6, t2);

    if (_winStack->cnt < _winStack->num) {
        _winStack->stack[_winStack->cnt].first = t2;
        _winStack->stack[_winStack->cnt].last = t5;
        _winStack->cnt++;
    }
    if (_winStack->cnt < _winStack->num) {
        _winStack->stack[_winStack->cnt].first = t3;
        _winStack->stack[_winStack->cnt].last = t7;
        _winStack->cnt++;
    }
    if (_winStack->cnt < _winStack->num) {
        _winStack->stack[_winStack->cnt].first = t5;
        _winStack->stack[_winStack->cnt].last = t6;
        _winStack->cnt++;
    }

    BigTurn(t1, 0);
    BigTurn(t2, 1);
    BigTurn(t3, 0);
    BigTurn(t4, 1);
    BigTurn(t5, 0);
    BigTurn(t6, 1);
    BigTurn(t7, 0);
    BigTurn(t8, 1);

    *delta =
            _eval->GetCost(t1, t6) + _eval->GetCost(t2, t5) +
            _eval->GetCost(t3, t8) + _eval->GetCost(t4, t7) -
            _eval->GetCost(t1, t2) - _eval->GetCost(t3, t4) -
            _eval->GetCost(t5, t6) - _eval->GetCost(t7, t8);

    return 0;
}

void ChainedLK::FirstKicker (int *t1, int *t2)
{
    const int n = _g->n;
    const int longCount = (int) ((double)(n) * HUNT_PORTION_LONG) + 10;
    int best, try1, len, next, prev, nextl, prevl;
    auto **goodList = _g->goodList;

    try1 = _eval->GetRand() % n;
    next = _flipper->Next(try1);
    prev = _flipper->Prev(try1);
    nextl = _eval->GetCost(try1, next);
    prevl = _eval->GetCost(try1, prev);
    if (nextl >= prevl) {
        *t1 = try1;
        *t2 = next;
        best = nextl - goodList[*t1][0].weight;
    } else {
        *t1 = prev;
        *t2 = try1;
        best = prevl - goodList[*t1][0].weight;
    }

    for (int i = 0; i < longCount; i++) {
        try1 = _eval->GetRand() % n;
        next = _flipper->Next(try1);
        prev = _flipper->Prev(try1);
        nextl = _eval->GetCost(try1, next);
        prevl = _eval->GetCost(try1, prev);
        if (nextl >= prevl) {
            len = nextl - goodList[try1][0].weight;
            if (len > best) {
                *t1 = try1;
                *t2 = next;
            }
        } else {
            len = prevl - goodList[try1][0].weight;
            if (len > best) {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
}

void ChainedLK::FindWalkTour (int *t1, int *t2, int *t3, int *t4,
        int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, c;

    FirstKicker(&s1, &s2);

    do {
        old = -1;
        c = s2;

        for (int i = 0;  i < WALK_STEPS; i++) {
            int j = _eval->GetRand() % (_g->degree[c]);
            if (old != _g->goodList[c][j].other) {
                old = c;
                c = _g->goodList[c][j].other;
            }
        }
        s3 = c;
        s4 = _flipper->Next(s3);

        c = s4;
        for (int i = 0; i < WALK_STEPS; i++) {
            int j = _eval->GetRand() % (_g->degree[c]);
            if (old != _g->goodList[c][j].other) {
                old = c;
                c = _g->goodList[c][j].other;
            }
        }
        s5 = c;
        s6 = _flipper->Next(s5);

        c = s6;
        for (int i = 0; i < WALK_STEPS; i++) {
            int j = _eval->GetRand() % (_g->degree[c]);
            if (old != _g->goodList[c][j].other) {
                old = c;
                c = _g->goodList[c][j].other;
            }
        }
        s7 = c;
        s8 = _flipper->Next(s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7 ||
             s1 == s8 ||
             s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6 || s2 == s7 ||
             s2 == s8 ||
             s3 == s5 || s3 == s6 || s3 == s7 || s3 == s8 ||
             s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8 ||
             s5 == s7 || s5 == s8 ||
             s6 == s7 || s6 == s8);

    *t1 = s1;  *t2 = s2;  *t3 = s3;  *t4 = s4;
    *t5 = s5;  *t6 = s6;  *t7 = s7;  *t8 = s8;
}

void ChainedLK::KickTurn (int c)
{
    int k;
    _q->AddToActiveQueue(c);
    k = _flipper->Next(c);
    _q->AddToActiveQueue(k);
    k = _flipper->Next(k);
    _q->AddToActiveQueue(k);
    k = _flipper->Prev(c);
    _q->AddToActiveQueue(k);
    k = _flipper->Prev(k);
    _q->AddToActiveQueue(k);
}

void ChainedLK::BigTurn (int c, int toNext)
{
    int k = c;

    _q->AddToActiveQueue(c);
    if (toNext) {
        for (int i = 0; i < MARK_LEVEL; i++) {
            k = _flipper->Next(k);
            _q->AddToActiveQueue(k);
        }
    } else {
        for (int i = 0; i < MARK_LEVEL; i++) {
            k = _flipper->Prev(k);
            _q->AddToActiveQueue(k);
        }
    }

    for (int i = 0; i < _g->degree[c]; i++) {
        _q->AddToActiveQueue(_g->goodList[c][i].other);
    }
}

ChainedLK::ChainedLK (const Evaluator *eval) :
    _eval(eval), _edgeList(nullptr), _edgeNum(0),
    _verbose(false)
{
    const int n = _eval->GetNumCity();
    _winCycle = new int[n];
    _flipper = new Flipper(n);

    InitEdgeList();

    const int hit = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    _winStack = new FlipStack(500 + n / 50, hit, _flipper);
    _stack = new FlipStack(hit, 0, _flipper);
    _g = new Graph(n, _edgeNum);
    _ad = new AddDel(n);
    _q = new Queue(n);
}

ChainedLK::~ChainedLK ()
{
    delete[] _winCycle;
    _winCycle = nullptr;

    delete _flipper;
    _flipper = nullptr;

    if (_edgeList) {
        delete[] _edgeList;
        _edgeList = nullptr;
        _edgeNum = 0;
    }

    delete _winStack;
    _winStack = nullptr;
    delete _stack;
    _stack = nullptr;
    delete _g;
    _g = nullptr;
    delete _ad;
    _ad = nullptr;
    delete _q;
    _q = nullptr;
}

void ChainedLK::InitEdgeList ()
{
    const int n = _eval->GetNumCity();
    const int maxNear = _eval->GetMaxNumNear();
    const int maxNum = maxNear * n * 2;

    _edgeList = new int[maxNum];
    _edgeNum = 0;

    for (int ci = 0; ci < n; ++ci) {
        int cnt = 0;
        for (int ni = 0; ni < maxNear; ++ni) {
            int cj = _eval->GetNear(ci, ni);
            if (cj > ci) {
                int idx = 2 * _edgeNum + 2 * cnt;
                _edgeList[idx] = ci;
                _edgeList[idx + 1] = cj;
                cnt++;
            }
        }
        _edgeNum += cnt;
    }
}

bool ChainedLK::DoIt (double &cost)
{
    _winStack->Reset();
    _stack->Reset();
    _g->Reset(_eval, _edgeNum, _edgeList);
    _ad->Reset();
    _q->Reset();

    const int n = _eval->GetNumCity();

    /* the max number of 4-swaps without progress */
    const int stallcount = 100000000;

    /* the number of 4-swap kicks */
    const int repeatcount = (n < 10) ? 0 : n;

    RUsage ru;
    int round = 0;
    int quitcount, delta;
    double best;
    const int *cycle;

    if (_verbose) {
        printf("linkern ...\n");
    }

    if (n < 10) {
        printf("Less than 10 nodes, setting repeatcount to 0\n");
    }

    cycle = _eval->MakeRand();
    _flipper->SetCycle(cycle);
    cost = _eval->ComputeCost(_flipper);
    best = cost;
    if (_verbose) {
        printf("Starting Cycle: %.0f\n", cost);
    }

    quitcount = stallcount;
    if (quitcount > repeatcount) {
        quitcount = repeatcount;
    }

    _stack->cnt = 0;
    _winStack->cnt = 0;
    _winCycle[0] = -1;

    /* init active_queue with random order */
    cycle = _eval->MakeRand();
    for (int i = 0; i < n; i++) {
        _q->AddToActiveQueue(cycle[i]);
    }

    LinKernighan(&best);

    _winStack->cnt = 0;
    _winCycle[0] = -1;

    if (_verbose) {
        if (quitcount > 0) {
            printf("%4d Steps   Best: %.0f\n", round, best);
            ru.Report("quitcount");
        } else {
            printf("LK Cycle: %.0f\n", best);
        }
    }

    while (round < quitcount) {
        int hit = 0;
        _stack->cnt = 0;

        if (IMPROVE_SWITCH == -1 || round < IMPROVE_SWITCH) {
            int rval = RandomFourSwap(&delta);
            if (rval) {
                fprintf(stderr, "RandomFourSwap failed\n");
                return false;
            }
        } else {
            delta = (int)KickImprove();
        }

        _stack->cnt = 0;
        double t = best + delta;
        LinKernighan(&t);

        if (t <= best) {
            _winStack->cnt = 0;
            _winCycle[0] = -1;
            if (t < best) {
                best = t;
                quitcount = round + stallcount;
                if (quitcount > repeatcount) {
                    quitcount = repeatcount;
                }
                hit++;
            }
        } else {
            if (_winCycle[0] == -1) {
                while (_winStack->cnt) {
                    _winStack->cnt--;
                    int c = _winStack->cnt;
                    _flipper->Flip(_winStack->stack[c].last, _winStack->stack[c].first);
                }
            } else {
                _flipper->SetCycle(_winCycle);
                while (_winStack->cnt) {
                    _winStack->cnt--;
                    int c = _winStack->cnt;
                    _flipper->Flip(_winStack->stack[c].last, _winStack->stack[c].first);
                }
                _winCycle[0] = -1;
            }
        }

        round++;
        if (_verbose && (hit || (round % 1000 == 999))) {
            printf("%4d Steps   Best: %.0f\n", round, best);
            ru.Report("hit");
        }
    }
    if (_verbose && round > 0) {
        printf("%4d Total Steps.\n", round);
    }

    double t = _eval->ComputeCost(_flipper);
    if (t != best) {
        printf("WARNING: LK incremental counter was off by %.0f\n", t-best);
        best = t;
    }
    cost = best;

    return true;
}

} /* namespace thu */
