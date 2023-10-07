#include "util.h"

namespace th { /* tsp heuristics */

RUsage::~RUsage ()
{
    if (!_tagRAII.empty()) {
        Report(_tagRAII.c_str());
    }
}

void RUsage::Report (const char* tag) const
{
    assert(tag);
    TimePoint end = Clock::now();
    std::chrono::duration<double> elapsed = end - _start;
    printf("[%s] Elapsed time: %.3f s.\n", tag, elapsed.count());
    fflush(stdout);
}

void RUsage::SetRAIIReport (const char* tagRAII)
{
    assert(tagRAII);
    _tagRAII = tagRAII;
}

TspLib::TspLib () :
    _dim(0), _x(nullptr), _y(nullptr), _edgeLen(nullptr)
{}

TspLib::~TspLib ()
{
    delete[] _x;
    delete[] _y;
}

int TspLib::EdgeLen_euclid (double xd, double yd)
{
    return (int)(sqrt(xd * xd + yd * yd) + 0.5);
}

int TspLib::EdgeLen_att (double xd, double yd)
{
    double rij = sqrt((xd * xd + yd * yd) / 10.0);

    int k = (int)rij;
    double tij = (double)k;

    int dij;
    if (tij < rij) {
        dij = (int)tij + 1;
    } else {
        dij = (int)tij;
    }
    return dij;
}

int TspLib::EdgeLen_euclid_ceiling (double xd, double yd)
{
    return (int)(ceil(sqrt(xd * xd + yd * yd)));
}

bool TspLib::Init (const char *fileName)
{
    if (_x || _y) {
        fprintf(stderr, "ERROR: init TspLib again\n");
        return false;
    }

    _dim = -1;

    char buf[256], key[256], field[256];
    int norm = -1;
    FILE *in = fopen(fileName, "r");

    if (!in) {
        fprintf(stderr, "ERROR: unable to open %s for input\n", fileName);
        return false;
    }

    while (fgets(buf, 254, in)) {
        char *p = buf;
        while (*p != '\0') {
            if (*p == ':') {
                *p = ' ';
            }
            p++;
        }
        p = buf;

        if (sscanf(p, "%s", key) != EOF) {
            p += strlen(key);
            while (*p == ' ') {
                p++;
            }

            if (!strcmp(key, "NAME")) {
                printf("Problem Name: %s", p);
            } else if (!strcmp(key, "TYPE")) {
                printf("Problem Type: %s", p);
                if (sscanf(p, "%s", field) == EOF || strcmp(field, "TSP")) {
                    fprintf(stderr, "ERROR: not a TSP problem\n");
                    goto CLEAN_UP;
                }
            } else if (!strcmp(key, "COMMENT")) {
                printf("%s", p);
            } else if (!strcmp(key, "DIMENSION")) {
                if (sscanf(p, "%s", field) == EOF) {
                    fprintf(stderr, "ERROR: in DIMENSION line\n");
                    goto CLEAN_UP;
                }
                _dim = atoi(field);
                printf("Number of Nodes: %d\n", GetDim());
            } else if (!strcmp(key, "EDGE_WEIGHT_TYPE")) {
                if (sscanf(p, "%s", field) == EOF) {
                    fprintf(stderr, "ERROR: in EDGE_WEIGHT_TYPE line\n");
                    goto CLEAN_UP;
                }
                if (!strcmp(field, "EUC_2D")) {
                    _edgeLen = EdgeLen_euclid;
                    printf("Rounded Euclidean Norm\n");
                } else if (!strcmp(field, "ATT")) {
                    _edgeLen = EdgeLen_att;
                    printf("ATT Norm\n");
                } else if (!strcmp(field, "CEIL_2D")) {
                    _edgeLen = EdgeLen_euclid_ceiling;
                    printf("Rounded Up Euclidean Norm\n");
                } else {
                    fprintf(stderr, "ERROR: not set up for norm %s\n", field);
                    goto CLEAN_UP;
                }
            } else if (!strcmp(key, "NODE_COORD_SECTION")) {
                if (GetDim() <= 0) {
                    fprintf(stderr, "ERROR: dimension not specified\n");
                    goto CLEAN_UP;
                }
                _x = new double[GetDim()];
                _y = new double[GetDim()];
                assert(_x && _y);
                for (int i = 0; i < GetDim(); i++) {
                    fscanf(in, "%*d %lf %lf", &(_x[i]), &(_y[i]));
                }
            }
        }
    }

CLEAN_UP:
    fclose(in);

    if (!_x || !_y) {
        fprintf(stderr, "ERROR: didn't find the data\n");
        return false;
    }

    return true;
}

Flipper::Flipper (int count)
{
    Init_0();
    Init_1(count);
}

Flipper::Flipper (int count, const int *cyc)
{
    Init_0();
    Init_1(count);
    SetCycle(count, cyc);
}

Flipper::~Flipper ()
{
    delete[] _parents;
    delete[] _children;
}

void Flipper::Init_0 ()
{
    _parents = nullptr;
    _children = nullptr;
    _reversed = 0;
    _groupSize = 0;
    _numSegments = 0;
    _splitCutoff = 0;
}

#define GROUPSIZE_FACTOR 0.50
#define SEGMENT_SPLIT_CUTOFF 0.30
void Flipper::Init_1 (int count)
{
    _reversed = 0;
    _groupSize = (int)(sqrt((double)count) * GROUPSIZE_FACTOR);
    _numSegments =  (count + _groupSize - 1) / _groupSize;
    _splitCutoff = _groupSize * SEGMENT_SPLIT_CUTOFF;
    _parents = new ParentNode[_numSegments];
    /* The +1 will stop a purify burp later */
    _children = new ChildNode[count + 1];
    assert(_parents && _children);

    int remain = count;
    int i = 0;
    int j = 2 * _groupSize;
    while (remain >= j) {
        _parents[i].size = _groupSize;
        remain -= _groupSize;
        i++;
    }
    if (remain > _groupSize) {
        _parents[i].size = remain / 2;
        remain -= (remain / 2);
        i++;
    }
    _parents[i].size = remain;
    i++;

    if (i != _numSegments) {
        fprintf (stderr, "ERROR: seg count is wrong\n");
        assert(0);
    }
}

void Flipper::SetCycle (int count, const int *cyc)
{
    ChildNode *c, *cPrev;
    ParentNode *p;

    c = &(_children[cyc[count - 1]]);
    p = _parents;
    for (int i = 0, cind = 0; i < _numSegments; p++, i++) {
        p->id = i;
        p->rev = 0;
        p->ends[0] = &(_children[cyc[cind]]);
        for (int j = p->size; j > 0; j--) {
            cPrev = c;
            c = &(_children[cyc[cind]]);
            c->id = cind;
            c->name = cyc[cind];
            c->parent = p;
            c->adj[0] = cPrev;
            cPrev->adj[1] = c;
            cind++;
        }
        p->ends[1] = c;
        p->adj[0] = p - 1;
        p->adj[1] = p + 1;
    }
    _parents[0].adj[0] = &(_parents[_numSegments - 1]);
    _parents[_numSegments - 1].adj[1] = &(_parents[0]);
}

void Flipper::GetCycle (int *x) const
{
    ChildNode *c, *start;
    int k = 0;

    start = &(_children[0]);
    c = start->adj[!((_reversed)^(start->parent->rev))];

    x[k++] = start->name;
    while (c != start) {
        x[k++] = c->name;
        c = c->adj[!((_reversed)^(c->parent->rev))];
    }
}

int Flipper::Next (int x) const
{
    return
      _children[x].adj[!((_reversed)^(_children[x].parent->rev))]->name;
}

int Flipper::Prev (int x) const
{
    return
      _children[x].adj[(_reversed)^(_children[x].parent->rev)]->name;
}

int Flipper::Sequence (int x, int y, int z) const
{
    ChildNode *a = &(_children[x]);
    ChildNode *b = &(_children[y]);
    ChildNode *c = &(_children[z]);
    ParentNode *pa = a->parent;
    ParentNode *pb = b->parent;
    ParentNode *pc = c->parent;

    if (pa == pb) {
        if (pa == pc) {
            if ((_reversed)^(pa->rev)) {
                if (a->id >= b->id) {
                    return (b->id >= c->id || c->id >= a->id);
                } else {
                    return (b->id >= c->id && c->id >= a->id);
                }
            } else {
                if (a->id <= b->id) {
                    return (b->id <= c->id || c->id <= a->id);
                } else {
                    return (b->id <= c->id && c->id <= a->id);
                }
            }
        } else {
            if ((_reversed)^(pa->rev)) {
                return (a->id >= b->id);
            } else {
                return (a->id <= b->id);
            }
        }
    } else if (pa == pc) {
            if ((_reversed)^(pa->rev)) {
                return (a->id <= c->id);
            } else {
                return (a->id >= c->id);
            }
    } else if (pb == pc) {
            if ((_reversed)^(pb->rev)) {
                return (b->id >= c->id);
            } else {
                return (b->id <= c->id);
            }
    } else {
        if (_reversed) {
            if (pa->id >= pb->id) {
                return (pb->id >= pc->id || pc->id >= pa->id);
            } else {
                return (pb->id >= pc->id && pc->id >= pa->id);
            }
        } else {
            if (pa->id <= pb->id) {
                return (pb->id <= pc->id || pc->id <= pa->id);
            } else {
                return (pb->id <= pc->id && pc->id <= pa->id);
            }
        }
    }
}

void Flipper::Flip (int x, int y)
{
    ChildNode *xc = &(_children[x]);
    ChildNode *yc = &(_children[y]);

    if (SameSegmant(xc, yc)) {
        if (xc != yc) {
            SameSegmentFlip(xc, yc);
        }
    } else {
        int xDir = ((_reversed)^(xc->parent->rev));
        int yDir = ((_reversed)^(yc->parent->rev));
        ChildNode *xPrev = xc->adj[xDir];
        ChildNode *yNext = yc->adj[!yDir];
        if (SameSegmant(yNext, xPrev)) {
            if (yNext != xPrev) {
                SameSegmentFlip(yNext, xPrev);
            }
            (_reversed) ^= 1;
        } else {
            int side;
            if (xc->parent->ends[xDir] == xc &&
                yc->parent->ends[!yDir] == yc) {
                if (_reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += _numSegments;
                if (side < _numSegments / 2) {
                    ConsecutiveSegmentFlip(xc->parent, yc->parent);
                } else {
                    ConsecutiveSegmentFlip(yc->parent->adj[!_reversed],
                                           xc->parent->adj[_reversed]);
                    (_reversed) ^= 1;
                }
            } else {
                if (xPrev->parent == xc->parent) {
                    SegmentSplit(xc->parent, xPrev, xc, 0);
                    if (SameSegmant(xc, yc)) {
                        if (xc != yc)
                            SameSegmentFlip(xc, yc);
                        return;
                    } else if (SameSegmant(yNext, xPrev)) {
                        if (yNext != xPrev) {
                            SameSegmentFlip(yNext, xPrev);
                        }
                        (_reversed) ^= 1;
                        return;
                    }
                }
                if (yNext->parent == yc->parent) {
                    SegmentSplit(yc->parent, yc, yNext, 0);
                    if (SameSegmant(xc, yc)) {
                        if (xc != yc)
                            SameSegmentFlip(xc, yc);
                        return;
                    } else if (SameSegmant(yNext, xPrev)) {
                        if (yNext != xPrev) {
                            SameSegmentFlip(yNext, xPrev);
                        }
                        (_reversed) ^= 1;
                        return;
                    }
                }
                if (_reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += _numSegments;
                if (side < _numSegments / 2) {
                    ConsecutiveSegmentFlip(xc->parent, yc->parent);
                } else {
                    ConsecutiveSegmentFlip(yc->parent->adj[!_reversed],
                                           xc->parent->adj[_reversed]);
                    (_reversed) ^= 1;
                }

            }
        }
    }
}

bool Flipper::SameSegmant (ChildNode *a, ChildNode *b) const
{
    return (
        a->parent == b->parent &&
        ( (!((_reversed)^(a->parent->rev)) && a->id <= b->id) ||
          (((_reversed)^(a->parent->rev)) && a->id >= b->id) )
    );
}

void Flipper::SameSegmentFlip (ChildNode *a, ChildNode *b)
{
    ParentNode *parent = a->parent;
    int dir = ((_reversed)^(parent->rev));
    ChildNode *aPrev = a->adj[dir];
    ChildNode *bNext = b->adj[!dir];
    ChildNode *c, *cNext;

    if ((dir && a->id - b->id > _splitCutoff) ||
       (!dir && b->id - a->id > _splitCutoff)) {
        if (aPrev->parent == parent)
            SegmentSplit(parent, aPrev, a, 1);
        if (bNext->parent == parent)
            SegmentSplit(parent, b, bNext, 2);
        aPrev->adj[!((_reversed)^(aPrev->parent->rev))] = b;
        bNext->adj[(_reversed)^(bNext->parent->rev)] = a;
        a->adj[dir] = bNext;
        b->adj[!dir] = aPrev;
        parent->rev ^= 1;
        return;
    }

    if (dir) {
        int id = a->id;
        aPrev->adj[!((_reversed)^(aPrev->parent->rev))] = b;
        bNext->adj[(_reversed)^(bNext->parent->rev)] = a;
        cNext = b->adj[1];
        b->adj[1] = aPrev;
        b->adj[0] = cNext;
        b->id = id--;
        c = cNext;
        while (c != a) {
            cNext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cNext;
            c->id = id--;
            c = cNext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bNext;
        a->id = id;
        if (parent->ends[1] == a)
            parent->ends[1] = b;
        if (parent->ends[0] == b)
            parent->ends[0] = a;
    } else {
        int id = a->id;
        aPrev->adj[!((_reversed)^(aPrev->parent->rev))] = b;
        bNext->adj[(_reversed)^(bNext->parent->rev)] = a;
        c = b->adj[0];
        b->adj[0] = aPrev;
        b->adj[1] = c;
        b->id = id++;
        while (c != a) {
            cNext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cNext;
            c->id = id++;
            c = cNext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bNext;
        a->id = id;
        if (parent->ends[0] == a)
            parent->ends[0] = b;
        if (parent->ends[1] == b)
            parent->ends[1] = a;
    }
}

void Flipper::ConsecutiveSegmentFlip (ParentNode *a, ParentNode *b)
{
    ParentNode *aPrev = a->adj[_reversed];
    ParentNode *bNext = b->adj[!_reversed];
    ParentNode *c, *cNext;
    ChildNode *aChild = a->ends[(_reversed)^(a->rev)];
    ChildNode *bChild = b->ends[!((_reversed)^(b->rev))];
    ChildNode *childPrev, *childNext;
    int id = a->id;

    if (_reversed) {
        childPrev = aChild->adj[!a->rev];
        childNext = bChild->adj[b->rev];
        childPrev->adj[childPrev->parent->rev] = bChild;
        childNext->adj[!childNext->parent->rev] = aChild;
        bChild->adj[b->rev] = childPrev;
        aChild->adj[!a->rev] = childNext;

        aPrev->adj[0] = b;
        bNext->adj[1] = a;
        c = b->adj[1];
        b->adj[1] = aPrev;
        b->adj[0] = c;
        b->id = id--;
        b->rev ^= 1;
        while (c != a) {
            cNext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cNext;
            c->id = id--;
            c->rev ^= 1;
            c = cNext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bNext;
        a->id = id;
        a->rev ^= 1;
    } else {
        childPrev = aChild->adj[a->rev];
        childNext = bChild->adj[!b->rev];
        childPrev->adj[!childPrev->parent->rev] = bChild;
        childNext->adj[childNext->parent->rev] = aChild;
        bChild->adj[!b->rev] = childPrev;
        aChild->adj[a->rev] = childNext;

        aPrev->adj[1] = b;
        bNext->adj[0] = a;
        c = b->adj[0];
        b->adj[0] = aPrev;
        b->adj[1] = c;
        b->id = id++;
        b->rev ^= 1;
        while (c != a) {
            cNext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cNext;
            c->id = id++;
            c->rev ^= 1;
            c = cNext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bNext;
        a->id = id;
        a->rev ^= 1;
    }
}

/* split between a and aPrev */
void Flipper::SegmentSplit (ParentNode *p, ChildNode *aPrev, ChildNode *a, int left_or_right)
{
    int side;
    int dir = ((_reversed)^(p->rev));
    int id;
    ParentNode *pNext;
    ChildNode *b, *bNext;

    if (dir) side = p->ends[1]->id - aPrev->id + 1;
    else     side = aPrev->id - p->ends[0]->id + 1;

    if ((left_or_right == 0 && side <= p->size / 2) || left_or_right == 1) {
        pNext = p->adj[_reversed];
        pNext->size += side;
        p->size -= side;
        if (pNext->rev == p->rev) {
            b = pNext->ends[!dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pNext;
                } while (b != aPrev);
            } else {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pNext;
                } while (b != aPrev);
            }
            pNext->ends[!dir] = aPrev;
            p->ends[dir] = a;
        } else {
            b = pNext->ends[dir];
            id = b->id;
            if (!dir) {
                bNext = b->adj[0];
                do {
                    b = bNext;
                    b->id = --id;
                    b->parent = pNext;
                    bNext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bNext;
                } while (b != aPrev);
            } else {
                bNext = b->adj[1];
                do {
                    b = bNext;
                    b->id = ++id;
                    b->parent = pNext;
                    bNext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bNext;
                } while (b != aPrev);
            }
            pNext->ends[dir] = aPrev;
            p->ends[dir] = a;
        }
    } else {
        pNext = p->adj[!_reversed];
        pNext->size += (p->size - side);
        p->size = side;
        if (pNext->rev == p->rev) {
            b = pNext->ends[dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pNext;
                } while (b != a);
            } else {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pNext;
                } while (b != a);
            }
            pNext->ends[dir] = a;
            p->ends[!dir] = aPrev;
        } else {
            b = pNext->ends[!dir];
            id = b->id;
            if (!dir) {
                bNext = b->adj[1];
                do {
                    b = bNext;
                    b->id = ++id;
                    b->parent = pNext;
                    bNext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bNext;
                } while (b != a);
            } else {
                bNext = b->adj[0];
                do {
                    b = bNext;
                    b->id = --id;
                    b->parent = pNext;
                    bNext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bNext;
                } while (b != a);
            }
            pNext->ends[!dir] = a;
            p->ends[!dir] = aPrev;
        }
    }
}

Evaluator::Evaluator () : _maxNumNear(50), _nearTbl(nullptr), _routeBuf(nullptr)
{
    std::random_device rd;
    _rand = new RandEngine(rd());
}

Evaluator::~Evaluator ()
{
    const int n = GetNumCity();
    for (int i = 0; i < n; ++i) {
        delete[] _nearTbl[i];
    }
    delete[] _nearTbl;

    delete[] _routeBuf;
    delete _rand;
}

bool Evaluator::Init (const char filename[])
{
    if (_nearTbl) {
        fprintf(stderr, "ERROR: init Evaluator again\n");
        return false;
    }

    if (!_tspLib.Init(filename)) {
        return false;
    }

    /* allocate memory */
    const int n = GetNumCity();
    const int maxNumNear = GetMaxNumNear();
    _nearTbl = new int*[n];
    for (int i = 0; i < n; ++i) {
        _nearTbl[i] = new int[maxNumNear];
    }
    _routeBuf = new int[n];

    BuildNeighborLists();

    return true;
}

/* times: O(n*n*logk)
 * space: O(n*k)
 */
void Evaluator::BuildNeighborLists ()
{
    const int n = GetNumCity();
    const int maxNumNear = GetMaxNumNear();
    RUsage ru;

    class LessCmp {
    public:
        explicit LessCmp (const Evaluator* eval, int ci) : _eval(eval), _ci(ci) {}
        ~LessCmp () = default;
        bool operator() (const int l, const int r) const {
            double costL = _eval->GetCost(l, _ci);
            double costR = _eval->GetCost(r, _ci);
            return costL < costR;
        }
    private:
        const Evaluator* _eval;
        int _ci;
    };

    for (int ci = 0; ci < n; ++ci) {
        LessCmp less(this, ci);
        std::priority_queue<int, std::vector<int>, LessCmp> pqu(less);
        for (int i = 0; i < maxNumNear + 1; ++i) {
            pqu.push(i);
        }
        for (int i = maxNumNear + 1; i < n; i++) {
            int top = pqu.top();
            if (less(i, top)) {
                pqu.pop();
                pqu.push(i);
            }
        }

        int cnt = 0;
        while (!pqu.empty()) {
            int cj = pqu.top();
            pqu.pop();
            if (ci != cj) {
                _nearTbl[ci][cnt++] = cj;
                if (cnt >= maxNumNear) {
                    break;
                }
            }
        }
        std::reverse(_nearTbl[ci], _nearTbl[ci] + maxNumNear);
        assert(cnt == maxNumNear);
    }
    ru.Report("neighbor-lists");
}

int Evaluator::DoIt (const int* route) const
{
    const int n = GetNumCity();
    int cost = 0;
    for (int i = 0; i < n; i++) {
        cost += GetCost(route[i], route[(i+1)%n]);
    }
    return cost;
}

int Evaluator::DoIt (const Flipper* f) const
{
    const int n = GetNumCity();
    int cost = 0;
    int c1, c2;

    c1 = 0;
    do {
        if (c1 == 0) {
            c2 = f->Prev(c1);
            cost += GetCost(c1, c2);
        }

        c2 = f->Next(c1);
        cost += GetCost(c1, c2);
        c1 = c2;
    } while (c1 != 0);

    return cost;
}

const int* Evaluator::MakeRand () const
{
    const int n = GetNumCity();
    for (int i = 0; i < n; ++i) {
        _routeBuf[i] = i;
    }
    std::shuffle(_routeBuf, _routeBuf + n, GetRandEnginge());
    return _routeBuf;
}

void TwoOpt::DoIt ()
{
    const int n = _eval->GetNumCity();
    const int* route = _eval->MakeRand();
    _flipper.SetCycle(n, route);
    TwoExchange();
}

void TwoOpt::TwoExchange ()
{
    const int n = _eval->GetNumCity();
    const int maxNumNear = _eval->GetMaxNumNear();
    int t1, t2, t3, t4, cost12, cost34, cost23, cost14;
    bool improved;

    /*******************************
     * from:
     *   t1 -> t2
     *    |     |
     *    |     |
     *   t3 <- t4
     * to:
     *   t1    t2         t1 -> t4
     *    | \ / |          |     |
     *    | / \ |    or:   |     |
     *   t3    t4         t3 <- t2
     *******************************/

    do {
        improved = false;
        for (t1 = 0; t1 < n; t1++) {
            t2 = _flipper.Next(t1);
            cost12 = _eval->GetCost(t1, t2);
            for (int i = 0; i < maxNumNear; i++) {
                t3 = _eval->GetNear(t2, i);
                cost23 = _eval->GetCost(t2, t3);

                /* speed-up: new edge-2-3 should be shorter than old edge-1-2 */
                if (cost23 >= cost12) {
                    break;
                } else {
                    t4 = _flipper.Prev(t3);
                    cost34 = _eval->GetCost(t3, t4);
                    cost14 = _eval->GetCost(t1, t4);
                    if (cost12 + cost34 > cost23 + cost14) {
                        _flipper.Flip(t2, t4);
                        improved = true;
                        break;
                    }
                }
            }
        }
    } while(improved);
}

} /* namespace th */
