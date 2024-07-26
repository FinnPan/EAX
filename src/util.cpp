#include "util.h"

#ifdef WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

namespace thu { /* tsp heuristics */

RUsage::RUsage ()
{
    Reset();
}

RUsage::~RUsage ()
{
}

void RUsage::Reset ()
{
    _dayTime0 = GetTimeOfDay();
    _cpuTime0 = GetTimeOfCPU();
}

void RUsage::Report (const char *tag) const
{
    std::chrono::duration<double> elapsedTime = GetTimeOfDay() - _dayTime0;
    double cpuTime = GetTimeOfCPU() - _cpuTime0;
    double physMem, virtMem;
    GetProcessMem(physMem, virtMem);

    printf("[%s] Elapsed = %.2fs; CPU = %.2fs; MEM = %.1fM\n",
            tag, elapsedTime.count(), cpuTime, physMem);
    fflush(stdout);
}

double RUsage::GetTimeOfCPU ()
{
    double cpu = 0.0;

#ifdef WIN32
    auto ConvertTimeFormat = [](const FILETIME *fTime) {
        SYSTEMTIME sTime;
        FileTimeToSystemTime(fTime, &sTime);
        double t = sTime.wHour * 3600.0 + sTime.wMinute * 60.0 +
                   sTime.wSecond + sTime.wMilliseconds / 1000.0;
        return t;
    };

    FILETIME createTime, exitTime, kernelTime, userTime;
    if (GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime)) {
        cpu = ConvertTimeFormat(&kernelTime) + ConvertTimeFormat(&userTime);
    }
#else
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
        cpu += ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1000000.0;
        cpu += ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1000000.0;
    }
#endif

    return cpu;
}

void RUsage::GetProcessMem (double &physMem, double &virtMem)
{
    physMem = 0.0;
    virtMem = 0.0;

#ifdef WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        physMem = pmc.WorkingSetSize / 1024.0 / 1024.0;
        virtMem = pmc.PagefileUsage / 1024.0 / 1024.0;
    }
#else
    FILE *fd = fopen("/proc/self/status", "r");
    if (fd) {
        char line[256];
        char key[128];
        char unit[32];
        unsigned long value;
        while (fgets(line, 256, fd)) {
            if (strstr(line, "Vm") && sscanf(line, "%s %lu %s", key, &value, unit) == 3) {
                if (strstr(key, "VmRSS")) {
                    physMem = value / 1024.0;
                } else if (strstr(key, "VmSize")) {
                    virtMem = value / 1024.0;
                }
            }
        }
        fclose(fd);
    }
#endif
}

TspLib::TspLib () :
    _dimension(0), _optimal(0), _x(nullptr), _y(nullptr), _edgeLen(nullptr)
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

    _dimension = -1;

    char buf[256], key[256], field[256];
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
                _dimension = atoi(field);
                printf("Number of Nodes: %d\n", GetDimension());
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
                if (GetDimension() <= 0) {
                    fprintf(stderr, "ERROR: dimension not specified\n");
                    goto CLEAN_UP;
                }
                _x = new double[GetDimension()];
                _y = new double[GetDimension()];
                for (int i = 0; i < GetDimension(); i++) {
                    fscanf(in, "%*d %lf %lf", &(_x[i]), &(_y[i]));
                }
            } else if (!strcmp(key, "OPTIMAL")) {
                if (sscanf(p, "%s", field) == EOF) {
                    fprintf(stderr, "ERROR: in OPTIMAL line\n");
                    goto CLEAN_UP;
                }
                _optimal = atoi(field);
                printf("Optimal cost: %d\n", GetOptimal());
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

Flipper::Flipper (int count) : _numCity(count)
{
    Init_0();
}

Flipper::Flipper (int count, const int *cyc) : _numCity(count)
{
    Init_0();
    SetCycle(cyc);
}

Flipper::~Flipper ()
{
    delete[] _parents;
    delete[] _children;
}

void Flipper::Init_0 ()
{
    const int n = _numCity;

    constexpr double GROUP_SIZE_FACTOR = 0.50;
    constexpr double SEGMENT_SPLIT_CUTOFF = 0.30;

    _groupSize = (int)(sqrt(n) * GROUP_SIZE_FACTOR);
    _numSegments =  (n + _groupSize - 1) / _groupSize;
    _splitCutoff = (int)(_groupSize * SEGMENT_SPLIT_CUTOFF);
    _reversed = 0;
    _parents = new ParentNode[_numSegments];
    _children = new ChildNode[n];
}

void Flipper::Init_1 ()
{
    _reversed = 0;

    int remain = _numCity;
    int cnt = 0;
    while (remain >= 2 * _groupSize) {
        int s = _groupSize;
        _parents[cnt++].size = s;
        remain -= s;
    }
    if (remain > _groupSize) {
        int s = remain / 2;
        _parents[cnt++].size = s;
        remain -= s;
    }
    _parents[cnt++].size = remain;

    assert(cnt == _numSegments);
}

void Flipper::SetCycle (const int *cyc)
{
    Init_1();

    ChildNode *c, *c0;
    ParentNode *p;

    c0 = _children + cyc[_numCity - 1];
    for (int i = 0, cnt = 0; i < _numSegments; i++) {
        p = _parents + i;
        p->id = i;
        p->rev = 0;
        p->ends[0] = _children + cyc[cnt];
        c = nullptr;
        for (int j = 0; j < p->size; j++) {
            c = _children + cyc[cnt];
            c->parent = p;
            c->id = cnt;
            c->name = cyc[cnt];
            c->adj[0] = c0;
            c0->adj[1] = c;
            c0 = c;
            cnt++;
        }
        p->ends[1] = c;
        p->adj[0] = p - 1;
        p->adj[1] = p + 1;
    }
    _parents[0].adj[0] = _parents + _numSegments - 1;
    _parents[_numSegments - 1].adj[1] = _parents;
}

void Flipper::GetCycle (int *cyc) const
{
    ChildNode *c0 = _children;
    ChildNode *c = c0;
    int cnt = 0;
    do {
        cyc[cnt++] = c->name;
        c = _children + Next(c->name);
    } while(c != c0);
}

int Flipper::Next (int x) const
{
    ChildNode *c = _children + x;
    return c->adj[IsForward(c->parent)]->name;
}

int Flipper::Prev (int x) const
{
    ChildNode *c = _children + x;
    return c->adj[IsBackward(c->parent)]->name;
}

bool Flipper::Sequence (int x, int y, int z) const
{
    ChildNode *a = _children + x;
    ChildNode *b = _children + y;
    ChildNode *c = _children + z;
    ParentNode *pa = a->parent;
    ParentNode *pb = b->parent;
    ParentNode *pc = c->parent;

    if (pa == pb) {
        if (pa == pc) {
            if (IsBackward(pa)) {
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
            if (IsBackward(pa)) {
                return (a->id >= b->id);
            } else {
                return (a->id <= b->id);
            }
        }
    } else if (pa == pc) {
        if (IsBackward(pa)) {
            return (a->id <= c->id);
        } else {
            return (a->id >= c->id);
        }
    } else if (pb == pc) {
        if (IsBackward(pb)) {
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
    ChildNode *xc = _children + x;
    ChildNode *yc = _children + y;
    const bool xBackward = IsBackward(xc->parent);
    const bool yForward = IsForward(yc->parent);
    ChildNode *xPrev = xc->adj[xBackward];
    ChildNode *yNext = yc->adj[yForward];

    if (SameSegment(xc, yc)) {
        SameSegmentFlip(xc, yc);
        return;
    }

    if (SameSegment(yNext, xPrev)) {
        SameSegmentFlip(yNext, xPrev);
        Reverse();
        return;
    }

    if (xc->parent->ends[xBackward] == xc &&
        yc->parent->ends[yForward] == yc) {
        Flip_0(xc, yc);
        return;
    }

    if (xPrev->parent == xc->parent) {
        SegmentSplit(xc->parent, xPrev, xc, 0);
        if (SameSegment(xc, yc)) {
            SameSegmentFlip(xc, yc);
            return;
        } else if (SameSegment(yNext, xPrev)) {
            SameSegmentFlip(yNext, xPrev);
            Reverse();
            return;
        }
    }

    if (yNext->parent == yc->parent) {
        SegmentSplit(yc->parent, yc, yNext, 0);
        if (SameSegment(xc, yc)) {
            SameSegmentFlip(xc, yc);
            return;
        } else if (SameSegment(yNext, xPrev)) {
            SameSegmentFlip(yNext, xPrev);
            Reverse();
            return;
        }
    }

    Flip_0(xc, yc);
}

void Flipper::Flip_0 (ChildNode *xc,  ChildNode *yc)
{
    const bool segBackward = _reversed;
    const bool segForward = !segBackward;
    int side = segBackward? (xc->parent->id - yc->parent->id)
                          : (yc->parent->id - xc->parent->id);
    if (side < 0) {
        side += _numSegments;
    }

    if (side < _numSegments / 2) {
        ConsecutiveSegmentFlip(xc->parent, yc->parent);
    } else {
        ConsecutiveSegmentFlip(yc->parent->adj[segForward],
                               xc->parent->adj[segBackward]);
        Reverse();
    }
}

bool Flipper::SameSegment (ChildNode *a, ChildNode *b) const
{
    if (a->parent == b->parent) {
        if (IsBackward(a->parent)) {
            return (a->id >= b->id);
        } else {
            return (a->id <= b->id);
        }
    }
    return false;
}

void Flipper::SameSegmentFlip (ChildNode *a, ChildNode *b) const
{
    ParentNode *parent = a->parent;
    const bool backward = IsBackward(parent);
    const bool forward = !backward;
    ChildNode *aPrev = a->adj[backward];
    ChildNode *bNext = b->adj[forward];

    if (a == b) {
        return;
    }

    if ((backward && a->id - b->id > _splitCutoff) ||
        (forward && b->id - a->id > _splitCutoff)) {
        if (aPrev->parent == parent) {
            SegmentSplit(parent, aPrev, a, 1);
        }
        if (bNext->parent == parent) {
            SegmentSplit(parent, b, bNext, 2);
        }
        aPrev->adj[IsForward(aPrev->parent)] = b;
        bNext->adj[IsBackward(bNext->parent)] = a;
        a->adj[backward] = bNext;
        b->adj[forward] = aPrev;
        parent->rev ^= 1;
    } else {
        int id_1 = backward? -1 : 1;
        int id = a->id;
        ChildNode *c = b;
        while (c != aPrev) {
            ChildNode *cNext = c->adj[backward];
            c->adj[backward] = c->adj[forward];
            c->adj[forward] = cNext;
            c->id = id;
            id += id_1;
            c = cNext;
        }
        a->adj[forward] = bNext;
        bNext->adj[IsBackward(bNext->parent)] = a;
        b->adj[backward] = aPrev;
        aPrev->adj[IsForward(aPrev->parent)] = b;
        if (parent->ends[backward] == a) {
            parent->ends[backward] = b;
        }
        if (parent->ends[forward] == b) {
            parent->ends[forward] = a;
        }
    }
}

void Flipper::ConsecutiveSegmentFlip (ParentNode *a, ParentNode *b) const
{
    const bool aBackward = IsBackward(a);
    const bool bForward = IsForward(b);
    ChildNode *aChild = a->ends[aBackward];
    ChildNode *bChild = b->ends[bForward];
    ChildNode *childPrev = aChild->adj[aBackward];
    ChildNode *childNext = bChild->adj[bForward];
    childPrev->adj[IsForward(childPrev->parent)] = bChild;
    bChild->adj[bForward] = childPrev;
    childNext->adj[IsBackward(childNext->parent)] = aChild;
    aChild->adj[aBackward] = childNext;

    const bool segBackward = _reversed;
    const bool segForward = !segBackward;
    ParentNode *aPrev = a->adj[segBackward];
    ParentNode *bNext = b->adj[segForward];
    int id_1 = segBackward? -1 : 1;
    int id = a->id;
    ParentNode *c = b;
    while (c != aPrev) {
        ParentNode *cNext = c->adj[segBackward];
        c->adj[segBackward] = c->adj[segForward];
        c->adj[segForward] = cNext;
        c->id = id;
        id += id_1;
        c->rev ^= 1;
        c = cNext;
    }
    a->adj[segForward] = bNext;
    bNext->adj[segBackward] = a;
    b->adj[segBackward] = aPrev;
    aPrev->adj[segForward] = b;
}

/* split between a and aPrev */
void Flipper::SegmentSplit (ParentNode *p, ChildNode *aPrev, ChildNode *a,
        int left_or_right) const
{
    const bool segBackward = _reversed;
    const bool segForward = !segBackward;
    const bool backward = IsBackward(p);
    const bool forward = !backward;
    const int side = backward? (p->ends[1]->id - aPrev->id + 1)
                             : (aPrev->id - p->ends[0]->id + 1);
    const bool do_left_side = ((left_or_right == 0 && side <= p->size / 2) ||
                               (left_or_right == 1));

    int id, id_1;
    ParentNode *pNext;
    ChildNode *b, *bNext;
    bool nextBackward, nextForward;

    if (do_left_side) {
        pNext = p->adj[segBackward];
        nextBackward = IsBackward(pNext);
        nextForward = !nextBackward;
        pNext->size += side;
        p->size -= side;

        b = pNext->ends[nextForward];
        bNext = b->adj[nextForward];
        id = b->id;
        id_1 = nextBackward? -1 : 1;
        do {
            b = bNext;
            id += id_1;
            b->id = id;
            b->parent = pNext;
            bNext = b->adj[forward];
            if (pNext->rev != p->rev) {
                b->adj[forward] = b->adj[backward];
                b->adj[backward] = bNext;
            }
        } while (b != aPrev);
        pNext->ends[nextForward] = aPrev;
        p->ends[backward] = a;
    } else {
        pNext = p->adj[segForward];
        nextBackward = IsBackward(pNext);
        nextForward = !nextBackward;
        pNext->size += (p->size - side);
        p->size = side;

        b = pNext->ends[nextBackward];
        bNext = b->adj[nextBackward];
        id = b->id;
        id_1 = nextForward? -1 : 1;
        do {
            b = bNext;
            id += id_1;
            b->id = id;
            b->parent = pNext;
            bNext = b->adj[backward];
            if (pNext->rev != p->rev) {
                b->adj[backward] = b->adj[forward];
                b->adj[forward] = bNext;
            }
        } while (b != a);
        pNext->ends[nextBackward] = a;
        p->ends[forward] = aPrev;
    }
}

Evaluator::Evaluator () :
    _maxNumNear(50), _randBuf(nullptr), _nearTbl(nullptr)
{
    std::random_device rd;
    _randEng = new RandEngine(rd());
}

Evaluator::~Evaluator ()
{
    delete _randEng;
    delete[] _randBuf;

    const int n = GetNumCity();
    for (int i = 0; i < n; ++i) {
        delete[] _nearTbl[i];
    }
    delete[] _nearTbl;
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

    const int n = GetNumCity();

    if (n <= 2) {
        fprintf(stderr, "ERROR: invalid city number (%d)\n", n);
        exit(1);
    }

    _randBuf = new int[n];
    for (int i = 0; i < n; ++i) {
        _randBuf[i] = i;
    }

    const int maxNumNear = GetMaxNumNear();
    _nearTbl = new int*[n];
    for (int i = 0; i < n; ++i) {
        _nearTbl[i] = new int[maxNumNear];
    }
    BuildNeighborLists();

    return true;
}

/* times: O(n*n*log_k)
 * space: O(n*k) */
void Evaluator::BuildNeighborLists ()
{
    const int n = GetNumCity();
    const int maxNumNear = GetMaxNumNear();
    RUsage ru;

    class LessCmp {
    public:
        explicit LessCmp (const Evaluator *eval, int ci) :
            _eval(eval), _ci(ci) {}
        ~LessCmp () = default;
        bool operator() (const int l, const int r) const {
            double costL = _eval->GetCost(l, _ci);
            double costR = _eval->GetCost(r, _ci);
            return costL < costR;
        }
    private:
        const Evaluator *_eval;
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

int Evaluator::ComputeCost (const int *arr) const
{
    const int n = GetNumCity();
    int cost = 0;
    for (int i = 0; i < n; i++) {
        cost += GetCost(arr[i], arr[(i+1)%n]);
    }
    return cost;
}

int Evaluator::ComputeCost (const Flipper *f) const
{
    int cost = 0;
    int c1 = 0;
    do {
        int c2 = f->Next(c1);
        cost += GetCost(c1, c2);
        c1 = c2;
    } while (c1 != 0);
    return cost;
}

const int *Evaluator::MakeRand () const
{
    const int n = GetNumCity();
    std::shuffle(_randBuf, _randBuf + n, GetRandEngine());
    return _randBuf;
}

double Evaluator::ComputeGap (int cost) const
{
    double gap = NAN;
    if (_tspLib.HasOptimal()) {
        double opt = _tspLib.GetOptimal();
        gap = (cost - opt) / opt;
    }
    return gap;
}

void TwoOpt::DoIt ()
{
    const int *arr = _eval->MakeRand();
    _flipper.SetCycle(arr);
    TwoExchange();
}

void TwoOpt::TwoExchange ()
{
    const int n = _eval->GetNumCity();
    const int maxNumNear = _eval->GetMaxNumNear();
    int t1, t2, t3, t4, cost12, cost34, cost23, cost14;
    bool improved;

    /*
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
     */

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

} /* namespace thu */
