#ifndef __THU_CLK_H
#define __THU_CLK_H

#include "util.h"

namespace thu { /* tsp heuristics */

class ChainedLK {
public:
    explicit ChainedLK (const Evaluator *e);
    ~ChainedLK ();
    void SetVerbose (bool v) { _verbose = v; }
    bool DoIt (double &cost);
private:
    void InitEdgeList ();
private:
    const Evaluator *_eval;
    int *_edgeList;
    int _edgeNum;
    bool _verbose;
};

} /* namespace thu */

#endif /* __THU_CLK_H */
