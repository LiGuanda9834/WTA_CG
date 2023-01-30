#ifndef _SUBPROBLEM_H__
#define _SUBPROBLEM_H__

#include "WTA.h"
#include "Scene.h"
#include "AlgorithmParameter.h"
#include <vector>

class SubProblem{
public:
    WTA*            wta;
    int             target_index;
    double*         weapon_dual_ui;
    double          target_dual_vj;
    AlgoParameter*  parameter;

public:
    SubProblem();
    ~SubProblem();

    vector<int>     cal_optimal_scene();
};

#endif