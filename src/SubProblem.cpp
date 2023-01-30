#include "SubProblem.h"

SubProblem::SubProblem(){

}

SubProblem::~SubProblem(){
    delete wta;
    if(weapon_dual_ui){
        delete[] weapon_dual_ui;
    }
}

vector<int> SubProblem::cal_optimal_scene(){
    int DEBUG_MODE = parameter->DEBUG_VERSION;
}