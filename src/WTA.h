#ifndef _WTA_H_
#define _WTA_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iterator>
#include <vector>
#include <algorithm>
#include <time.h>
#include <cassert>
#include <functional>
#include <unordered_map>
#include "Scene.h"

using std::vector;
using std::string;
using std::unordered_map;

#define EPS     1e-6
#define BUFSIZE 1000

class WTA
{
public:
    int         weapon_num_m;
    int         target_num_n;
    // P = p[weapon][target]
    double**    probability_matrix;
    
    // V = V[t]
    double*     target_value;
    int         DEBUG_MODE; //
public:
    WTA();
    WTA(int n_, int m_, int seed);
    ~WTA();

    // Initialize the WTA by Random Generator or ReadFile
    void Init(int n_, int m_, int seed);
    void Init_sparse(int n_, int m_, int sparsity, int seed);

    // calculate the x_is, q_js
    double cal_linear_coef(int _j, vector<int> _scene);
    bool is_weapon_in_scene(int _i, vector<int> _scene);
    

    double set_scene_qjs(Scene& scene);

    // print the prob matrix and value of the WTA
    void print_model();
};

#endif