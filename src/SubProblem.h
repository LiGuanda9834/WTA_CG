#ifndef _SUBPROBLEM_H__
#define _SUBPROBLEM_H__

#include "WTA.h"
#include "Scene.h"
#include "AlgorithmParameter.h"
#include <vector>

#include <ilcplex/cplex.h>

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

#define WTA_EPS (1e-6)

class SubProblem{
public:
   /**
    * 1 for solved and optval >= 0
    * 0 for solved and optval < 0
   */
    int             sub_status;
    WTA*            wta;

    // variables to describe WTA model
    int             weapon_num;
    int             target_index;
    double*         weapon_dual_ui;
    double          target_dual_vj;

    double          target_value;
    double*         j_prob_vector;

    // varible used in LP solve, cutind, cutval and rhs is the linear constraint created by x
    double*         x;
    int*            cutind;
    double*         cutval;
    double          rhs;


    double          sub_optval;

    /**
     * @brief   to set algorithm param
     * @todo    Add these params to the parameter of total algorithm
    */ 
    int             seed;
    int             separate_fractional_solutions;
    /**
     * @brief   Is any bug here
     * 0 for no bug
    */

    AlgoParameter*  parameter;

public:
    SubProblem() = default;
    ~SubProblem() = default;

    void Init(WTA* wta, int t_index, double* w_dual_ui, double t_dual_vj, int seed, int is_fractional, AlgoParameter* param);

    double destroy_prob();

    // cal fx and cutval is depend on the member variable x. Only the x is set, this two function could be correct. 
    double cal_fx();
    double cal_cutval(int i);
    double cal_rhs();

    double cal_reduced_cost(Scene& temp_scene);
    
    void    cal_constraint_by_x(vector<int> weapon_sets);

    double cal_eta_lower_bound();
    double cal_eta_upper_bound();

// 之后将该返回值改为return status 确保能够得到问题的求解状态
   vector<int> cal_optimal_scene(Scene& scene);
   vector<int> cal_optimal_scene_by_enum();

   // A bad name, need change
    bool is_optval_geq_zero();

   // Use this to check weather a solution is frac or integer
    bool print_solution();


    void print_LP_info();
    void print_model();
    void print_debug();
   
   /**
    * @param cut_num:    Init how many cut
    * @param wt_ratio:   In one cut how many nonzero weapon
    * @return cut_point: take down the activated weapon indices in each vector<int>
   */
    void find_init_cut(vector<vector<int>>& cut_point, int cut_num, int wt_ratio);

    void Delete();
};

/* Declaration of the data structure for the function benders_callback */

typedef struct {
   /* Parameter to decide when Benders' cuts are going to be separated:
0: only when an integer solution if found
(i.e., wherefrom == CPX_CALLBACK_MIP_CUT_FEAS )
1: even to cut-off fractional solutions, 
at the end of the cplex cut-loop
(i.e., wherefrom == CPX_CALLBACK_MIP_CUT_LAST || 
wherefrom == CPX_CALLBACK_MIP_CUT_FEAS ) */
   int separate_fractional_solutions; 
   /* Environment for the worker LP used to generate Benders' cuts */
   CPXENVptr env; 
   /* Worker LP used to generate Benders' cuts */
   CPXLPptr  lp;  
   /* number of nodes in the ATSP instance */
   int num_nodes; 
   /* Number of columns in the master ILP */
   int num_x_cols;
   /* Number of columns in the worker LP */
   int num_v_cols, num_u_cols;
   /* Data structure to 
      -- read the solution from the master ILP 
      -- update the objective function of the worker LP */
   double *x;
   int *indices;
   /* Data structure to read an unbounded direction (ray) from the worker LP */
   double *ray;
   /* Data structure to add a Benders' cut to the master ILP */
   double *cutval;
   int *cutind;

} USER_CBHANDLE;


/* Declarations for functions in this program */

static int CPXPUBLIC 
benders_callback  (CPXCENVptr env, void *cbdata, int wherefrom, 
      void *cbhandle, int *useraction_p);

static int
set_benders_callback  (CPXENVptr env, SubProblem* sub_prob),
create_master_ILP     (CPXENVptr env, CPXLPptr lp, SubProblem* sub_prob),
init_user_cbhandle    (USER_CBHANDLE *user_cbhandle, int num_nodes, int separate_fractional_solutions),
free_user_cbhandle    (USER_CBHANDLE *user_cbhandle);


static void
free_and_null  (char **ptr),
usage          (char *progname);

#endif