#include "Pricing.h"
#include "Node.h"
#include "Scene.h"
#include "WTA.h"
#include "SubProblem.h"

Pricing::Pricing(WTA* _wta, AlgoParameter _parameter) 
{
   wta = _wta;
   ptrNode = NULL;
   
   scenes = {};
   states = {};
   int n = wta->target_num_n;
   int m = wta->weapon_num_m;
   int dual_num = n + m;
   dual_value = new double[dual_num];
   parameter = _parameter;
}


/**
 * @brief Use this to get information from the node and Initialize a feasible solution
 * @todo Use a greedy method or more flexible solution
*/

void Pricing::Set(Node &node) {
   ptrNode = &node;
    int weapon_num = wta->weapon_num_m;
    int target_num = wta->target_num_n;
   bool is_test = 1;

   // current method, just random assign a feasible solution
   if(is_test){
      for(int j = 0; j < target_num; j++){
         Scene scene_j(j, vector<int>());
         scenes.push_back(scene_j);
      }
      for(int i = 0; i < weapon_num; i++){
         int assigned_target = i % target_num;
         scenes[assigned_target].weapon_indices.push_back(i);
      }
   }
    //find one scene for each target by greedy
   else{
   }
}
/**
 * @brief Solve the subproblems of target t
 * @param : dual values, conflict information
 * @return column index 
 */
void Pricing::Solve(vector<double> &dual) 
{
   scenes = {};
   vector<Scene> temp_scenes = {};
   bool all_constraint_satisfied = true;
   for(int i = 0; i < dual.size(); i++){
      dual_value[i] = dual[i];
   }
   int m = wta->weapon_num_m;
   int n = wta->target_num_n;
   double* weapon_dual = new double[m];
   for(int i = 0; i < m; i++){
      weapon_dual[i] = dual_value[i];
   }
   vector<int> Is_correct(n, 1);
   // @todo change those params to the real params
   int _seed = 0; 
   int _is_frac = 1;
   AlgoParameter temp_parameter = parameter;
   // generate a scene for each target from subproblem
   
   // if This Flag == 1, test correct by enumeration
   int DEBUG_CORRECT = 0;

   for(int t = 0; t < n ; t++)
   {  
      printf("Now calculate the scene %d\n", t);
      SubProblem sub_prob;
      sub_prob.Init(wta, t, weapon_dual, dual_value[m + t], _seed, _is_frac, &parameter);
      Scene temp_scene;
      vector<int> solve_by_OA = sub_prob.cal_optimal_scene(temp_scene);


      if(DEBUG_CORRECT){
         vector<int> solve_by_Enum = sub_prob.cal_optimal_scene_by_enum();
         printf("weapon_num : %d\n", sub_prob.weapon_num); 
         for(int i = 0; i < sub_prob.weapon_num; i++){
            if(solve_by_OA[i] != solve_by_Enum[i]){
               Is_correct[i] = 0;
            }
         }
      }
      temp_scenes.push_back(temp_scene);
      // @todo decide weather do not push scene when one target optval > 0 or all  
      // 子问题最优解大于0 说明此时对偶约束已经满足 即不需要添加该场景
      bool target_j_satisfied = sub_prob.is_optval_geq_zero();

      double real_redcost = sub_prob.cal_reduced_cost(temp_scene);

      printf("opt_val : %e, real_redcost : %e\n", sub_prob.sub_optval, real_redcost);
      if(!(real_redcost > -1e-8)){
         all_constraint_satisfied = false;
      }
      sub_prob.Delete();
   }


   if(DEBUG_CORRECT){
      for(int i = 0; i < n; i++){
         if(Is_correct[i]){
            printf("target %d is Correct!\n", i);
         }
      }
   }
   if(!all_constraint_satisfied){
      scenes = temp_scenes;
   }
   printf("Gain %d scenes\n", scenes.size());
   delete[] weapon_dual;
}

Pricing::~Pricing(){
   if(dual_value){
      delete[] dual_value;
      dual_value = NULL;
   }
}