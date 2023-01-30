#include "Pricing.h"
#include "Node.h"
#include "Scene.h"
#include "WTA.h"

Pricing::Pricing(WTA* _wta) 
{
   wta = _wta;
   ptrNode = NULL;
   
   scenes = {};
   states = {};
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
   for(int t = 0; t < wta->target_num_n ; t++)
   {  
      Scene scene;
      
      // @todo generate a scene from the subproblem
      scenes.push_back(scene);
   }
}