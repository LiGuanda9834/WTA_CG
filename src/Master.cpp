#include "Master.h"

Master::Master(
   WTA* _wta, 
   const AlgoParameter &_param
): wta(_wta), param(_param)
{
   scenepool = ScenePool();
   numLp = 0;
   time = 0;
   lp = new LP_ALL_IN_ONE(wta->target_num_n, wta->weapon_num_m);
}

// Use this to get information from the node
void Master::Set(const Node &node) 
{
   // @todo Update this function
   for( int i = 0; i < (int)scenepool.size(); i++ )
   {
      const Scene &scene = scenepool[i];
      if( node.valid(path) )
      {
         if( y[i].getLB() != 0 || y[i].getUB() != 1 )
            y[i].setBounds(0, 1);
      }
      else
      {
         if( y[i].getLB() != 0 || y[i].getUB() != 0 )
            y[i].setBounds(0, 0);
      }
   }
   printf("Master Problem has been Initialized \n");
}

/**
 * @brief Use this to solve the LP by solver like HIGHS
 * @todo need to update the result
*/
bool Master::Solve() 
{
   numLp++;
   bool result = true;
   return result;
}

void Master::GetDualValues(vector<double> &dual) 
{
}


/**
 * @brief Use this to add scenes to the current master RP's ScenePool
 * @note  Do not delete scenes from the scenes Pool
*/
void Master::AddCol(vector<Scene> &scenes) 
{
   int scene_num = scenes.size();
   for(int i = 0; i < scene_num; i++){
      scenepool.add_scene(scenes[i]);
   }

   // @Todo : Add Information to the LP
   for( Path &path: paths )
   {
      auto pathId = pathpool.size();
      IloNumColumn column(env);  
      column += obj(path.cost);
      for( int i = 1; i < (int)path.size(); i++ )
      {
         int p = path[i];
         // column +=   
      }
      y.add(IloNumVar(column, 0, param.maxValue, ILOFLOAT, ("y"+std::to_string(pathId)).c_str()));
   }

}


/**
 * @brief Use this function to pass LP solution to the Node
 * @todo callback the objective value from the solver
 * @todo callback weather the ith Scene has been used in this solution, if so pass the value to val
*/
void Master::GetSol(Node &node) 
{
   node.lpsol = new LpSol();
   node.lpsol->obj;
   for( int i = 0; i < scenepool.size(); i++ )
   {
      double val = 0; // @todo callback the 
      if( val > param.eps )
      {
         node.lpsol->push_back(std::make_pair(&scenepool[i], val));
      }
   }
}


LP_ALL_IN_ONE::LP_ALL_IN_ONE(){}

LP_ALL_IN_ONE::LP_ALL_IN_ONE(int _N_rows, int _N_cols){
   N_rows = _N_rows;
   N_cols = _N_cols;
   A = new vector<double>[N_rows];
   rhs = new double[N_rows];
   opt_dual_sol = new double[N_rows];


   cost = vector<double>(N_cols, 0);
   lb = vector<double>(N_cols, 0);
   up = vector<double>(N_cols, 0);

   for(int i = 0; i < N_rows; i++){
      A[i] = vector<double>(N_cols, 0);
      rhs[i] = 0;
      opt_dual_sol = 0;
   }
   opt_value = 0;
   LP_STATE = 0;
}

LP_ALL_IN_ONE::~LP_ALL_IN_ONE(){

   if(N_rows > 0){
      delete[] rhs;
      delete[] A;
      delete[] opt_dual_sol;
   }
}

double* LP_ALL_IN_ONE::GET_DUAL_VALUE(){
   if(opt_dual_sol){
      return opt_dual_sol;
   }
   else{
      printf("LP haven't been solved \n");
      return NULL;
   }
}