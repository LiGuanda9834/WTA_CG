#include "Master.h"

#include "Highs.h"
#include <cassert>

using std::cout;
using std::endl;

void LP_ALL_IN_ONE::INITIALIZE_LP_BY_MASTER(WTA* wta, ScenePool scene_pool){
   
   printf("Now Initlized the LP my master problem\n");
   int m = wta->weapon_num_m;
   int n = wta->target_num_n;
   N_rows = wta->weapon_num_m + wta->target_num_n;
   N_cols = scene_pool.size();

   if(rhs){
      delete[] rhs;
      rhs = NULL;
   }
   if(A){
      delete[] A;
      A = NULL;
   }
   if(opt_dual_sol){
      delete[] opt_dual_sol;
      opt_dual_sol = NULL;
   }

   A = new vector<double>[N_rows];
   rhs = new double[N_rows];
   opt_dual_sol = new double[N_rows];
   opt_sol = vector<double>(N_cols, 0);
   

   cost = vector<double>(N_cols, 0);

   // Init cost
   for(int i = 0; i < N_cols; i++){
      Scene temp_scene = scene_pool[i];
      cost[i] = wta->set_scene_qjs(temp_scene);
   }

   // Init A, weapon row coef, change from column to row
   for(int i = 0; i < N_cols; i++){
      Scene temp_scene = scene_pool[i];
      vector<int> temp_weapon_indices = temp_scene.weapon_indices;
      int nz_col_num = temp_weapon_indices.size();
      for(int j = 0; j < nz_col_num; j++){
         int temp_weapon_index = temp_weapon_indices[j];
         A[temp_weapon_index].push_back(i);
      }
      // Init A, target row coef
      int target_index = temp_scene.target_index;
      int row_index = target_index + m;
      A[row_index].push_back(i);
   }

   

   lb = vector<double>(N_cols, 0);
   ub = vector<double>(N_cols, 0);

   for(int i = 0; i < N_rows; i++){
      rhs[i] = 1;
      opt_dual_sol[i] = 0;
   }
   opt_value = 0;
   LP_STATE = 0;
}

void LP_ALL_IN_ONE::PRINT_LP_INFO(){
   printf("Now print the LP info\n");
   printf("LP has %d rows and %d columns\n", N_rows, N_cols);
   printf("c:\n");
   for(int i = 0; i < N_cols; i++){
      printf("%.2f\t", cost[i]);
   }
   // print Ax \leq b
   printf("\n A: \n");
   for(int i = 0; i < N_rows; i++){
      int nz_num = A[i].size();
      printf("r%d:\t", i);
      if(nz_num == 0){
         for(int j = 0; j < N_cols; j++){
               printf("0\t");
         } 
      }
      else{
         for(int j = 0; j < N_cols; j++){
            int temp_nz_id = 0;
            int temp_nz_weapon = A[i][temp_nz_id];
            if(temp_nz_weapon == j){
               printf("1\t");
               temp_nz_id++;
            }
            else{
               printf("0\t");
            }
         }
      }
      
      printf("<=\t %.2f\n", rhs[i]);
   }

   // print lb and ub
   printf("lb & ub\n");
   printf("size: %d - %d\n", lb.size(), ub.size());

   for(int i = 0; i < N_cols; i++){
      printf("%d'th scene:\t%.2f - %.2f\n", i, lb[i], ub[i]);

   }

   printf("Opt Sol:");
   printf("size: %d\n", opt_sol.size());
   for(int i = 0; i < N_cols; i++){
      printf(" %.2f", opt_sol[i]);
   }
   printf("\n");
   printf("dual variable\n");
   for(int i = 0; i < N_rows; i++){
      printf("%.2f\n", opt_dual_sol[i]);
   }
   printf("\n");
}

void LP_ALL_IN_ONE::Delete(){
   if(A){
      delete[] A;
      A = NULL;
   }
   if(rhs){
      delete[] rhs;
      rhs = NULL;
   }
   if(opt_dual_sol){
      delete[] opt_dual_sol;
      opt_dual_sol = NULL;
   }
}

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
   /**
    * @brief Set the upper bound and the lower bound of all the scenes
    * @todo Is this step neccessary
   */



   for( int i = 0; i < (int)scenepool.size(); i++ )
   {
      const Scene &scene = scenepool[i];
      if( node.valid(scene) )
      {
         if( lp->lb[i] != 0 || lp->ub[i] != 1 )
         {
            lp->lb[i] = 0;
            lp->ub[i] = 1;
         }
      }
      else
      {
         if( lp->lb[i] != 0 || lp->ub[i] != 0 )
         {
            lp->lb[i] = 0;
            lp->ub[i] = 0;
         }
      }
   }

}

void Master::Initialize_HIGHS(){
  printf("\n----- Now Initialize master's HIGHS --- \n");

   printf("Now the Scene Pool has : %d\n", scenepool.size());
  int DEBUG = 0;

  // Although the first constraint could be expressed as an upper
  // bound on x_1, it serves to illustrate a non-trivial packed
  // column-wise matrix.
  //

   int m = wta->weapon_num_m;
   int n = wta->target_num_n;
   int N_cols = scenepool.size();
   int N_rows = m + n;

  model.lp_.num_col_ = N_cols;
  model.lp_.num_row_ = m + n;
  model.lp_.sense_ = ObjSense::kMinimize;
  model.lp_.offset_ = 0;
  
   vector<double> cost(N_cols, 0);

   // Init cost
   for(int i = 0; i < N_cols; i++){
      Scene temp_scene = scenepool[i];
      cost[i] = wta->set_scene_qjs(temp_scene);
   }

model.lp_.col_cost_ = cost;


// Init lb and ub of LP
  model.lp_.col_lower_ = vector<double>(N_cols, 0.0);
  model.lp_.col_upper_ = vector<double>(N_cols, 1.0);
  model.lp_.row_lower_ = vector<double>(N_rows, 1.0);
  model.lp_.row_upper_ = vector<double>(N_rows, 1.0);

for(int i = 0; i < m; i++){
   model.lp_.row_lower_[i] = 0.0;
}
// Init A, weapon row coef, change from column to row
// Init A, pass each cols to column_compress data structure
vector<int>       start_index;
vector<int>       nz_rows;
vector<double>    nz_vals;

int               nz_counter = 0;

   for(int i = 0; i < N_cols; i++){
      start_index.push_back(nz_counter);
      Scene temp_scene = scenepool[i];
      vector<int> temp_weapon_indices = temp_scene.weapon_indices;
      int nz_row_num = temp_weapon_indices.size();
      for(int j = 0; j < nz_row_num; j++){
         nz_rows.push_back(temp_weapon_indices[j]);
         nz_vals.push_back(1.0);
      }
      nz_counter += nz_row_num;

      // Init A, target row coef
      int target_index = temp_scene.target_index;
      int target_row_index = target_index + m;
      nz_rows.push_back(target_row_index);
      nz_vals.push_back(1.0);
      nz_counter++;
   }
   start_index.push_back(nz_counter);
  // Here the orientation of the matrix is column-wise
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  // a_start_ has num_col_1 entries, and the last entry is the number
  // of nonzeros in A, allowing the number of nonzeros in the last
  // column to be defined
  model.lp_.a_matrix_.start_ = start_index;
  model.lp_.a_matrix_.index_ = nz_rows;
  model.lp_.a_matrix_.value_ = nz_vals;

   printf("N_rows : %d", m + n);
   printf("N_cols : %d",  N_cols);
   int real_nz_num = nz_vals.size();
   if(real_nz_num == nz_counter){
      printf("column index : ");
      for(int i = 0; i < start_index.size(); i++){
         printf(" %d ", start_index[i]);
      }
      printf("\nnz index : ");
      for(int i = 0; i < nz_rows.size(); i++){
         printf("\t%d", nz_rows[i]);
      }
      printf("\nnz val : ");
      for(int i = 0; i < nz_vals.size(); i++){
         printf("\t%.2f", nz_vals[i]);
      }
      printf("\n");
   }
   else{
      printf("!!!!! counter has problem\n");
   }
   return_status = highs.passModel(model);
   assert(return_status==HighsStatus::kOk);
   printf("\n---- Highs of Master Problem has been Initialized ---- \n");
}



Master::~Master(){
   if(lp){
      lp->Delete();
      delete lp;
      lp = NULL;
   }
}
void LP_ALL_IN_ONE::SOLVE_LP(){
   printf("Now solve LP \n");
   
   //simplest 
   int is_simplest = 1;
   if(is_simplest == 1){
      printf("use this to calculate the simplest version\n");
      printf("Alert!!! Only correct when w_num = t_num\n");
      int m = N_rows / 2;
      int n = m;
      for(int i = 0; i < m; i++){
         opt_dual_sol[i] = 0;
      }
      for(int i = 0; i < n; i++){
         opt_dual_sol[i + m] = cost[i];
      }

   }
}
/**
 * @brief Use this to solve the LP by solver like HIGHS
 * @todo need to update the result
*/
bool Master::Solve() 
{
   opt_scene_indices = {};
   numLp++;
   printf("Now solve Mater Problem\n");
   printf("scene pool size : %d\n", scenepool.size());
   //lp->INITIALIZE_LP_BY_MASTER(this->wta, this->scenepool);
   //lp->SOLVE_LP();
   //lp->PRINT_LP_INFO();
   bool result = true;
   printf("------------- Start Highs ---------------------!!!!!!!!!!!!!!!!\n\n");
   /*
  HighsModel model;
  model.lp_.num_col_ = 2;
  model.lp_.num_row_ = 3;
  model.lp_.sense_ = ObjSense::kMinimize;
  model.lp_.offset_ = 3;
  model.lp_.col_cost_ = {1.0, 1.0};
  model.lp_.col_lower_ = {0.0, 1.0};
  model.lp_.col_upper_ = {4.0, 1.0e30};
  model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
  model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};
  //
  // Here the orientation of the matrix is column-wise
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  // a_start_ has num_col_1 entries, and the last entry is the number
  // of nonzeros in A, allowing the number of nonzeros in the last
  // column to be defined
  model.lp_.a_matrix_.start_ = {0, 2, 5};
  model.lp_.a_matrix_.index_ = {1, 2, 0, 1, 2};
  model.lp_.a_matrix_.value_ = {1.0, 3.0, 1.0, 2.0, 2.0};
  //

  // Create a Highs instance
  Highs highs;
  HighsStatus return_status;
  //
  // Pass the model to HiGHS
  */


   // This flag is used to check weather the LP Solution is integer, if not, throw an alert 
   bool is_integer = 1;


  printf("\n --- now write the incumbent LP  --- \n");
  return_status = highs.writeModel("wta.lp");

  //
  // Get a const reference to the LP data in HiGHS
  const HighsLp& lp = highs.getLp();
  //
  // Solve the model
  return_status = highs.run();
  assert(return_status==HighsStatus::kOk);
  //
  // Get the model status
  const HighsModelStatus& model_status = highs.getModelStatus();
  assert(model_status==HighsModelStatus::kOptimal);
  cout << "Model status: " << highs.modelStatusToString(model_status) << endl;
  //
  // Get the solution information
  const HighsInfo& info = highs.getInfo();
  cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
  cout << "Objective function value: " << info.objective_function_value << endl;
  cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << endl;
  cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << endl;
  cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;
  const bool has_values = info.primal_solution_status;
  const bool has_duals = info.dual_solution_status;
  const bool has_basis = info.basis_validity;
  //
  // Get the solution values and basis
  const HighsSolution& solution = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  //
  // Report the primal and solution values and basis
  for (int col=0; col < lp.num_col_; col++) {
    cout << "Column " << col;
    if (has_values) cout << "; value = " << solution.col_value[col];
    if(solution.col_value[col] > 1e-6 && solution.col_value[col] < 1 - 1e-6){
      printf("Alert!!!!!! We find a non-integer LP solution in Master Problem !!!!!!!\n");
    }
    if( solution.col_value[col] > 1 - 1e-6){
      opt_scene_indices.push_back(col);
    }
    if (has_duals) cout << "; dual = " << solution.col_dual[col];
    if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.col_status[col]);
    cout << endl;
  }
  for (int row=0; row < lp.num_row_; row++) {
    cout << "Row    " << row;
    if (has_values) cout << "; value = " << solution.row_value[row];
    if (has_duals) cout << "; dual = " << solution.row_dual[row];
    if (has_basis) cout << "; status: " << highs.basisStatusToString(basis.row_status[row]);
    cout << endl;
  }
  



   /*
---------- This Blocl change the LP to MIP, but I think this is not Neccessary and make dual value to 0,
           So I Note this block 

  // Now indicate that all the variables must take integer values
  model.lp_.integrality_.resize(lp.num_col_);
  for (int col=0; col < lp.num_col_; col++)
    model.lp_.integrality_[col] = HighsVarType::kInteger;

  highs.passModel(model);
  // Solve the model
  return_status = highs.run();
  assert(return_status==HighsStatus::kOk);
  // Report the primal solution values
  for (int col=0; col < lp.num_col_; col++) {
    cout << "Column " << col;
    if (info.primal_solution_status) cout << "; value = " << solution.col_value[col];
    cout << endl;
  }
  for (int row=0; row < lp.num_row_; row++) {
    cout << "Row    " << row;
    if (info.primal_solution_status) cout << "; value = " << solution.row_value[row];
    cout << endl;
  }
   const HighsLp& lp_ = highs.getLp();
   const HighsSolution& solution_ = highs.getSolution();
   int N_Dual_sol = lp_.num_row_;
   printf("Print Dual after MIP\n");
   for(int i = 0; i < N_Dual_sol; i++){
      printf("dual %d : %.2f\n", i, solution.row_dual[i]);
   }
   */
   printf("------------- End Highs ---------------------\n\n");
   return result;
}

// get all the Dual_Values from the lp, m + n
void Master::GetDualValues(vector<double> &dual) 
{
   printf("Now get dual sols from lp to master \n");
   const HighsLp& lp = highs.getLp();
   const HighsSolution& solution = highs.getSolution();

   int N_Dual_sol =  lp.num_row_;
   printf("size : %d\n", solution.row_dual.size());
   printf("dual size : %d\n", dual.size());
   dual = vector<double>(N_Dual_sol, 0);
   for(int i = 0; i < N_Dual_sol; i++){
      printf("dual %d : %.2f\n", i, solution.row_dual[i]);
      dual[i] = solution.row_dual[i];
   }
   return;
}


/**
 * @brief Use this to add scenes to the current master RP's ScenePool
 * @note  Do not delete scenes from the scenes Pool
*/
void Master::AddCol(vector<Scene> &scenes) 
{
   int target_num = wta->target_num_n;
   int scene_num = scenes.size();
   for(int i = 0; i < scene_num; i++){
      double lower_ = 0.0;
      double upper_ = 1.0;
      double cost_ = wta->set_scene_qjs(scenes[i]);
      int nz_num_ = 1 + scenes[i].weapon_indices.size();
      int* indices_ = new int[nz_num_];
      double* values_ = new double[nz_num_];
      for(int w = 0; w < nz_num_ - 1; w++){
         indices_[w] = scenes[i].weapon_indices[w];
         values_[w] = 1;
      }
      indices_[nz_num_ - 1] = wta->weapon_num_m + scenes[i].target_index;
      values_[nz_num_ - 1] = 1;
      
      //int Scene_id = scenepool.size();// Is this neccessary?
      //lp->ADD_NEW_COLUMN(scenes[i], Scene_id);
      highs.addCol(cost_, lower_, upper_, nz_num_, indices_, values_);

      scenepool.add_scene(scenes[i]);
      delete[] indices_;
      delete[] values_;
   }
   return_status = highs.writeModel("New_Master.lp");
   printf("Now Print all the scenes in the scene pool \n");
   //scenepool.print_all_scene();
   for(int line = 0; line <50; line++){
      printf("\n");
   }
   scenepool.print_scene_by_target(target_num);
   for(int line = 0; line <50; line++){
      printf("\n");
   }
   // @Todo : Add Information to the LP
}

bool Master::Check_is_scenes_new(std::vector<Scene> &scenes){
   if(scenes.size() == 1){
      printf("No New Scenes \n");
      return false;
   }
   bool scenes_is_new = false;
   int target_num = scenes.size();
   int last_scene_index = scenepool.size();
   int first_scene_index = last_scene_index - target_num;
   for(int i = 0; i < target_num; i++){
      Scene old_scene = scenepool[first_scene_index + i];
      Scene new_scene = scenes[i];
      bool is_same_i = old_scene.is_same_w_t(new_scene);
      if(!is_same_i){
         scenes_is_new = true;
      }
   }
   if(!scenes_is_new){
      printf("We don't find new scene\n");
   }
   return scenes_is_new;
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
      double val = 0; // @todo callback the ??
      if( val > param.eps )
      {
         node.lpsol->push_back(std::make_pair(&scenepool[i], val));
      }
   }
}



LP_ALL_IN_ONE::LP_ALL_IN_ONE(int _N_rows, int _N_cols){
   N_rows = _N_rows;
   N_cols = _N_cols;
   A = new vector<double>[N_rows];
   rhs = new double[N_rows];
   opt_dual_sol = new double[N_rows + N_cols];


   cost = vector<double>(N_cols, 0);
   lb = vector<double>(N_cols, 0);
   ub = vector<double>(N_cols, 0);

   for(int i = 0; i < N_rows; i++){
      rhs[i] = 0;
      opt_dual_sol[i] = 0;
   }
   for(int i = 0; i < N_cols; i++){
      opt_dual_sol[N_rows + i] = 0;
   }
   opt_value = 0;
   LP_STATE = 0;
}

LP_ALL_IN_ONE::~LP_ALL_IN_ONE(){
   if(rhs){
      delete[] rhs;
      rhs = NULL;
   }
   if(A){
      delete[] A;
      A = NULL;
   }
   if(opt_dual_sol){
      delete[] opt_dual_sol;
      opt_dual_sol = NULL;
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

double Master::cal_obj_val(){
   double obj_value = 0;
   printf("the opt_sol is : \n");
   for(int i = 0; i < opt_scene_indices.size(); i++){
      double target_value = wta->cal_target_value(scenepool[opt_scene_indices[i]]);
      printf("value of target %d : %.2f\n", i, target_value);
      obj_value += target_value;
      scenepool[opt_scene_indices[i]].PrintScene();
   }
   printf("optimal value is %.4f\n", obj_value);
   return obj_value;
}