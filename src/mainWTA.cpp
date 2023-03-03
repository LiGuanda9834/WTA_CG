#include "BranchAndCut.h"
#include "WTA.h"
#include "print.h"
#include "Pricing.h"
#include "Scene.h"
#include "SubProblem.h"
#include <string>
#include <ctime>


int main(int argc, char *argv[]) {
   //info::print(true, "Start main function");
   printf("Start main function\n");
   // 保存数据
   // 设置参数
   int _temp_testsize = 0;
   string Name_1 = "/share/home/liguanda/WTA-Problem/Code/column_generation/Column_generation_structure/data/wta";
   string Name_2 = "5";
   string Name_3 = ".txt";
   if(argc > 1){
      _temp_testsize = atoi(argv[1]);
      Name_2 = argv[1];
   }
   Name_1 += Name_2;
   Name_1 += Name_3;

   /**
    * 0  for branch and cut
    * 1  for subproblem
    * -1 for short test
    * */ 

   clock_t start_time, finish_time;
   start_time = clock();

   AlgoParameter parameter(
      "WTA",              // problem name
      "BranchAndPrice",   // algorithm name
      3600,               // time limit = 3600 (seconds)
      "AA", "",           // teston_prefix, teston_extension
      "./src/mcfdr/data", // path_data
      true,               // rootOnly
      false,              // enableCuts
      false,              // depthFirst if true, otherwise it is primalBoundFirst
      true,               // branchOnSum
      false,              // enableCplexLog
      true                // debug
   );


   // 列生成
   
   // algo.Run();


   int problem_size = _temp_testsize;

   WTA* wta = new WTA(problem_size, problem_size);
   FILE* fin = NULL;
   printf("\n ----- start anylisis the structure ----- \n");
   
   /* use branch-and-price algorithm to solve mcfdr problem */


   char* test_file_name = (char*)Name_1.c_str();
   //char* test_file_name = (char*)"/share/home/liguanda/WTA-Problem/Code/column_generation/Column_generation_structure/data/wta5.txt";
   
   fin = fopen(test_file_name, "r");
   int test_size = 0;
   if(fin == NULL){
      printf("Fail to read file\n");
   }
   else{
      wta->Init_readfile(fin, 0);
      printf("Initialized finished\n");
      wta->print_model();
   }

   //printf("start read file \n");

   
   // @date 2023/02/21 进行内存泄露检测

   LP_ALL_IN_ONE test_lp(5,5);
   //test_lp.Delete();
   Master test_master(wta, parameter);
   Pricing test_pricing(wta, parameter);
   Tree test_tree(1);

   BranchAndCut test_BAC(wta, parameter);
   
   std::vector<int> _scene = {1,2,3,4};
   Scene test_scene(4, _scene);
   
   ScenePool test_scene_pool;

   //Node test_node;
   test_BAC.Run();
   /*
   test_pricing.Set(test_node);
   int num_scene = test_pricing.scenes.size();

   LP_ALL_IN_ONE lp(5,5);
   
   //lp.PRINT_LP_INFO();

*/
  // 2023/02/21 进行内存泄露检测 终止

   // ------------------ Use this to test Branch and Cut method ------------
   //

   // Use this to show the Branch and cut search strategy
   /*
   int _size1 = test_BAC.globalLb;
   int isPrimal = test_BAC.tree.searchStrategy;
   printf("%s \t: %d\n ", "Is Primal?", isPrimal);
   */

   printf("\n ----- anylisis the structure finished -----\n");
   finish_time = clock();

   double sec = double(finish_time - start_time) / double(CLOCKS_PER_SEC);
   std::cout << finish_time - start_time   << "/" << CLOCKS_PER_SEC << " = " << sec << " (s) "<< std::endl;
   int seed = 0;
   srand(seed);

   /**
    * @brief   Use this part to test the subproblem
    * @todo    figure out how to combine the subproblem and cplex, note three lines
   */
  
   printf("------ Now test the subproblem -------\n");
   int temp_m = 5;
   int temp_n = 5;
   double* temp_dual_u = new double[5];
   int target_id = 0;
   //Set simplest dual value to subproblem
   int trial_size = 0;
   if(!trial_size){
      printf("Do not test subproblem in this test \n");
   }
   int* Is_correct = NULL;
   if(trial_size){
      Is_correct = new int[trial_size];
   }
   for(int trial_num = 0; trial_num < trial_size; trial_num++){
      Is_correct[trial_num] = 1;
      for(int i = 0; i < temp_m; i++){
         int random_dual_value = rand()%5;
         temp_dual_u[i] = -random_dual_value;
      }
      double temp_dual_v = 7.9;
      SubProblem temp_sub;
      temp_sub.Init(wta, 0, temp_dual_u, temp_dual_v, 0, 1, &parameter);
      Scene temp_scene;
      vector<int> solve_by_OA = temp_sub.cal_optimal_scene(temp_scene);
      vector<int> solve_by_Enum = temp_sub.cal_optimal_scene_by_enum();
      printf("weapon_num : %d\n", temp_sub.weapon_num); 
      for(int i = 0; i < temp_sub.weapon_num; i++){
         if(solve_by_OA[i] != solve_by_Enum[i]){
            Is_correct[trial_num] = 0;
         }
      }
      temp_sub.Delete();
   }

   //temp_sub.print_debug();

   /* This block is used to test cal_constraint_by_x(int weapon_set) and find_init_cut() function
      Test is finished

   vector<vector<int>> temp_init_cut;
   temp_sub.find_init_cut(temp_init_cut, 1, 1);
   printf("size: %d \n", temp_init_cut[0].size());

   temp_sub.print_LP_info();


    
   vector<int> temp_weapon_set = {0,3,4};
   printf("weapon set: \n");
   for(int i = 0; i < temp_weapon_set.size(); i++){
      printf("%d ", temp_weapon_set[i]);
   }
   printf("\n");
   
   temp_sub.cal_constraint_by_x(temp_weapon_set);
   temp_sub.print_LP_info();
   */



  for(int i = 0; i < trial_size; i++){
      if(Is_correct[i]){
         printf("The OA method find problem %dth\tcorrect answer !!!! \n", i);
      }
      else{
         printf("The OA method find the wrong answer ????? \n");
      }
  }



   
   printf("------ subproblem test finished -------\n");
   if(Is_correct){
      delete[] Is_correct;
   }
   printf("\n%s \t: %.2f\n ", "test_size", 1.2);
   


   /**
    * @todo learn what those parameter means and I can use whitch parameter
    * @todo need node.h
   */
   //BranchAndCut algorithm(mcfdr, parameter);
   //algorithm.Run();
   if(wta){
      wta->Delete();
      delete wta;
   }
   

   if(temp_dual_u){
      delete[] temp_dual_u;
   }
   
   printf("Size : %d \n", _temp_testsize );
   std::cout << Name_1 << std::endl;
   return 0;
}