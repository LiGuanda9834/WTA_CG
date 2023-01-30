#include "BranchAndCut.h"
#include "WTA.h"
#include "print.h"
#include "Pricing.h"
#include "Scene.h"
#include <string>



int main(int argc, char *argv[]) {
   info::print(true, "AAAA");
   
   // 保存数据
   //设置参数

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
   
   printf("\n ----- start anylisis the structure ----- \n \n");
   /* use branch-and-price algorithm to solve mcfdr problem */
   WTA* wta = new WTA(4,6,0);
   wta->print_model();
   
   Master test_master(wta, parameter);
   Pricing test_pricing(wta);
   Tree test_tree(1);

   BranchAndCut test_BAC(wta, parameter);
   
   std::vector<int> _scene = {1,2,3,4};
   Scene test_scene(4, _scene);
   
   ScenePool test_scene_pool;

   Node test_node;
   test_pricing.Set(test_node);
   int num_scene = test_pricing.scenes.size();
   

   test_BAC.Run();
   /*
   int _size1 = test_BAC.globalLb;
   int isPrimal = test_BAC.tree.searchStrategy;
   printf("%s \t: %d\n ", "Is Primal?", isPrimal);
   */

   printf("\n%s \t: %.2f\n ", "test_size", 1.2);

   
   

   printf("\n ----- anylisis the structure finished -----\n");
   /**
    * @todo learn what those parameter means and I can use whitch parameter
    * @todo need node.h
   */
   //BranchAndCut algorithm(mcfdr, parameter);
   //algorithm.Run();
  
   delete wta;
   return 0;
}