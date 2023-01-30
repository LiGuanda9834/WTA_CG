/**
 * @file BranchAndCut.cpp
 * @author Zhou Xiao 
 * @brief 
 * @date 2022-11-07
 * @copyright Copyright (c) 2022
 */

#include <cassert>
#include <iostream>
#include "BranchAndCut.h"
#include "print.h"


BranchAndCut::BranchAndCut(
   WTA* _wta,
   const AlgoParameter &_parameter
) : wta(_wta), parameter(_parameter), master(wta, parameter), pricing(wta), tree(parameter.depthFirst)
{
   globalLb = -parameter.maxValue;
   globalUb = parameter.maxValue;
   nodeCnt = 0;
   cutsCnt = 0;
   timeOnHeur = 0.0;

   dual = {};
   
   col = NULL;
   row = NULL;
   incumbent = NULL;

   timer = compute::Timer();
   rootLb = -parameter.maxValue;
   rootUb = parameter.maxValue;
}

// 执行函数
void BranchAndCut::Run() {
   printf(" --- Start Initialize the Root Note ---\n");
   /* create and Initialize the root node*/
   Node root;
   /* the initial status of root is unsolved */
   root.status = NodeStatus::Unsolved; 
   root.lpsol = NULL;
   root.lb = 0;
   
   printf(" --- Start solve the Root Note ---\n");
   // Solve the LP of the root Node
   SolveRootNode(root);

   return;
   /* solve the leaves of branch and bound tree */
   if( !parameter.rootOnly )
   {
      while( !tree.empty() && !TimeLimit() && !Optimal() )
      {
         Node node = tree.next();
         tree.pop();
         Solve(node);
         if( node.status == NodeStatus::Integral || node.status == NodeStatus::Fractional )
         {
            globalLb = tree.updateGlobalLb();
         }
         if( node.status == NodeStatus::Fractional )
         {
            Branch(node);
         }
      }
      if( tree.empty() )
      {
         globalLb = globalUb;
      }
   }
   // 输出信息函数
   OutputPerformanceMetrics();
}

/**
 * @brief Solve root node 
 * @param node
 * @return node status
 */
void BranchAndCut::SolveRootNode(Node &node)
{
   info::print(debuginfo, "enter SolveRootNode");
   master.Set(node);
   pricing.Set(node);
   bool feasible = InitialColumns(node); 
   //infeasible: 
   master.scenepool.print_all_scene();
   return;
   if( feasible )
   {
      /* if root is feasible, then get the optimal solution of LP relaxtion of root node */
      ColumnGeneration(node);
      master.GetSol(node);
      node.lb = node.lpsol->obj;
   }
   else //infeasible
   {
      /* otherwise, if the root node is infeasible, then ... */
   }

   if( node.status == NodeStatus::Fractional )
   {
      // Branch();
   }

}

void BranchAndCut::Solve(Node &node)
{

}

/**
 * @brief master lp initial columns
 * @param : root node 
 * @return index of columns
 * @todo add cols, add rows((5b) and part of conflict cons[lazy constraint])
 */
bool BranchAndCut::InitialColumns(
   Node          &node
)
{  
   // pricing.Solve(/*no dual variables*/);
   master.AddCol(pricing.scenes);
   return true;
}


bool BranchAndCut::ColumnGeneration(
   Node          &node
) 
{
   if( !master.Solve() )
      return false;

   master.GetDualValues(dual);
   // pricing.Solve(mcfdr, conflict, dual);

   /* add column(s) gradually to master problem */
   while (!pricing.scenes.empty())  /* while routes empty, stop*/
   {
      /* Add path(s) from the candidate pool */
      master.AddCol(pricing.scenes);

      master.Solve();

      /* Get the value of dual variables of master problem */
      master.GetDualValues(dual);
      
      pricing.Solve(dual);
   
   }

   return true;
}

bool BranchAndCut::BigM(Node &node)
{
   master.AddSlackToCardinality();
   bool feasible = master.Solve();
   if( !feasible )
   {}
   else
   {
      double slackValue = master.GetSlackValue();
      master.GetDualValues(dual);
      pricing.Solve(dual);
      while( !pricing.scenes.empty() && slackValue > parameter.eps )
      {
         master.AddCol(pricing.scenes);
         master.Solve();
         slackValue = master.GetSlackValue();
         master.GetDualValues(dual);
         pricing.Solve(dual);
      }
   }
   master.DelSlack();
   return feasible;
}

void BranchAndCut::Branch(Node &parent)
{
   if( parameter.enableBranchOnSum )
   {
      double total = parent.lpsol->accumulate();
      if( !compute::integral(total, parameter.eps) )
      {
         int floor = int(std::floor(total));
         // tree.emplace(NULL, ++nodeCnt, true, )
      }
      
   }
}

void BranchAndCut::Print(Node &node, bool diving) const 
{
   string tag = diving ? "dive" : "Node";
   info::print_tab(tag, tree.size(), nodeCnt, node.statusStr(), node.lb, globalLb, globalUb, Gap(),
                     timer.elapsed(), master.time, master.numLp, master.numCol());
}

void BranchAndCut::OutputPerformanceMetrics() {
   // cout << "Porblem Name         : " << mcfdr->probname << endl;
   std::cout << "Parameter Infomation : " << std::endl;
   std::cout << "   Algorithm Name    : " << parameter.algo_name << std::endl;
   std::cout << "   Time Limit        : " << parameter.time_limit << std::endl;
   std::cout << "   Is Root Only      : " << parameter.rootOnly << std::endl;
   std::cout << "   Is Enable Cuts    : " << parameter.enableCuts << std::endl;
   std::cout << "Root Gap(%)          : " << 100 * (rootUb - rootLb) / rootLb << std::endl;
   // cout << "Gap Closed(%)        : " << 100 * (rootLb - FirstLb) / (objValue - FirstLb) << std::endl;
   std::cout << "Global Gap(%)        : " << Gap() << std::endl;
   std::cout << "Time                 : " << timer.elapsed() << std::endl;
   std::cout << "Time On Master       : " << master.time << std::endl;
   std::cout << "Node Count           : " << nodeCnt - tree.size() << std::endl;
   std::cout << "Number of Columns In Master : " << master.numCol() << std::endl;
   std::cout << "Current Date Time    : " << timer.getCurrentDateTime() << std::endl;
}