#ifndef ORA_BRANCHANDCUT_H
#define ORA_BRANCHANDCUT_H

#include "Master.h"
#include "Pricing.h"
#include "Solution.h"
#include "Tree.h"
#include "WTA.h"
#include "AlgorithmParameter.h"


class COL
{
   public:
      void PrintCol(){ printf("%.2f<%s>[%d]\n", this->obj, this->name.c_str(), this->index); };

   public:
      int index;   /* index of a column */
      string name;
      double obj;
      double redcost;
      vector<int> edge;
};

class ROW
{
   public:
      void PrintRow(){ printf("<%s>[%d]\n", this->name.c_str(), this->index); };

   public:
      int index;
      string name;
      double rhs;
      double dualsol;
      // vector<int> 
};

class BranchAndCut 
{
   public:
      BranchAndCut(WTA *wta, const AlgoParameter &parameter);

      ~BranchAndCut() = default;

      void Run();

   private:
      void SolveRootNode(Node &node);
      void Solve(Node &node);

      // @todo What is this? do I need this?
      bool InitialColumns(Node &node);
      bool ColumnGeneration(Node &node);

      // @todo what is this? In strongbanrching?
      bool BigM(Node &node);

      void Branch(Node &parent);
      
      void IncumbentHeuristic(Node &node);

      inline bool TimeLimit() const
      {
         return timer.elapsed() > parameter.time_limit;
      }

      inline bool Optimal() const
      {
         return globalLb + parameter.optimalityGap > globalUb;
      }

      inline double Gap() const
      {
         return 100.0 * (globalUb - globalLb) / globalLb;
      }

      void Print(Node &node, bool diving = false) const;

      //Print the parameters and the result of the algorithm
      void OutputPerformanceMetrics();

   public:
      WTA* wta;
      const AlgoParameter &parameter;
      Master master;
      Pricing pricing;
      Tree tree;

      // @todo Current dual sol?
      vector<double> dual; 

      /* index, edges, order of edges */
      COL *col; // path
      ROW *row;
      
      // a sulution, 0-1matrix, not suitable for WTA, need change
      int** incumbent;

      double globalUb;
      double globalLb;

      // statistics
      long nodeCnt;
      int cutsCnt;
      double timeOnHeur;

   private:
      
      compute::Timer timer;
      double rootUb;
      double rootLb;
};


#endif //ORA_BRANCHANDCUT_H
