// LP relaxation of subproblem at each node in the branch-and-bound tree

#ifndef _MASTER_H_
#define _MASTER_H_

// @todo Need import HIGHS instead of cplex
#include "WTA.h"
#include "AlgorithmParameter.h"
#include "Node.h"

/**
 * @brief This class save the information of the LP
 * @todo It will be repalced by the LP Solver like HIGHS

*/
class LP_ALL_IN_ONE{
public:
   int               N_rows;
   int               N_cols;

   vector<double>*   A;
   double*           rhs;
   vector<double>    cost;
   vector<double>    lb;
   vector<double>    up;

   vector<double>    opt_sol;
   double*           opt_dual_sol;
   double            opt_value;

   int               LP_STATE;

   LP_ALL_IN_ONE();
   LP_ALL_IN_ONE(int _N_rows, int _N_cols);

   void INITIALIZE_LP_BY_MASTER(){}
   void SOLVE_LP(){}
   double* GET_DUAL_VALUE();

   ~LP_ALL_IN_ONE();
};



class Master {
   public:
      Master();
      Master(WTA* wta, const AlgoParameter &param);

      ~Master()
      {      }


      void Set(const Node &node);
      bool Solve();
      void GetDualValues(std::vector<double> &dual);
      void AddCol(std::vector<Scene> &scenes);
      void AddSlackToCardinality(){}
      void GetSol(Node &node);

      // Delete a variable ??   
      void DelSlack(){}


      // Get a variable value ??
      double GetSlackValue(){}

      // Writes the active model to the file specified by filename.
      void writeModel(const std::string &name = "master.lp"){}

       int numCol() const {
        return scenepool.size();
       }

   public:
      WTA*                 wta;        // @brief The wta model
      const AlgoParameter  &param;     // @brief the Param of the algorithm
      ScenePool            scenepool;  // @brief Use this to takedown all the weapon scnen
      long                 numLp;      // @brief Use this to calculate the number of LP
      double               time;       // @brief take down the time of the master problem

   public:
   /**
    * @brief Use this member varialbs to solve LP
    * @todo Replace the following member varible to HIGHS
   */
      LP_ALL_IN_ONE*       lp;
  /*
      IloEnv env;
      IloModel model;
      IloCplex cplex;
      IloObjective obj;
      IloRangeArray rng;
      IloRange *ptrNumVehicle;
      IloNumVarArray y;
      IloNumVar slack;
  */
};


#endif //_MASTER_H_
