#ifndef ORA_SOLUTION_H
#define ORA_SOLUTION_H

#include "Scene.h"

class Solution : public vector<Scene> 
{
   public:
      int cost;
   public:
      Solution() : cost(0) {

    }

      void calculate() 
      {
         cost = 0;
         // for (Path &path: *this)
            // cost += path.cost;
      }

};

#endif //ORA_SOLUTION_H
