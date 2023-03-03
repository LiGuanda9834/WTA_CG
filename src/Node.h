#ifndef MCFDR_NODE_H
#define MCFDR_NODE_H

#include <queue>
#include "Scene.h"


enum NodeStatus {
    Unsolved, Infeasible, PrunedByBound, Fractional, Integral
};

// This class to take down a spcific lp solution. 
// Still need to update to the WTA form
class LpSol : public std::vector<std::pair<Scene *, double>>
{
   public:
      double obj;

   public:
      LpSol() = default;
      ~LpSol(){};

      double accumulate()
      {
         double total = 0;
         for (auto &p : *this)
         {
            total += p.second;
         }
         return total;
      }

      double cost(int x, int y)
      {
         
      }
      //  bool integerFeasible() {
      //      for (auto v: *this) {
      //          if (!huq::integral(v.second)) return false;
      //      }
      //      return true;
      //  }
};

class Node
{
   public:
      Node() = default;

      Node(const Node *_parent, long _nodeId, bool branchOnNum, int a, int b, int c = -1)
          : parent(_parent), nodeId(_nodeId),
            lbNumVehicle(branchOnNum ? a : _parent->lbNumVehicle),
            ubNumVehicle(branchOnNum ? b : _parent->ubNumVehicle),
            forbidden(_parent->forbidden),
            status(NodeStatus::Unsolved), lpsol(nullptr),
            lb(_parent == nullptr ? 0 : _parent->lb)
      {
         if (!branchOnNum)
         {
            if (c == 0)
            {
               forbidden[a][b] = true;
            }
            else if (c == 1)
            {
               auto nVertex = forbidden.size();
               for (int i = 0; i < nVertex; ++i)
               {
                  if (i != a && !forbidden[i][b])
                  {
                     forbidden[i][b] = true;
                  }
                  if (i != b && !forbidden[a][i])
                  {
                     forbidden[a][i] = true;
                  }
               }
            }
            else
            {
               std::cerr << "invalid branch\r\n";
            }
         }
      }

      ~Node()
      {
         parent = nullptr;
         if (lpsol != nullptr)
         {
            delete lpsol;
            lpsol = NULL;
         }
      }

      std::string statusStr() const
      {
         if (status == NodeStatus::Unsolved)
            return "Unsolved";
         else if (status == NodeStatus::Infeasible)
            return "Infeasible";
         else if (status == NodeStatus::PrunedByBound)
            return "Pruned";
         else if (status == NodeStatus::Fractional)
            return "Fractional";
         else if (status == NodeStatus::Integral)
            return "Integral";
         else
            return "Unknown";
      }
      // Use this tho check 
      bool valid(const Scene &temp_scene) const
      {
         for (int i = 0; i < temp_scene.weapon_indices.size(); ++i)
         {
            int target_id = temp_scene.target_index, weapon_id = temp_scene.weapon_indices[i];
            if (forbidden[target_id][weapon_id]) {
                return false;
            }
        }
        return true;
      }

   public:
      const Node *parent;
      long nodeId;

      int lbNumVehicle;
      int ubNumVehicle;

      // Use this to checkout weather a weapon could attack a target
      // forbiddeb[target][weapon] == 0 means weapon could not reach target
      std::vector<std::vector<bool>> forbidden;

      NodeStatus status;
      LpSol *lpsol;
      double lb;
};

#endif //MCFDR_NODE_H
