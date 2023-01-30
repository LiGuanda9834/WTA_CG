#ifndef _PRICING_H__
#define _PRICING_H__

#include "WTA.h"
#include "Node.h"


/**
 * @brief use this class to takedown the state of the Pricing 
*/
class Label {
public:
    Label() = default;
    Label(const Label & _parent, int j);
    ~Label() = default;
public:
    Label *parent;
    int cur;
    int serviceStartTime;
    int load;
    double cost;
    std::vector<int> onboardOrders;
    std::vector<bool> unreachableOrders; // inst.n+1

    bool pruned; // prune by bound or is dominated
};

// auto CompareLabelAscending = [](const Label *a, const Label *b) 
// {
//     return a->cost > b->cost;
// };

class Pricing {
public:
   Pricing();
   Pricing(WTA* _wta);
   ~Pricing() = default;

   void Set(Node &node);
   void Solve(vector<double> &dual);
   // void extend(const Label &parent);
   // bool dominance(int j);

public:
   WTA* wta;

/**
 * @category 1 : initalize, find a feasible solution
 * @category 2 : in CG
 * genereted scenes with negative reduced cost, 
 * if empty, means the Subproblem has reach the optimal 
*/
   std::vector<Scene> scenes;

   Node *ptrNode;

   std::vector<std::vector<Label>> states; // @todo What is this?
};


#endif //ORA_PRICING_H
