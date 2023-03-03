#ifndef _ALGO_PARAMETER_H__
#define _ALGO_PARAMETER_H__

#include <string>
#include <random>
#include <iostream>
#include <limits>


class AlgoParameter
{
public:
    // 0 for release without debug
    // 5 for test root node
    int DEBUG_VERSION;

    int seed;
    std::default_random_engine rnd;
    double eps;

    std::string problem_name;
    std::string algo_name;
    int time_limit;
    std::string test_prefix;
    std::string test_extension;
    std::string path_data;
    std::string path_result_sol;
    std::string path_result_csv;

    bool rootOnly;
    bool enableCuts;
    bool depthFirst;
    bool enableBranchOnSum;
    bool enableCplexLog;

    long maxValue = std::numeric_limits<long>::max();
    double optimalityGap = 1.0e-6;

public:
   AlgoParameter(){}
    AlgoParameter(std::string _prob, std::string _algo, int _tlim, std::string _prefix, std::string _ext,
                    std::string _data,
                    bool _rootOnly, bool _enableCuts, bool _depthFirst,
                    bool _enableBranchOnSum = false, bool _enableCplexLog = false, bool _DEBUG_VERSION = 0)
          : seed(3), eps(1.0e-6),
            problem_name(_prob), algo_name(_algo), time_limit(_tlim),
            test_prefix(_prefix), test_extension(_ext), path_data(_data),
            rootOnly(_rootOnly), enableCuts(_enableCuts), depthFirst(_depthFirst),
            enableBranchOnSum(_enableBranchOnSum),
            enableCplexLog(_enableCplexLog), DEBUG_VERSION(_DEBUG_VERSION)
      {
         std::srand(seed); // ensure the same random sequence for each run
         rnd = std::default_random_engine(seed);
      }

      ~AlgoParameter() = default;

      std::string csv_name()
      {
         return "result_" + algo_name + ".csv";
      }
};


#endif