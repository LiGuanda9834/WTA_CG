#ifndef _SCENE_H__
#define _SCENE_H__

#include "print.h"
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>

/**
 * @brief Each Scene meaning a set of weapons
 * 
*/
class Scene{
   public:
                                 // name of Scene. @todo Is that neccessary
      int index;                 // index of a column @todo how to give a id to a columns
      int target_index;          // y           : the index of the target
      double obj_qjs;            // q_js        : objective coff q_js
      double redcost;            // sum(xu)+v-q : redust cost of the specific scene
      std::vector<int> weapon_indices;  // s           : all the weapons in the Scene
   public:
      Scene() = default;
      Scene(int _j, std::vector<int> _scene);
      ~Scene() = default;
      

      inline bool emptyScene() const
      {
         return weapon_indices.size() < 1;
      }

      void PrintScene();
};



class ScenePool : public std::vector<Scene> {
public:
   int Scene_num;

   ScenePool(){Scene_num = 0;}
   ~ScenePool() = default;

   void print_all_scene();
   void add_scene(Scene scene);


   // @todo print scenes by target
   void print_scene_by_target(){}
};

#endif //ORA_PATH_H
