#include "Scene.h"
#include <functional>

Scene::Scene(int _j, std::vector<int> _scene){
    target_index = _j;
    weapon_indices = std::vector<int>();
    for(int i = 0; i < _scene.size(); i++){
        weapon_indices.push_back(_scene[i]);
    }
}

void Scene::PrintScene(){
    printf("Target %d : ", target_index);
    for(int i = 0; i < weapon_indices.size(); i++){
        printf("%d ", weapon_indices[i]);
    }
    printf("\n");
}

void ScenePool::print_all_scene(){
    printf("There are \t %d Scenes in the Scene Pool. \n", Scene_num);

    for_each(this->begin(), this->end(), std::mem_fun_ref(&Scene::PrintScene));
}

void ScenePool::add_scene(Scene scene){
      Scene_num++;
      push_back(scene);
}