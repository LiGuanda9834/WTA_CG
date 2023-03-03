#include "Scene.h"
#include <functional>

Scene::Scene(int _j, std::vector<int> _scene){
    target_index = _j;
    weapon_indices = std::vector<int>();
    for(int i = 0; i < _scene.size(); i++){
        weapon_indices.push_back(_scene[i]);
    }
}

void Scene::Set_Scene(int _j, std::vector<int> _scene){
    target_index = _j;
    weapon_indices = std::vector<int>();
    for(int i = 0; i < _scene.size(); i++){
        if(_scene[i] == 1){
            weapon_indices.push_back(i);
        }
    }
}

void Scene::PrintScene(){
    printf("Target %d : ", target_index);
    for(int i = 0; i < weapon_indices.size(); i++){
        printf("%d ", weapon_indices[i]);
    }
    printf("\n");
}

bool Scene::is_same_w_t(Scene compare_scene){
    bool is_same = 1;
    if(compare_scene.target_index != target_index){
        is_same = 0;
        return is_same;
    }
    if(compare_scene.weapon_indices.size() != weapon_indices.size()){
        is_same = 0;
        return is_same;  
    }
    for(int i = 0; i < weapon_indices.size(); i++){
        if(compare_scene.weapon_indices[i] != weapon_indices[i]){
            is_same = 0;
            return is_same;
        }
    }
    return is_same;
}

void ScenePool::print_all_scene(){
    printf("There are \t %d Scenes in the Scene Pool. \n", Scene_num);

    for_each(this->begin(), this->end(), std::mem_fun_ref(&Scene::PrintScene));
}

void ScenePool::add_scene(Scene scene){
      Scene_num++;
      push_back(scene);
}

void ScenePool::print_scene_by_target(int target_num){
    int scene_num = size();
    vector<Scene>* target_scenes = new vector<Scene>[target_num];
    for(int i = 0; i < scene_num; i++){
        Scene temp_scene = at(i);
        target_scenes[temp_scene.target_index].push_back(temp_scene);
    }
    for(int t = 0; t < target_num; t++){
        printf("Target : %d, Size : %d\n", t, target_scenes[t].size() );
        for(int s = 0; s < target_scenes[t].size(); s++){
            target_scenes[t][s].PrintScene();
        }
    }
    if(target_scenes){
        delete[] target_scenes;
        target_scenes = NULL;
    }

}