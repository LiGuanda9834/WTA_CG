#include "WTA.h"

/**
 * @brief Generate a instance by random
 * @details P[i][j] in [0.6, 0.9], Value of targets W in [25, 100]
*/
WTA::WTA(int n_, int m_, int seed)
{
    printf("Create WTA\n");
    target_num_n = n_;
    weapon_num_m = m_;
    target_value = new double[target_num_n];
    //P[i][j]   i : weapon, j : target
    probability_matrix = new double*[weapon_num_m]; 
    srand(seed);
    for(int i = 0; i < target_num_n; i++)
    {
    target_value[i] = rand()%75 + 25;
    }
    DEBUG_MODE = 0;
    
    
    for(int i = 0; i < weapon_num_m; i++)
    {  
        probability_matrix[i] = new double[target_num_n];
        for(int j = 0; j < target_num_n; j++){
            probability_matrix[i][j] = double(rand())/double(RAND_MAX) * 0.3 + 0.6;
            if(DEBUG_MODE == 5){
                printf("%d -> %d : %.4f \n", i, j, probability_matrix[i][j]);
            }
        }
    }
}

// n for target num, m for weapon num
WTA::WTA(int n_, int m_){
    printf("Create empty WTA with size W: %d, T:%d \n Need Initialized\n", m_, n_);
    target_num_n = n_;
    weapon_num_m = m_;
    target_value = new double[target_num_n];
    //P[i][j]   i : weapon, j : target
    probability_matrix = new double*[weapon_num_m]; 
    for(int i = 0; i < target_num_n; i++)
    {
        target_value[i] = 0;
    }
    for(int i = 0; i < weapon_num_m; i++)
    {  
        probability_matrix[i] = new double[target_num_n];
        for(int j = 0; j < target_num_n; j++){
            probability_matrix[i][j] = 0;
        }
    }
    DEBUG_MODE = 0;
}


void WTA::Init(int n_, int m_, int seed = 0)
{
    printf("Initialize by random\n");
    target_num_n = n_;
    weapon_num_m = m_;
    target_value = new double[target_num_n];
    //P[i][j]   i : weapon, j : target
    probability_matrix = new double*[weapon_num_m]; 
    srand(seed);
    for(int i = 0; i < target_num_n; i++)
    {
        target_value[i] = rand()%75 + 25;
    }
    for(int i = 0; i < weapon_num_m; i++)
    {  
        probability_matrix[i] = new double[target_num_n];
        for(int j = 0; j < target_num_n; j++){
            probability_matrix[i][j] = rand()/RAND_MAX * 0.3 + 0.6;
            printf("%d -> %d : %.4f \n", i, j, probability_matrix[i][j]);
        }
    }
}

void WTA::Init_sparse(int n_, int m_, int sparsity_prob_, int seed = 0)
{
    target_num_n = n_;
    weapon_num_m = m_;

    int sparsity_prob = sparsity_prob_;

    int nZero_counter = 0;

    for(int i = 0; i < target_num_n; i++)
    {
        target_value[i] = rand()%75 + 25;
    }
    for(int i = 0; i < weapon_num_m; i++)
    {  
        for(int j = 0; j < target_num_n; j++){
            if (rand()%100 > sparsity_prob){
                probability_matrix[i][j] = double(rand())/double(RAND_MAX) * 0.3 + 0.6;
                nZero_counter++;
            } 
            else{
                probability_matrix[i][j] = 0;
            }
        }
    }
    int probMatrix_size = weapon_num_m * target_num_n;
    double real_sparsity = 0;
    real_sparsity = (double)(probMatrix_size - nZero_counter) / (double)probMatrix_size;
}


void WTA::Init_readfile(FILE* fin, int seed_)
{
    printf("Initialzed by read file \n");
    seed = seed_;
    int status = 0;
    int f = fscanf (fin, "%d", &target_num_n);
    printf("num : %d\n", weapon_num_m);
    if(!target_value){
        printf("new Again !!!!!!!!!!!!!!!!\n");
        target_value = new double[target_num_n];
    }
    for(int i = 0; i < target_num_n; i++)
    {
        fscanf (fin, "%lf", &target_value[i]);
        printf("%d:%.0f\n", i, target_value[i]);
    }
    for(int i = 0; i < target_num_n; i++)
    {
        if(!probability_matrix[i]){
            probability_matrix[i] = new double[target_num_n];
        }
        for(int j = 0; j < weapon_num_m; j++)
        {
            fscanf(fin, "%lf", &(probability_matrix[i][j]));
            //printf("%.2f\n", p[i][j]);
        }
    }
}




void WTA::print_model(){
    printf("->\t");
    for(int i = 0; i < weapon_num_m; i++){
        printf("W%d \t", i);
    }
    printf("value \n");
    for(int j = 0; j < target_num_n; j++){
        printf("T%d\t", j);
        for(int i = 0; i < weapon_num_m; i++){
            printf("%.2f \t", probability_matrix[i][j]);
        }
        printf("%.2f \n", target_value[j]);
    }
}

void WTA::Delete(){
    if(this){
        if(target_value){
            delete[] target_value;
        }
        if(probability_matrix){
            for(int t = 0; t < target_num_n; t++){
                if(probability_matrix[t]){
                    delete[] probability_matrix[t];
                }
            }
            delete[] probability_matrix;
        }
    }
}

double WTA::cal_linear_coef(int _j, vector<int> _scene){
    double q_js = target_value[_j];
    int activated_weapon_num = _scene.size();
    for(int i = 0; i < activated_weapon_num; i++){
        int _weapon = _scene[i];
        double q_js_temp = q_js;
        q_js = q_js * (1 - probability_matrix[_weapon][_j]);
        int DEBUG_qjs = 0;
        if(DEBUG_qjs){
            printf("%s : %.2f * %.2f -> %.2f \n", "q_js", q_js_temp, 1 - probability_matrix[_weapon][_j], q_js);
        }
    }
    return q_js;
}


bool WTA::is_weapon_in_scene(int _i, vector<int> _scene){
    if (std::find(_scene.begin(), _scene.end(), _i) != _scene.end()){
        return 1;
    }
    else{
        return 0;
    }
}
double WTA::set_scene_qjs(Scene& scene){
    double q_js = cal_linear_coef(scene.target_index, scene.weapon_indices);
    scene.obj_qjs = q_js;
    return q_js;
}

double WTA::cal_target_value(Scene temp_scene){
    int t_index = temp_scene.target_index;
    double return_value = target_value[t_index];
    int activated_weapon_num = temp_scene.weapon_indices.size();
    for(int i = 0; i < activated_weapon_num; i++){
        return_value = return_value * (1 - probability_matrix[temp_scene.weapon_indices[i]][t_index]);
    }
    return return_value;
}