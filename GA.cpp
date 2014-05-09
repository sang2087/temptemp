#define MAX_GENERATION 500000
#define ELITISM true
#define WORST_DELETE false
#define TWO_OPT true

#include "GA.h"

void GA::Initialize(){
  //random init
  srand(time(NULL));

  //set values

  parent = new Chromosome[population_number];
  offspring = new Chromosome[population_number];
  total_best_offspring.fitness = 1000000;

  Reset();

  Analysis();
}

GA::GA(){
  //  SetValue(string filename, int population_number, int rullet_constant, int mutation_ratio, int worst_delete, int opt_number)
  SetValue("cycle.in.318", 50, 10, 0.003, 10, 1, 0.7);

  Initialize();
  time_t start_t = time(0);
  time_t end_t = time(0);
  int i = 0;
  while(i < MAX_GENERATION){
//    cout << "\nGENERATION " << i << endl;
//     cout << "crossover" << endl;

    Crossover(); //all offsprings are crossovered, include selection.
//      cout << "Mutation" << endl;
    Mutate();

//     cout << "Analysis" << endl;
    if(end_t - start_t > limit_time * time_ratio){
      start_optimization = true;
    }

    Analysis(); //analyze several feature. ex)best, worst, etc.

    if(TWO_OPT && start_optimization){

      cout << "Optimization" << endl;
      //cout << "!!!" << endl;
      start_optimization = true;
      Optimization();

  //  cout << "Analy2" << endl;
      Analysis();
      //cout << "\n" << endl;
    }
   // cout << "replace" << endl;
    Replace();
    i++;

    end_t = time(0);
    if(end_t - start_t > limit_time * 0.90){
      SetOutput();
      cout << "BEST  " << total_best_offspring.fitness << endl;
      // for(int j=0;j<gene_number;j++){
      //   cout << total_best_offspring.gene[j].number <<" ";
      // }
      //cout << endl;
      //cout << "opt_counter :" << opt_counter<< endl;

      break;
    }
  }
  delete[] parent;
  delete[] input_gene;
}

void GA::Crossover(){
  //ConstructRulletwheel();
  for(int i=0;i<population_number;i++){
    if(ELITISM){
      if(i == best_offspring_index)
       continue;
    }
    //int* point= Select();
    //int parent1_index = point[0];
    //int parent2_index = point[1];
    int parent1_index = Tournament();
    int parent2_index = Tournament();

    // cout << parent1_index << " " << parent2_index << endl;
    // cout << "before" << endl;
    // for(int j=0;j<gene_number;j++){
    //   cout << parent[parent1_index].gene[j].number << " ";
    // }
    // cout << endl;
    // for(int j=0;j<gene_number;j++){
    //   cout << parent[parent2_index].gene[j].number << " ";
    // }
    // cout << endl;
    //CyclicCrossover(i, parent1_index, parent2_index);
    MultipointCrossover(i, parent1_index, parent2_index,(int)(gene_number / 25));
    //MultipointCrossover(i, parent1_index, parent2_index,5);
    // cout << "after" << endl;
    // for(int j=0;j<gene_number;j++){
    //   cout << offspring[i].gene[j].number << " ";
    // }
    // cout << endl;
  }
}
void GA::Mutate(){
  for(int i = 0;i<population_number;i++){
    if(ELITISM){
      if(i == best_offspring_index)
        continue;
    }
    int mutation_number = mutation_distribution.get_binary_distribution();
    for(int j=0;j<mutation_number;j++){
      Inversion(i);
    }
  }
}

void GA::Inversion(int parent_index){
  int position1 = rand() % (gene_number - 1) + 1;
  int position2 = rand() % (gene_number - 1) + 1;

  swap(offspring[parent_index].gene[position1], offspring[parent_index].gene[position2]);
}

int* GA::Select(){ //point rulletwheel
  if((int)rullet[population_number] != 0){ //if completly converge(100%), rullet[population_number] == 0
    int* point = new int[2];
    point[0] = rand() % (int)rullet[population_number];
    point[1] = rand() % (int)rullet[population_number];
    if(point[0] > point[1]){
      swap(point[0], point[1]);
    }
    int k=0;
    //cout << (int)rullet[population_number] << " " << point[0] << " " << point[1] << endl;
    for(int i=1;i<=population_number;i++){
      if(point[k] < rullet[i]){
        point[k] = i-1;
        if(k==1){
          return point;
        }
        k++;
        if(point[k] < rullet[i]){
          point[k] = i-1;
          return point;
        }
      }
    }
  }else{
    //cout << "COMPLETLY CONVERGED!!" << endl;
    int* point=new int[2];
    point[0]=rand() % population_number; 
    point[1]=rand() % population_number; 
    return point;
  }
  //cout << "SELECTION ERROR" << endl;
  return NULL;
}

void GA::Replace(){
  if(WORST_DELETE){
    for(int i=0;i<worst_delete_number;i++){
      //offspring[worst_delete_index[i]] = offspring[(int)rand()%population_number];
      offspring[worst_delete_index[i]].fitness = -1;
      //cout << offspring[worst_delete_index[i]].fitness << "(" << worst_delete_index[i] << ") ";
      // int gene[gene_number];
      // for(int k=0;k<gene_number;k++){
      //   gene[k] = k + 1;
      // }

      // Gene* new_gene = new Gene[gene_number];
      // new_gene[0] = input_gene[0]; //set first gene to 1

      // for(int j=1;j<gene_number;j++){
      //   int r = (rand() % (gene_number - j)) + 1;
      //   new_gene[j] = input_gene[gene[r] - 1];
      //   swap(gene[r], gene[gene_number-j]);
      // }

      // offspring[worst_delete_index[i]].gene = new_gene;
      // offspring[worst_delete_index[i]].fitness = CalculateFitness(offspring[worst_delete_index[i]].gene);

    }
    //cout << endl;
  }
  parent = offspring;
}

void GA::CyclicCrossover(int offspring_index, int parent1_index, int parent2_index){
  Gene* offspring_gene = new Gene[gene_number];

  for(int i=0;i<gene_number;i++){
    offspring_gene[i].number = 0;
  }

  int gene_hash[2][gene_number+1];
  for(int i=0;i<gene_number;i++){
    gene_hash[0][parent[parent1_index].gene[i].number] = i;
    gene_hash[1][parent[parent2_index].gene[i].number] = i;
  }

  int cyclic_offspring_index = rand() % gene_number;
  bool loop_change = false;
  for(int i=0;i<gene_number;i++){
    if(!loop_change){
      offspring_gene[cyclic_offspring_index] = parent[parent1_index].gene[cyclic_offspring_index];
    }else{
      offspring_gene[cyclic_offspring_index] = parent[parent2_index].gene[cyclic_offspring_index];
    }

    if(!loop_change){
      cyclic_offspring_index = gene_hash[loop_change][parent[parent2_index].gene[cyclic_offspring_index].number];
    }else{
      cyclic_offspring_index = gene_hash[loop_change][parent[parent1_index].gene[cyclic_offspring_index].number];
    }

    int start_cyclic_offspring_index = cyclic_offspring_index;
    while(offspring_gene[cyclic_offspring_index].number != 0){
      cyclic_offspring_index++;
      if(cyclic_offspring_index == gene_number){
        cyclic_offspring_index = 0;
      }
      if(start_cyclic_offspring_index == cyclic_offspring_index){
        i = gene_number;
        break;
      }
      loop_change = !loop_change;
    }
  }
  offspring[offspring_index].gene = offspring_gene;
  // offspring[offspring_index].fitness = CalculateFitness(offspring[offspring_index].gene);

  // if(parent[parent1_index].fitness <= parent[parent2_index].fitness){
  //   parent[parent2_index] = offspring[offspring_index];
  // }else{
  //   parent[parent1_index] = offspring[offspring_index];
  // }
  
}

void GA::ConstructRulletwheel(){
  rullet = NULL;
  delete[] rullet;

  rullet = new float[population_number + 1];
  rullet[0] = 0;
  for(int i=0;i<population_number;i++){
    float ratio;
    if(parent[i].fitness == -1){
//      cout << "??" << endl;
      ratio = 0.0;
    }else{
      ratio = (worst_parent.fitness - parent[i].fitness) + (worst_parent.fitness - best_parent.fitness)/(rullet_constant - 1);
    }
//    cout << parent[i].fitness << " " << ratio << endl;
    rullet[i+1] = rullet[i] + ratio;
    //cout << rullet[i] << " ";
  }
  //cout << endl;
}

void GA::Analysis(){
  if(WORST_DELETE){
    worst_delete_index = new int[worst_delete_number];
    for(int i=0;i<worst_delete_number;i++){
      worst_delete_index[i] = i;
    }
    for(int i=0;i<worst_delete_number;i++){
      for(int j=0;j<i;j++){
        if(offspring[worst_delete_index[i]].fitness > offspring[worst_delete_index[j]].fitness){
          swap(worst_delete_index[i], worst_delete_index[j]);
        }
      }
    }
  }

  if(TWO_OPT && start_optimization){
    opt_index = new int[opt_number];
    for(int i=0;i<opt_number;i++){
      opt_index[i] = i;
    }
    for(int i=0;i<opt_number;i++){
      for(int j=0;j<i;j++){
        if(offspring[opt_index[i]].fitness < offspring[opt_index[j]].fitness){
          swap(opt_index[i], opt_index[j]);
        }
      }
    }
  }

  offspring[0].fitness = CalculateFitness(offspring[0].gene); //Calculate Fitness

  best_offspring = offspring[0];
  best_offspring_index = 0;

  for(int i=0;i<population_number;i++){
    offspring[i].fitness = CalculateFitness(offspring[i].gene); //Calculate Fitness
    if(offspring[i].fitness < best_offspring.fitness){
      best_offspring = offspring[i];
      best_offspring_index = i;
    }
    if(WORST_DELETE && i > worst_delete_number){
      if(offspring[i].fitness >= offspring[worst_delete_index[worst_delete_number-1]].fitness){
        int k = worst_delete_number;
        while(k >=1 && offspring[i].fitness >= offspring[worst_delete_index[k-1]].fitness){
          k--;
        }

        for(int j = worst_delete_number - 1; j > k;j--){
          worst_delete_index[j] = worst_delete_index[j-1];
        }
        worst_delete_index[k] = i;
      }
    }

    if(TWO_OPT && start_optimization && i > opt_number){
      if(offspring[i].fitness <= offspring[opt_index[opt_number-1]].fitness){
        int k = opt_number;
        while(k >=1 && offspring[i].fitness <= offspring[opt_index[k-1]].fitness){
          k--;
        }

        for(int j = opt_number - 1; j > k;j--){
          opt_index[j] = opt_index[j-1];
        }
        opt_index[k] = i;
      }
    }
  }

  // cout << best_offspring.fitness<<" " << sum/population_number;
  // if(start_optimization)
  //   cout << "!!!!!!" << endl;
  // else
  //   cout << endl;

  if(best_offspring.fitness < total_best_offspring.fitness){
    cout << "BEST  " << best_offspring.fitness << endl;
    //cout << "WORST " << worst_offspring.fitness << endl;
    total_best_offspring = best_offspring;
  }
  best_parent = best_offspring;
  worst_parent = worst_offspring;
}

void GA::Reset(){
  int gene[gene_number];
  for(int i=0;i<gene_number;i++){
    gene[i] = i + 1;
  }

  for(int i=0;i<population_number;i++){
    parent[i].gene = new Gene[gene_number];
    parent[i].gene[0] = input_gene[0]; //set first gene to 1

    for(int j=1;j<gene_number;j++){
      int r = (rand() % (gene_number - j)) + 1;
      parent[i].gene[j] = input_gene[gene[r] - 1];
      swap(gene[r], gene[gene_number-j]);
    }
  }
  offspring = parent;
  start_optimization = false;
}

float GA::CalculateFitness(Gene* gene){
  float fitness = 0;
  for(int i = 0;i < gene_number;i++){
    int j = i + 1;
    if(i == gene_number - 1){
      j = 0;
    }
    float dx = pow(gene[i].x - gene[j].x, 2);
    float dy = pow(gene[i].y - gene[j].y, 2);
    fitness += sqrt(dx + dy);
  }
  return fitness;
}

void GA::SetValue(string filename, int p, int r, float m, int pre, int o, float t){
  GetInput(filename);
  population_number = p;

  rullet_constant = r;
  mutation_ratio = m;
  worst_delete_number = pre;
  opt_number = o;
  tournament_ratio = t;
  on_breaker = true;

  if(gene_number<150){
    time_ratio = 0.05;
  }else if(gene_number >=150 && gene_number<300){
    time_ratio = 0.2;
  }else if(gene_number >=300 && gene_number<500){
    time_ratio = 0.3;
  }else if(gene_number>=500){
    population_number = p/5;
    mutation_ratio = m/3;
    time_ratio = 0.5;
  }
  //cout << "population_number : " << population_number << endl;
  //cout << "rullet_constant : " << rullet_constant << endl;
  //cout << "mutation_ratio : " << mutation_ratio << endl;
  //cout << "worst_delete_number : " << worst_delete_number << endl;
  //cout << "tournament_ratio : " << tournament_ratio<< endl;
  //cout << "limit_time : " << limit_time << endl;

  mutation_distribution.set_value(gene_number, mutation_ratio);
}

void GA::GetInput(string filename){
  ifstream fin(filename.c_str());
  fin >> gene_number;
  input_gene = new Gene[gene_number];

  for(int i = 0; i < gene_number; i++){
    input_gene[i].number = i + 1;
    fin >> input_gene[i].x;
    fin >> input_gene[i].y;
  }
  fin >> limit_time;

  fin.close();
}


void GA::Optimization(){
  for(int opt_iterator = 0;opt_iterator<opt_number;opt_iterator++){
    opt_counter++;
  //for(int iterator = 0;iterator<population_number;iterator++){
    int offspring_index = opt_index[opt_iterator];
    //int offspring_index = best_offspring_index;
    //int offspring_index = iterator;
    Chromosome best_offspring = ChromosomeClone(offspring[offspring_index]);
    int r = rand() % gene_number;
    for(int k=1;k<gene_number;k++){
      bool breaker = false;
      int i = (r+k) % gene_number;
      for(int j=i+2;j<gene_number;j++){
        Gene* candidate_gene = new Gene[gene_number];
        int iterator = 0;
        while(iterator < gene_number){
          if(iterator < i){
            candidate_gene[iterator] = offspring[offspring_index].gene[iterator];
          }else if(iterator < j){
            int index = i + j - iterator -1;
            candidate_gene[iterator] = offspring[offspring_index].gene[index];
          }else{
            candidate_gene[iterator] = offspring[offspring_index].gene[iterator];
          }
          iterator++;
        }

        float candidate_fitness = CalculateFitness(candidate_gene);
        if(candidate_fitness < best_offspring.fitness){
          for(int i=0;i<gene_number;i++){
            best_offspring.gene[i] = candidate_gene[i];
          }
          best_offspring.fitness = candidate_fitness;
          
          breaker= true;
          break;
          delete[] candidate_gene;
        }
        delete[] candidate_gene;
      }
      if(on_breaker && breaker)
        break;
    }
    //cout << offspring[offspring_index].fitness << " " << best_offspring.fitness << endl;
    //
    offspring[offspring_index].fitness = best_offspring.fitness;
    for(int i=0;i<gene_number;i++){
      offspring[offspring_index].gene[i] = best_offspring.gene[i];
    }
  }
}
int GA::Tournament(){
  int n = 8;
  int player_index[n];
  float best_fitness=0;
  for(int i=0;i<n;i++){
    do{
      player_index[i] = (rand() % population_number);
    }while(parent[player_index[i]].fitness == -1);
  }

  for(int i=n/2;i>=1;i=i/2){
    for(int j=0;j<i;j++){
      int large;
      int small;
      if(parent[player_index[j*2]].fitness > parent[player_index[j*2+1]].fitness){
        large = player_index[j*2];
        small = player_index[j*2+1];
      }else{
        large = player_index[j*2+1];
        small = player_index[j*2];
      }

      if(rand() % 100 > tournament_ratio*100){
        player_index[j] = large;
      }else{
        player_index[j] = small;
      }
    }
  }

  // for(int i=0;i<n;i++){
  //   for(int j=0;j<i;j++){
  //     if(parent[player_index[i]].fitness_score < parent[player_index[j]].fitness_score){
  //       int temp = player_index[i];
  //       player_index[i] = player_index[j];
  //       player_index[j] = temp;
  //     }
  //   }
  // }

  //int k = player_distribution.get_binary_distribution();
  return player_index[0];
}


Chromosome GA::ChromosomeClone(Chromosome src){
  Chromosome dst;
  dst.gene = new Gene[gene_number];
  dst.fitness = src.fitness;
  for(int i=0;i<gene_number;i++){
    dst.gene[i] = src.gene[i];
  }
  return dst;
}

void GA::MultipointCrossover(int offspring_index,int parent1_index, int parent2_index, int cut_number){
  Gene* offspring_gene = new Gene[gene_number];
  int gene[gene_number];
  for(int i=0;i<gene_number-1;i++){
    gene[i] = i+1;
  }
  int cut_point[cut_number];

  for(int i=0; i<cut_number;i++){
    int r = (rand() % (gene_number - i - 2)) ;
    cut_point[i] = gene[r];
    //swap
    swap(gene[r],gene[gene_number-i-2]);
  }

  sort(cut_point, cut_point+cut_number);

  int k = 0;
  bool parent_change = true;
  for(int i=0;i<gene_number;i++){
    if(i == cut_point[k]){
      parent_change = !parent_change;
      k++;
    }
    if(parent_change){
      offspring_gene[i] = parent[parent1_index].gene[i];
    }else{
      offspring_gene[i] = parent[parent2_index].gene[i];
    }
  }

  offspring_gene = repair(offspring_gene);

  for(int i=0;i<gene_number;i++){
    offspring[offspring_index].gene[i] = offspring_gene[i];
  }
  delete[] offspring_gene;
}

Gene* GA::repair(Gene* damaged_gene){
  bool check_gene[gene_number+1];
  damaged_gene[0] = input_gene[0];
  for(int i=1;i<=gene_number;i++){
    check_gene[i] = false;
  }

  int disappeared_gene[gene_number];
  int disappeared_index = 0;
  int duplicated_gene[gene_number];
  int duplicated_index = 0;

  for(int i=0;i<gene_number;i++){
    if(check_gene[damaged_gene[i].number] == false){
      check_gene[damaged_gene[i].number] = true;
    }else{ //중복이 일어나면
      duplicated_gene[duplicated_index] = i; //index i인 곳이 중복되어 있음
      duplicated_index++;
    }
  }
  for(int i=1;i<=gene_number;i++){
    if(check_gene[i] == false){
      disappeared_gene[disappeared_index] = i; //number가 i인 애가 비어있음
      disappeared_index++;
    }
  }

  for(int i=0;i<duplicated_index;i++){
    damaged_gene[duplicated_gene[i]] = input_gene[disappeared_gene[i]-1];
  }

  return damaged_gene;
}
void GA::SetOutput(){
  ofstream fout("cycle.out");
  for(int i=0;i<gene_number;i++){
    fout << total_best_offspring.gene[i].number << " ";
  }
  fout << endl;
  fout.close();
}


