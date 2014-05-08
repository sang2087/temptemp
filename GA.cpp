#define MAX_GENERATION 500000
#define ELITISM true
#define PRESELECT true
#define TWO_OPT false

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
  int mm = 1;
  for(int k=0;k<50;k++){
  //  SetValue(string filename, int population_number, int rullet_constant, int mutation_ratio, int preselect, int opt_number)
    if(k%10 == 0){
      mm=mm+1;
    }
    SetValue("grading1/cycle.in.5", 100, mm, 0.03, 5, 10);

    Initialize();
    time_t start_t = time(0);
    time_t end_t = time(0);
    int i = 0;
    while(i < MAX_GENERATION){
  //    cout << "\nGENERATION " << i << endl;
  //    cout << "crossover" << endl;
      Crossover(); //all offsprings are crossovered, include selection.
  //    cout << "Mutation" << endl;
      Mutate();
  //    cout << "Analysis" << endl;
      if(end_t - start_t > limit_time * 0.10){
        start_optimization = true;
      }
      Analysis(); //analyze several feature. ex)best, worst, etc.
      if(TWO_OPT && start_optimization){
        //어느정도 수렴이 되면 최적화를 시킬 것!
        cout << "!!!" << endl;
        start_optimization = true;
        Optimization();
        Analysis();
        //cout << "\n" << endl;
      }
        //    cout << "replace" << endl;
      Replace();
      i++;

      end_t = time(0);
      if(end_t - start_t > limit_time * 0.90){
        cout << "BEST  " << total_best_offspring.fitness << endl;
        for(int j=0;j<gene_number;j++){
          cout << total_best_offspring.gene[j].number <<" ";
        }
        cout << endl;

        break;
      }
    }
  }
}

void GA::Crossover(){
  ConstructRulletwheel();

  for(int i=0;i<population_number;i++){
    if(ELITISM){
      if(i == best_offspring_index)
       continue;
    }
    int parent1_index = Select();
    int parent2_index = Select();

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
    CyclicCrossover(i, parent1_index, parent2_index);
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

int GA::Select(){ //point rulletwheel
  if((int)rullet[population_number] != 0){ //if completly converge(100%), rullet[population_number] == 0
    int point = rand() % (int)rullet[population_number];
    for(int i=1;i<=population_number;i++){
      if(point < rullet[i]){
        return i-1;
      }
    }
  }else{
    //cout << "COMPLETLY CONVERGED!!" << endl;
    return rand() % population_number;
  }
  cout << "SELECTION ERROR" << endl;
  return -1;
}

void GA::Replace(){
  if(PRESELECT){
    for(int i=0;i<preselect_number;i++){
      offspring[preselect_index[i]] = offspring[best_offspring_index];
    }
  }
  preselect_index = NULL;
  delete preselect_index;

  parent = NULL;
  delete parent;
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

  int cyclic_offspring_index = 0;
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
}

void GA::ConstructRulletwheel(){
  rullet = NULL;
  delete rullet;

  rullet = new float[population_number + 1];
  rullet[0] = 0;
  for(int i=0;i<population_number;i++){
    float ratio = (worst_parent.fitness - parent[i].fitness) + (worst_parent.fitness - best_parent.fitness)/(rullet_constant - 1);
    rullet[i+1] = rullet[i] + ratio;
    //cout << rullet[i] << " ";
  }
  //cout << endl;
}

void GA::Analysis(){
  if(PRESELECT){
    preselect_index = new int[preselect_number];
    for(int i=0;i<preselect_number;i++){
      preselect_index[i] = i;
    }
    for(int i=0;i<preselect_number;i++){
      for(int j=0;j<i;j++){
        if(offspring[preselect_index[i]].fitness > offspring[preselect_index[j]].fitness){
          swap(preselect_index[i], preselect_index[j]);
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
    for(int i=0;i<opt_number;i++){
      cout << offspring[opt_index[i]].fitness << " ";
    }
    cout << endl;
  }

  offspring[0].fitness = CalculateFitness(offspring[0].gene); //Calculate Fitness

  best_offspring = offspring[0];
  worst_offspring = offspring[0];
  best_offspring_index = 0;
  worst_offspring_index = 0;

  for(int i=0;i<population_number;i++){
    offspring[i].fitness = CalculateFitness(offspring[i].gene); //Calculate Fitness
    if(offspring[i].fitness < best_offspring.fitness){
      best_offspring = offspring[i];
      best_offspring_index = i;
    }
    if(offspring[i].fitness > worst_offspring.fitness){
      worst_offspring = offspring[i];
      worst_offspring_index = i;
    }

    if(PRESELECT){
      if(offspring[i].fitness >= offspring[preselect_index[preselect_number-1]].fitness){
        int k = preselect_number;
        while(k >=1 && offspring[i].fitness >= offspring[preselect_index[k-1]].fitness){
          k--;
        }

        for(int j = preselect_number - 1; j > k;j--){
          preselect_index[j] = preselect_index[j-1];
        }
        preselect_index[k] = i;
      }
    }
    if(TWO_OPT && start_optimization){
      if(offspring[i].fitness <= offspring[opt_index[preselect_number-1]].fitness){
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
  if(best_offspring.fitness < total_best_offspring.fitness){
    //cout << "BEST  " << best_offspring.fitness << endl;
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
  offspring = NULL;
  delete offspring;
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

void GA::SetValue(string filename, int p, int r, float m, int pre, int o){
  GetInput(filename);
  population_number = p;
  rullet_constant = r;
  mutation_ratio = m;
  preselect_number = pre;
  opt_number = o;
  cout << "population_number : " << population_number << endl;
  cout << "rullet_constant : " << rullet_constant << endl;
  cout << "mutation_ratio : " << mutation_ratio << endl;
  cout << "preselect_number : " << preselect_number << endl;
  cout << "limit_time : " << limit_time << endl;

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
    int offspring_index = opt_index[opt_iterator];
    //int offspring_index = best_offspring_index;
    Chromosome best_offspring = ChromosomeClone(offspring[offspring_index]);
    for(int i=1;i<gene_number;i++){
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
          best_offspring.gene = candidate_gene;
          best_offspring.fitness = candidate_fitness;
        }
      }
    }
    offspring[offspring_index] = best_offspring;
  }
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
  offspring[offspring_index].gene = offspring_gene;
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

