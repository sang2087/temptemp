#ifndef GA_H_
#define GA_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <algorithm>
#include "BD.h"

using namespace std;

struct Gene {
  int number;
  float x;
  float y;
};

struct Chromosome {
  Gene* gene;
  float fitness;
};

class GA {
  public:
    GA();
    void Initialize();
    void Crossover();
    void Mutate();
    void Replace();
    void Optimization();
  private:

    void GetInput(string filename);
    int* Select();
    void SetValue(string filename, int p, int r, float m, int pre, int o, float t);
    void Reset();

    void Inversion(int parent_index);

    Gene* repair(Gene* damaged_gene);

    void MultipointCrossover(int offspring_index, int parent1_index, int parent2_index, int cut_number);
    void CyclicCrossover(int offspring_index, int parent1_index, int parent2_index);
    float CalculateFitness(Gene* gene);
    void ConstructRulletwheel();
    void Analysis();

    int Tournament();
    float tournament_ratio;



    Chromosome ChromosomeClone(Chromosome src);

    //chromosomes
    Chromosome* parent;
    Chromosome* offspring;

    //default values
    Gene* input_gene; //inputed genes
    int population_number;
    int gene_number;
    float rullet_constant;
    int limit_time;
    float mutation_ratio;

    //analysis values
    Chromosome total_best_offspring;
    Chromosome best_offspring;
    Chromosome worst_offspring;
    int best_offspring_index;
    int worst_offspring_index;
    Chromosome best_parent;
    Chromosome worst_parent;

    //rullet value
    float* rullet;

    BD mutation_distribution;
    BD tournament_distribution;

    int* worst_delete_index;
    int worst_delete_number;
    bool start_optimization;
    int opt_number;
    int* opt_index;
};

#endif
