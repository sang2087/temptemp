#include "BD.h"
#include <iostream>
#include "math.h"
#include <cstdlib>
#include <ctime>

using namespace std;

void BD::set_value(int input_n, float input_p){
  n_multiplier = 1;
  if(n > 100){
    n_multiplier = 1.0 * input_n / 100;
    input_n = 100;
  }
  n = input_n;
  p = input_p;

  dist_arr = new double[n+1];
  for(int k=1;k<=n;k++){
    double c=10000.0;
    for(int i=n-k+1;i<=n;i++){
      c *=i;
    }
    for(int i=1;i<k;i++){
      c/=i;
    }
    dist_arr[k] = c * pow(p,k) * pow(1-p, n-k);
  }
}

int BD::get_binary_distribution(){
  double sum = 0;
  double r = 1.0 * (rand() % 10000 * n * p);
  for(int i=0;i<n;i++){
   sum += dist_arr[i];
   if(sum>r)
     return (int)((i-1)*n_multiplier);
  }
  return 0;
}

