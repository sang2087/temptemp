class BD{
  public:
    void set_value(int n, float p);
    int get_binary_distribution();
  private:
    double* dist_arr;
    int n;
    float n_multiplier;
    float p;
};
