#ifndef FOUR_HPP_
#define FOUR_HPP_
#include "Matrix.hpp"

using namespace std;

class Fourierseries {
private:
  int nmax;
  int w;
  Matrix Pa;

public:
  // constructor
  Fourierseries();
  Fourierseries(int n, double w);
  Fourierseries(int n, double w, const Matrix pa);
  Fourierseries(const Fourierseries &F);
  virtual ~Fourierseries(); // deconstructor

  // member function
  void getFunctionalValues(Matrix &t, Matrix &l);
  void getDesignMatrix(Matrix &t, Matrix &A);
};
#endif