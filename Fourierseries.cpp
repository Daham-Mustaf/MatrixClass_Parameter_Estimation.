#include "Fourierseries.hpp"
#include <cmath>
#include <fstream>

Fourierseries::Fourierseries(){}

Fourierseries::Fourierseries(int n0, double w0) : nmax(n0), w(w0) {}
Fourierseries::Fourierseries(int n0, double w0, const Matrix pa)
    : nmax(n0), w(w0), Pa(pa) {}

Fourierseries::Fourierseries(const Fourierseries &F) {
  nmax = F.nmax;
  w = F.w;
  Pa = F.Pa;
}

Fourierseries::~Fourierseries(){}

void Fourierseries::getFunctionalValues(Matrix &t, Matrix &l) {

  int n0 = nmax;
  double w0 = w;
  int N = t.getrow();
  for (int i = 1; i <= N; i++) {
    double sum = Pa.get(1, 1) * 0.5;
    for (int j = 1; j <= n0; j++) {
      sum += Pa.get(2 * j, 1) * cos(j * w0 * t.get(i, 1)) +
             Pa.get(2 * j + 1, 1) * sin(j * w0 * t.get(i, 1));
    }
    l.set(i, 1, sum);
  }
}

void Fourierseries::getDesignMatrix(Matrix &t, Matrix &A) {
  int n0 = nmax;
  double w0 = w;
  int N = t.getrow();
  for (int i = 1; i <= N; i++) {
    A.set(i, 1, 0.5);
    for (int j = 1; j <= n0; j++) {
      A.set(i, 2 * j, cos(j * w0 * t.get(i, 1)));
      A.set(i, 2 * j + 1, sin(j * w0 * t.get(i, 1)));
    }
  }
}