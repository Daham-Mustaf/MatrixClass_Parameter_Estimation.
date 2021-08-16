#ifndef MATRIX_HPP_
#define MATRIX_HPP_
#include "f77blas.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class Matrix {
private:
  int rows;
  int columns;
  vector<double> values;
  static int m_total; // total memory

public:
  // constructor
  Matrix();
  Matrix(int r, int c);
  Matrix(int r, int c, const vector<double> val);
  Matrix(const Matrix &copy_M);
  virtual ~Matrix();

  // operator
  Matrix &operator=(const Matrix &m);
  Matrix operator+(const Matrix &m2) const;
  Matrix operator+(double x) const;
  Matrix operator-(const Matrix &m2) const;
  Matrix operator-(double x) const;
  Matrix operator*(double x) const;
  Matrix &operator+=(const Matrix &m2);
  Matrix &operator+=(double x);
  Matrix &operator-=(const Matrix &m2);
  Matrix &operator-=(double x);
  Matrix &operator*=(double x);
  double operator()(int r, int c) const;

  // member fucntion
  void reshape(int r, int c);
  double get(int r, int c);
  int getrow();
  int getcolumn();
  void set(int r, int c, double v);
  double max();
  double min();
  void max_P();
  void min_P();
  static int c_totalMemory();
  void chol(Matrix &L, int n);

  // IO
  void AsciiRead(const string filename);
  void AsciiWrite(const string filename);
  void BinaryRead(const string filename);
  void BinaryWrite(const string filename);
  void print();

  // OpenBLAS
  void isSymmProductOf(const Matrix &A, char TransA);
  void plusSymmProductOf(const Matrix &A, char TransA, double w = 1.0);
  void isProductOf(const Matrix &A, const Matrix &B, char transA, char transB);
  void plusProductOf(const Matrix &A, const Matrix &B, char transA, char transB,
                     double w = 1.0);
  Matrix operator*(const Matrix &m2) const;
  Matrix &operator*=(const Matrix &m2);

  // LAPACK
  void chol();
  void solveWithCholReduced(Matrix &R);
  void invCholReduced();
};
#endif