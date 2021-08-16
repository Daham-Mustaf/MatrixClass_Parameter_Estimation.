#include "Matrix.hpp"

// using namespace std;

Matrix::Matrix() {
  m_total = m_total + 2 * 4 + values.size() * sizeof(values[0]);
  c_totalMemory();
}

Matrix::Matrix(int r, int c) : rows(r), columns(c) {
  rows = r;
  columns = c;
  for (int i = 0; i < rows * columns; i++)
    values.push_back(0);
  m_total = m_total + 2 * 4 + values.size() * sizeof(values[0]);
  c_totalMemory();
}

Matrix::Matrix(int r, int c, const vector<double> val)
    : rows(r), columns(c), values(val) {
  rows = r;
  columns = c;
  for (int i = 0; i < rows * columns; i++) {
    values.push_back(val[i]);
  }
  m_total = m_total + 2 * 4 + values.size() * sizeof(values[0]);
  c_totalMemory();
}

Matrix::~Matrix() { values.clear(); }

// copy constructor
Matrix::Matrix(const Matrix &copy_M) {
  rows = copy_M.rows;
  columns = copy_M.columns;
  values = copy_M.values;
  m_total = m_total + 2 * 4 + values.size() * sizeof(values[0]);
  c_totalMemory();
}

// operator=
Matrix &Matrix::operator=(const Matrix &m) {
  m_total = m_total - values.size() * sizeof(values[0]);
  rows = m.rows;
  columns = m.columns;
  values = m.values;
  m_total = m_total + values.size() * sizeof(values[0]);
  c_totalMemory();

  return (*this);
}

// operator+=
Matrix &Matrix::operator+=(const Matrix &m2) {
  if (rows != m2.rows || columns != m2.columns) {
    cout << "Wrong! Only Matrices with same row and column could plus!\n"
         << endl;
    return (*this);
  }
  for (int i = 0; i < rows * columns; i++) {
    values[i] += m2.values[i];
  }
  return (*this);
}
Matrix &Matrix::operator+=(double x) {
  for (int i = 0; i < rows * columns; i++) {
    values[i] += x;
  }
  return (*this);
}

// operator+
Matrix Matrix::operator+(const Matrix &m2) const {
  if (rows != m2.rows || columns != m2.columns) {
    cout << "Wrong! Only Matrices with same row and column could plus!\n"
         << endl;
    return (*this);
  }
  Matrix m(*this);
  m += m2;
  return (m);
}
Matrix Matrix::operator+(double x) const {
  Matrix m(*this);
  m += x;
  return (m);
}

// operator-=
Matrix &Matrix::operator-=(const Matrix &m2) {
  if (rows != m2.rows || columns != m2.columns) {
    cout << "Wrong! Only Matrices with same row and column could minus!\n";
    return (*this);
  }
  for (int i = 0; i < rows * columns; i++) {
    values[i] -= m2.values[i];
  }
  return (*this);
}
Matrix &Matrix::operator-=(double x) {
  for (int i = 0; i < rows * columns; i++) {
    values[i] -= x;
  }
  return (*this);
}

// operator-
Matrix Matrix::operator-(const Matrix &m2) const {
  if (rows != m2.rows || columns != m2.columns) {
    cout << "Wrong! Only Matrices with same row and column could minus!\n";
    return (*this);
  }
  Matrix m(*this);
  m -= m2;
  return (m);
}
Matrix Matrix::operator-(double x) const {
  Matrix m(*this);
  m -= x;
  return (m);
}

// operator*=
Matrix &Matrix::operator*=(double x) {
  for (int i = 0; i < rows * columns; i++) {
    values[i] *= x;
  }
  return (*this);
}

// operator*
Matrix Matrix::operator*(double x) const {
  Matrix m(*this);
  m *= x;
  return (m);
}

// operator()
double Matrix::operator()(int r, int c) const {
  return values[(c - 1) * rows + r - 1];
}

// member function
// resize
void Matrix::reshape(int r, int c) {
  values.resize(r * c);
  rows = r;
  columns = c;
}

// get entries
double Matrix::get(int r, int c) { return values[(c - 1) * rows + r - 1]; }

int Matrix::getrow() { return rows; }
int Matrix::getcolumn() { return columns; }

// set entries
void Matrix::set(int r, int c, double v) { values[(c - 1) * rows + r - 1] = v; }

// largest entry
double Matrix::max() { return *max_element(values.begin(), values.end()); }

// smallest entry
double Matrix::min() { return *min_element(values.begin(), values.end()); }

// largest entry position
void Matrix::max_P() {
  int r, c;
  int maxPosition = max_element(values.begin(), values.end()) - values.begin();
  if ((maxPosition + 1) % columns != 0) {
    r = floor((maxPosition + 1) / columns) + 1;
    c = maxPosition + 1 - (r - 1) * columns;
  } else {
    r = floor((maxPosition + 1) / columns);
    c = columns;
  }
  cout << "the largest entry is at row " << r << " column " << c << endl;
}

// smallest entry position
void Matrix::min_P() {
  int r, c;
  int minPosition = min_element(values.begin(), values.end()) - values.begin();
  if ((minPosition + 1) % columns != 0) {
    r = floor((minPosition + 1) / columns) + 1;
    c = minPosition + 1 - (r - 1) * columns;
  } else {
    r = floor((minPosition + 1) / columns);
    c = columns;
  }
  cout << "the smallest entry is at row " << r << " column " << c << endl;
}

// cholesky factorization
void Matrix::chol(Matrix &L, int n) {
  L.set(1, 1, sqrt(values.at(0)));
  for (int i = 2; i <= n; i++) {
    L.set(i, 1, values.at((i - 1) * n) / L.get(1, 1));
  }
  for (int j = 2; j <= n; j++) {
    double temp = 0;
    for (int k = 1; k < j; k++) {
      temp += L.get(j, k) * L.get(j, k);
    }
    L.set(j, j, sqrt(values.at((j - 1) * n + j - 1) - temp));
    for (int i = j + 1; i <= n; i++) {
      temp = 0;
      for (int k = 1; k < j; k++) {
        temp += L.get(i, k) * L.get(j, k);
      }
      L.set(i, j, (values.at((i - 1) * n + j - 1) - temp) / L.get(j, j));
    }
  }
}

// IO
void Matrix::AsciiRead(const string filename) {
  ifstream in(filename);
  if (!in.is_open()) {
    cerr << "Error! Can't open file:" << filename << '\n';
    return;
  }
  m_total = m_total - values.size() * sizeof(values[0]); // calculate m_total
  in >> rows;
  in >> columns;
  values.resize(rows * columns);
  int i = 0;
  while (!in.eof() && i < rows * columns) {
    in >> values.at(i);
    i++;
  }
  // cout << "INFO::AsciiRead. Read a Matrix(" << rows << ',' << columns<< ")
  // from file.\n";
  m_total = m_total + values.size() * sizeof(values[0]);
  c_totalMemory();
}

void Matrix::AsciiWrite(const string filename) {
  ofstream out(filename);
  if (!out.is_open()) {
    cerr << "Error! Can't open file:" << filename << '\n';
    return;
  }
  out << rows << endl;
  out << columns << endl;
  out << setprecision(10);
  copy(values.begin(), values.end(), ostream_iterator<double>(out, "\n"));
  cout << "INFO::AsciiWrite. Write Matrix in file '" << filename << "'\n";
}

void Matrix::BinaryRead(const string filename) {
  fstream in(filename, ios_base::in | ios_base::binary);
  if (!in.is_open()) {
    cerr << "Error! Can't open file:" << filename << '\n';
    return;
  }
  m_total = m_total - values.size() * sizeof(values[0]); // calculate m_total
  in.read(reinterpret_cast<char *>(&rows), sizeof(rows));
  in.read(reinterpret_cast<char *>(&columns), sizeof(columns));
  values.resize(rows * columns);
  in.read(reinterpret_cast<char *>(&values[0]),
          values.size() * sizeof(values[0]));

  cout << "INFO::BinaryRead. Read a Matrix(" << rows << ',' << columns
       << ") from file.\n";
  m_total = m_total + values.size() * sizeof(values[0]);
  c_totalMemory();
}

void Matrix::BinaryWrite(const string filename) {
  ofstream out(filename, ios_base::binary);
  if (!out.is_open()) {
    cerr << "Error! Can't open file:" << filename << '\n';
    return;
  }
  out.write(reinterpret_cast<char *>(&rows), sizeof(rows));
  out.write(reinterpret_cast<char *>(&columns), sizeof(columns));
  out.write(reinterpret_cast<char *>(&values.front()),
            values.size() * sizeof(values.front()));
  cout << "INFO::BinaryWrite. Write Matrix in file '" << filename << "'\n";
}

// standard output print
void Matrix::print() {
  cout << "[";
  for (int i = 0; i < rows - 1; i++) {
    for (int j = 0; j < columns; j++) {
      cout << " ";
      cout << values[rows * j + i];
    }
    cout << ";" << endl;
  }
  for (int j = 0; j < columns ; j++) {
    cout << " ";
    cout << values[(j + 1) * rows - 1];
  }
  cout << "] ; % dim " << rows << "x" << columns << endl << endl;
}

int Matrix::m_total = 0;
int Matrix::c_totalMemory() {
  // cout << "Current memory usage: " << m_total << " Byte\n";
  return m_total;
}

// OpenBLAS
void Matrix::isSymmProductOf(const Matrix &A, char TransA) {
  // calculate size
  if (TransA == 'N' || TransA == 'n') // A*A'
  {
    columns = A.rows;
    rows = A.rows;
    values.resize(columns * rows);
    f77blas_dgemm('N', 'T', A.rows, A.rows, A.columns, 1.0, A.values.data(),
                  A.rows, A.values.data(), A.rows, 0.0, values.data(), rows);
  } else // TransA == 'T', A'*A
  {
    columns = A.columns;
    rows = A.columns;
    values.resize(columns * rows);
    f77blas_dgemm('T', 'N', A.columns, A.columns, A.rows, 1.0, A.values.data(),
                  A.rows, A.values.data(), A.rows, 0.0, values.data(), rows);
  }
}

void Matrix::plusSymmProductOf(const Matrix &A, char TransA, double w) {
  if (TransA == 'N' || TransA == 'n') // A*A'
  {
    columns = A.rows;
    rows = A.rows;
    values.resize(columns * rows);
    f77blas_dgemm('N', 'T', A.rows, A.rows, A.columns, w, A.values.data(),
                  A.rows, A.values.data(), A.rows, 1.0, values.data(), rows);
  } else // TransA == 'T', A'*A
  {
    columns = A.columns;
    rows = A.columns;
    values.resize(columns * rows);
    f77blas_dgemm('T', 'N', A.columns, A.columns, A.rows, w, A.values.data(),
                  A.rows, A.values.data(), A.rows, 1.0, values.data(), rows);
  }
}

void Matrix::isProductOf(const Matrix &A, const Matrix &B, char transA,
                         char transB) {
  int k;
  if (transA == 'N' || transA == 'n') {
    rows = A.rows;
    k = A.columns;
  } else {
    rows = A.columns;
    k = A.rows;
  }

  if (transB == 'N' || transB == 'n')
    columns = B.columns;
  else
    columns = B.rows;

  values.resize(columns * rows);

  f77blas_dgemm(transA, transB, rows, columns, k, 1.0, A.values.data(), A.rows,
                B.values.data(), B.rows, 0.0, values.data(), rows);
}

void Matrix::plusProductOf(const Matrix &A, const Matrix &B, char transA,
                           char transB, double w) {
  int k;
  if (transA == 'N' || transA == 'n') {
    rows = A.rows;
    k = A.columns;
  } else {
    rows = A.columns;
    k = A.rows;
  }

  if (transB == 'N' || transB == 'n')
    columns = B.columns;
  else
    columns = B.rows;

  values.resize(columns * rows);

  f77blas_dgemm(transA, transB, rows, columns, k, w, A.values.data(), A.rows,
                B.values.data(), B.rows, 1.0, values.data(), rows);
}

Matrix Matrix::operator*(const Matrix &B) const {
  Matrix C(rows, B.columns);
  C.isProductOf(*this, B, 'N', 'N');
  return (C);
}

Matrix &Matrix::operator*=(const Matrix &B) {
  Matrix C(rows, B.columns);
  C.isProductOf(*this, B, 'N', 'N');
  *this = C;
  return (*this);
}

// LAPACK
void Matrix::chol() {
  int info1;
  f77blas_dpotrf('U', rows, values.data(), rows, info1);
  if (info1 != 0)
    throw std::logic_error(
        "Error using chol. Matrix must be positive definite.");
}
void Matrix::solveWithCholReduced(Matrix &R) {
  f77blas_dtrsm('L', 'U', 'T', 'N', rows, columns, 1.0, R.values.data(), R.rows,
                values.data(), rows);
  f77blas_dtrsm('L', 'U', 'N', 'N', rows, columns, 1.0, R.values.data(), R.rows,
                values.data(), rows);
}
void Matrix::invCholReduced() {
  Matrix R(*this);
  solveWithCholReduced(R);
}