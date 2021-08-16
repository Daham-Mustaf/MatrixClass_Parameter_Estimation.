#include "Matrix.hpp"
#include <random>
#include <vector>
using namespace std;

int main() {

  // generates N uniformly distributed control points in a range xmin to xmax
  cout << "Please input the number of control points:" << endl;
  int N;
  cin >> N;
  cout << "Please input the lower and upper boundary of control points:"
       << endl;
  double xmin, xmax;
  cin >> xmin >> xmax;

  // generate control points in a matrix
  Matrix controlpoint(N, 1);
  random_device r; // seed
  mt19937 myBasicGenerator(r());
  uniform_real_distribution<double> aUniDist(xmin, xmax);
  for (int i = 1; i <= N; i++)
    controlpoint.set(i, 1, aUniDist(myBasicGenerator));
  // write control points in a file
  controlpoint.AsciiWrite("ControlPoints.txt");

  cout << "Please input the maximum order of Fourier Series:" << endl;
  int nmax; // maximum order of Fourier Series
  cin >> nmax;

  // generate random parameters ai and bi, store in maxtix and write in a file
  uniform_real_distribution<double> paraUniDist(-2.0,
                                                2.0); 
  double n = 2 * nmax + 1;
  Matrix Pa(n, 1);
  for (int i = 1; i <= n; i++)
    Pa.set(i, 1, paraUniDist(myBasicGenerator));
  Pa.AsciiWrite("FourierParameters.txt");

  return 0;
}