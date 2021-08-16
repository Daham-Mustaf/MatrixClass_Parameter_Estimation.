#include "Fourierseries.hpp"
#include "Matrix.hpp"
#include <cmath>
#include <random>
#include <vector>
using namespace std;

int main() {

  // read the control points and true parameters
  Matrix controlpoint, Pa;
  controlpoint.AsciiRead("ControlPoints.txt");
  Pa.AsciiRead("FourierParameters.txt");

  // specify Fourierseries variables
  double w0, sigma0;
  cout << "Please input the fundamental frequency of Fourier Series:" << endl;
  cin >> w0;
  cout << "Please input the variance factor:" << endl;
  cin >> sigma0;

  // generate the N 'true' observation for control points
  int nmax;
  nmax = (Pa.getrow() - 1) / 2; // maximum order
  Fourierseries F(nmax, w0, Pa);
  Matrix ob(controlpoint.getrow(), 1);
  F.getFunctionalValues(controlpoint, ob);
  ob.AsciiWrite("TrueObservation.txt");

  // simulated observation
  random_device r; // seed
  mt19937 myBasicGenerator(r());
  normal_distribution<double> aNormalDistribution(0, sqrt(sigma0));

  // add noise
  Matrix Noise(controlpoint.getrow(), 1);
  for (int i = 1; i <= controlpoint.getrow(); i++)
    Noise.set(i, 1, aNormalDistribution(myBasicGenerator));
  Matrix S(controlpoint.getrow(), 1);
  S = ob + Noise;
  S.AsciiWrite("SimulatedObservation.txt");

  // store w0 and sigma0
  vector<double> set{w0, sigma0};
  Matrix FourierSetting(2, 1, set);
  FourierSetting.AsciiWrite("FourierSetting.txt");

  return 0;
}