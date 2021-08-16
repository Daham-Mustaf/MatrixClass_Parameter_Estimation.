#include "Matrix.hpp"
#include <iostream>
using namespace std;

int main() {

  Matrix diff1;
  Matrix Pa, x;
  Pa.AsciiRead("FourierParameters.txt");
  x.AsciiRead("EstimatedFourierParameters.txt");
  diff1 = Pa - x;
  cout << "Difference between true and estimated parameters:" << endl;
  diff1.print();

  Matrix diff2;
  Matrix ob, E;
  ob.AsciiRead("TrueObservation.txt");
  E.AsciiRead("EstimatedObservation.txt");
  diff2 = ob - E;
  cout << "Difference between true and adjusted Observation:" << endl;
  diff2.print();

  Matrix FourierSetting, EVF;
  FourierSetting.AsciiRead("FourierSetting.txt");
  EVF.AsciiRead("EstimatedVarianceFactor.txt");
  cout << " Variance0:" << FourierSetting(2, 1) << endl;
  cout << " Estimated Variance Factor:" << EVF(1, 1) << endl;

  return 0;
}