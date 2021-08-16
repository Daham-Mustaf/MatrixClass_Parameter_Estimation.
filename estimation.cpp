#include "Fourierseries.hpp"
#include "Matrix.hpp"

int main() {

  // read neccessary data
  Matrix controlpoint, Pa, S, FourierSetting;
  controlpoint.AsciiRead("ControlPoints.txt");
  Pa.AsciiRead("FourierParameters.txt");
  S.AsciiRead("SimulatedObservation.txt");
  FourierSetting.AsciiRead("FourierSetting.txt");

  // form DesignMatrix A
  double w0 = FourierSetting(1, 1);
  int nmax = (Pa.getrow() - 1) / 2;
  Fourierseries F(nmax, w0);
  Matrix A(controlpoint.getrow(), Pa.getrow());
  F.getDesignMatrix(controlpoint, A);

  Matrix N, n;
  Matrix R(Pa.getrow(),Pa.getrow());
  N.isSymmProductOf(A, 'T');
  n.isProductOf(A, S, 'T', 'N');
  R = N;
  R.chol();

  // get parameter x
  Matrix x;
  x = n;
  x.solveWithCholReduced(R);
  x.print();

  Matrix E;
  E = A * x;

  // get residuals v, variance factor
  Matrix v = E - S;
  // determine the covariance of the parameters
  Matrix coVar(R);
  coVar.invCholReduced();

  Matrix EVF;
  EVF.isProductOf(v, v, 'T', 'N');
  EVF.set(1, 1, EVF(1, 1) / (controlpoint.getrow() - x.getrow()));

  // save in file:  parameters, adjusted observations, variance factor,
  // residuals, covariance matrix
  x.AsciiWrite("EstimatedFourierParameters.txt");
  E.AsciiWrite("EstimatedObservation.txt");
  EVF.AsciiWrite("EstimatedVarianceFactor.txt");
  v.AsciiWrite("Residual.txt");
  coVar.AsciiWrite("Covariance.txt");

  return 0;
}