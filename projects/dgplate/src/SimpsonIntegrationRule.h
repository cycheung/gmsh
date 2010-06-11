//
// C++ Interface: terms
//
// Description: Simpson integration Rule (integration of thickness for thin bodies)
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
# ifndef SIMPSONINTEGRATIONRULE_H_
# define SIMPSONINTEGRATIONRULE_H_
double SimpsonIntegration(const std::vector<double> &y, const std::vector<double> &zk){
  // init
  int nsimpminus1 = zk.size()-1;
  // Value of step
  double twotimeshdiv3 = 2.*(zk[nsimpminus1]-zk[0])/(3.*nsimpminus1);

  // integration
  double val = 0.5*(y[0]*zk[0]+y[nsimpminus1]*zk[nsimpminus1]);
  for(int i=2;i<nsimpminus1;i+=2){
    val+=2*y[i-1]*zk[i-1]+y[i]*zk[i];
  }
  val+=2*y[nsimpminus1-1]*zk[nsimpminus1-1];
  return twotimeshdiv3*val;
}

double SimpsonIntegration(const std::vector<double> &y, const double bminusa){
  // init
  int nsimpminus1 = y.size()-1;
  // Value of step
  double twotimeshdiv3 = 2.*(bminusa)/(3.*nsimpminus1);

  // integration
  double val = 0.5*(y[0]+y[nsimpminus1]);
  for(int i=2;i<nsimpminus1;i+=2){
    val+=2*y[i-1]+y[i];
  }
  val+=2*y[nsimpminus1-1];
  return twotimeshdiv3*val;
}
#endif // SIMPSONINTEGRATIONRULE_H_
