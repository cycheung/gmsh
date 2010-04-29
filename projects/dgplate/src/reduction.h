//
// C++ Interface: terms
//
// Description: Function needed to compute element of reduction n^alpha and m^alpha for thin bodies
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef REDUCTION_H_
#define REDUCTION_H_
#include "LinearElasticShellHookeTensor.h"
#include "LocalBasis.h"
#include "SVector3.h"
inline double scaldot(const SVector3 &a, const SVector3 &b);

static inline void diff(const SVector3 &a,const SVector3 &b,SVector3 &c);

double epsilongd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads, const int nbdof, const fullMatrix<double> &disp);

double rhogd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::HessType> &hess, const int nbdof, const fullMatrix<double> &disp);

void stressReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::GradType> &Grads,const LocalBasis *lb,const fullMatrix<double> &disp, const int nbdof, std::vector<SVector3> &n);

void momentReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::HessType> &hess,const LocalBasis *lb,const fullMatrix<double> &disp, const int nbdof, std::vector<SVector3> &m);

void stressReductionHat(const std::vector<SVector3> &n,const LocalBasis *lb,std::vector<SVector3> &nhat);
// Should be somewhere else ??
void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,const int n_p,const fullMatrix<double> & disp,SVector3 &ujump);
 #endif // reduction.h
