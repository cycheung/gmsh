//
// Description : Parametric representation of a beam
//
//
// Author:   <Boris Sedji>,  04/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

BeamTerm::BeamTerm(FunctionSpace<SVector3> *LagSpace,MElement *e, BeamParam *Beam)
{
  _Beam = Beam;
  _e = e;
  _LagSpace = LagSpace;
}

void BeamTerm::get(fullMatrix<double> &m)
{

     bool sym  = true;
     double gp[3],uvw[3];
     double EA = 1e10;
     gp[2] = 0;
     SPoint2 gpxy = _Beam->getPoint(0.5);
     gp[0] = gpxy.x();
     gp[1] = gpxy.y();
     _e->xyz2uvw(gp,uvw);
     int npts = 1;
     IntPt GP[1];
     GP[0].pt[0] = uvw[0];GP[0].pt[1] = uvw[1] ;GP[0].pt[1] = uvw[2];
     GP[0].weight = 2;
     double cos = _Beam->getCos();
     double sin = _Beam->getSin();
     if (sym)
      {
        int nbFF = _LagSpace->getNumKeys(_e);
        double jac[3][3];
        fullMatrix<double> B(1, nbFF);
        fullMatrix<double> BT(nbFF, 1);
        printf("npts : %d\n",npts);
        m.resize(nbFF, nbFF);
        m.setAll(0.);
        std::cout << m.size1() << "  " << m.size2() << std::endl;
        for (int i = 0; i < npts; i++)
        {
          const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
          if (_e->getParent()) _e=_e->getParent();
          const double weight = GP[i].weight;
          //const double detJ = ele->getJacobian(u, v, w, jac);
          const double detJ = _Beam->getLength();
          std::vector<TensorialTraits<SVector3>::GradType> Grads;
          _LagSpace->gradf(_e,u, v, w, Grads); // a optimiser ??
          for (int j = 0; j < nbFF; j++)
          {
            BT(j, 0) = B(0, j) = Grads[j](0,0)*cos*cos + Grads[j](0,1)*sin*cos + Grads[j](1,0)*sin*cos + Grads[j](1,1)*sin*sin;
          }
          m.gemm(BT, B, EA*weight * detJ, 1.);
        }
    }
    m.print();
}
