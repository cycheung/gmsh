//
// Description : Parametric representation of a beam
//
//
// Author:   <Boris Sedji>,  04/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _BEAMLIN_H_
#define _BEAMLIN_H_

#include "BeamParam.h"
#include "MElement.h"

class SPoint2;

class BeamLin : public BeamParam
{

  protected :

    std::vector <SPoint2> _Points;
    double _L;
    double _cos;
    double _sin;
    double _c;

  public :

    BeamLin(std::vector <SPoint2> Points)
    {
      _Points = Points;
      _L = (_Points[1].x() - _Points[0].x())*(_Points[1].x() - _Points[0].x()) +  (_Points[1].y() - _Points[0].y() )*(_Points[1].y() - _Points[0].y() ) ;
      _L = sqrt(_L);
      _cos = (_Points[1].x() - _Points[0].x() ) / _L;
      _sin = (_Points[1].y() - _Points[0].y() ) / _L;
      _c = -_sin*_Points[0].x() + _cos*_Points[0].y();
    }
    virtual SPoint2 getPoint(double t);
    virtual BeamParam* cut(MElement *e);
    virtual double getSin() {return _sin;};
    virtual double getCos() {return _cos;};
    virtual double getC(){return _c;};
    virtual double getLength(){return _L;};
    virtual int getType(){return 1;};
    virtual std::vector <SPoint2>* getPoints(){return &_Points;};

};



SPoint2 BeamLin::getPoint(double t)

{

  double x;
  double y;

  x = _Points[0].x() + t*_L*_cos ;
  y = _Points[0].y() + t*_L*_sin ;

  SPoint2 X(x,y);

  return X;

}

BeamParam *BeamLin::cut(MElement *e)
{

  int j = 0;
  bool intersect=false;
  // take the vertices and verify if it s crossed by the line
  // cartesian equation of the line : sin * x - cos * y + c = 0
  for (int i=0; i<e->getNumVertices(); i++)
  {
    if (_sin*e->getVertex(i)->x() - _cos*e->getVertex(i)->y() + _c >=0) j++;
    else j--;
  }

  // if all same sign then same side
  if (j == e->getNumVertices() || j == -e->getNumVertices()) return NULL;
  else // if cut by BeamLin direction line
  {

  std::vector <SPoint2> InterpolPoints;
  double coef;  // intersection with element side ratio
  double Lc;  // element side length
  double cn,sn;  // unit director vector of the side

  // we ll take each side delimited by vertices (i,i+1)
  for (int i=0;i<e->getNumVertices();i++)
  {
    double sin = _sin;
    double cos = _cos;
    double c = _c;
    double a ;
    int j ;
    double b;

    // ditance from line of vertex i
    a = sin * e->getVertex(i)->x() - cos * e->getVertex(i)->y() + c;
    j = (i+1) % e->getNumVertices();
    // ditance from line of vertex j = i+1
    b = sin * e->getVertex(j)->x() - cos * e->getVertex(j)->y() + c;

    // if it s a crossed side
    if ((a > 0 && b < 0) || (a < 0 && b > 0))
    {
      // side length
      Lc = (e->getVertex(j)->x() - e->getVertex(i)->x())*(e->getVertex(j)->x() - e->getVertex(i)->x()) +  (e->getVertex(j)->y() - e->getVertex(i)->y() )*(e->getVertex(j)->y() - e->getVertex(i)->y() ) ;
      Lc = sqrt(Lc);

      cn = (e->getVertex(j)->x() - e->getVertex(i)->x() ) / Lc;
      sn = (e->getVertex(j)->y() - e->getVertex(i)->y() ) / Lc;

      // thales law of proportional sides
      coef = std::abs(a/(std::abs(a)+std::abs(b)));

      double x;
      double y;
      // intersection point
      x = e->getVertex(i)->x() + coef * cn * Lc;  // intersection point
      y = e->getVertex(i)->y() + coef * sn * Lc;

      std::vector <SPoint2>* BeamPoints = &_Points;

      double L1;
      double L2;

      // distance from beam vertices
      L1 = sqrt(( (*BeamPoints)[0].x() - x)*((*BeamPoints)[0].x()  - x) +  ((*BeamPoints)[0].y()  - y )*((*BeamPoints)[0].y()  - y ));
      L2 = sqrt(( (*BeamPoints)[1].x() - x)*((*BeamPoints)[1].x()  - x) +  ((*BeamPoints)[1].y()  - y )*((*BeamPoints)[1].y()  - y ));

      // if length is greater than beam length then it s not in the beam
      if (L1 + L2 > _L)
      {
        // and we take the nearest point as interpolation point
        if (L1 < L2) InterpolPoints.push_back((*BeamPoints)[0]);
        else InterpolPoints.push_back((*BeamPoints)[1]);
      }
      else
      {
        // at least one point has to be different of beam vertices
        intersect = true;
        SPoint2 p(x,y);
        InterpolPoints.push_back(p);
      }
    }
   }
   if (intersect)
   {
      BeamLin *Beam = new BeamLin(InterpolPoints);
      return Beam ;
   }
   else return NULL;
  }
}

#endif

