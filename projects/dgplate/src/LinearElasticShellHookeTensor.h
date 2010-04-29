//
// C++ Interface: terms
//
// Description: Hooke Tensor for shell (linear elastic)
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

# ifndef LINEARELASTICSHELLHOOKETENSOR_H_
# define LINEARELASTICSHELLHOOKETENSOR_H_
#include <stdio.h>
#include "LocalBasis.h"
// tensor of dimension 4 in 2D (used for Hooke tensor)
class Tensor4dim2{
  protected :
    double tensor[2][2][2][2];
  public :
    Tensor4dim2(){
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++)
              tensor[a][b][c][d]=0.;
    }

    Tensor4dim2(const Tensor4dim2 &_in){
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
          for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++){
              this->set(i,j,ii,jj,_in(i,j,ii,jj));
            }
    }
    ~Tensor4dim2(){}

    void set(const int a,const int b,const int c,const int d,const double dd){tensor[a][b][c][d]=dd;}
    void setAll(const double dd){
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++)
              tensor[a][b][c][d]=dd;

    }

    Tensor4dim2 & operator = (const Tensor4dim2 &_in){
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
          for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++){
              this->set(i,j,ii,jj,_in(i,j,ii,jj));
            }
      return *this;
    }

    double operator() (const int a,const int b,const int c,const int d) const {return tensor[a][b][c][d];}
    double get(const int a,const int b,const int c,const int d) const {return tensor[a][b][c][d];}

    Tensor4dim2 operator+ (Tensor4dim2 &H){
      Tensor4dim2 Hsum;
      double temp;
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++){
              temp=tensor[a][b][c][d]+H(a,b,c,d);
              Hsum.set(a,b,c,d,temp);
            }
    return Hsum;
    }

    Tensor4dim2  operator* (const double dd){
      Tensor4dim2 dH;
      double temp;
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++){
              temp=dd*tensor[a][b][c][d];
              dH.set(a,b,c,d,temp);
            }
      return dH;
    }

    void operator+= (const Tensor4dim2 &H){
      double temp;
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++){
              temp=tensor[a][b][c][d]+H(a,b,c,d);
              this->set(a,b,c,d,temp);
            }
    }
    void operator*= (const double dd){
      double temp;
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++){
              temp=dd*tensor[a][b][c][d];
              this->set(a,b,c,d,temp);
            }
    }
    void print(){
      for(int a=0;a<2;a++)
        for(int b=0;b<2;b++)
          for(int c=0;c<2;c++)
            for(int d=0;d<2;d++)
              printf("(%d,%d,%d,%d)= %f\n",a,b,c,d,tensor[a][b][c][d]);

    }
};

class LinearElasticShellHookeTensor : public Tensor4dim2{

  public :
    LinearElasticShellHookeTensor(){};
    ~LinearElasticShellHookeTensor(){};
    void set(const LocalBasis *lb,const double &C11,const double &nu);
    void set(const LocalBasis *lb,const double &nu);
    void hat(const LocalBasis *lb, const double C11, const double nu);
    void hat(const LocalBasis *lb,const LinearElasticShellHookeTensor *H);
    void mean(const LinearElasticShellHookeTensor *Hp, const LinearElasticShellHookeTensor *Hm);
};
#endif // LinearElasticShellHooketensor
