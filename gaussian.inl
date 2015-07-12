/***************************************************************************************************
 * File:              $URL: svn://bolitho.dyndns.org/Research/Poisson/Branches/Weighting/Src/PPolynomial.inl $
 * Last Changed By:   $Author: misha $
 * Revision:				  $Revision: 81 $
 * Last Updated:      $Date: 2006-04-14 02:23:51 -0400 (Fri, 14 Apr 2006) $
 *
 ***************************************************************************************************/
#include <math.h>
#include "Gaussian.h"

#define PI     3.1415926535897932384
#define SQRT_2 1.4142135623730950488

////////////////////////
// GaussianPolynomial //
////////////////////////
template<int Degree>
GaussianPolynomial<Degree>::GaussianPolynomial(void){
	ctr=0;
	scl=0;
}
template<int Degree>
template<int Degree2>
GaussianPolynomial<Degree>::GaussianPolynomial(const GaussianPolynomial<Degree2>& G){
	ctr=G.ctr;
	scl=G.scl;
	P=G.P;
}

template<int Degree>
template<int Degree2>
GaussianPolynomial<Degree>& GaussianPolynomial<Degree>::operator  = (const GaussianPolynomial<Degree2> &G){
	this->ctr=G.ctr;
	this->scl=G.scl;
	this->P=G.P;
	return (*this);
}

template<int Degree>
double GaussianPolynomial<Degree>::operator()(const double& t) const{
	double tt=t-ctr;
	return P(t)*exp(-scl*tt*tt);
}
GaussianPolynomial<1> GaussianPolynomial<0>::derivative(void) const{
	GaussianPolynomial<1> gp;
	Polynomial<1> poly;
	gp.ctr=ctr;
	gp.scl=scl;
	poly.coefficients[0]=-ctr;
	poly.coefficients[1]=1.0;
	gp.P=-poly*P*2.0*scl;
	return gp;
}
template<int Degree>
GaussianPolynomial<Degree+1> GaussianPolynomial<Degree>::derivative(void) const{
	GaussianPolynomial<Degree+1> gp;
	Polynomial<1> poly;
	gp.ctr=ctr;
	gp.scl=scl;
	poly.coefficients[0]= ctr*2.0*scl;
	poly.coefficients[1]=-2.0*scl;
	gp.P =Polynomial<Degree+1>(P.derivative())+poly*P;
	return gp;
}
template<int Degree>
GaussianPolynomial<Degree> GaussianPolynomial<Degree>::operator * (const double& t) const{
	GaussianPolynomial gp=(*this);
	gp.P*=t;
	return gp;
}
template<int Degree>
GaussianPolynomial<Degree>& GaussianPolynomial<Degree>::operator *= (const double& t){
	P*=t;
	return (*this);
}
template<int Degree>
GaussianPolynomial<Degree> GaussianPolynomial<Degree>::operator / (const double& t) const{
	GaussianPolynomial gp=(*this);
	gp.P/=t;
	return gp;
}
template<int Degree>
GaussianPolynomial<Degree>& GaussianPolynomial<Degree>::operator /= (const double& t){
	P/=t;
	return (*this);
}
template<int Degree>
GaussianPolynomial<Degree> GaussianPolynomial<Degree>::scale(const double& t) const{
	GaussianPolynomial<Degree> G=(*this);
	G.P.scale(t);
	G.scl/=t*t;
	return G;
}
template<int Degree>
GaussianPolynomial<Degree> GaussianPolynomial<Degree>::shift(const double& t) const{
	GaussianPolynomial<Degree> G=(*this);
	G.P=G.P.shift(t);
	G.ctr+=t;
	return G;
}
template<int Degree>
template<int Degree2>
GaussianPolynomial<Degree+Degree2> GaussianPolynomial<Degree>::operator * (const Polynomial<Degree2>& P) const{
	GaussianPolynomial<Degree+Degree2> G;
	G.ctr=ctr;
	G.scl=scl;
	G.P=this->P*P;
	return G;
}
template<int Degree>
template<int Degree2>
GaussianPolynomial<Degree+Degree2> GaussianPolynomial<Degree>::operator * (const GaussianPolynomial<Degree2>& G) const{
	GaussianPolynomial<Degree+Degree2> GP;
	GP.ctr=(scl*ctr+G.ctr*G.scl)/(scl+G.scl);
	GP.scl=scl+G.scl;
	double s=exp((ctr*scl+G.ctr*G.scl)*(ctr*scl+G.ctr*G.scl)/(scl+G.scl)-scl*ctr*ctr-G.scl*G.ctr*G.ctr);
	GP.P=P*G.P*s;
	return GP;
}
double GaussianPolynomial<0>::Integral(void) const{return sqrt(PI/scl)*P.coefficients[0];}
double GaussianPolynomial<1>::Integral(void) const{
	GaussianPolynomial<0> G;
	G.scl=scl;
	G.P.coefficients[0]=P.coefficients[0]+P.coefficients[1]*ctr;
	return G.Integral();
}
template<int Degree>
double GaussianPolynomial<Degree>::Integral(void) const{
	if(Degree%2){
		GaussianPolynomial<Degree-1> G=shift(-ctr);
		return G.Integral();
	}
	else{
		GaussianPolynomial<Degree-2> G=shift(-ctr);
		G.P.coefficients[Degree-2]+=(Degree-1)/(2.0*scl)*P.coefficients[Degree];
		return G.Integral();
	}
}
template<int Degree>
void GaussianPolynomial<Degree>::printnl(void) const{
	printf("[ ");
	for(int j=0;j<=Degree;j++){
		printf("%6.4f x^%d ",P.coefficients[j],j);
		if(j<Degree && P.coefficients[j+1]>=0){printf("+");}
	}
	printf(" ] * exp[ %6.4f (x",-scl);
	if(ctr<0){printf("+%6.4f",-ctr);}
	else{printf("%6.4f",-ctr);}
	printf(")^2 ]\n");
}

GaussianPolynomial<0> UnitGaussian(const double& deviation){
	GaussianPolynomial<0> G;
	G.P.coefficients[0]=1;
	G.ctr=0;
	G.scl=1.0/(2.0*deviation*deviation);
	return G/G.Integral();
}
