/***************************************************************************************************
 * File:              $URL: svn://bolitho.dyndns.org/Research/Poisson/Branches/Weighting/Src/PPolynomial.h $
 * Last Changed By:   $Author: matthew $
 * Revision:				  $Revision: 23 $
 * Last Updated:      $Date: 2006-03-03 15:31:48 -0500 (Fri, 03 Mar 2006) $
 *
 ***************************************************************************************************/

#ifndef GAUSSIAN_INCLUDED
#define GAUSSIAN_INCLUDED

#include "Polynomial.h"

template<int Degree>
class GaussianPolynomial{
public:
	double ctr;
	double scl;
	Polynomial<Degree> P;

	GaussianPolynomial(void);
	template<int Degree2>
	GaussianPolynomial(const GaussianPolynomial<Degree2>& G);
	template<int Degree2>
	GaussianPolynomial& operator  = (const GaussianPolynomial<Degree2> &p);

	double operator() (const double& t) const;
	GaussianPolynomial<Degree+1> derivative(void) const;
	double Integral(void) const;

	GaussianPolynomial scale(const double& t) const;
	GaussianPolynomial shift(const double& t) const;

	GaussianPolynomial operator * (const double& scl) const;
	GaussianPolynomial& operator *= (const double& scl);
	GaussianPolynomial operator / (const double& scl) const;
	GaussianPolynomial& operator /= (const double& scl);

	template<int Degree2>
	GaussianPolynomial<Degree+Degree2>  operator *  (const Polynomial<Degree2>& P) const;

	template<int Degree2>
	GaussianPolynomial<Degree+Degree2>  operator *  (const GaussianPolynomial<Degree2>& G) const;

	void printnl(void) const;

};
GaussianPolynomial<0> UnitGaussian(const double& standardDeviation);
#include "Gaussian.inl"
#endif // GAUSSIAN_INCLUDED
