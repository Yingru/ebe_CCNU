#include"Int.h"
//#define TEST_INT
/******************************************************
*
*	gauss
*
*
* Gauss-Legendre Quadrature w/ switchable no of points 
* 4 Jun 90, es
********************************************************/
#ifdef TEST_INT
double testfunc(double x){
//   return 1/sqrt(100*100-x*x);
     return x*x;
}

double testgalafunc(double x){
   return sin(2*x)*exp(-x);
   //the analytical solution is 0.4
   //The numerical intergral gives 0.399917
}


int main()
{
  double xlo=0, xhi=1.0;
  double invslope=0.25;
  void *a;
  cout<<"test gauss np=20 int="<<gauss(20,testfunc, xlo, xhi, a)<<endl;
  cout<<"test gaula np=15 int="<<gala(15,testgalafunc, xlo, invslope, a)<<endl;
  return 0;
}
#endif



double gauss(int n, double (*f)(double), double xlo, double xhi, void *optvec)
{
        int ndiv= 10;
	double	xoffs, xdiff; 
	int	ix, idiv;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gaulep4; w=gaulew4; break;
		case 8:		p= gaulep8; w=gaulew8; break;
		case 10:	p=gaulep10; w=gaulew10; break;
		case 12:	p=gaulep12; w=gaulew12; break;
		case 16:	p=gaulep16; w=gaulew16; break;
		case 20:	p=gaulep20; w=gaulew20; break;
		case 48:	p=gaulep48; w=gaulew48; break;
		default:	printf("\ngauss():%d points not in list\n",n);
				exit(0);
		}
	s = 0;
	xdiff = 0.5 * ( xhi - xlo ) / ndiv;
	for( idiv=1; idiv<=ndiv; idiv++ ) {
		xoffs = xlo + (2*idiv-1) * xdiff;
		for( ix=0; ix<n/2; ix++ ) 	/* n is even */
			s += w[ix] * ( f(xoffs+xdiff*p[ix]) 
				     + f(xoffs-xdiff*p[ix]) );
		}
	return( s * xdiff );
}


/***************************************************
*
*	gala()
*
* Gauss-Laguerre quadrature extends to +oo
* adapted from old FORTRAN, 24 JUL 90, es, 17 aug 90, 27 Mar 91
****************************************************/

double	gala(int n, double (*f)(double), double xlo, double invslope, void *optvec )
{
	double	*x, *w;
	int	i;
	double	sum	=0;

	if( n == 4 )		{ x = gala4x; w = gala4w; }
	else if( n == 8 )	{ x = gala8x; w = gala8w; }
	else if( n == 12 )	{ x = gala12x;w = gala12w;}
	else if( n == 15 )	{ x = gala15x;w = gala15w;}
	else {
		printf("\ngala():n=%d not in list\n", n );
		return 0;
		}
	for( i=0; i<n; i++ ) 
		sum += w[i] * f( invslope*x[i] + xlo);
	return invslope * sum;	/* make up for transformation */
}
