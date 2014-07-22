#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

using std::string;

void nrerror(string error_text)
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text.data());
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}



void gcf(float *gammcf, float a, float x, float *gln)
{
  float gammln(float xx);
  void nrerror(string error_text);
  int i;
  float an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}



void gser(float *gamser, float a, float x, float *gln)
{
  float gammln(float xx);
  void nrerror(string error_text);
  int n;
  float sum,del,ap;
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}



float gammp(float a, float x)
{
  void gcf(float *gammcf, float a, float x, float *gln);
  void gser(float *gamser, float a, float x, float *gln);
  void nrerror(string error_text);
  float gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}


/*#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double SpecialGammln(double xx)
{
assert( xx == 2.0);
return (0.0); // cause (2-1)! = 1 and log 1 = 0
}
void gcf(double *gammcf, double a, double x, double *gln)
{
// returns P(a,x) using continue fracs, also rtn ln $\Gamma(a)$ as *gln 
int i;
double an,b,c,d,del,h;
*gln = SpecialGammln(a);
b = x + 1.0 - a;
c=1.0/FPMIN;
d=1.0/b;
h = d;
assert( x >= a + 1.0);
for ( i=1;i<=ITMAX;i++)
{
an = -i * (i-a);
b += 2.0;
d = an* d + b;
if ( fabs(d) < FPMIN) d=FPMIN;
c = b + an/c;
if ( fabs(c) < FPMIN) c=FPMIN;
d = 1.0/d;
del = d * c;
h *= del;
if ( fabs(del-1.0) < EPS) break;
}
assert ( i <= ITMAX);
*gammcf = exp( -x + a * log(x) -(*gln)) * h; // put factors in front 
fprintf(stdout,"debug b=%le\ndebug c=%le\ndebug d=%le\ndebug h=%le\n",b,c,d,h);
}
main(int argc,char **argv)
{
double gln;
double fx;
double a=2.0;
double x=3;
assert( x >= a+1);
gcf(&fx, 2.0, 3.0, &gln);
fprintf(stdout," Answer P(a=%lf x=%lf)=%lf\n",a,x, 1.0 - fx);
}*/
