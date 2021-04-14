#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <nicklib.h>  
#include <getpars.h>  

char *iname = NULL, *oname = NULL ; 

int ispos = YES ;
int symmode = NO ;
int nsmooth = 1 ;

void fsmooth(double *a, int n, int nsmooth)  ;
void readcommands(int argc, char **argv)  ;

#define WVERSION 110 

int main(int argc, char **argv) 
{
  double *x , *y, ya, yb ; 
  double **xx ; 
  double y1, ym ; 
    
 /* Note: y[0] == y[3] for periodic data */

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  const gsl_interp_type *tgsl = gsl_interp_cspline_periodic;
  gsl_spline *spline ; 

  int i, n, jlo, jhi, z, zz, yvlen, t, toff ; 
  double xi, yi, lo, hi ;
  double *yval, *zval, *ww  ;  
  int nsnmooth = 1 ; 


  FILE *ofile ;  

  readcommands(argc, argv) ;

  if (symmode) printf("## symmode set!\n") ;

  n = numlines(iname) ;  
  xx = initarray_2Ddouble(2, n, 0) ; 
  n = getxx(xx, n, 2, iname) ; 
  x = xx[0] ; 
  y = xx[1] ; 
  y1 = y[0] + y[n-1]   ; 
  y[0] = y[n-1] = 0.5*y1 ;  // force periodic boundary
  lo = x[0]; hi = x[n-1] ; 

  jlo = nnint(ceil(lo)) ; 
  jhi = nnint(floor(hi)) ; 
  if (nsmooth > 1) fsmooth(y, n, nsmooth) ; 


  t = MAX(abs(jlo), abs(jhi)) ; 
  toff = 4*t + 100 ; 
  yvlen = jhi-jlo + toff + 50  ; 
  ZALLOC(yval, yvlen, double) ;

  spline = gsl_spline_alloc(tgsl, n) ; 
  gsl_spline_init (spline, x, y, n);

  ofile = stdout ; 
  if (oname != NULL) openit(oname, &ofile, "w") ; 
  fprintf(ofile, "## smooth -i %s\n", iname) ;

  vzero(yval, yvlen) ; 
  zval = yval + toff ;  


  for (z=jlo; z<=jhi; ++z) { 
   xi = z ; 
   yi = gsl_spline_eval (spline, xi, acc);
   if  (ispos) yi = MAX (yi, 0) ;
   zval[z] = yi ;  // smooth but force positive  
  }

  if (symmode) { 
   jlo = MIN(jlo, -jhi) ; 
   jhi = MAX(jhi, -jlo) ; 
  }

  jlo -= 10 ;
  jhi += 10 ;

  for (z=jlo; z<=jhi; ++z) { 
   if (symmode == NO) break ;
   zz  = toff - z  ;
   if (zz<0) fatalx("bad z %d %d\n", z, zz)  ;
   if (zz>=yvlen) fatalx("bad z %d %d\n", z, zz)  ;  
   ya = zval[z] ;  
   yb = zval[-z] ;
   zval[z] = zval[-z] = 0.5*(ya + yb) ; 
  }
   
   

  bal1(yval, yvlen) ; 
  ym = 0 ;
  for (z=jlo; z<=jhi; ++z) { 
   yi = zval[z] ; 
   fprintf(ofile, "%6d %12.6f\n", z, yi) ;
   ym += (double) z * yi ; 
  }
  printf("## mean: %9.3f\n", ym) ;

  fclose(ofile) ;
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(yval) ;

  return 0;
}
void readcommands(int argc, char **argv) 

{
  int i;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "i:o:f:vnis")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'o':
	oname = strdup(optarg) ;
	break;

      case 'f':
	nsmooth = atoi(optarg) ;
	break;

      case 'v':
        printf("version: %s\n", WVERSION) ;
	break;

      case 'n':
        ispos = NO ;
	break;

      case 's':
        symmode = YES ;
	break;

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

}

void fsmooth(double *x, int n, int nsmooth)  
{

 double *filt,  *ww, *ff, *w1 ; 
 int t, k, off ;
  
 t = nsmooth % 2 ; 
 if (t==0) fatalx("-f parameter myst be odd\n") ;

 t = nsmooth/2 ; 
 ZALLOC(ff, 2*nsmooth, double) ; 
 filt = ff+nsmooth ; 
 vclear(filt-t, 1.0, nsmooth) ;
 bal1(filt, nsmooth) ; 

 ZALLOC(ww, n + 2*nsmooth, double) ;
 w1 = ww + nsmooth ; 
 for (k=0; k<n; ++k) { 
  off = k-t ; 
  w1[k] = vdot(filt-t, x+off, nsmooth) ; 
 }
 copyarr(w1, x, n) ;

 free(ww) ;
 free(ff) ;

}
