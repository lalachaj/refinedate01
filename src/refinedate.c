#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 

#define N  10*1000000   
int verbose = NO ; 

#define WVERSION 100

typedef struct  {
 int lo ; 
 int hi ; 
 int len ;
 int *xtime ;
 double *xprob ; 
 double *xprob2 ; 
} TDATA  ; 

TDATA *qa, *qb, *qdis ; 

char *aname, *bname, *disname, *zname = NULL ; 
char *outaname=NULL, *outbname=NULL ; 

void readcommands(int argc, char **argv)  ; 
int readxx(TDATA *qqq,  double **xx, char *name)  ;
double disprob(int dis)  ;
void printq(TDATA *qqq, char *tag)  ;
void printqq(char *zname, double *x, TDATA *qa, TDATA *qb,  int n1, int n2)  ;
void printqf(char *oname, TDATA *qqq, char *hdr, int abswitch) ;

int main(int argc, char **argv) 
{
  double **xx ;  
  int lenxx ; 
  double *w1, *w2, *w1p, *w2p ; 
  double ya, ydis, y ; 
  int i, j, k, n1, n2 ;
  int ta, tb ; 
  double *x2 ; 
  char hdr[1024] ;

  aname = bname = disname = NULL ; 
  ZALLOC(qa, 1, TDATA) ; 
  ZALLOC(qb, 1, TDATA) ; 
  ZALLOC(qdis, 1, TDATA) ; 

  xx = initarray_2Ddouble(2, 10000, -999) ; 
   
  readcommands(argc, argv) ; 
  readxx(qa, xx, aname) ;
  readxx(qb, xx, bname) ;
  readxx(qdis, xx, disname) ;

  ZALLOC(w1, 10000, double) ; 
  ZALLOC(w2, 10000, double) ; 
  ZALLOC(w1p, 10000, double) ; 
  ZALLOC(w2p, 10000, double) ; 

  n1 = qa -> len ; 
  n2 = qb -> len ; 
  ZALLOC(x2, n1*n2, double) ; 
  for (i=0; i<n1; ++i) { 
   ta = qa -> xtime[i] ; 
   ya = qa -> xprob[i] ;  
   vst(w2, qb -> xprob, ya, n2) ;
   for (j=0; j<n2; ++j) { 
    tb = qb -> xtime[j] ; 
    k = ta-tb ;  // a - b date 
    ydis = disprob(k) ; 
// dangerous bend.  make sure we have the sign right here
    w2[j] *= ydis ;    
    y = qa->xprob[i] * qb->xprob[j] * ydis ; 
    x2[i*n2+j] += y ; 
   }
/** 
    vvp(w2p, w2p, w2, n2) ; 
    w1p[i] = asum(w2, n2) ; 
*/
  } 
  vzero(w1p, n1) ;
  vzero(w2p, n2) ;
  for (i=0; i<n1; ++i) { 
   for (j=0; j<n2; ++j) { 
    y  = x2[i*n2+j] ;
    w1p[i] += y ; 
    w2p[j] += y ; 
  }} 

  bal1(w1p, n1) ; 
  bal1(w2p, n2) ;

  copyarr(w1p, qa -> xprob2, n1) ; 
  copyarr(w2p, qb -> xprob2, n2) ; 
  sprintf(hdr, "## %s\n", aname) ; 
//printf("aaa ##time prior posterior\n") ; 
  printqf(outaname, qa, hdr, 1) ; 
  sprintf(hdr, "## %s\n", bname) ; 
  printqf(outbname, qb, hdr, 2) ; 
  bal1(x2, n1*n2) ;  
  printqq(zname, x2, qa, qb, n1, n2) ;
   
  printf("## end of refinedate\n") ; 
  return 0 ; 
}
void printq(TDATA *qqq, char *tag)  
{
 int k, n = qqq -> len ;
 for (k=0; k<n; ++k) { 
  printf("%s ", tag) ; 
  printf("%5d ", qqq -> xtime[k]) ; 
  printf("%12.6f ", qqq -> xprob[k]) ; 
  printf("%12.6f ", qqq -> xprob2[k]) ; 
  printnl() ;
 }
} 
void printqq(char *zname, double *x, TDATA *qa, TDATA *qb,  int n1, int n2) 
{

 int i, j ; 
 double y, thresh ;
 FILE *fff ; 


if (zname == NULL) return ; 
 openit(zname, &fff, "w") ;  
 fprintf(fff, "### joint posterior %s %s\n", aname, bname) ; 

 thresh = 1.0e-6/(double) (n1*n2) ;

for (i=0; i<n1; ++i) { 
 for (j=0; j<n2; ++j) { 
   y = x[i*n2+j] ; 
   if (y<thresh) continue ;
   fprintf(fff, "%5d ", qa -> xtime[i]) ;
   fprintf(fff, "%5d ", qb -> xtime[j]) ;
   fprintf(fff, "%15.9f", y) ; 
   fprintf(fff, "\n") ;
 }} 

 fclose(fff) ;

}
double disprob(int dis) 
{
  int t ; 
  if (dis < qdis -> lo) return 0 ; 
  if (dis > qdis -> hi) return 0 ; 
  t = dis - qdis -> lo ; 
  return qdis -> xprob[t]  ; 
}

int readxx(TDATA *qqq,  double **xx, char *name) 
{
  int n ;  

  n = getxx(xx, 10000, 2, name) ;
  ZALLOC(qqq  -> xtime, n, int) ; 
  ZALLOC(qqq  -> xprob, n, double) ; 
  ZALLOC(qqq  -> xprob2, n, double) ; 

  copyarr(xx[1], qqq -> xprob, n) ; 
  bal1(qqq -> xprob, n) ;
// make prob. sum to 1
  fixit(qqq -> xtime, xx[0], n) ; 

  qqq -> len = n ; 
  qqq -> lo = qqq -> xtime[0] ; 
  qqq -> hi = qqq -> xtime[n-1] ; 

  return n ; 
  
}

void readcommands(int argc, char **argv) 

{
  int i;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "a:b:d:x:y:z:v")) != -1) {

    switch (i)
      {

      case 'a':
	aname = strdup(optarg) ;
	break;

      case 'b':
	bname = strdup(optarg) ;
	break;

      case 'd':
	disname = strdup(optarg) ;
	break;

      case 'x':
	outaname = strdup(optarg) ;
	break;

      case 'y':
	outbname = strdup(optarg) ;
	break;

      case 'z':
	zname = strdup(optarg) ;
        break ; 

      case 'v':
        printf("version: %s\n", WVERSION) ;
	break;

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }
}

void calcmsq(double *x, double *prob, int n, double *pmean, double *psdev) 
{

 double *ww ; 
 double ym, yvar ; 

 ZALLOC(ww, n, double) ; 

 *pmean = ym = vdot(x, prob, n) / asum(prob, n) ; 
 

 vsp(ww, x, -ym, n) ; 
 vvt(ww, ww, ww, n) ;
 yvar = vdot(ww, prob, n) / asum(prob, n) ;

 *psdev = sqrt(yvar + 1.0e-20) ; 

 free(ww) ;
}
void printqf(char *oname, TDATA *qqq, char *hdr, int ab) 
{ 
 char xtag[10] ; 
 FILE *fff ; 
 int k, n = qqq -> len ;
 double xm, xsdev  ; 
 char cc[3] = "AB" ; 
 double *tim ; 
 
 if (ab==1) strcpy(xtag, "aaa") ;
 if (ab==2) strcpy(xtag, "bbb") ;
 if (oname != NULL) strcpy(xtag, "") ;

 ZALLOC(tim, n, double) ; 
 floatit(tim, qqq -> xtime, n) ; 
 calcmsq(tim, qqq -> xprob, n, &xm, &xsdev) ;
 printf("mean and std.dev. for prior (sample %c) :: ", cc[ab-1]) ; 
 printf("%9.1f %9.1f", xm, xsdev) ; 
 printnl() ; 
 calcmsq(tim, qqq -> xprob2, n, &xm, &xsdev) ;
 printf("mean and std.dev. for posterior (sample %c) :: ", cc[ab-1]) ; 
 printf("%9.1f %9.1f", xm, xsdev) ; 
 printnl() ; 


 if (oname == NULL) { 
  printf("%s %s\n", xtag, hdr) ; 
  printq(qqq, xtag) ;
  return ; 
 }

 openit(oname, &fff, "w") ; 
 strcpy(xtag, "") ;
 for (k=0; k<n; ++k) {
  fprintf(fff, "%s ", xtag) ;
  fprintf(fff, "%5d ", qqq -> xtime[k]) ;
  fprintf(fff, "%12.6f ", qqq -> xprob[k]) ;
  fprintf(fff, "%12.6f ", qqq -> xprob2[k]) ;
  fprintf(fff, "\n") ;
 }

 free(tim) ;

 fclose(fff) ;
 return ;
}

